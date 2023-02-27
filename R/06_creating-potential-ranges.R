## script that creates potential range shapefiles for each species based on their thermal tolerance limits
## note: to run this script successfully, run functions at the end first!
library(tidyverse)
library(sp)
library(sf)
library(raster)
library(ncdf4)
library(rnaturalearth)
library(smoothr)
library(lwgeom)
library(rmapshaper)
library(RColorBrewer)
select <- dplyr::select
sf_use_s2(FALSE)
## read in function to rasterize realized ranges: 
source("R/shp2rast.R") ## function is shp2rast(rast, shp)
source("R/potential-niche-functions.R")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Read in temperature rasters           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read temperature data from script 01
raster_terr_low <- stack("data-processed/temperature-data/terrestrial/raster_terr_low.grd")
raster_terr_high <- stack("data-processed/temperature-data/terrestrial/raster_terr_high.grd")
raster_marine_low <- stack("data-processed/temperature-data/marine/raster_marine_low.grd")
raster_marine_high <- stack("data-processed/temperature-data/marine/raster_marine_high.grd")
raster_intertidal_low <- stack("data-processed/temperature-data/intertidal/raster_intertidal_low.grd")
raster_intertidal_high <- stack("data-processed/temperature-data/intertidal/raster_intertidal_high.grd")
## read in operative temperature from script 05
Te_sun_min <- stack("large-files/operative-temperatures/Te_sun_min.grd")
Te_sun_max <- stack("large-files/operative-temperatures/Te_sun_max.grd")
Te_shade_min <- stack("large-files/operative-temperatures/Te_shade_min.grd")
Te_shade_max <- stack("large-files/operative-temperatures/Te_shade_max.grd")

## read in realm mask layers:
t_mask <- raster("data-processed/masks/raster_terr_mask.grd") 
t_mask_nichemapr <- raster("data-processed/masks/raster_terr_mask_nichemapr.grd") 
m_mask <- raster("data-processed/masks/raster_marine_mask.grd")
i_mask <- raster("data-processed/masks/raster_intertidal_mask.grd")

## read in depth and elevation data:
elevs <- stack("data-processed/masks/elev-and-bath-layers/elevs_minmax.grd") 
depths <-  stack("data-processed/masks/elev-and-bath-layers/depths_minmax.grd") 
shelf_mask <- raster("data-processed/masks/elev-and-bath-layers/raster_200mdepth_mask.grd")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####    Read in thermal tolerance and trait data        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in thermal limits:
thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = "_"))

## read in acclimation data:
acc_data <- read.csv("data-processed/traits/acclimation-data.csv")

## read in species traits:
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = "_"))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####           Rasterize realized ranges                #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in realized ranges:
realized_ranges <- st_read("large-files/realized-ranges/realized-ranges.shp") %>%
  mutate(range_id = paste(species, source, sep = "_")) 

r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1, vals = 1)

i=1
while (i < nrow(realized_ranges)+1) {
  
  ## get speices' realized range:
  range <- realized_ranges[i, ]
  
  ## rasterize realized range:
  sp_range <- as_Spatial(range) 
  rr_raster <- shp2rast(shp = sp_range, rast = r) 
  rr_og <- rr_raster # store unrestricted version of rasterized range for comparison later
  
  ## get the species' traits:
  sp <- str_replace_all(traits$genus_species, "_", " ")
  sp_traits <- traits[which(sp == range$species), ]
  
  ## constrain realized range rasters by:
  ## - raster cells that were permitted to be in the potential range (in the realm)
  ## - elevational range for terrestrial species
  ## - depth distribution for marine and intertidal species
  if(range$realm == "Terrestrial") {
    rr_raster <- mask(rr_raster, t_mask) ## get rid of cells in realized range outside of the realm
    
    ## get elevational range limits 
    elev_max <- sp_traits$upper_elevation_limit_meters
    elev_min <- sp_traits$lower_elevation_limit_meters
    
    ## if both limits are not NA, restrict by both limits
    if (!is.na(elev_min) & !is.na(elev_max)) {
      rr_raster[elevs$elev_max < elev_min] <- NA 
      rr_raster[elevs$elev_min > elev_max] <- NA 
      ## exclude cells with max elev less than species min elevation
      ## exclude cells with min elev greater than species max elevation
    }
    else if (!is.na(elev_min)) {
      rr_raster[elevs$elev_max < elev_min] <- NA 
      ## exclude cells with max elev less than species min elevation
    }
    else if (!is.na(elev_max)) {
      rr_raster[elevs$elev_min > elev_max] <- NA 
      ## exclude cells with min elev greater than species max elevation
    }
  }
  else if(range$realm %in% c("Marine", "Intertidal")) {
    
    ## get rid of cells in realized range outside of the realm
    if (range$realm == "Marine") {
      rr_raster <- mask(rr_raster, m_mask) 
    }
    else if (range$realm == "Intertidal") {
      rr_raster <- mask(rr_raster, m_mask) 
    }
    
    ## get depth distribution limits 
    depth_upper <- sp_traits$upper_depth_limit
    depth_lower <- sp_traits$lower_depth_limit
    
    ## if both limits are not NA, restrict by both depth limits
    if (!is.na(depth_lower) & !is.na(depth_upper)) {
      rr_raster[depths$depth_min > depth_upper] <- NA 
      rr_raster[depths$depth_max < depth_lower] <- NA 
      ## exclude cells with deepest depth above species upper limit  
      ## exclude cells with shallowest depth less than species lower limit 
    }
    else if (!is.na(depth_upper)) {
      rr_raster[depths$depth_min > depth_upper] <- NA 
      ## exclude cells with deepest depth above species upper limit  
    }
    else if (!is.na(depth_lower)) {
      rr_raster[depths$depth_max < depth_lower] <- NA 
      ## exclude cells with shallowest depth less than species lower limit 
    }
    ## if lower depth limit is unknown, restrict to continental shelf if it is a coastal spp
    else if (is.na(depth_lower) & !is.na(sp_traits$coastal_or_oceanic)) {
      is_coastal <- sp_traits$coastal_or_oceanic == "coastal"
      
      if(is_coastal) {
        rr_raster[is.na(shelf_mask)] <- NA
      }
    }
  }
  
  ## add to list of rasters
  if (i == 1) {
    rasterized_rrs <- rr_raster
    og_rasterized_rrs <- rr_og
  }
  else {
    rasterized_rrs <- addLayer(rasterized_rrs, rr_raster, updatevalue = NA)
    og_rasterized_rrs <- addLayer(og_rasterized_rrs, rr_og, updatevalue = NA)
  }
  
  print(paste("Finished range number:", i))
  i = i + 1
}

names(rasterized_rrs) <- realized_ranges$range_id
names(og_rasterized_rrs) <- realized_ranges$range_id

##saveRDS(rasterized_rrs, "data-processed/realized-ranges/rasterized_rrs.rds")
rasterized_rrs <- readRDS("data-processed/realized-ranges/rasterized_rrs.rds")
##saveRDS(og_rasterized_rrs, "data-processed/realized-ranges/og_rasterized_rrs.rds")
og_rasterized_rrs <- readRDS("data-processed/realized-ranges/og_rasterized_rrs.rds")

rasterized_rrs_nichemapr <- mask(rasterized_rrs, t_mask_nichemapr) 
##saveRDS(rasterized_rrs_nichemapr, "data-processed/realized-ranges/rasterized_rrs_nichemapr.rds")
rasterized_rrs_nichemapr <- readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Create potential range shapefiles          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## partition thermal limits into upper and lower limits
upper_limits <- thermal_limits %>%
  filter(type == "max")
lower_limits <- thermal_limits %>%
  filter(type == "min") 

## keep limits of species with both upper and lower thermal limits
both_upper_all <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower_all <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]

## keep limits of species with only an upper or lower thermal limit
only_upper_all <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower_all <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Terrestrial species          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## get rid of thermal limits of species in other realms:
both_upper <- filter(both_upper_all, realm == "Terrestrial")
both_lower <- filter(both_lower_all, realm == "Terrestrial")
only_upper <- filter(only_upper_all, realm == "Terrestrial")
only_lower <-  filter(only_lower_all, realm == "Terrestrial")

## prepare zoogeographic realms for use in restricting potential ranges:
zoo <- read_sf("data-raw/polygons/CMEC regions & realms/newRealms.shp")

# simplify the object to make it faster to use
sf::sf_use_s2(FALSE)
zoo <- zoo %>% 
  st_simplify(dTolerance = 0.01) %>% 
  group_by(Realm) %>% 
  dplyr::summarise()
plot(zoo)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## SPECIES WITH BOTH THERMAL LIMITS: ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## CREATING POTENTIAL THERMAL RANGES:
## filter out habitat with temperatures outside of the fundamental thermal niche 
Te <- filter_by_tolerance_Te(both_upper = both_upper,
                             both_lower = both_lower,
                             Te_max = Te_shade_max, 
                             Te_min = Te_shade_min)
Te_in_sun <- filter_by_tolerance_Te(both_upper = both_upper, 
                                    both_lower = both_lower,
                                    Te_max = Te_sun_max, 
                                    Te_min = Te_shade_min)
acclimatized <- filter_by_tolerance_acclimatized(both_upper = both_upper,
                                             both_lower = both_lower,
                                             raster_high = Te_shade_max, 
                                             raster_low = Te_shade_min,
                                             acc_data = acc_data, 
                                             realm = "Terrestrial")

## restrict range to:
## - zoogeographic realms that overlap species' realized range
## - habitat within the species elevational range 
prs_terrestrial_Te <- create_potential_ranges(clumped_temps = zoo, 
                                              realized_ranges = realized_ranges,
                                              filtered = Te,
                                              type = 'Te')
prs_terrestrial_Te_in_sun <- create_potential_ranges(clumped_temps = zoo, 
                                                     realized_ranges = realized_ranges,
                                                     filtered = Te_in_sun,
                                                     type = 'Te_in_sun')
prs_terrestrial_acclimatized <- create_potential_ranges(clumped_temps = zoo, 
                                                      realized_ranges = realized_ranges,
                                                      filtered = acclimatized,
                                                      type = 'acclimatized')

# saveRDS(prs_terrestrial_Te, "data-processed/potential-ranges/terrestrial/prs_terrestrial_Te.rds")
# saveRDS(prs_terrestrial_Te_in_sun, "data-processed/potential-ranges/terrestrial/prs_terrestrial_Te_in_sun.rds")
# saveRDS(prs_terrestrial_acclimatized, "data-processed/potential-ranges/terrestrial/prs_terrestrial_acclimatized.rds")
prs_terrestrial_Te  <- readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_Te.rds")
prs_terrestrial_Te_in_sun <- readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_Te_in_sun.rds")
prs_terrestrial_acclimatized <- readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_acclimatized.rds")


## EXTRACTING TEMPS IN THERMAL NICHES:
niche_Te_terr <- extract_thermal_niche(type = "Te", 
                                            prs = prs_terrestrial_Te,
                                            realm = 'Terrestrial')
niche_Te_in_sun_terr <- extract_thermal_niche(type = "Te_in_sun", 
                                              prs = prs_terrestrial_Te_in_sun,
                                              realm = 'Terrestrial')
niche_acclimatized_terr <- extract_thermal_niche(type = "acclimatized", 
                                            prs = prs_terrestrial_acclimatized,
                                            realm = 'Terrestrial')

## save:
write.csv(niche_Te_terr, "data-processed/thermal-niches/terrestrial/niche_Te_terr.csv", 
          row.names = FALSE)
write.csv(niche_Te_in_sun_terr, "data-processed/thermal-niches/terrestrial/niche_Te_in_sun_terr.csv", 
          row.names = FALSE)
write.csv(niche_acclimatized_terr, "data-processed/thermal-niches/terrestrial/niche_acclimatized_terr.csv", 
          row.names = FALSE)

niche_Te_terr <- read.csv("data-processed/thermal-niches/terrestrial/niche_Te_terr.csv")
niche_Te_in_sun_terr <- read.csv("data-processed/thermal-niches/terrestrial/niche_Te_in_sun_terr.csv")
niche_acclimatized_terr <- read.csv("data-processed/thermal-niches/terrestrial/niche_acclimatized_terr.csv")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  SPECIES WITH ONE THERMAL LIMIT:  ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## CREATING POTENTIAL THERMAL RANGES:
## filter out habitat with temperatures that exceed the species' single fundamental thermal niche limit
upper_Te <- filter_by_tolerance_Te_one_limit(lims = only_upper,
                                                      limit_type = "tmax",
                                                      Te_max = Te_shade_max, 
                                                      Te_min = Te_shade_min)

lower_Te <- filter_by_tolerance_Te_one_limit(lims = only_lower,
                                                      limit_type = "tmin",
                                                      Te_max = Te_shade_max, 
                                                      Te_min = Te_shade_min)

upper_Te_in_sun <- filter_by_tolerance_Te_one_limit(lims = only_upper,
                                                    limit_type = "tmax",
                                                    Te_max = Te_sun_max, 
                                                    Te_min = Te_shade_min)

upper_acclimatized <- filter_by_tolerance_acclimatized_one_limit(lims = only_upper,
                                                     limit_type = "tmax",
                                                     raster_high = Te_shade_max, 
                                                     raster_low = Te_shade_min,
                                                     acc_data = acc_data, 
                                                     realm = "Terrestrial")

lower_acclimatized <- filter_by_tolerance_acclimatized_one_limit(lims = only_lower,
                                                     limit_type = "tmin",
                                                     raster_high = Te_shade_max, 
                                                     raster_low = Te_shade_min,
                                                     acc_data = acc_data, 
                                                     realm = "Terrestrial")

## restrict range to:
## - zoogeographic realms that overlap species' realized range
## - habitat within the species elevational range
prs_terrestrial_upper_Te <- create_potential_ranges_one_limit(clumped_temps = zoo,
                                           realized_ranges = realized_ranges,
                                           filtered = upper_Te,
                                           type = 'Te')
prs_terrestrial_lower_Te <- create_potential_ranges_one_limit(clumped_temps = zoo,
                                                       realized_ranges = realized_ranges,
                                                       filtered = lower_Te,
                                                       type = 'Te')
prs_terrestrial_upper_Te_in_sun <- create_potential_ranges_one_limit(clumped_temps = zoo,
                                                                     realized_ranges = realized_ranges,
                                                                     filtered = upper_Te_in_sun,
                                                                     type = 'Te_in_sun')
prs_terrestrial_upper_acclimatized <- create_potential_ranges_one_limit(clumped_temps = zoo,
                                                                 realized_ranges = realized_ranges,
                                                                 filtered = upper_acclimatized,
                                                                 type = 'acclimatized')
prs_terrestrial_lower_acclimatized <- create_potential_ranges_one_limit(clumped_temps = zoo,
                                                                 realized_ranges = realized_ranges,
                                                                 filtered = lower_acclimatized,
                                                                 type = 'acclimatized')

# saveRDS(prs_terrestrial_upper_Te, "data-processed/potential-ranges/terrestrial/prs_terrestrial_upper_Te.rds")
# saveRDS(prs_terrestrial_lower_Te, "data-processed/potential-ranges/terrestrial/prs_terrestrial_lower_Te.rds")
# saveRDS(prs_terrestrial_upper_Te_in_sun, "data-processed/potential-ranges/prs_terrestrial_upper_Te_in_sun.rds")
# saveRDS(prs_terrestrial_upper_acclimatized, "data-processed/potential-ranges/terrestrial/prs_terrestrial_upper_acclimatized.rds")
# saveRDS(prs_terrestrial_lower_acclimatized, "data-processed/potential-ranges/terrestrial/prs_terrestrial_lower_acclimatized.rds")
prs_terrestrial_upper_Te <- readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_upper_Te.rds")
prs_terrestrial_lower_Te <- readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_lower_Te.rds")
prs_terrestrial_upper_Te_in_sun <- readRDS("data-processed/potential-ranges/prs_terrestrial_upper_Te_in_sun.rds")
prs_terrestrial_upper_acclimatized <- readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_upper_acclimatized.rds")
prs_terrestrial_lower_acclimatized <- readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_lower_acclimatized.rds")

## EXTRACTING TEMPS IN THERMAL NICHES:
niche_upper_Te_terr <- extract_thermal_niche_one_limit(type = "Te", 
                                       prs = prs_terrestrial_upper_Te,
                                       realm = 'Terrestrial', 
                                       limit_type = "tmax")
niche_lower_Te_terr <- extract_thermal_niche_one_limit(type = "Te", 
                                                 prs = prs_terrestrial_lower_Te,
                                                 realm = 'Terrestrial',
                                             limit_type = "tmin")
niche_upper_Te_in_sun_terr <- extract_thermal_niche_one_limit(type = "Te_in_sun", 
                                                              prs = prs_terrestrial_upper_Te_in_sun,
                                                              realm = 'Terrestrial',
                                                              limit_type = "tmax")
niche_upper_acclimatized_terr <- extract_thermal_niche_one_limit(type = "acclimatized", 
                                                 prs = prs_terrestrial_upper_acclimatized,
                                                 realm = 'Terrestrial', 
                                                 limit_type = "tmax")
niche_lower_acclimatized_terr <- extract_thermal_niche_one_limit(type = "acclimatized", 
                                                       prs = prs_terrestrial_lower_acclimatized,
                                                       realm = 'Terrestrial', 
                                                       limit_type = "tmin")

## save:
write.csv(niche_upper_Te_terr, "data-processed/thermal-niches/terrestrial/niche_upper_Te_terr.csv", 
          row.names = FALSE)
write.csv(niche_lower_Te_terr, "data-processed/thermal-niches/terrestrial/niche_lower_Te_terr.csv", 
          row.names = FALSE)
write.csv(niche_upper_Te_in_sun_terr, "data-processed/thermal-niches/terrestrial/niche_upper_Te_in_sun_terr.csv", 
          row.names = FALSE)
write.csv(niche_upper_acclimatized_terr, "data-processed/thermal-niches/terrestrial/niche_upper_acclimatized_terr.csv", 
          row.names = FALSE)
write.csv(niche_lower_acclimatized_terr, "data-processed/thermal-niches/terrestrial/niche_lower_acclimatized_terr.csv", 
          row.names = FALSE)

niche_upper_Te_terr <- read.csv("data-processed/thermal-niches/terrestrial/niche_upper_Te_terr.csv")
niche_lower_Te_terr <- read.csv("data-processed/thermal-niches/terrestrial/niche_lower_Te_terr.csv")
niche_upper_Te_in_sun_terr <- read.csv("data-processed/thermal-niches/terrestrial/niche_upper_Te_in_sun_terr.csv")
niche_upper_acclimatized_terr <- read.csv("data-processed/thermal-niches/terrestrial/niche_upper_acclimatized_terr.csv")
niche_lower_acclimatized_terr <- read.csv("data-processed/thermal-niches/terrestrial/niche_lower_acclimatized_terr.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Marine species          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## get rid of thermal limits of species in other realms:
both_upper <- filter(both_upper_all, realm == "Marine")
both_lower <- filter(both_lower_all, realm == "Marine")
only_upper <- filter(only_upper_all, realm == "Marine")
only_lower <-  filter(only_lower_all, realm == "Marine")

## clump contiguous habitat together:
clumped_temps <- raster_marine_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## SPECIES WITH BOTH THERMAL LIMITS: ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## CREATING POTENTIAL THERMAL RANGES:
## filter out habitat with temperatures outside of the fundamental thermal niche 
reg <- filter_by_tolerance(both_upper = both_upper,
                                    both_lower = both_lower,
                                    raster_high = raster_marine_high, 
                                    raster_low = raster_marine_low)
acclimatized <- filter_by_tolerance_acclimatized(both_upper = both_upper,
                                                      both_lower = both_lower,
                                                      raster_high = raster_marine_high, 
                                                      raster_low = raster_marine_low,
                                                      acc_data = acc_data, 
                                                      realm = "Marine")

## restrict range to:
## - contiguous marine habitat that overlaps species' realized range
## - habitat within the species' depth distribution
prs_marine <- create_potential_ranges(clumped_temps = clumped_temps, 
                                      realized_ranges = realized_ranges, filtered = reg,
                                      type = 'reg')
prs_marine_acclimatized <- create_potential_ranges(clumped_temps = clumped_temps, 
                                                 realized_ranges = realized_ranges, 
                                                 filtered = acclimatized,
                                                 type = 'acclimatized')

# saveRDS(prs_marine, "data-processed/potential-ranges/marine/prs_marine.rds")
# saveRDS(prs_marine_acclimatized, "data-processed/potential-ranges/marine/prs_marine_acclimatized.rds")
prs_marine <- readRDS("data-processed/potential-ranges/marine/prs_marine.rds")
prs_marine_acclimatized <- readRDS("data-processed/potential-ranges/marine/prs_marine_acclimatized.rds")

## EXTRACTING TEMPS IN THERMAL NICHES:
niche_reg_marine <- extract_thermal_niche(type = "reg", prs_marine, realm = 'Marine')
niche_acclimatized_marine <- extract_thermal_niche(type = "acclimatized", prs_marine_acclimatized, realm = 'Marine')

write.csv(niche_reg_marine, "data-processed/thermal-niches/marine/niche_reg_marine.csv", row.names = FALSE)
write.csv(niche_acclimatized_marine, "data-processed/thermal-niches/marine/niche_acclimatized_marine.csv", row.names = FALSE)

niche_reg_marine = read.csv("data-processed/thermal-niches/marine/niche_reg_marine.csv")
niche_acclimatized_marine = read.csv("data-processed/thermal-niches/marine/niche_acclimatized_marine.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  SPECIES WITH ONE THERMAL LIMIT:  ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## CREATING POTENTIAL THERMAL RANGES:
## filter out habitat with temperatures outside of the fundamental thermal niche 
upper_reg <- filter_by_tolerance_one_limit(limit_type = 'tmax', 
                                           lims = only_upper,
                           raster_high = raster_marine_high, 
                           raster_low = raster_marine_low)
lower_reg <- filter_by_tolerance_one_limit(limit_type = 'tmin',
                                           lims = only_lower,                                  
                                           raster_high = raster_marine_high, 
                                           raster_low = raster_marine_low)
upper_acclimatized <- filter_by_tolerance_acclimatized_one_limit(limit_type = 'tmax', 
                                           lims = only_upper,
                                           raster_high = raster_marine_high, 
                                           raster_low = raster_marine_low,
                                           acc_data = acc_data, 
                                           realm = "Marine")
lower_acclimatized <- filter_by_tolerance_acclimatized_one_limit(limit_type = 'tmin',
                                           lims = only_lower,                                  
                                           raster_high = raster_marine_high, 
                                           raster_low = raster_marine_low,
                                           acc_data = acc_data, 
                                           realm = "Marine")

## restrict range to:
## - contiguous marine habitat that overlaps species' realized range
## - habitat within the species' depth distribution
prs_marine_upper_reg <- create_potential_ranges_one_limit(clumped_temps = clumped_temps,
                                                                 realized_ranges = realized_ranges,
                                                                 filtered = upper_reg,
                                                                 type = 'reg')
prs_marine_lower_reg <- create_potential_ranges_one_limit(clumped_temps = clumped_temps,
                                                                 realized_ranges = realized_ranges,
                                                                 filtered = lower_reg,
                                                                 type = 'reg')
prs_marine_upper_acclimatized <- create_potential_ranges_one_limit(clumped_temps = clumped_temps,
                                                                   realized_ranges = realized_ranges,
                                                                   filtered = upper_acclimatized,
                                                                   type = 'acclimatized')
prs_marine_lower_acclimatized <- create_potential_ranges_one_limit(clumped_temps = clumped_temps,
                                                                   realized_ranges = realized_ranges,
                                                                   filtered = lower_acclimatized,
                                                                   type = 'acclimatized')

# saveRDS(prs_marine_upper_reg, "data-processed/potential-ranges/marine/prs_marine_upper_reg.rds")
# saveRDS(prs_marine_lower_reg, "data-processed/potential-ranges/marine/prs_marine_lower_reg.rds")
# saveRDS(prs_marine_upper_acclimatized, "data-processed/potential-ranges/marine/prs_marine_upper_reg_acclimatized.rds")
# saveRDS(prs_marine_lower_acclimatized, "data-processed/potential-ranges/marine/prs_marine_lower_reg_acclimatized.rds")
prs_marine_upper_reg <- readRDS("data-processed/potential-ranges/marine/prs_marine_upper_reg.rds")
prs_marine_lower_reg <- readRDS("data-processed/potential-ranges/marine/prs_marine_lower_reg.rds")
prs_marine_upper_acclimatized <- readRDS("data-processed/potential-ranges/marine/prs_marine_upper_reg_acclimatized.rds")
prs_marine_lower_acclimatized <- readRDS("data-processed/potential-ranges/marine/prs_marine_lower_reg_acclimatized.rds")

## EXTRACTING TEMPS IN THERMAL NICHES:
niche_upper_reg_marine <- extract_thermal_niche_one_limit(type = "reg", 
                                             prs = prs_marine_upper_reg,
                                             realm = 'Marine',
                                             limit_type = "tmax")
niche_lower_reg_marine <- extract_thermal_niche_one_limit(type = "reg", 
                                             prs = prs_marine_lower_reg,
                                             realm = 'Marine',
                                             limit_type = "tmin")
niche_upper_acclimatized_marine <- extract_thermal_niche_one_limit(type = "acclimatized", 
                                                       prs = prs_marine_upper_acclimatized,
                                                       realm = 'Marine',
                                                       limit_type = "tmax")
niche_lower_acclimatized_marine <- extract_thermal_niche_one_limit(type = "acclimatized", 
                                                       prs = prs_marine_lower_acclimatized,
                                                       realm = 'Marine',
                                                       limit_type = "tmin")

## save:
write.csv(niche_upper_reg_marine, "data-processed/thermal-niches/marine/niche_upper_reg_marine.csv", 
          row.names = FALSE)
write.csv(niche_lower_reg_marine, "data-processed/thermal-niches/marine/niche_lower_reg_marine.csv", 
          row.names = FALSE)
write.csv(niche_upper_acclimatized_marine, "large-files/thermal-niches/marine/niche_upper_acclimatized_marine.csv", 
          row.names = FALSE)
write.csv(niche_lower_acclimatized_marine, "data-processed/thermal-niches/marine/niche_lower_acclimatized_marine.csv", 
          row.names = FALSE)

niche_upper_reg_marine <- read.csv("data-processed/thermal-niches/marine/niche_upper_reg_marine.csv")
niche_lower_reg_marine <- read.csv("data-processed/thermal-niches/marine/niche_lower_reg_marine.csv")
niche_upper_acclimatized_marine <- read.csv("large-files/thermal-niches/marine/niche_upper_acclimatized_marine.csv")
niche_lower_acclimatized_marine <- read.csv("data-processed/thermal-niches/marine/niche_lower_acclimatized_marine.csv")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Intertidal species           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## get rid of thermal limits of species in other realms:
both_upper <- filter(both_upper_all, realm == "Intertidal")
both_lower <- filter(both_lower_all, realm == "Intertidal")
only_upper <- filter(only_upper_all, realm == "Intertidal")
only_lower <-  filter(only_lower_all, realm == "Intertidal")

## clump contiguous habitat together:
clumped_temps <- raster_intertidal_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps$geometry)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## SPECIES WITH BOTH THERMAL LIMITS: ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## CREATING POTENTIAL THERMAL RANGES:
## filter out habitat with temperatures outside of the fundamental thermal niche 
reg <- filter_by_tolerance(both_upper = both_upper,
                                    both_lower = both_lower,
                                    raster_high = raster_intertidal_high,
                                    raster_low = raster_intertidal_low)
acclimatized <- filter_by_tolerance_acclimatized(both_upper = both_upper,
                                                      both_lower = both_lower,
                                                      raster_high = raster_intertidal_high, 
                                                      raster_low = raster_intertidal_low,
                                                      acc_data = acc_data, 
                                                      realm = "Intertidal")

## restrict range to:
## - contiguous intertidal habitat (continental shelf beside a land mass) that overlaps species' realized range
## - habitat within the species' depth distribution
prs_intertidal <- create_potential_ranges(clumped_temps = clumped_temps, 
                                          realized_ranges = realized_ranges, filtered = reg,
                                          type = 'reg')
prs_intertidal_acclimatized <- create_potential_ranges(clumped_temps = clumped_temps, 
                                                     realized_ranges = realized_ranges, 
                                                     filtered = acclimatized,
                                                     type = 'acclimatized')

# saveRDS(prs_intertidal, "data-processed/potential-ranges/intertidal/prs_intertidal.rds")
# saveRDS(prs_intertidal_acclimatized, "data-processed/potential-ranges/intertidal/prs_intertidal_acclimatized.rds")
prs_intertidal <- readRDS("data-processed/potential-ranges/intertidal/prs_intertidal.rds")
prs_intertidal_acclimatized <- readRDS("data-processed/potential-ranges/intertidal/prs_intertidal_acclimatized.rds")


## EXTRACTING TEMPS IN THERMAL NICHES:
niche_reg_intertidal <- extract_thermal_niche(type = "reg", prs_intertidal, realm = 'Intertidal')
niche_acclimatized_intertidal <- extract_thermal_niche(type = "acclimatized", prs_intertidal_acclimatized, realm = 'Intertidal')

write.csv(niche_reg_intertidal, "data-processed/thermal-niches/intertidal/niche_reg_intertidal.csv", row.names = FALSE)
write.csv(niche_acclimatized_intertidal, "data-processed/thermal-niches/intertidal/niche_acclimatized_intertidal.csv", row.names = FALSE)

niche_reg_intertidal <- read.csv("data-processed/thermal-niches/intertidal/niche_reg_intertidal.csv")
niche_acclimatized_intertidal <- read.csv("data-processed/thermal-niches/intertidal/niche_acclimatized_intertidal.csv")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  SPECIES WITH ONE THERMAL LIMIT:  ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## CREATING POTENTIAL THERMAL RANGES:
## filter out habitat with temperatures outside of the fundamental thermal niche 
## note: no species had only lower thermal limit, hence functions not applied to lower limit
upper_reg <- filter_by_tolerance_one_limit(limit_type = 'tmax', 
                                           lims = only_upper,
                                           raster_high = raster_marine_high, 
                                           raster_low = raster_marine_low)
upper_acclimatized <- filter_by_tolerance_acclimatized_one_limit(limit_type = 'tmax', 
                                                                 lims = only_upper,
                                                                 raster_high = raster_marine_high, 
                                                                 raster_low = raster_marine_low,
                                                                 acc_data = acc_data, 
                                                                 realm = "Intertidal")

## restrict range to:
## - contiguous intertidal habitat (continental shelf beside a land mass) that overlaps species' realized range
## - habitat within the species' depth distribution
prs_intertidal_upper_reg <- create_potential_ranges_one_limit(clumped_temps = clumped_temps,
                                                          realized_ranges = realized_ranges,
                                                          filtered = upper_reg,
                                                          type = 'reg')
prs_intertidal_upper_acclimatized <- create_potential_ranges_one_limit(clumped_temps =
                                                                                  clumped_temps,
                                                                                realized_ranges = realized_ranges,
                                                                                filtered = upper_acclimatized,
                                                                                type = 'acclimatized')

# saveRDS(prs_intertidal_upper_reg, "data-processed/potential-ranges/intertidal/prs_intertidal_upper_reg.rds")
# saveRDS(prs_intertidal_upper_acclimatized, "data-processed/potential-ranges/intertidal/prs_intertidal_upper_reg_acclimatized.rds")
prs_intertidal_upper_reg <- readRDS("data-processed/potential-ranges/intertidal/prs_intertidal_upper_reg.rds")
prs_intertidal_upper_acclimatized <- readRDS("data-processed/potential-ranges/intertidal/prs_intertidal_upper_reg_acclimatized.rds")

## EXTRACTING TEMPS IN THERMAL NICHES:
niche_upper_reg_intertidal <- extract_thermal_niche_one_limit(type = "reg", 
                                                prs = prs_intertidal_upper_reg,
                                                realm = 'Intertidal', limit_type = "tmax")
niche_upper_acclimatized_intertidal <- extract_thermal_niche_one_limit(type = "acclimatized", 
                                                         prs = prs_intertidal_upper_acclimatized,
                                                         realm = 'Intertidal',
                                                         limit_type = "tmax")

## save:
write.csv(niche_upper_reg_intertidal, "data-processed/thermal-niches/intertidal/niche_upper_reg_intertidal.csv", 
          row.names = FALSE)
write.csv(niche_upper_acclimatized_intertidal, "data-processed/thermal-niches/intertidal/niche_upper_acclimatized_intertidal.csv", 
          row.names = FALSE)

niche_upper_reg_intertidal <- read.csv("data-processed/thermal-niches/intertidal/niche_upper_reg_intertidal.csv")
niche_upper_acclimatized_intertidal <- read.csv("data-processed/thermal-niches/intertidal/niche_upper_acclimatized_intertidal.csv")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####     Collate potential ranges across realms      #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## regular:
potential_ranges_both_limits_Te <- stack(prs_terrestrial_Te[[1]],
                                                 prs_marine[[1]],
                                                 prs_intertidal[[1]])
potential_ranges_upper_limits_Te <- stack(prs_terrestrial_upper_Te[[1]],
                                         prs_marine_upper_reg[[1]],
                                         prs_intertidal_upper_reg[[1]])
potential_ranges_lower_limits_Te <- stack(prs_terrestrial_lower_Te[[1]],
                                         prs_marine_lower_reg[[1]])

potential_ranges_filling_script <- list(stack(prs_terrestrial_Te[[1]],
                                              prs_marine[[1]],
                                              prs_intertidal[[1]]), 
                                        stack(prs_terrestrial_Te[[2]],
                                              prs_marine[[2]],
                                              prs_intertidal[[2]]),
                                        stack(prs_terrestrial_Te[[3]],
                                              prs_marine[[3]],
                                              prs_intertidal[[3]]),
                                        stack(prs_terrestrial_Te[[4]],
                                              prs_marine[[4]],
                                              prs_intertidal[[4]]),
                                        stack(prs_terrestrial_Te[[5]],
                                              prs_marine[[5]],
                                              prs_intertidal[[5]]))

saveRDS(potential_ranges_filling_script, 'data-processed/potential-ranges/all-realms/potential_ranges_filling_script.rds')

## sun:
potential_ranges_both_limits_Te_in_sun <- stack(prs_terrestrial_Te_in_sun[[1]],
                                                prs_marine[[1]],
                                                prs_intertidal[[1]])

## acclimatized:
potential_ranges_both_limits_acclimatized <- stack(prs_terrestrial_acclimatized[[1]],
                                                 prs_marine_acclimatized[[1]],
                                                 prs_intertidal_acclimatized[[1]])
potential_ranges_upper_limits_acclimatized <- stack(prs_terrestrial_upper_acclimatized[[1]],
                                                   prs_marine_upper_acclimatized[[1]],
                                                   prs_intertidal_upper_acclimatized[[1]])

potential_ranges_lower_limits_acclimatized <- stack(prs_terrestrial_lower_acclimatized[[1]],
                                                   prs_marine_lower_acclimatized[[1]])

potential_ranges_filling_script_acc <- list(stack(prs_terrestrial_acclimatized[[1]],
                                              prs_marine_acclimatized[[1]],
                                              prs_intertidal_acclimatized[[1]]), 
                                        stack(prs_terrestrial_acclimatized[[2]],
                                              prs_marine_acclimatized[[2]],
                                              prs_intertidal_acclimatized[[2]]),
                                        stack(prs_terrestrial_acclimatized[[3]],
                                              prs_marine_acclimatized[[3]],
                                              prs_intertidal_acclimatized[[3]]),
                                        stack(prs_terrestrial_acclimatized[[4]],
                                              prs_marine_acclimatized[[4]],
                                              prs_intertidal_acclimatized[[4]]),
                                        stack(prs_terrestrial_acclimatized[[5]],
                                              prs_marine_acclimatized[[5]],
                                              prs_intertidal_acclimatized[[5]]))

saveRDS(potential_ranges_filling_script_acc, 'data-processed/potential-ranges/all-realms/potential_ranges_filling_script_acc.rds')
saveRDS(potential_ranges_both_limits_Te, "data-processed/potential-ranges/all-realms/potential_ranges_Te.rds")
saveRDS(potential_ranges_lower_limits_Te, "data-processed/potential-ranges/all-realms/potential_ranges_lower_Te.rds")
saveRDS(potential_ranges_upper_limits_Te, "data-processed/potential-ranges/all-realms/potential_ranges_upper_Te.rds")
saveRDS(potential_ranges_both_limits_acclimatized, "data-processed/potential-ranges/all-realms/potential_ranges_acclimatized.rds")
saveRDS(potential_ranges_both_limits_Te_in_sun, "data-processed/potential-ranges/all-realms/potential_ranges_Te_in_sun.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### Collate temps in thermal niche across species with one + both limits #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## make so that niche_Te, niche_Te_in_sun, etc. have data from species with one limit and both limits
## add column "limit_type"
niche_Te_terr$limit_type  <- "tmax, tmin"
niche_Te_in_sun_terr$limit_type  <- "tmax, tmin"
niche_acclimatized_terr$limit_type  <- "tmax, tmin"
niche_upper_Te_terr$limit_type  <- "tmax"
niche_lower_Te_terr$limit_type  <- "tmin"
niche_upper_acclimatized_terr$limit_type  <- "tmax"
niche_lower_acclimatized_terr$limit_type <- "tmin"
niche_Te_in_sun_terr$limit_type <- 'tmax, tmin'
niche_upper_Te_in_sun_terr$limit_type <- 'tmax'

niche_upper_acclimatized_marine$limit_type = "tmax"
niche_lower_acclimatized_marine$limit_type = "tmin"
niche_upper_reg_marine$limit_type = "tmax"
niche_lower_reg_marine$limit_type = "tmin"
niche_reg_marine$limit_type = "tmax, tmin"
niche_acclimatized_marine$limit_type = "tmax, tmin"

niche_upper_acclimatized_intertidal$limit_type = "tmax"
niche_upper_reg_intertidal$limit_type = "tmax"
niche_reg_intertidal$limit_type = "tmax, tmin"
niche_acclimatized_intertidal$limit_type = "tmax, tmin"

## combine
niche_reg_intertidal <- rbind(niche_reg_intertidal, niche_upper_reg_intertidal)
niche_acclimatized_intertidal <- rbind(niche_acclimatized_intertidal, niche_upper_acclimatized_intertidal)
niche_reg_marine <- rbind(niche_reg_marine, niche_upper_reg_marine, niche_lower_reg_marine)
niche_acclimatized_marine <- rbind(niche_acclimatized_marine, niche_upper_acclimatized_marine,
                                   niche_lower_acclimatized_marine)
niche_Te_terr <- rbind(niche_Te_terr, niche_upper_Te_terr, niche_lower_Te_terr)
niche_Te_in_sun_terr <- rbind(niche_Te_in_sun_terr, niche_upper_Te_in_sun_terr, niche_lower_Te_terr)
niche_acclimatized_terr <- rbind(niche_acclimatized_terr, niche_upper_acclimatized_terr,
                                 niche_lower_acclimatized_terr)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### Collate temps in thermal niche across realms  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## do some rearranging to make niches across realms have the same number of types
niche_Te_intertidal <- niche_reg_intertidal %>%
  mutate(type = ifelse(type == "realized_reg", "realized_Te", 'potential_Te'))
niche_Te_in_sun_intertidal <- niche_reg_intertidal %>%
  mutate(type = ifelse(type == "realized_reg", "realized_Te_in_sun", 'potential_Te_in_sun'))

niche_Te_marine <- niche_reg_marine %>%
  mutate(type = ifelse(type == "realized_reg", "realized_Te", 'potential_Te'))
niche_Te_in_sun_marine <- niche_reg_marine %>%
  mutate(type = ifelse(type == "realized_reg", "realized_Te_in_sun", 'potential_Te_in_sun'))

niche_Te <- rbind(niche_Te_terr, niche_Te_marine, niche_Te_intertidal)
niche_Te_in_sun <- rbind(niche_Te_in_sun_terr, niche_Te_in_sun_marine, 
                         niche_Te_in_sun_intertidal)
niche_acclimatized <- rbind(niche_acclimatized_terr, niche_acclimatized_marine, 
                          niche_acclimatized_intertidal)

write.csv(niche_Te, "large-files/thermal-niches/all-realms/niche_Te.csv", row.names = FALSE)
write.csv(niche_Te_in_sun, "large-files/thermal-niches/all-realms/niche_Te_in_sun.csv", row.names = FALSE)
write.csv(niche_acclimatized, "large-files/thermal-niches/all-realms/niche_acclimatized.csv", row.names = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~##
#####      Plots       #####
##~~~~~~~~~~~~~~~~~~~~~~~~##
## function that plots potential ranges and corresponding realized ranges into a folder 
plot_ranges_overlap  <- function (folder, potential_ranges, type) {
  
  land <- as.data.frame(raster("data-processed/masks/raster_terr_mask.grd"), xy=TRUE)
  colnames(land)[1:3] <- c("longitude", "latitude", "mask")
  ocean <- as.data.frame(raster("data-processed/masks/raster_marine_mask.grd"), xy=TRUE)
  colnames(ocean)[1:3] <- c("longitude", "latitude", "mask")
  intertidal <- as.data.frame(raster("data-processed/masks/raster_intertidal_mask.grd"), xy=TRUE)
  colnames(intertidal)[1:3] <- c("longitude", "latitude", "mask")
  
  if (type %in% c("Te", "Te_in_sun", "acclimatized")) {
    land <- as.data.frame(raster("data-processed/masks/raster_terr_mask_nichemapr.grd"), xy=TRUE)
    colnames(land)[1:3] <- c("longitude", "latitude", "mask")
    
    rasterized_rrs <- readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds")
  }
  
  rrs <- as.data.frame(rasterized_rrs, xy=TRUE)
  colnames(rrs)[1:2] <- c("longitude", "latitude")
  
  prs <- as.data.frame(potential_ranges, xy=TRUE)
  colnames(prs)[1:2] <- c("longitude", "latitude")
  
  
  ## ggplots of all:
  i = 1
  while (i < ncol(prs) - 1) {
    range <- colnames(prs)[i+2]
    
    ## get realized range:
    species <- range %>%
      str_replace_all("_", ".")
   
    split <- str_split_fixed(range, "\\_", n = 3)
    source <- split[1,3]
    
    ## get species name and realm
    species <- paste(split[1,1], split[1,2], sep = ".")
    realm <- traits$Realm[which(traits$genus_species == paste(split[1,1], split[1,2], sep = "_"))]
    
    if(realm == "Terrestrial") {
      backdrop = land
    }
    else if (realm == "Marine") {
      backdrop = ocean
    }
    else {
      backdrop = intertidal
    }
    
    rr_index <- which(str_detect(colnames(rrs), species) & str_detect(colnames(rrs), source))
    
    r <- rrs[,c(1:2, rr_index)] 
    colnames(r)[3] <- "rr"
    r <- left_join(r, prs[,c(1:2, which(colnames(prs) == range))]) 
    colnames(r)[4] <- "pr"
    r$rr <- as.factor(r$rr)
    
    r_gg <- r %>%
      ggplot(., aes(x = longitude, y = latitude)) +
      xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
      geom_raster(aes(fill = rr)) + 
      scale_fill_manual(values = "yellow", na.value = "transparent") +
      annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.4,
               fill = r$pr) +
      annotate(geom="raster", x=land$longitude, y=land$latitude, alpha=.1,
               fill = backdrop$mask) +
      labs(title = range,
           y = "Latitude",
           x = "Longitude") +
      scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
      scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) +
      theme(legend.position = "none")
    
    
    ## write to file:
    ggsave(r_gg, path = paste("figures/additional-figures/range-plots/potential/", folder, sep = ""), 
           filename = paste(range, ".png", sep = ""), 
           height = 6, width = 10, units = "in", device = "png")
    
    i = i + 1
  }
  
}

isna <-c()
for (i in 1:nlayers(rasterized_rrs)) {
  isna <- append(isna, all(is.na(values(rasterized_rrs[[i]]))))
}

isna <- which(isna)
names(rasterized_rrs)[isna]

## plot the ranges:
plot_ranges_overlap(folder = "main-analysis/both", potential_ranges = potential_ranges_both_limits_Te,
                    type = "Te")
plot_ranges_overlap(folder = "main-analysis/upper", potential_ranges = potential_ranges_upper_limits_Te,
                    type = "Te")
plot_ranges_overlap(folder = "main-analysis/lower", potential_ranges = potential_ranges_lower_limits_Te,
                    type = "Te")

plot_ranges_overlap(folder = "acclimatised/both", 
                    potential_ranges = potential_ranges_both_limits_acclimatized,
                    type = "acclimatized")
plot_ranges_overlap(folder = "acclimatised/upper", 
                    potential_ranges = potential_ranges_upper_limits_acclimatized,
                    type = "acclimatized")
plot_ranges_overlap(folder = "acclimatised/lower", 
                    potential_ranges = potential_ranges_lower_limits_acclimatized,
                    type = "acclimatized")

plot_rrs <- function () {
  rasterized_rrs <- readRDS("data-processed/realized-ranges/rasterized_rrs.rds")
  og_rasterized_rrs <- readRDS("data-processed/realized-ranges/og_rasterized_rrs.rds")
  
  rrs <- as.data.frame(rasterized_rrs, xy=TRUE)
  colnames(rrs)[1:2] <- c("longitude", "latitude")
  
  og_rrs <- as.data.frame(og_rasterized_rrs, xy=TRUE)
  colnames(og_rrs)[1:2] <- c("longitude", "latitude")
  
  land <- as.data.frame(raster("data-processed/masks/raster_terr_mask.grd"), xy=TRUE)
  colnames(land)[1:3] <- c("longitude", "latitude", "mask")
  ocean <- as.data.frame(raster("data-processed/masks/raster_marine_mask.grd"), xy=TRUE)
  colnames(ocean)[1:3] <- c("longitude", "latitude", "mask")
  intertidal <- as.data.frame(raster("data-processed/masks/raster_intertidal_mask.grd"), xy=TRUE)
  colnames(intertidal)[1:3] <- c("longitude", "latitude", "mask")
  
  ## ggplots of all:
  i = 1
  while (i < ncol(og_rrs) - 1) {
    range <- colnames(og_rrs)[i+2]
    
    ## get realized range:
    species <- range %>%
      str_replace_all("\\_", ".")
    
    split <- str_split_fixed(species, "\\.", n = 3)
    source <- split[1,3]
    
    ## get species name and realm
    species <- paste(split[1,1], split[1,2], sep = ".")
    realm <- traits$Realm[which(paste(traits$Genus, traits$Species, sep = "_") == 
                                  paste(split[1,1], split[1,2], sep = "_"))]
    
    if(realm == "Terrestrial" | realm == "Freshwater") {
      backdrop = land
    }
    else if (realm == "Marine") {
      backdrop = ocean
    }
    else {
      backdrop = intertidal
    }
    
    r <- rrs[,c(1:2, i+2)] 
    colnames(r)[3] <- "rr"
    r <- left_join(r, og_rrs[,c(1:2, i+2)]) 
    colnames(r)[4] <- "og_rr"
    
    r_gg <- r %>%
      ggplot(., aes(x = longitude, y = latitude)) +
      xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
      geom_raster(aes(fill=as.factor(og_rr))) + 
      scale_fill_manual(values = c("red"), aesthetics = 'fill', labels = ) +
      annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.5,
               fill = r$rr) +
      annotate(geom="raster", x=land$longitude, y=land$latitude, alpha=.1,
               fill = backdrop$mask) +
      labs(title = range,
           y = "Latitude",
           x = "Longitude") +
      scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
      scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) +
      theme(legend.position = "none")
    
    
    ## write to file:
    ggsave(r_gg, path = paste("figures/additional-figures/range-plots/realized/", sep = ""), 
           filename = paste(range, ".png", sep = ""), 
           height = 6, width = 10, units = "in", device = "png")
    
    
    i = i + 1
  }
  
}

plot_rrs()

