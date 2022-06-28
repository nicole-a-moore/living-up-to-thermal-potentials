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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
#####         Functions           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## functions to filter rasterized temperature data (raster_high, raster_low) by thermal tolerance limits (either upper, lower, or both) 

## functions for filtering by both upper and lower thermal limits:
## terrestrial:
filter_by_tolerance_Te <- function(both_upper, 
                                   both_lower, 
                                   Te_max, 
                                   Te_min) {
  ## create an individual raster layer of difference between thermal limit and Te
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_high <- high
  rr_low <- low
  
  ## subset to species that we have a Te for:
  upper <- filter(both_upper, both_upper$genus_species %in% names(Te_max))
  lower <- filter(both_lower, both_lower$genus_species %in% names(Te_min))
  
  species = 1
  while (species < nrow(upper) + 1) {
    
    ## get species name and extract its Te layers:
    sp <- upper$genus_species[species] 
    
    h <- Te_max[[which(names(Te_max) == sp)]]
    l <- Te_min[[which(names(Te_min) == sp)]]
    
    rr_high <- addLayer(rr_high, h) 
    rr_low <- addLayer(rr_low, l) 
    
    ## exclude raster cells outside of the thermal tolerance (where Tmax - hottest Te > 0 and where Tmin - coldest Te< 0)
    below_ctmax <- upper$thermal_limit[species] - h
    h[below_ctmax < 0] <- NA
    above_ctmin <- lower$thermal_limit[species] - l
    l[above_ctmin > 0] <- NA
    
    high <- addLayer(high, h) 
    low <- addLayer(low, l) 
    
    species = species + 1
  }
  
  names(high) <- upper$genus_species
  names(low) <- lower$genus_species
  
  ## combine to find cells where Te_sun is less than CTmax AND Te_shade is greater than CTmin 
  combined_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                          crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  combined_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                         crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  i = 1  
  while (i < nlayers(high) + 1) {
    combined_high <- addLayer(combined_high, mask(high[[i]], low[[i]]), updatevalue = NA)
    combined_low <- addLayer(combined_low, mask(low[[i]], high[[i]]), updatevalue = NA)
    
    i = i + 1
  }
  names(combined_high) <- names(high)
  names(combined_low) <- names(low)
  names(rr_high) <- names(high)
  names(rr_low) <- names(low)
  #plot(combined_high[[1]])
  #plot(combined_low[[1]])
  
  return (list(combined_high, combined_low, rr_high, rr_low))
}
## marine and intertidal:
filter_by_tolerance <- function(both_upper, 
                                both_lower, 
                                raster_high, 
                                raster_low) {
    ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species
    high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                  crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    
    ## create set of temps to use to extract realized thermal niches later 
    rr_high <- high
    rr_low <- low
    
    species = 1
    high_names <- c()
    low_names <- c()
    while (species < nrow(both_upper) + 1) {
      species_traits <- traits[which(traits$genus_species == both_upper$genus_species[species]),]
      
      if (species_traits$cold_season_dormancy_ == "Yes" & species_traits$hot_season_dormancy_ == "Yes") {
        h <- max(stack(raster_high$hot_dormancy_6mo,  raster_high$cold_dormancy_6mo), na.rm = TRUE)
        l <- min(stack(raster_low$cold_dormancy_6mo,  raster_low$hot_dormancy_6mo), na.rm = TRUE)
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 
      }
      else if (species_traits$cold_season_dormancy_ == "Yes") {
        h <- raster_high$cold_dormancy_6mo
        l <-  raster_low$cold_dormancy_6mo
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 
      }
      else if (species_traits$hot_season_dormancy_ == "Yes") {
        h <- raster_high$hot_dormancy_6mo
        l <- raster_low$hot_dormancy_6mo
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l)
      }
      else {
        h <- raster_high$seasonal_high_temp
        l <- raster_low$seasonal_low_temp
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 
      }
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      species = species + 1
    }
    
    names(high) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
    names(low) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
    ##plot(high)
    ##plot(low)
    
    ## combine to find cells where seasonal high temp is less than CTmax AND seasonal low temp is greater than CTmin 
    combined_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    combined_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    
    i = 1  
    while (i < nlayers(high) + 1) {
      combined_high <- addLayer(combined_high, mask(high[[i]], low[[i]]), updatevalue = NA)
      combined_low <- addLayer(combined_low, mask(low[[i]], high[[i]]), updatevalue = NA)
      
      i = i + 1
    }
    names(combined_high) <- names(high)
    names(combined_low) <- names(low)
    names(rr_high) <- names(high)
    names(rr_low) <- names(low)
    #plot(combined_high[[1]])
    #plot(combined_low[[1]])
    
    return (list(combined_high, combined_low, rr_high, rr_low))
}

filter_by_tolerance_acclimatized <- function(both_upper,
                                           both_lower,
                                           raster_high, 
                                           raster_low,
                                           acc_data, realm) {
  
  ## create an individual raster layer of difference between thermal limit and Te
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_high <- high
  rr_low <- low
  
  tlims_upper_all <- high
  tlims_lower_all <- low
  
  ## read in acclimation temperatures:
  if (realm == "Terrestrial") {
    acc_cold <- stack("large-files/operative-temperatures/terr_acc_cold.grd")
    acc_hot <- stack("large-files/operative-temperatures/terr_acc_hot.grd")
    
    ## subset to species that we have a Te for:
    upper <- filter(both_upper, both_upper$genus_species %in% names(raster_high))
    lower <- filter(both_lower, both_lower$genus_species %in% names(raster_low))
    
    species = 1
    while (species < nrow(upper) + 1) {
      
      ## get species name and extract its Te layers:
      sp <- upper$genus_species[species] 
      
      h <- raster_high[[which(names(raster_high) == sp)]]
      l <- raster_low[[which(names(raster_low) == sp)]]
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## create raster layers of thermal limits if species was acclimatized to the cell
      tlims_upper <- acc_hot[[which(names(acc_hot) == sp)]]
      tlims_lower <- acc_cold[[which(names(acc_cold) == sp)]]
      
      equ <- acc_data[which(acc_data$genus_species == sp),]
      hot_equ <- equ[which(equ$type == "max"),]
      cold_equ <- equ[which(equ$type == "min"),]
      
      if (nrow(hot_equ) != 0) {
        if (!is.na(hot_equ$acclimation_temperature)){
          slope = hot_equ$ARR_equ_slope
          int = hot_equ$ARR_equ_int
          
          tlims_upper <- slope*tlims_upper + int
          
          # ## make sure acclimation temps are beyond tested acclimation temperatures:
          # tlims_upper[tlims_upper > slope*hot_equ$max_upper_acc_temp + int] <-
          #   slope*hot_equ$max_upper_acc_temp + int
          # tlims_upper[tlims_upper < slope*hot_equ$min_upper_acc_temp + int] <-
          #   slope*hot_equ$min_upper_acc_temp + int
        }
        else {
          tlims_upper[!is.na(tlims_upper)] <- upper$thermal_limit[species]
        }
      }
      else {
        tlims_upper[!is.na(tlims_upper)] <- upper$thermal_limit[species]
      }
      if (nrow(cold_equ) != 0) {
        if (!is.na(cold_equ$acclimation_temperature)){
          slope = cold_equ$ARR_equ_slope
          int = cold_equ$ARR_equ_int
          
          tlims_lower <- slope*tlims_lower + int
          
          # ## make sure acclimation temps are beyond tested acclimation temperatures:
          # tlims_lower[tlims_lower < slope*cold_equ$min_lower_acc_temp + int] <-
          #   slope*cold_equ$min_lower_acc_temp + int
          # tlims_lower[tlims_lower > slope*cold_equ$max_lower_acc_temp + int] <-
          #   slope*cold_equ$max_lower_acc_temp + int
        }
        else {
          tlims_lower[!is.na(tlims_lower)] <- lower$thermal_limit[which(lower$genus_species == sp)]
        }
      }
      else {
        tlims_lower[!is.na(tlims_lower)] <- lower$thermal_limit[which(lower$genus_species == sp)]
      }
      
      ## exclude raster cells outside of the thermal tolerance (where acclimatized Tmax - hottest Te > 0 and where acclimatized Tmin - coldest Te< 0)
      h[tlims_upper < h]  <- NA
      l[tlims_lower > l] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      tlims_upper_all <- addLayer(tlims_upper_all, tlims_upper) 
      tlims_lower_all <- addLayer(tlims_lower_all, tlims_lower) 
      
      species = species + 1
    }
  }
  else {
    upper <- both_upper
    lower <- both_lower
    
    species = 1
    high_names <- c()
    low_names <- c()
    while (species < nrow(both_upper) + 1) {
      species_traits <- traits[which(traits$genus_species == both_upper$genus_species[species]),]
      sp <- species_traits$genus_species[1] 
      
      if (species_traits$cold_season_dormancy_ == "Yes" & species_traits$hot_season_dormancy_ == "Yes") {
        h <- max(stack(raster_high$hot_dormancy_6mo,  raster_high$cold_dormancy_6mo), na.rm = TRUE)
        l <- min(stack(raster_low$cold_dormancy_6mo,  raster_low$hot_dormancy_6mo), na.rm = TRUE)
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 

      }
      else if (species_traits$cold_season_dormancy_ == "Yes") {
        h <- raster_high$cold_dormancy_6mo
        l <-  raster_low$cold_dormancy_6mo
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 

      }
      else if (species_traits$hot_season_dormancy_ == "Yes") {
        h <- raster_high$hot_dormancy_6mo
        l <- raster_low$hot_dormancy_6mo
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 
  
      }
      else {
        h <- raster_high$seasonal_high_temp
        l <- raster_low$seasonal_low_temp
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 
      }
      
      ## create raster layer of thermal limits if species was acclimatized to the cell
      equ <- acc_data[which(acc_data$genus_species == species_traits$genus_species),]
      hot_equ <- equ[which(equ$type == "max"),]
      cold_equ <- equ[which(equ$type == "min"),]
      
      tlims_lower <- raster_low[[4]]
      tlims_upper <- raster_high[[4]]
      
      if (nrow(hot_equ) != 0) {
        if (!is.na(hot_equ$acclimation_temperature)){
          slope = hot_equ$ARR_equ_slope
          int = hot_equ$ARR_equ_int
          
          tlims_upper <- slope*tlims_upper + int
          
          # ## make sure acclimation temps are beyond tested acclimation temperatures:
          # tlims_upper[tlims_upper > slope*hot_equ$max_upper_acc_temp + int] <-
          #   slope*hot_equ$max_upper_acc_temp + int
          # tlims_upper[tlims_upper < slope*hot_equ$min_upper_acc_temp + int] <-
          #   slope*hot_equ$min_upper_acc_temp + int
        }
        else {
          tlims_upper[!is.na(tlims_upper)] <- upper$thermal_limit[species]
        }
      }
      else {
        tlims_upper[!is.na(tlims_upper)] <- upper$thermal_limit[species]
      }
      if (nrow(cold_equ) != 0) {
        if (!is.na(cold_equ$acclimation_temperature)){
          slope = cold_equ$ARR_equ_slope
          int = cold_equ$ARR_equ_int
          
          tlims_lower <- slope*tlims_lower + int
          
          # ## make sure acclimation temps are beyond tested acclimation temperatures:
          # tlims_lower[tlims_lower < slope*cold_equ$min_lower_acc_temp + int] <-
          #   slope*cold_equ$min_lower_acc_temp + int
          # tlims_lower[tlims_lower > slope*cold_equ$max_lower_acc_temp + int] <-
          #   slope*cold_equ$max_lower_acc_temp + int
        }
        else {
          tlims_lower[!is.na(tlims_lower)] <- lower$thermal_limit[which(lower$genus_species == sp)]
        }
      }
      else {
        tlims_lower[!is.na(tlims_lower)] <- lower$thermal_limit[which(lower$genus_species == sp)]
      }
      
      ## exclude raster cells outside of the acclimatized thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      h[tlims_upper < h]  <- NA
      l[tlims_lower > l] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
    
      tlims_upper_all <- addLayer(tlims_upper_all, tlims_upper) 
      tlims_lower_all <- addLayer(tlims_lower_all, tlims_lower) 
      
      species = species + 1
    }
    
    names(high) <- paste(upper$Genus, upper$Species, sep = "_")
    names(low) <- paste(upper$Genus, upper$Species, sep = "_")
    names(tlims_upper_all) <- paste(upper$Genus, upper$Species, sep = "_")
    names(tlims_lower_all) <- paste(upper$Genus, upper$Species, sep = "_")
    ##plot(high)
    ##plot(low)
  }
  
  
  ## combine to find cells where Te_sun is less than CTmax AND Te_shade is greater than CTmin 
  combined_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                          crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  combined_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                         crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  i = 1  
  while (i < nlayers(high) + 1) {
    combined_high <- addLayer(combined_high, mask(high[[i]], low[[i]]), updatevalue = NA)
    combined_low <- addLayer(combined_low, mask(low[[i]], high[[i]]), updatevalue = NA)
    
    i = i + 1
  }
  names(combined_high) <- names(high)
  names(combined_low) <- names(low)
  names(rr_high) <- names(high)
  names(rr_low) <- names(low)
  #plot(combined_high[[1]])
  #plot(combined_low[[1]])
  
  return (list(combined_high, combined_low, rr_high, rr_low, tlims_upper_all, tlims_lower_all))
}

## functions for filtering by only one thermal limit:
filter_by_tolerance_one_limit <- function(lims, 
                                raster_high, 
                                raster_low, 
                                limit_type) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species
  temps <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_temps <- temps
  
  species = 1
  while (species < nrow(lims) + 1) {
    species_traits <- traits[which(traits$genus_species == lims$genus_species[species]),]
    
    if (species_traits$cold_season_dormancy_ == "Yes" & species_traits$hot_season_dormancy_ == "Yes") {
      if (limit_type == "tmax") {
        t <- max(stack(raster_high$hot_dormancy_6mo,  raster_high$cold_dormancy_6mo), na.rm = TRUE)
      }
      else if (limit_type == 'tmin') {
        t <- min(stack(raster_low$cold_dormancy_6mo,  raster_low$hot_dormancy_6mo), na.rm = TRUE)
      }
      rr_temps <- addLayer(rr_temps, t) 
    }
    else if (species_traits$cold_season_dormancy_ == "Yes") {
      if (limit_type == "tmax") {
        t <- raster_high$cold_dormancy_6mo
      }
      else if (limit_type == 'tmin') {
        t <- raster_low$cold_dormancy_6mo
      }
      rr_temps <- addLayer(rr_temps, t) 
    }
    else if (species_traits$hot_season_dormancy_ == "Yes") {
      if (limit_type == "tmax") {
        t <- raster_high$hot_dormancy_6mo
      }
      else if (limit_type == 'tmin') {
        t <- raster_low$hot_dormancy_6mo
      }
      rr_temps <- addLayer(rr_temps, t) 
    }
    else {
      if (limit_type == "tmax") {
        t <- raster_high$seasonal_high_temp
      }
      else if (limit_type == 'tmin') {
        t <- raster_low$seasonal_low_temp
      }
      rr_temps <- addLayer(rr_temps, t) 
    }
    
    ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 or where Tmin - seasonal_low< 0)
    if (limit_type == "tmax") {
      below_ctmax <- lims$thermal_limit[species] - t
      t[below_ctmax < 0] <- NA
    }
    else if (limit_type == 'tmin') {
      above_ctmin <- lims$thermal_limit[species] - t
      t[above_ctmin > 0] <- NA
    }
    temps <- addLayer(temps, t) 
    
    species = species + 1
  }
  
  names(temps) <- paste(lims$Genus, lims$Species, sep = "_")
  ##plot(temps)
  names(rr_temps) <- names(temps)
  
  return (list(temps, rr_temps))
}

filter_by_tolerance_Te_one_limit <- function(limit_type, lims, 
                                             Te_max, 
                                             Te_min) {
  ## create an individual raster layer of difference between thermal limit and Te
  temps <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_temps <- temps
  
  ## subset to species that we have a Te for:
  lims <- filter(lims, lims$genus_species %in% names(Te_max))
  
  species = 1
  while (species < nrow(lims) + 1) {
    
    ## get species name and extract its Te layers:
    sp <- lims$genus_species[species] 
    
    if (limit_type == "tmax") {
      t = Te_max[[which(names(Te_max) == sp)]]
    }
    else if (limit_type == "tmin") {
      t <- Te_min[[which(names(Te_min) == sp)]]
    }
    
    rr_temps <- addLayer(rr_temps, t) 
    
    ## exclude raster cells outside of the thermal tolerance limit (where Tmax - hottest Te > 0 or where Tmin - coldest Te< 0)
    if (limit_type == "tmax") {
      below_ctmax <- lims$thermal_limit[species] - t
      t[below_ctmax < 0] <- NA
    }
    else if (limit_type == "tmin") {
      above_ctmin <- lims$thermal_limit[species] - t
      t[above_ctmin > 0] <- NA
    }
    temps <- addLayer(temps, t) 
    
    species = species + 1
  }
  
  names(temps) <- lims$genus_species
  names(rr_temps) <- names(temps)
  
  return (list(temps, rr_temps))
}

filter_by_tolerance_acclimatized_one_limit <- function(lims, limit_type,
                                                       raster_high, 
                                                       raster_low,
                                                       acc_data, realm) {
  ## create an individual raster layer of difference between thermal limit and Te
  temps <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_temps <- temps
  
  tlims_all <- temps
  
  ## read in acclimation temperatures:
  if (realm == "Terrestrial") {
    if (limit_type == 'tmax') {
      acc <- stack("large-files/operative-temperatures/terr_acc_hot.grd")
    }
    else if (limit_type == "tmin") {
      acc <- stack("large-files/operative-temperatures/terr_acc_cold.grd")
    }
  
    ## subset to species that we have a Te for:
    lims <- filter(lims, lims$genus_species %in% names(raster_high))
    
    species = 1
    while (species < nrow(lims) + 1) {
      
      ## get species name and extract its Te layers:
      sp <- lims$genus_species[species] 
      
      if (limit_type == 'tmax') {
        t <- raster_high[[which(names(raster_high) == sp)]]
      }
      else if (limit_type == "tmin") {
        t <- raster_low[[which(names(raster_low) == sp)]]
      }
      rr_temps <- addLayer(rr_temps, t) 
      
      ## create raster layers of thermal limits if species was acclimatized to the cell
      tlims <- acc[[which(names(acc) == sp)]]
      
      equ <- acc_data[which(acc_data$genus_species == sp),]
      hot_equ <- equ[which(equ$type == "max"),]
      cold_equ <- equ[which(equ$type == "min"),]
      
      if (nrow(hot_equ) != 0) {
        if (!is.na(hot_equ$acclimation_temperature) & limit_type == "tmax"){
          slope = hot_equ$ARR_equ_slope
          int = hot_equ$ARR_equ_int
          
          tlims <- slope*tlims + int
        }
        else {
          tlims[!is.na(tlims)] <- lims$thermal_limit[species]
        }
      }
      else if (nrow(cold_equ) != 0 & limit_type == "tmin") {
        if (!is.na(cold_equ$acclimation_temperature)){
          slope = cold_equ$ARR_equ_slope
          int = cold_equ$ARR_equ_int
          
          tlims <- slope*tlims + int
        }
        else {
          tlims[!is.na(tlims)] <- lims$thermal_limit[which(lims$genus_species == sp)]
        }
      }
      else {
        tlims[!is.na(tlims)] <- lims$thermal_limit[which(lims$genus_species == sp)]
      }
      
      ## exclude raster cells outside of the thermal tolerance (where acclimatized Tmax - hottest Te > 0 or where acclimatized Tmin - coldest Te< 0)
      if (limit_type == 'tmax') {
        t[tlims < t]  <- NA
      }
      else if (limit_type == "tmin") {
        t[tlims > t] <- NA
      }
      
      temps <- addLayer(temps, t)
      tlims_all <- addLayer(tlims_all, tlims) 
      
      species = species + 1
    }
  }
  else {
    species = 1
    while (species < nrow(lims) + 1) {
      species_traits <- traits[which(traits$genus_species == lims$genus_species[species]),]
      sp <- species_traits$genus_species[1] 
      
      if (species_traits$cold_season_dormancy_ == "Yes" & species_traits$hot_season_dormancy_ == "Yes") {
        if (limit_type == "tmax") {
          t <- max(stack(raster_high$hot_dormancy_6mo,  raster_high$cold_dormancy_6mo), na.rm = TRUE)
        }
        else if (limit_type == 'tmin') {
          t <- min(stack(raster_low$cold_dormancy_6mo,  raster_low$hot_dormancy_6mo), na.rm = TRUE)
        }
        rr_temps <- addLayer(rr_temps, t) 
      }
      else if (species_traits$cold_season_dormancy_ == "Yes") {
        if (limit_type == "tmax") {
          t <- raster_high$cold_dormancy_6mo
        }
        else if (limit_type == 'tmin') {
          t <- raster_low$cold_dormancy_6mo
        }
        rr_temps <- addLayer(rr_temps, t) 
      }
      else if (species_traits$hot_season_dormancy_ == "Yes") {
        if (limit_type == "tmax") {
          t <- raster_high$hot_dormancy_6mo
        }
        else if (limit_type == 'tmin') {
          t <- raster_low$hot_dormancy_6mo
        }
        rr_temps <- addLayer(rr_temps, t) 
      }
      else {
        if (limit_type == "tmax") {
          t <- raster_high$seasonal_high_temp
        }
        else if (limit_type == 'tmin') {
          t <- raster_low$seasonal_low_temp
        }
        rr_temps <- addLayer(rr_temps, t) 
        
      }
      
      ## create raster layer of thermal limits if species was acclimatized to the cell
      equ <- acc_data[which(acc_data$genus_species == species_traits$genus_species),]
      hot_equ <- equ[which(equ$type == "max"),]
      cold_equ <- equ[which(equ$type == "min"),]
      
      if (limit_type == "tmax") {
        tlims <- raster_high[[4]]
      }
      else if (limit_type == 'tmin') {
        tlims <- raster_low[[4]]
      }
      
      if (nrow(hot_equ) != 0 & limit_type == "tmax") {
        if (!is.na(hot_equ$acclimation_temperature)){
          slope = hot_equ$ARR_equ_slope
          int = hot_equ$ARR_equ_int
          
          tlims <- slope*tlims + int
        }
        else {
          tlims[!is.na(tlims)] <- lims$thermal_limit[species]
        }
      }
      else if (nrow(cold_equ) != 0 & limit_type == "tmin") {
        if (!is.na(cold_equ$acclimation_temperature)){
          slope = cold_equ$ARR_equ_slope
          int = cold_equ$ARR_equ_int
          
          tlims <- slope*tlims + int
        }
        else {
          tlims[!is.na(tlims)] <- lims$thermal_limit[which(lims$genus_species == sp)]
        }
      }
      else {
        tlims[!is.na(tlims)] <- lims$thermal_limit[which(lims$genus_species == sp)]
      }
      
      ## exclude raster cells outside of the acclimatized thermal tolerance (where Tmax - seasonal_high > 0 or where Tmin - seasonal_low< 0)
      if (limit_type == "tmax") {
        t[tlims < t]  <- NA
      }
      else if (limit_type == "tmin") {
        t[tlims > t]  <- NA
      }
      temps <- addLayer(temps, t) 
      tlims_all <- addLayer(tlims_all, tlims) 
      
      species = species + 1
    }
  }
  
  names(temps) <- paste(lims$Genus, lims$Species, sep = "_")
  names(tlims_all) <- paste(lims$Genus, lims$Species, sep = "_")
  names(rr_temps) <- names(temps)
  
  return (list(temps, rr_temps, tlims_all))
}


## creates potential range rasterstack containing a layer for each realized range that represents the species potential range containing only contiguous habitat in clumps that their realized range touches 
## returns 3 versions of potential ranges (boolean values, high temperature values and low temperature values) 
create_potential_ranges <- function (clumped_temps, 
                                     filtered, 
                                     realized_ranges, 
                                     type) {
  
  filtered_high <- filtered[[1]]
  filtered_low <- filtered[[2]]
  
  rr_high <- filtered[[3]]
  rr_low <- filtered[[4]]
  
  rr_high_new <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                        crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  rr_low_new <- rr_high_new
  
  if (type == "acclimatized") {
    tlims_upper_all <- filtered[[5]]
    tlims_lower_all <- filtered[[6]]
    tlims_upper <- tlims_lower <- rr_high_new
  }
  
  pr_restricted_bool <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  pr_restricted_high <- pr_restricted_bool
  pr_restricted_low <- pr_restricted_bool
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1, vals = 1)
  
  rasterized_rrs <- readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds")
  
  names <- c()
  i = 1
  while (i < length(names(filtered_high)) + 1) {
    ## get unrestricted potential range 
    potential_raster <- filtered_high[[i]] 
    
    ## get realized range and realized range raster
    if(type == 'reg') {
      split <- str_split_fixed(names(potential_raster), "_", n = 3)
      species <- paste(split[1,1], split[1,2], sep = " ")
    }
    else if (type == "acclimatized") {
      species <- str_replace_all(names(potential_raster), "_", " ")
    }
    else {
      species <- names(potential_raster) %>%
        str_replace_all("_", " ") 
    }
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0") 
    
    ## get traits 
    sp_traits <- traits[which(str_replace_all(traits$genus_species, "_", " ") == species), ]
    
    ## if terrestrial, restrict by elevational range limits:
    if (sp_traits$Realm == "Terrestrial") {
      ## get elevational range limits 
      elev_max <- sp_traits$upper_elevation_limit_meters
      elev_min <- sp_traits$lower_elevation_limit_meters
      
      ## if both limits are not NA, restrict by both limits
      if (!is.na(elev_min) & !is.na(elev_max)) {
        potential_raster[elevs$elev_max < elev_min] <- NA 
        potential_raster[elevs$elev_min > elev_max] <- NA 
        ## exclude cells with max elev less than species min elevation
        ## exclude cells with min elev greater than species max elevation
      }
      else if (!is.na(elev_min)) {
        potential_raster[elevs$elev_max < elev_min] <- NA 
        ## exclude cells with max elev less than species min elevation
      }
      else if (!is.na(elev_max)) {
        potential_raster[elevs$elev_min > elev_max] <- NA 
        ## exclude cells with min elev greater than species max elevation
      }
    }
    ## if marine, restrict by depth distribution limits:
    else if (sp_traits$Realm %in% c("Marine", "Intertidal")) {
      ## get depth distribution limits 
      depth_upper <- sp_traits$upper_depth_limit
      depth_lower <- sp_traits$lower_depth_limit
      
      ## if both limits are not NA, restrict by both depth limits
      if (!is.na(depth_lower) & !is.na(depth_upper)) {
        potential_raster[depths$depth_min > depth_upper] <- NA 
        potential_raster[depths$depth_max < depth_lower] <- NA 
        ## exclude cells with deepest depth above species upper limit  
        ## exclude cells with shallowest depth less than species lower limit 
      }
      else if (!is.na(depth_upper)) {
        potential_raster[depths$depth_min > depth_upper] <- NA 
        ## exclude cells with deepest depth above species upper limit  
      }
      else if (!is.na(depth_lower)) {
        potential_raster[depths$depth_max < depth_lower] <- NA 
        ## exclude cells with shallowest depth less than species lower limit 
      }
      ## if lower depth limit is unknown, restrict to continental shelf if it is a coastal spp
      else if (is.na(depth_lower) & !is.na(sp_traits$coastal_or_oceanic)) {
        is_coastal <- sp_traits$coastal_or_oceanic == "coastal"
        
        if(is_coastal) {
          potential_raster[is.na(shelf_mask)] <- NA
        }
      }
    }
    
    ## if species has multiple realized ranges, loop through them:
    num = 1
    while (num < nrow(realized_range)+1) {
      
      ## find areas of habitat that realized range intersects
      intersects <- st_intersects(clumped_temps, realized_range[num,], sparse = FALSE)[,]
      intersects <- filter(clumped_temps, intersects == TRUE) 
      intersects$geometry <- st_union(intersects) # combine all habitat into one for rasterizing
      
      ## rasterize them:
      if (is.na(intersects[[1,1]])) {
        num = num + 1
      }
      else {
        r_intersects <- shp2rast(r, as_Spatial(intersects[1,]))
        
        if (type == "acclimatized") {
          ## save acclimatized thermal limits across area allowed to be in potential range
          tlims_upper <- addLayer(tlims_upper, mask(tlims_upper_all[[i]], r_intersects))
          tlims_lower <- addLayer(tlims_lower, mask(tlims_lower_all[[i]], r_intersects))
        }
        
        pr_restricted <- mask(potential_raster, intersects) ## get rid of potential range outside of areas that realized range intersects
        pr_restricted[!is.na(pr_restricted)] <- 1 ## make it boolean
        
        ## create boolean version of potential range 
        pr_restricted_bool <- addLayer(pr_restricted_bool, pr_restricted, updatevalue = NA)
        
        ## create version of prs with temperature values for extracting thermal niche later 
        pr_restricted_high <- addLayer(pr_restricted_high, 
                                       mask(filtered_high[[i]], pr_restricted, updatevalue = NA), 
                                       updatevalue = NA)
        pr_restricted_low <- addLayer(pr_restricted_low, 
                                      mask(filtered_low[[i]], pr_restricted, updatevalue = NA), 
                                      updatevalue = NA)
        
        ## add layer to new rr:
        rr_high_new <- addLayer(rr_high_new, rr_high[[i]])
        rr_low_new <- addLayer(rr_low_new, rr_low[[i]])
        
        names <- append(names, paste(names(potential_raster), realized_range$source[num], sep = "_"))
        
        num = num + 1
      }
    }
    
    i = i + 1
  }
  names(pr_restricted_bool) <- names
  names(pr_restricted_high) <- names
  names(pr_restricted_low) <- names
  names(rr_high_new) <- names
  names(rr_low_new) <- names
  #plot(pr_restricted_bool)
  
  if (type == "acclimatized") {
    names(tlims_upper) <- names(tlims_lower) <- names
    return(list(pr_restricted_bool, pr_restricted_high, pr_restricted_low, rr_high_new, rr_low_new,
                tlims_upper, tlims_lower))
  }
  else {
    return(list(pr_restricted_bool, pr_restricted_high, pr_restricted_low, rr_high_new, rr_low_new))
  }
}

create_potential_ranges_one_limit <- function (clumped_temps,  
                                               filtered, 
                                               realized_ranges, 
                                               type) {
  
  filtered_temp <- filtered[[1]]
  rr_temp <- filtered[[2]]
  
  rr_new <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1, vals = 1)
  
  if (type == "acclimatized") {
    tlims_all <- filtered[[3]]
    tlims <- rr_new
  }
  
  pr_restricted_bool <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  pr_restricted_temps <- pr_restricted_bool

  rasterized_rrs <- readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds")
  
  names <- c()
  i = 1
  while (i < length(names(filtered_temp)) + 1) {
    ## get unrestricted potential range 
    potential_raster <- filtered_temp[[i]] 
    
    ## get realized range and realized range raster
    if(type == 'reg') {
      split <- str_split_fixed(names(potential_raster), "_", n = 3)
      species <- paste(split[1,1], split[1,2], sep = " ")
    }
    else if (type == "acclimatized") {
      species <- str_replace_all(names(potential_raster), "_", " ")
    }
    else {
      species <- names(potential_raster) %>%
        str_replace_all("_", " ") 
    }
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0") 
    
    ## get traits 
    sp_traits <- traits[which(str_replace_all(traits$genus_species, "_", " ") == species), ]
    
    ## if terrestrial, restrict by elevational range limits:
    if (sp_traits$Realm == "Terrestrial") {
      ## get elevational range limits 
      elev_max <- sp_traits$upper_elevation_limit_meters
      elev_min <- sp_traits$lower_elevation_limit_meters
      
      ## if both limits are not NA, restrict by both limits
      if (!is.na(elev_min) & !is.na(elev_max)) {
        potential_raster[elevs$elev_max < elev_min] <- NA 
        potential_raster[elevs$elev_min > elev_max] <- NA 
        ## exclude cells with max elev less than species min elevation
        ## exclude cells with min elev greater than species max elevation
      }
      else if (!is.na(elev_min)) {
        potential_raster[elevs$elev_max < elev_min] <- NA 
        ## exclude cells with max elev less than species min elevation
      }
      else if (!is.na(elev_max)) {
        potential_raster[elevs$elev_min > elev_max] <- NA 
        ## exclude cells with min elev greater than species max elevation
      }
    }
    ## if marine, restrict by depth distribution limits:
    else if (sp_traits$Realm %in% c("Marine", "Intertidal")) {
      ## get depth distribution limits 
      depth_upper <- sp_traits$upper_depth_limit
      depth_lower <- sp_traits$lower_depth_limit
      
      ## if both limits are not NA, restrict by both depth limits
      if (!is.na(depth_lower) & !is.na(depth_upper)) {
        potential_raster[depths$depth_min > depth_upper] <- NA 
        potential_raster[depths$depth_max < depth_lower] <- NA 
        ## exclude cells with deepest depth above species upper limit  
        ## exclude cells with shallowest depth less than species lower limit 
      }
      else if (!is.na(depth_upper)) {
        potential_raster[depths$depth_min > depth_upper] <- NA 
        ## exclude cells with deepest depth above species upper limit  
      }
      else if (!is.na(depth_lower)) {
        potential_raster[depths$depth_max < depth_lower] <- NA 
        ## exclude cells with shallowest depth less than species lower limit 
      }
      ## if lower depth limit is unknown, restrict to continental shelf if it is a coastal spp
      else if (is.na(depth_lower) & !is.na(sp_traits$coastal_or_oceanic)) {
        is_coastal <- sp_traits$coastal_or_oceanic == "coastal"
        
        if(is_coastal) {
          potential_raster[is.na(shelf_mask)] <- NA
        }
      }
    }
    
    ## if species has multiple realized ranges, loop through them:
    num = 1
    while (num < nrow(realized_range)+1) {
      
      ## find areas of habitat that realized range intersects
      intersects <- st_intersects(clumped_temps, realized_range[num,], sparse = FALSE)[,]
      intersects <- filter(clumped_temps, intersects == TRUE) 
      intersects$geometry <- st_union(intersects) # combine all habitat into one for rasterizing
      
      ## rasterize them:
      if (nrow(intersects) == 0) {
        num = num + 1
      }
      else if(is.na(intersects[[1,1]])) {
        num = num + 1
      }
      else {
        r_intersects <- shp2rast(r, as_Spatial(intersects[1,]))
        
        if (type == "acclimatized") {
          ## save acclimatized thermal limits across area allowed to be in potential range
          tlims <- addLayer(tlims, mask(tlims_all[[i]], r_intersects))
        }
        
        pr_restricted <- mask(potential_raster, intersects) ## get rid of potential range outside of areas that realized range intersects
        pr_restricted[!is.na(pr_restricted)] <- 1 ## make it boolean
        
        ## create boolean version of potential range 
        pr_restricted_bool <- addLayer(pr_restricted_bool, pr_restricted, updatevalue = NA)
        
        ## create version of prs with temperature values for extracting thermal niche later 
        pr_restricted_temps <- addLayer(pr_restricted_temps, 
                                       mask(filtered_temp[[i]], pr_restricted, updatevalue = NA), 
                                       updatevalue = NA)
        
        ## add layer to new rr:
        rr_new <- addLayer(rr_new, rr_temp[[i]])
        names <- append(names, paste(names(potential_raster), realized_range$source[num], sep = "_"))
        
        num = num + 1
      }
    }
    
    i = i + 1
  }
  names(pr_restricted_bool) <- names(pr_restricted_temps) <- names(rr_new) <- names
  # plot(pr_restricted_bool)
  
  if (type == "acclimatized") {
    names(tlims) <- names
    return(list(pr_restricted_bool, pr_restricted_temps, rr_new, tlims))
  }
  else {
    return(list(pr_restricted_bool, pr_restricted_temps, rr_new))
  }
}

extract_thermal_niche <- function(type, prs, realm) {
  
  rrs <- rasterized_rrs_nichemapr
  names(rrs) <- str_replace_all(names(rrs), '\\.', '_')
  
  prs_high <- prs[[2]]
  prs_low <- prs[[3]]
  
  rrs_high <- prs[[4]]
  rrs_low <- prs[[5]]
  
  if (type == "acclimatized") {
    tlims_high <- prs[[6]]
    tlims_low <- prs[[7]]
  }
  
  range = 1
  while (range < nlayers(prs_high) + 1) {
    
    ## get species info
    ## get realized range and realized range raster
    split <- str_split_fixed(names(prs_high)[[range]], '_', n = 3) 
    species <- paste(split[1,1], split[1,2], sep = '_')
    range_id <- paste(species, split[1,3], sep = '_')
    
    species_traits <- traits[which(traits$genus_species == species),]
    
    cold_dormancy <- species_traits$cold_season_dormancy_ == "Yes" & !is.na(species_traits$cold_season_dormancy_)
    hot_dormancy <- species_traits$hot_season_dormancy_ == "Yes" & !is.na(species_traits$hot_season_dormancy_)
    dormancy <- cold_dormancy || hot_dormancy
    dormancy_type <- ifelse(!dormancy, "none", ifelse(cold_dormancy & hot_dormancy, "both", 
                                                      ifelse(cold_dormancy, "cold", "hot")))
    ## make sure pr is not empty:
    pr_high_empty <- length(which(is.na(values(prs_high[[range]])) == FALSE)) <= 0
    pr_low_empty <- length(which(is.na(values(prs_low[[range]])) == FALSE)) <= 0
    
    ## get temps within potential range 
    if (!pr_high_empty) {
      high_vals <- data.frame(high_or_low = "high", 
                              temps = values(prs_high[[range]])[which(is.na(values(prs_high[[range]])) == FALSE)],
                              type = paste("potential", type, sep = "_"),
                              dormancy = dormancy_type)
    }
    else {
      high_vals <- data.frame(high_or_low = "high", 
                              temps = NA,
                              type = paste("potential", type, sep = "_"),
                              dormancy = dormancy_type)
    }
    if (!pr_low_empty) {
      low_vals <- data.frame(high_or_low = "low", 
                             temps = values(prs_low[[range]])[which(is.na(values(prs_low[[range]])) == FALSE)],
                             type = paste("potential", type, sep = "_"),
                             dormancy = dormancy_type)
    }
    else {
      low_vals <- data.frame(high_or_low = "low", 
                             temps = NA,
                             type = paste("potential", type, sep = "_"),
                             dormancy = dormancy_type)
    }
   
    ## if acclimatized version, extract acclimatized thermal limits:
    if (type == "acclimatized") {
      ## make sure tlims is not empty:
      tlims_high_empty <- length(which(is.na(values(tlims_high[[range]])) == FALSE)) <= 0
      tlims_low_empty <- length(which(is.na(values(tlims_low[[range]])) == FALSE)) <= 0
     
      if (!tlims_high_empty) {
        high_vals <- rbind(high_vals, data.frame(high_or_low = "high", 
                                               temps = values(tlims_high[[range]])[which(is.na(values(tlims_high[[range]]))
                                                                                        == FALSE)],
                                               type = "tlims",
                                               dormancy = dormancy_type))
      }
      else {
        high_vals <- rbind(high_vals, data.frame(high_or_low = "high", 
                                               temps = NA,
                                               type = "tlims",
                                               dormancy = dormancy_type))
      }
      if (!tlims_low_empty) {
        low_vals <- rbind(low_vals, data.frame(high_or_low = "low", 
                               temps = values(tlims_low[[range]])[which(is.na(values(tlims_low[[range]]))
                                                                        == FALSE)],
                               type = "tlims",
                               dormancy = dormancy_type))
      }
      else {
        low_vals <- rbind(low_vals, data.frame(high_or_low = "low", 
                                               temps = NA,
                                               type = "tlims",
                                               dormancy = dormancy_type))
      }
    }
    
    ## get temps within corresponding rasterized realized range
    rr <- rrs[[which(names(rrs) == range_id)]]
    
    rrs_high[[range]] <- mask(rrs_high[[range]], rr, updatevalue = NA)
    rrs_low[[range]] <- mask(rrs_low[[range]], rr, updatevalue = NA)
    
    ## make sure rr is not empty:
    rr_high_empty <- length(which(is.na(values(rrs_high[[range]])) == FALSE)) <= 0
    rr_low_empty <- length(which(is.na(values(rrs_low[[range]])) == FALSE)) <= 0
    
    ## get temps within realized range 
    if (!rr_high_empty) {
      high_vals <- rbind(high_vals, data.frame(high_or_low = "high", 
                                               temps = values(rrs_high[[range]])[which(is.na(values(rrs_high[[range]])) == FALSE)],
                                               type = paste("realized", type, sep = "_"),
                                               dormancy = dormancy_type))
    }
    else {
      high_vals <- rbind(high_vals, data.frame(high_or_low = "high", 
                                               temps = NA,
                                               type = paste("realized", type, sep = "_"),
                                               dormancy = dormancy_type))
    }
    if (!rr_low_empty) {
      low_vals <- rbind(low_vals, data.frame(high_or_low = "low", 
                                             temps = values(rrs_low[[range]])[which(is.na(values(rrs_low[[range]])) == FALSE)],
                                             type = paste("realized", type, sep = "_"),
                                             dormancy = dormancy_type))
    }
    else {
      low_vals <- rbind(low_vals, data.frame(high_or_low = "low", 
                                             temps = NA,
                                             type = paste("realized", type, sep = "_"),
                                             dormancy = dormancy_type))
    }
    
    ## bind dfs together:
    if(range == 1) {
      thermal_niche <- rbind(high_vals, low_vals) %>%
        mutate(range = range_id) %>%
        select(range, everything()) 
    }
    else {
      thermal_niche <- rbind(high_vals, low_vals) %>%
        mutate(range = range_id) %>%
        select(range, everything()) %>%
        rbind(thermal_niche, .)
    }
    
    range = range + 1
  }
  
  thermal_niche$realm <- realm
  
  ## if type is Te, get rid of temperatures exceeding 5th percentile of minimum operative temperatures in the shade:
  if (type %in% c("Te", "Te_in_sun", "acclimatized") & realm == "Terrestrial") {
    niche_type = paste("realized", type, sep = "_")
    
    realized_thermal_niche <- thermal_niche %>%
      filter(type == niche_type) %>%
      filter((high_or_low == "low" & temps > -30) | (high_or_low == "high")) 
    
    thermal_niche <- thermal_niche %>%
      filter(type != niche_type) %>%
      rbind(., realized_thermal_niche)
  }
  
  return(thermal_niche)
}

extract_thermal_niche_one_limit <- function(type, prs, realm, limit_type) {
  
  rrs <- rasterized_rrs_nichemapr
  names(rrs) <- str_replace_all(names(rrs), '\\.', '_')
  
  prs_temps <- prs[[2]]
  
  rrs_temps <- prs[[3]]

  if (type == "acclimatized") {
    tlims_temps <- prs[[4]]
  }
  
  range = 1
  while (range < nlayers(prs_temps) + 1) {
    
    ## get species info
    ## get realized range and realized range raster
    
    split <- str_split_fixed(names(prs_temps)[[range]], '_', n = 3) 
    species <- paste(split[1,1], split[1,2], sep = '_')
    range_id <- paste(species, split[1,3], sep = '_')
  
    species_traits <- traits[which(traits$genus_species == species),]
    
    cold_dormancy <- species_traits$cold_season_dormancy_ == "Yes" & !is.na(species_traits$cold_season_dormancy_)
    hot_dormancy <- species_traits$hot_season_dormancy_ == "Yes" & !is.na(species_traits$hot_season_dormancy_)
    dormancy <- cold_dormancy || hot_dormancy
    dormancy_type <- ifelse(!dormancy, "none", ifelse(cold_dormancy & hot_dormancy, "both", 
                                                      ifelse(cold_dormancy, "cold", "hot")))
    ## make sure pr is not empty:
    pr_empty <- length(which(is.na(values(prs_temps[[range]])) == FALSE)) <= 0
    
    ## get temps within potential range 
    if (!pr_empty) {
      vals <- data.frame(high_or_low = NA, 
                         temps = values(prs_temps[[range]])[which(is.na(values(prs_temps[[range]])) == FALSE)],
                              type = paste("potential", type, sep = "_"),
                              dormancy = dormancy_type)
    }
    else {
      vals <- data.frame(high_or_low = NA, 
                         temps = NA,
                         type = paste("potential", type, sep = "_"),
                         dormancy = dormancy_type)
    }
    
    ## if acclimatized version, extract acclimatized thermal limits:
    if (type == "acclimatized") {
      ## make sure tlims is not empty:
      tlims_empty <- length(which(is.na(values(tlims_temps[[range]])) == FALSE)) <= 0
      
      if (!tlims_empty) {
        vals <- rbind(vals, data.frame(high_or_low = NA, 
                                      temps = values(tlims_temps[[range]])
                                      [which(is.na(values(tlims_temps[[range]])) == FALSE)],
                                                 type = "tlims",
                                                 dormancy = dormancy_type))
      }
      else {
        vals <- rbind(vals, data.frame(high_or_low = NA, 
                                               temps = NA,
                                               type = "tlims",
                                               dormancy = dormancy_type))
      }
    }
    
    ## get temps within corresponding rasterized realized range
    rr <- rrs[[which(names(rrs) == range_id)]]
    
    rrs_temps[[range]] <- mask(rrs_temps[[range]], rr, updatevalue = NA)
    
    ## make sure rr is not empty:
    rr_empty <- length(which(is.na(values(rrs_temps[[range]])) == FALSE)) <= 0
    
    ## get temps within realized range 
    if (!rr_empty) {
      vals <- rbind(vals, data.frame(high_or_low = NA, 
                                          temps = values(rrs_temps[[range]])[which(is.na(values(rrs_temps[[range]])) == FALSE)],
                                               type = paste("realized", type, sep = "_"),
                                               dormancy = dormancy_type))
    }
    else {
      vals <- rbind(vals, data.frame(high_or_low = NA, 
                                             temps = NA,
                                             type = paste("realized", type, sep = "_"),
                                             dormancy = dormancy_type))
    }
    
    if (limit_type == "tmax") {
      vals$high_or_low = "high"
    }
    else if(limit_type == 'tmin') {
      vals$high_or_low = "low"
    }
    
    if(range == 1) {
      thermal_niche <- vals %>%
        mutate(range = range_id) %>%
        select(range, everything()) 
    }
    else {
      thermal_niche <- vals %>%
        mutate(range = range_id) %>%
        select(range, everything()) %>%
        rbind(thermal_niche, .)
    }
    
    range = range + 1
  }
  
  thermal_niche$realm <- realm
  
  ## if type is Te, get rid of temperatures exceeding 5th percentile of minimum operative temperatures in the shade:
  if (type %in% c("Te", "Te_in_sun", "acclimatized") & realm == "Terrestrial" & limit_type == "tmin") {
    niche_type = paste("realized", type, sep = "_")
    
    realized_thermal_niche <- thermal_niche %>%
      filter(type == niche_type) %>%
      filter(high_or_low == "low" & temps > -30)
    
    thermal_niche <- thermal_niche %>%
      filter(type != niche_type) %>%
      rbind(., realized_thermal_niche)
  }
  
  return(thermal_niche)
}

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

