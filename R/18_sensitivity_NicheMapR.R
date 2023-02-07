## script to:
# 1) validate operative temperature estimates with field body temperature 
# 2) test the sensitivity of results to parameters chosen in biophysical models
library(tidyverse)
library(raster)
library(readxl)
library(NicheMapR)
select = dplyr::select

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######    1) validate operative temperatures       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in mean Te for each species across sun and shade 
meanTe <- stack("large-files/operative-temperatures/Te_mean.grd")

## read in field body temperature estimates from Algar et al 2018
algar <- read_excel("data-raw/traits/Algar_Morley_Boyd_Tb_microclimate_data.xlsx")

## see which species are in algar
sp <- str_replace_all(names(meanTe), "\\_", " ")

length(which(sp %in% algar$species)) ## 76 species 

al_sub <- filter(algar, species %in% sp)

## for these species, see if max Te at collection location during months of collection matches mean Tb during months of collection 
## make column for start month of collection, end month of collection:
al_sub = al_sub %>%
  mutate(start_month = str_split_fixed(months, "\\-", 2)[,1],
         end_month = str_split_fixed(months, "\\-", 2)[,2]) %>%
  mutate(start_month = ifelse(start_month == "U", NA, start_month),
         end_month = ifelse(start_month == "U", NA, 
                            ifelse(end_month == "", NA, end_month))) %>%
  mutate(genus_species = str_replace_all(species, " ", "_")) %>%
  mutate(start_month = ifelse(start_month == "Sept", "Sep", start_month)) %>%
  mutate(end_month = ifelse(end_month == "Jul", "July", end_month)) %>%
  mutate(num = 1:nrow(.))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######          OPERATIVE TEMPERATURE MODEL                #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# This function uses microclimatic data to estimate operative (=equilibrium) temperature of an animal in the sun or in the shade
Te_function <- function(S,    # Solar radiation (Wm-2)
                        Ta,   # Air temperature (?C)
                        Tg,   # Soil temperature (?C)
                        v,    # Wind speed (m/s)
                        RH,   # Relative humidity (%)
                        d,    # Body length (m)
                        r,    # Total resistance to water loss (s m-1)
                        alpha_lw=0.965, # Skin absorbance (long wave) (Buckley 2007)
                        alpha_s=0.9,    # Skin absorbance (short wave) (Buckley 2007)
                        eps=0.965,      # Skin IR emissivity (Buckley 2007)
                        Fa=0.5,         # View factor from sky (Algar et al. 2018)
                        Fg=0.5){        # View factor from the ground (Algar et al. 2018)
  
  ## Emission and absorption of long-wave radiation
  sigma = 5.67e-8 # W m-2 K-4, Stefan-Boltzmann constant
  eps_sky = 9.2e-6 * (Ta + 273)^2 # Clear sky emissivity of long wave radiation (Buckley 2007)
  La = eps_sky * sigma * (Ta + 273)^4 # long wave radiation from the sky
  
  eps_g = 0.965 # emisivity of the soil surface (Algar et al. 2018)
  Lg = eps_g * sigma * (Tg + 273)^4 # long wave radiation from the soil surface
  
  Rlw = alpha_lw * (Fa * La + Fg * Lg) # absorbed long-wave radiation
  
  ## Absorption of short-wave solar radiation
  # This is a simplification of Lauren's Buckley's (2007) model for the geometry of a lizard. 
  # Here we  assume that the the upper half of the body receives both direct and scattered solar radiation
  Rsol = alpha_s * Fa * S 
  
  ## Evaporative cooling
  ps_a = exp(77.3450 + 0.0057 * (Ta+273) - 7235 / (Ta+273)) / (Ta+273)^8.2  # Air vapor pressure (Pa)
  rho_saturated = 2.2 * ps_a / (Ta+273) # Trasform into density (g m-3)
  RH_prop <- RH * 1e-2
  
  EWL = (rho_saturated - RH_prop * rho_saturated) / r  # Evaporative water loss (g s-1) (Spotila and Berman 1976)
  Qewl = 2257 * EWL # Evaporative cooling (W) = latent heat of vaporization of water (J g-1) x EWL (g s-1)
  
  ## Operative temperature
  cp = 29.3 # J mol-1 K-1, specific heat of the air
  Te = Ta + (Rsol + Rlw - eps * sigma * (Ta + 273)^4 - Qewl) / (4 * sigma * (Ta + 273)^3 + cp*(1.4 + 0.135*sqrt(v/d)))
  
  return(Te)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####      1. Prepare map      #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## map of terrestrial areas:
map <- raster("data-processed/masks/raster_terr_mask.gri")
regions0 <- shapefile("data-raw/polygons/Shapes/newRealms.shp") # this is a shapefile to crop terrestrial regions
regions <- aggregate(rbind(regions0))
map <- mask(map, regions)

plot(map)

xy <- rasterToPoints(map)
cells <- which(values(map) == 1)
xy.values <- xy[,1:2]
nrow(xy.values) # n cells - 14912

## read in microclimate data 
microclim_data <- readRDS("large-files/operative-temperatures/NicheMapR-microclimate.rds")

## read in species traits and subset to terrestrial species with known body length
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv") %>% 
  filter(Realm == "Terrestrial") %>%
  filter(!is.na(.$maximum_body_size_SVL_HBL_cm_)) 

## convert body length to m
traits$maximum_body_size_SVL_HBL_cm_ <- traits$maximum_body_size_SVL_HBL_cm_/100

Te_allspp <- c()

## loop through species and calculate monthly max Te in sun in cells it was collected in
x=1
while (x < (nrow(al_sub)) + 1) {
  
  cur_obs <- al_sub[x,]
  sp_traits <- filter(traits, genus_species == cur_obs$genus_species)
  
  ## set species-specific Te model parameters 
  # body length (m)
  d <- sp_traits$maximum_body_size_SVL_HBL_cm_
  
  # skin resistance to water loss 
  if (sp_traits$Class == "Amphibia") {
    # approx. 300 sm-1 for amphibians
    r <- 300
  }
  else {
    # 6e5 sm-1 for lizards and other dry-skinned ectotherms
    r = 6e5
  }
  
  Te_data <- data.frame(xy.values, cells)
  
  ## filter to cell where species was sampled 
  lat = Te_data$y[which(abs(Te_data$y - cur_obs$Y) == min(abs(Te_data$y - cur_obs$Y)))]
  lon = Te_data$x[which(abs(Te_data$x - cur_obs$X) == min(abs(Te_data$x - cur_obs$X)))]
  
  Te_data <- filter(Te_data, x == unique(lon) & y == unique(lat))
  
  ## if no cells found, skip
  if(nrow(Te_data) == 0) {
    x = x + 1
  }
  else {
    i = which(xy.values[,1] == unique(lon) & xy.values[,2] == unique(lat))
    
    # open microclimatic data of this cell
    if(length(i) == 0) {
      ## if no climate data for a cell, skip
      x = x + 1
    }
    else {
      microclim_data_cell <- microclim_data[[i]]
      
      # Operative temperature in the sun
      S_sun <- microclim_data_cell$S_sun # solar radiation (Wm-2)
      Ta_sun <- microclim_data_cell$Ta_sun # air temperature (?C)
      Tg_sun <- microclim_data_cell$Tg_sun  # Soil surface temperature (?C)
      v_sun <- microclim_data_cell$v_sun   # Wind velocity (m/s)
      RH_sun <- microclim_data_cell$RH_sun # relative humidity (%)
      
      microclim_data_cell$Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)
      
      # Operative temperature in the shade
      S_shade <- microclim_data_cell$S_shade # solar radiation (Wm-2)
      Ta_shade <- microclim_data_cell$Ta_shade # air temperature (?C)
      Tg_shade <- microclim_data_cell$Tg_shade  # Soil surface temperature (?C)
      v_shade <- microclim_data_cell$v_shade   # Wind velocity (m/s)
      RH_shade <- microclim_data_cell$RH_shade # relative humidity (%)
      
      microclim_data_cell$Te_shade <- Te_function(S_shade, Ta_shade, Tg_shade, v_shade, RH_shade, d, r)
      
      # Extract Te across each month 
      months_key = data.frame(DOY = 1:365, month = factor(c(rep("Jan", 31),
                                                            rep("Feb", 28),
                                                            rep("March", 31),
                                                            rep("April", 30),
                                                            rep("May", 31),
                                                            rep("June", 30),
                                                            rep("July", 31),
                                                            rep("Aug", 31),
                                                            rep("Sep", 30),
                                                            rep("Oct", 31),
                                                            rep("Nov", 30),
                                                            rep("Dec", 31)), ordered = T,
                                                          levels = c("Jan", "Feb", "March", "April", "May", "June", 
                                                                     "July", "Aug",
                                                                     "Sep", "Oct", "Nov", "Dec")))
      
      microclim_data_cell <- left_join(microclim_data_cell, months_key, by = "DOY") 
      
      ## subset data to months of collection
      # mos = unique(months_key$month)
      # s_i = which(mos == cur_obs$start_month)
      # e_i = which(mos == cur_obs$end_month)
      
      ## if no collection dates, skip 
      # if(length(s_i) == 0 & length(e_i) == 0) {
      #   x = x + 1
      # }
      # else {
        #  if(length(s_i) != 0 & length(e_i) == 0) {
        #   months = mos[s_i]
        # }
        # else {
        #   months = mos[s_i:e_i]
        # }
        # 
        microclim_data_cell <- microclim_data_cell %>%
          # filter(month %in% months) %>%
          select(Te_sun, Te_shade)
        
        ## save  
        Te_data <- cbind(cur_obs, microclim_data_cell)
        
        ## add to dataframe
        Te_allspp <- rbind(Te_allspp, Te_data)
        
        x = x + 1
        # }
    }
    
    print(paste("Done species number: ", x, sep = ""))
  }
  
}

saveRDS(Te_allspp, "data-processed/intermediate-files/Tb_Te_comparison.rds")
Te_allspp <- readRDS("data-processed/intermediate-files/Tb_Te_comparison.rds")

## plot mean Tb field across latitude, showing the distribution of Te values at each collection location
te_dots <- select(Te_allspp, -c(Te_sun, Te_shade)) %>%
  unique()

te_tb_fig <- Te_allspp %>%
  gather(key = "Te_type", value = "Te", c(Te_sun, Te_shade)) %>%
  ggplot(aes(x = abs(Y), y = Te, group = num)) + 
  geom_point(shape = 95, alpha = 0.2) + 
  geom_point(data = te_dots, aes(y = meanTb), colour = "red", size = 0.9) + ## add dot for mean Tb
  theme_classic() +
  labs(y = "Temperature (°C)", x = "Absolute latitude of collection location (°N)")

ggsave(te_tb_fig, path = "figures/extended-data/", filename = "Te-Tb-fig.png", 
       width = 6, height = 3, device = "png")

te_tb_fig_diff <- Te_allspp %>%
  gather(key = "Te_type", value = "Te", c(Te_sun, Te_shade)) %>%
  group_by(num) %>%
  mutate(max_Te = max(Te)) %>%
  mutate(min_Te = min(Te)) %>%
  select(species, max_Te, min_Te, num, X, Y, meanTb) %>%
  gather(key = "Te_type", value = "ext_Te", c(max_Te, min_Te)) %>%
  distinct() %>%
  ggplot(aes(x = abs(Y), y = ext_Te - meanTb, colour = Te_type)) + 
  geom_point() + 
  scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Hottest Te", "Coldest Te")) +
  theme_classic() +
  labs(y = "Difference between extreme Te\nand mean Tb at collection location (°C)", 
       x = "Absolute latitude of collection location (°N)", 
       colour = "Extreme Te:")

ggsave(te_tb_fig_diff, path = "figures/extended-data/", filename = "Te-Tb-fig-diff.png", 
       width = 6, height = 3, device = "png")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######     2) test sensitivity to parameters       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## for one species, compare the distribution of extreme operative temperatures used in our analysis to those generated when habitat characteristics, skin absorbance, and burrow depth are allowed to vary
## vary type of soil (rock, sand, loam), soil reflectiveness (25%, 50%, 75%), and burrow depth (10cm, 20cm, 30cm)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####      2. Extract microclimatic conditions in each pixel         #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Run NicheMapR across the world map (~4 hours)
minshade <- 0 # minimum shade level (full sun: 0% of shade)
maxshade <- 90 # maximum shade (full shade: 90%)
Usrhyt <- 0.01 # animal's height (1 cm)

soiltype <- c(0, 1, 4)
s <- c("rock", "sand", "loam")
REFL <- c(0.25, 0.5, 0.75)
DEP <- c(0, 2.5,  5,  10,  15,  20,  30,  40,  50,  60)

# run NicheMapR to compute hourly estimations of microclimatic conditions in each grid cell
# DOY (day of the year), TIME (minutes), TALOC (Air temperature), RHLOC (relative humidity), VLOC (wind speed), and SOLR (solar radiation)
for (soil in soiltype) {
  for (refl in REFL) {
    
    microclim_data <- list()
    i=1
    while (i <= nrow(xy.values)) { 
      start <- Sys.time()
      
      loc <- c(x=xy.values[i,1], y=xy.values[i,2]) # set lon / lat
      
      
        micro <- micro_global(loc = loc, minshade = minshade, maxshade = maxshade, Usrhyt = Usrhyt, 
                              timeinterval = 365, 
                              soiltype = soiltype[soil], 
                              REFL = REFL[refl], 
                              DEP = DEP)
        
        micro_sun <- as.data.frame(micro$metout) # meteorological conditions in the sun 
        micro_shade <- as.data.frame(micro$shadmet) # and in the shade
        micro_soil_sun <- as.data.frame(micro$soil) # soil temperature in the sun
        micro_soil_shade <- as.data.frame(micro$shadsoil) # soil temperature in the sun
        
        # Make a dataset with all varables we need
        microclim_data_cell <- data.frame("DOY" = micro_sun$DOY,
                                          "TIME" = micro_sun$TIME,
                                          "S_sun" = micro_sun$SOLR, 
                                          "Ta_sun" = micro_sun$TALOC, 
                                          "Tg_sun_0cm" = micro_soil_sun$D0cm,  
                                          "Tg_sun_10cm" = micro_soil_sun$D10cm,
                                          "Tg_sun_20cm" = micro_soil_sun$D20cm,
                                          "Tg_sun_30cm" = micro_soil_sun$D30cm,
                                          "Tg_sun_40cm" = micro_soil_sun$D40cm,
                                          "Tg_sun_50cm" = micro_soil_sun$D50cm,
                                          "v_sun" = micro_sun$VLOC,   
                                          "RH_sun" = micro_sun$RH,
                                          "S_shade" = micro_shade$SOLR * 0.1, 
                                          "Ta_shade" = micro_shade$TALOC, 
                                          "Tg_shade" = micro_soil_shade$D0cm,  
                                          "v_shade" = micro_shade$VLOC,   
                                          "RH_shade" = micro_shade$RH)
        
        # Store microclimatic data of each cell
        microclim_data[[i]] <- microclim_data_cell
        
      
      
      end <- Sys.time()
      elapsed <- end - start
      expected <- elapsed[[1]] * (nrow(xy.values)-i) / 60
      print(paste0("remaining: ", round(expected,3), " minutes = ", round(expected/60,3), " hours")) # not super precise... but will give you an idea of the remaining time
      
      print(paste("Done cell number: ", i, sep = ""))
      i=i+1
    }
    ## save
    saveRDS(microclim_data, paste("large-files/operative-temperatures/NicheMapR-microclimate_sensitivity_refl-", REFL[refl], 
                                  "_soil-", s[soil], ".rds", sep = ""))
    # each entry of the list "microclim_data" contains a dataframe with all the microclimatic variables we need to compute Te at each cell
  }
}









## run the model many times and generate not just one body temperature but a mean and a range
# discuss whether the conclusions are the same for other temperatures within that range of uncertainty













