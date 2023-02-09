# script to estimate extreme operative body temperatures of terrestrial species

######################## Install NicheMapR from github ########################
# require(devtools)
# install_github("mrke/NicheMapR")
library(NicheMapR)
library(evobiR)
library(tidyverse)
library(raster)
# download global climate database
# get.global.climate(folder = "large-files/air-temperatures")

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

writeRaster(map, "data-processed/masks/raster_terr_mask_nichemapr.grd", 
            overwrite = TRUE, format = "raster")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####      2. Extract microclimatic conditions in each pixel         #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Run NicheMapR across the world map (~4 hours)
minshade <- 0 # minimum shade level (full sun: 0% of shade)
maxshade <- 90 # maximum shade (full shade: 90%)
Usrhyt <- 0.01 # animal's height (1 cm)

# run NicheMapR to compute hourly estimations of microclimatic conditions in each grid cell
# DOY (day of the year), TIME (minutes), TALOC (Air temperature), RHLOC (relative humidity), VLOC (wind speed), and SOLR (solar radiation)
microclim_data <- list()
for (i in 1:nrow(xy.values)) { 
  start <- Sys.time()
  
  loc <- c(x=xy.values[i,1], y=xy.values[i,2]) # set lon / lat
  tryCatch({
    micro <- micro_global(loc = loc, minshade = minshade, maxshade = maxshade, Usrhyt = Usrhyt, 
                          timeinterval = 365)
    
    micro_sun <- as.data.frame(micro$metout) # meteorological conditions in the sun 
    micro_shade <- as.data.frame(micro$shadmet) # and in the shade
    micro_soil_sun <- as.data.frame(micro$soil) # soil temperature in the sun
    micro_soil_shade <- as.data.frame(micro$shadsoil) # soil temperature in the sun
    
    # Make a dataset with all varables we need
    microclim_data_cell <- data.frame("DOY" = micro_sun$DOY,
                                      "TIME" = micro_sun$TIME,
                                      "S_sun" = micro_sun$SOLR, 
                                      "Ta_sun" = micro_sun$TALOC, 
                                      "Tg_sun" = micro_soil_sun$D0cm,  
                                      "v_sun" = micro_sun$VLOC,   
                                      "RH_sun" = micro_sun$RH,
                                      
                                      "S_shade" = micro_shade$SOLR * 0.1, 
                                      "Ta_shade" = micro_shade$TALOC, 
                                      "Tg_shade" = micro_soil_shade$D0cm,  
                                      "v_shade" = micro_shade$VLOC,   
                                      "RH_shade" = micro_shade$RH)
    
    # Store microclimatic data of each cell
    microclim_data[[i]] <- microclim_data_cell
    
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
  
  end <- Sys.time()
  elapsed <- end - start
  expected <- elapsed[[1]] * (nrow(xy.values)-i) / 60
  print(paste0("remaining: ", round(expected,3), " minutes = ", round(expected/60,3), " hours")) # not super precise... but will give you an idea of the remaining time
  
  print(paste("Done cell number: ", i, sep = ""))
}

#saveRDS(microclim_data, "large-files/operative-temperatures/NicheMapR-microclimate.rds")
microclim_data <- readRDS("large-files/operative-temperatures/NicheMapR-microclimate.rds")
# each entry of the list "microclim_data" contains a dataframe with all the microclimatic variables we need to compute Te at each cell



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        3. Calculate when periods of dormancy start/end         #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Figure out which days to remove for species that are dormant during the hot and cold season
i=1
while(i < length(microclim_data) + 1) {
  
  # open microclimatic data of this cell
  microclim_data_cell <- microclim_data[[i]]
  
  if (is.null(microclim_data_cell)) {
    microclim_data[[i]] <- microclim_data_cell
  }
  else {
    ## find which 6 consecutive months are hottest, which are coldest:
    daily_highs <- microclim_data_cell %>%
      group_by(DOY) %>%
      do(mutate(., Ta_sun = max(.$Ta_sun, na.rm = TRUE))) %>%
      ungroup() %>%
      dplyr::select(DOY, Ta_sun) %>%
      unique(.) %>%
      rbind(., .[1:183,])
    
    daily_lows <- microclim_data_cell %>%
      group_by(DOY) %>%
      do(mutate(., Ta_sun = min(.$Ta_sun, na.rm = TRUE))) %>%
      ungroup() %>%
      dplyr::select(DOY, Ta_sun) %>%
      unique(.) %>%
      rbind(., .[1:183,])
    
    sw_high <- SlidingWindow(daily_highs$Ta_sun, window = 180, FUN = sum, step = 1)
    sw_low <- SlidingWindow(daily_lows$Ta_sun, window = 180, FUN = sum, step = 1)
    
    ## get index of days to block out
    hot_start <- which(sw_high == max(sw_high))[1]
    cold_start <- which(sw_low == min(sw_low))[1]
    
    hot_end <- hot_start + 179
    cold_end <- cold_start + 179
    
    if(hot_end > 365) {
      hot_end = hot_end - 365
    }
    if(cold_end > 365) {
      cold_end = cold_end - 365
    }
    
    microclim_data[[i]] <- microclim_data_cell %>%
      mutate(hot_start = hot_start, hot_end = hot_end, cold_start = cold_start, cold_end = cold_end)
  }
  
  print(paste("Done cell number: ", i, sep = ""))
  i=i+1
}

#saveRDS(microclim_data, "large-files/operative-temperatures/NicheMapR-microclimate_dormancy-indecies.rds")
microclim_data <- readRDS("large-files/operative-temperatures/NicheMapR-microclimate_dormancy-indecies.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        4. Compute Te of each species         #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Compute Te of each species at each hour of the day, every day of the year 
# In each pixel, select highest and lowest Te in the sun and shade and compute acclimation temperatures 
# For species that are dormant, don't include Te during times they are dormant 

## read in species traits and subset to terrestrial species with known body length
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv") %>% 
  filter(Realm == "Terrestrial") %>%
  filter(!is.na(.$maximum_body_size_SVL_HBL_cm_)) 

## convert body length to m
traits$maximum_body_size_SVL_HBL_cm_ <- traits$maximum_body_size_SVL_HBL_cm_/100

Te_allspp <- list()
## for each species:
x=1
while (x < (nrow(traits) + 1)) {
  
  sp_traits <- traits[x,]
  
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
  
  Te_data <- data.frame(xy.values, cells, "maxTe_sun"=NA, "maxTe_shade"=NA, "minTe_sun"=NA,
                        "minTe_shade"=NA,"doy_hot" = NA, "doy_cold" = NA, 
                        "meanTe_sun" = NA, "meanTe_shade" = NA)
  
  for(i in 1:nrow(Te_data)){
    
    # open microclimatic data of this cell
    microclim_data_cell <- microclim_data[[i]]
    
    ## if no climate data for a cell, set all values of Te to NA
    if (is.null(microclim_data_cell)) {
      Te_data$maxTe_sun[i] <- NA
      Te_data$maxTe_shade[i] <- NA
      Te_data$minTe_sun[i] <- NA
      Te_data$minTe_shade[i] <- NA
    }
    else {
      
      ## get dormancy indecies
      hot_start <- microclim_data_cell$hot_start[1]
      cold_start <- microclim_data_cell$cold_start[1]
      
      hot_end <- microclim_data_cell$hot_end[1]
      cold_end <- microclim_data_cell$cold_end[1]
      
      if(hot_end < hot_start) {
        hot_indecies <- c(hot_start:365, 1:hot_end)
      }
      else {
        hot_indecies <- c(hot_start:hot_end)
      }
      if(cold_end < cold_start) {
        cold_indecies <- c(cold_start:365, 1:cold_end)
      }
      else {
        cold_indecies <- c(cold_start:cold_end)
      }
      
      # if seasonal dormancy, block out temps in 6 consecutive warmest/coldest months:
      if(sp_traits$cold_season_dormancy_ == "Yes" & sp_traits$hot_season_dormancy_ == "Yes") {
        
        ## block out temperatures in 6 hottest and coldest consecutive months
        both_d <- microclim_data_cell %>%
          filter(!DOY %in% hot_indecies & !DOY %in% cold_indecies) 
        
        microclim_data_cell <- both_d
      }
      else if (sp_traits$cold_season_dormancy_ == "Yes") {
        
        ## block out temperatures in 6 coldest consecutive months
        c_d <- microclim_data_cell %>%
          filter(!DOY %in% cold_indecies) 
        
        microclim_data_cell <- c_d
      }
      else if (sp_traits$hot_season_dormancy_ == "Yes") {
       
        ## block out temperatures in 6 hottest consecutive months
        h_d <- microclim_data_cell %>%
          filter(!DOY %in% hot_indecies) 
        
        microclim_data_cell <- h_d
      }
      
      # Operative temperature in the sun
      S_sun <- microclim_data_cell$S_sun # solar radiation (Wm-2)
      Ta_sun <- microclim_data_cell$Ta_sun # air temperature (?C)
      Tg_sun <- microclim_data_cell$Tg_sun  # Soil surface temperature (?C)
      v_sun <- microclim_data_cell$v_sun   # Wind velocity (m/s)
      RH_sun <- microclim_data_cell$RH_sun # relative humidity (%)
      
      Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)
      
      # Operative temperature in the shade
      S_shade <- microclim_data_cell$S_shade # solar radiation (Wm-2)
      Ta_shade <- microclim_data_cell$Ta_shade # air temperature (?C)
      Tg_shade <- microclim_data_cell$Tg_shade  # Soil surface temperature (?C)
      v_shade <- microclim_data_cell$v_shade   # Wind velocity (m/s)
      RH_shade <- microclim_data_cell$RH_shade # relative humidity (%)
      
      Te_shade <- Te_function(S_shade, Ta_shade, Tg_shade, v_shade, RH_shade, d, r)
      
      # Extract maximum and minimum Te, mean Te
      Te_data$maxTe_sun[i] <- max(Te_sun)
      Te_data$maxTe_shade[i] <- max(Te_shade)
      Te_data$minTe_sun[i] <- min(Te_sun)
      Te_data$minTe_shade[i] <- min(Te_shade)
      Te_data$meanTe_sun[i] <- mean(Te_sun)
      Te_data$meanTe_shade[i] <- mean(Te_shade)
      
      ## acclimatization temperature:
      doy_hot_extreme <- microclim_data_cell$DOY[which(Te_shade == max(Te_shade))]
      doy_cold_extreme <- microclim_data_cell$DOY[which(Te_shade == min(Te_shade))]
      
      acclim_data = microclim_data[[i]]
      days <- unique(acclim_data$DOY)
      
      ## compute the max daily air temperature 7 days before hottest day
      if (which(days == doy_hot_extreme) <= 7) {
        if (doy_hot_extreme == 1) {
          doys_before_hot = days[365-(7-doy_hot_extreme)]:days[365]
        }
        else {
          doys_before_hot = append(days[365-(7-doy_hot_extreme)]:days[365],days[1]:days[doy_hot_extreme-1])
        }
      }
      else {
        doys_before_hot <- days[which(days == doy_hot_extreme) - 7]:days[which(days == doy_hot_extreme) - 1]
      }
      ## compute the min daily air temperature 7 weeks before the coldest day 
      if (which(days == doy_cold_extreme) <= 49) {
        if (doy_cold_extreme == 1) {
          doys_before_cold = days[365-(49-doy_cold_extreme)]:days[365]
        }
        else {
          doys_before_cold = append(days[365-(49-doy_cold_extreme)]:days[365],
                                    days[1]:days[doy_cold_extreme-1])
        }
      }
      else {
        doys_before_cold <- days[which(days == doy_cold_extreme) - 49]:days[which(days == doy_cold_extreme) - 1]
      }
      
      Te_data$hot_acc_temp[i] <- filter(acclim_data, DOY %in% doys_before_hot) %>%
        summarize(max(Ta_sun, na.rm = TRUE)) %>%
        as.numeric(as.character(.))
      Te_data$cold_acc_temp[i] <- filter(acclim_data, DOY %in% doys_before_cold) %>%
        summarize(min(Ta_sun, na.rm = TRUE))%>%
        as.numeric(as.character(.))
      
      Te_data$doy_hot[i] = doy_hot_extreme
      Te_data$doy_cold[i] = doy_cold_extreme
      
      print(paste("Done cell number: ", i, sep = "", " for species number ", x))
    }
    
  }
  
  Te_allspp[[x]] <- Te_data
  names(Te_allspp)[x] <- as.character(sp_traits$genus_species)
  
  print(paste("Done species number: ", x, sep = ""))
  
  x = x+1
}

#saveRDS(Te_allspp, "large-files/operative-temperatures/Te_allspp.rds")
Te_allspp <- readRDS("large-files/operative-temperatures/Te_allspp.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        5. Map Te of each species             #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Make maps of max and min Te in the sun and shade and acclimation temperatures for each species 
x=1
while (x < (length(Te_allspp) + 1)) {
  Te_data <- Te_allspp[[x]]
  
  map_minTe_sun <- map_maxTe_sun <- map_minTe_shade <- map_maxTe_shade <- map
  map_acc_cold <- map_acc_hot <-map_meanTe <- map
  for(i in 1:nrow(Te_data)){
    cell <- Te_data$cells[i]
    map_minTe_sun[cell] <- Te_data$minTe_sun[i]
    map_maxTe_sun[cell] <- Te_data$maxTe_sun[i]
    map_minTe_shade[cell] <- Te_data$minTe_shade[i]
    map_maxTe_shade[cell] <- Te_data$maxTe_shade[i]
    map_acc_cold[cell] <- Te_data$cold_acc_temp[i]
    map_acc_hot[cell] <- Te_data$hot_acc_temp[i]
    map_meanTe[cell] <- (Te_data$meanTe_sun[i] + Te_data$meanTe_shade[i])/2
  }
  
  #plot(map_maxTe_sun)
  #plot(map_minTe_shade)
  
  if (x == 1) {
    Te_sun_min <- map_minTe_sun
    Te_sun_max <- map_maxTe_sun
    Te_shade_min <- map_minTe_shade
    Te_shade_max <- map_maxTe_shade
    terr_acc_hot <- map_acc_hot
    terr_acc_cold <- map_acc_cold
    Te_mean <- map_meanTe
  }
  else {
    Te_sun_min <- addLayer(Te_sun_min, map_minTe_sun)
    Te_sun_max <- addLayer(Te_sun_max, map_maxTe_sun)
    Te_shade_min <- addLayer(Te_shade_min, map_minTe_shade)
    Te_shade_max <- addLayer(Te_shade_max, map_maxTe_shade)
    terr_acc_hot <- addLayer(terr_acc_hot, map_acc_hot)
    terr_acc_cold <- addLayer(terr_acc_cold, map_acc_cold)
    Te_mean <- addLayer(Te_mean, map_meanTe)
  }
  print(paste("Done species number: ", x, sep = ""))
x=x+1
}

names(Te_sun_min) <- names(Te_sun_max) <- names(Te_shade_min) <- names(Te_shade_max) <-
  names(terr_acc_cold) <- names(terr_acc_hot) <- names(Te_mean) <- names(Te_allspp)

## save the maps: 
writeRaster(Te_sun_min, "large-files/operative-temperatures/Te_sun_min.grd", overwrite = TRUE, format = "raster")
writeRaster(Te_sun_max, "large-files/operative-temperatures/Te_sun_max.grd", overwrite = TRUE, format = "raster")
writeRaster(Te_shade_min, "large-files/operative-temperatures/Te_shade_min.grd", overwrite = TRUE, format = "raster")
writeRaster(Te_shade_max, "large-files/operative-temperatures/Te_shade_max.grd", overwrite = TRUE, format = "raster")
writeRaster(terr_acc_cold, "large-files/operative-temperatures/terr_acc_cold.grd", overwrite = TRUE, format = "raster")
writeRaster(terr_acc_hot, "large-files/operative-temperatures/terr_acc_hot.grd", overwrite = TRUE, format = "raster")
writeRaster(Te_mean, "large-files/operative-temperatures/Te_mean.grd", overwrite = TRUE, format = "raster")
