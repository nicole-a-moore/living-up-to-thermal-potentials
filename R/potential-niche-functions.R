## script containing all functions needed to create potential ranges & niches 
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



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       Metric calculating functions          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## function that takes all temperatures in a species realized and potential thermal niche, ctmax and ctmin,
## and extracts niche limits
## collapses potential_niche to a single minimum and maximum available temperature 
## collapses realized_niche to a single minimum and maximum occupied temperature 
extract_niche_limits <- function(r_niches, p_niches, thermal_limits, type) {
  
  ## loop through each range:
  range = 1
  while (range < length(unique(r_niches$range)) + 1) {
    
    range_id <- unique(r_niches$range)[range]
    
    ## extract the p and r niche temps for range:
    rniche <- r_niches[which(r_niches$range == range_id), ]
    pniche <- p_niches[which(p_niches$range == range_id), ]
    
    ## get species thermal limits
    limit_type = unique(rniche$limit_type)
    
    split <- str_split_fixed(range_id, '_', n = 3) 
    species <- paste(split[1,1], split[1,2], sep = ' ')
    
    ## if acclimatized thermal limits, get fundamental niche limits by taking max and min value
    if (type == "acclimatized") {
      tlims_upper <- thermal_limits %>%
        filter(high_or_low == "high", range == range_id) 
      tlims_lower <- thermal_limits %>%
        filter(high_or_low == "low", range == range_id)
      
      ctmax <- max(tlims_upper$temps, na.rm = TRUE)
      ctmin <- min(tlims_lower$temps, na.rm = TRUE)
    }
    else {
      lims <- thermal_limits[which(thermal_limits$genus_species == species),]
      ctmax <- filter(lims, type == 'max') %>%
        .$thermal_limit
      ctmin <- filter(lims, type == 'min') %>%
        .$thermal_limit
    }
    
    ## extract niche limits:
    r_upper <- max(rniche$temps, na.rm = TRUE)
    r_lower <- min(rniche$temps, na.rm = TRUE)
    p_upper <- max(pniche$temps, na.rm = TRUE)
    p_lower <- min(pniche$temps, na.rm = TRUE)
    
    ## add to df:
    if (range == 1) {
      if (limit_type == 'tmax, tmin') {
        niche_lims <- data.frame(range = range_id, 
                                 type = type,
                                 dormancy = unique(rniche$dormancy),
                                 realm = unique(rniche$realm),
                                 r_niche_upper = r_upper,
                                 p_niche_upper = p_upper,
                                 ctmax = ctmax,
                                 r_niche_lower = r_lower,
                                 p_niche_lower = p_lower, 
                                 ctmin = ctmin, 
                                 limit_type = limit_type)
      }
      else if (limit_type == 'tmax') {
        niche_lims <- data.frame(range = range_id, 
                                 type = type,
                                 dormancy = unique(rniche$dormancy),
                                 realm = unique(rniche$realm),
                                 r_niche_upper = r_upper,
                                 p_niche_upper = p_upper,
                                 ctmax = ctmax,
                                 r_niche_lower = NA,
                                 p_niche_lower = NA, 
                                 ctmin = NA, 
                                 limit_type = limit_type)
      }
      else if (limit_type == 'tmin') {
        niche_lims <- data.frame(range = range_id, 
                                 type = type,
                                 dormancy = unique(rniche$dormancy),
                                 realm = unique(rniche$realm),
                                 r_niche_upper = NA,
                                 p_niche_upper = NA,
                                 ctmax = NA,
                                 r_niche_lower = r_lower,
                                 p_niche_lower = p_lower, 
                                 ctmin = ctmin,
                                 limit_type = limit_type)
      }
    }
    else {
      if (limit_type == 'tmax, tmin') {
        niche_lims <- rbind(niche_lims, data.frame(range = range_id, 
                                                   type = type,
                                                   dormancy = unique(rniche$dormancy),
                                                   realm = unique(rniche$realm),
                                                   r_niche_upper = r_upper,
                                                   p_niche_upper = p_upper,
                                                   ctmax = ctmax,
                                                   r_niche_lower = r_lower,
                                                   p_niche_lower = p_lower, 
                                                   ctmin = ctmin, 
                                                   limit_type = 'tmax, tmin'))
      }
      else if (limit_type == 'tmax') {
        niche_lims <-  rbind(niche_lims, data.frame(range = range_id, 
                                                    type = type,
                                                    dormancy = unique(rniche$dormancy),
                                                    realm = unique(rniche$realm),
                                                    r_niche_upper = r_upper,
                                                    p_niche_upper = p_upper,
                                                    ctmax = ctmax,
                                                    r_niche_lower = NA,
                                                    p_niche_lower = NA, 
                                                    ctmin = NA, 
                                                    limit_type = 'tmax'))
      }
      else if (limit_type == 'tmin') {
        niche_lims <-  rbind(niche_lims, data.frame(range = range_id, 
                                                    type = type,
                                                    dormancy = unique(rniche$dormancy),
                                                    realm = unique(rniche$realm),
                                                    r_niche_upper = NA,
                                                    p_niche_upper = NA,
                                                    ctmax = NA,
                                                    r_niche_lower = r_lower,
                                                    p_niche_lower = p_lower, 
                                                    ctmin = ctmin,
                                                    limit_type = 'tmin'))
      }
    }
    range = range + 1
  }
  
  return(niche_lims)
}


## function that takes a species realized + potential thermal niche limits and calculates:
##    1. warm niche underfilling: °C of available thermal niche above the maximum realized temperature but 
##       below CTmax, multiplied by -1
##    2. warm niche overfilling: °C of occupied thermal niche above the maximum available temperature 
##    3. cold niche underfilling: °C of available thermal niche below the minimum realized temperature but 
##       above CTmin, multiplied by -1 
##    4. cold niche overfilling: °C of occupied thermal niche below the minimum available temperature 
##    5. warm fundamental niche overfilling: °C of occupied thermal niche above the ctmax
##    6. cold fundamental niche overfilling: °C of occupied thermal niche below the ctmin
calculate_filling <- function(niche_lims) {
  ## warm under: (p niche upper - r niche upper) - if negative, assign 0 and if positive, multiply by -1
  ## warm over: (r niche upper - p niche upper) - if negative, assign 0 and if positive, leave
  ## warm fundamental under: (r niche upper - ctmax) - if positive, assign 0 and if negative, multiply by -1
  ## warm fundamental over: (r niche upper - ctmax) - if negative, assign 0 and if positive, leave
  ## cold under: (p niche lower - r niche lower) - if positive, assign 0 and if negative, leave
  ## cold over: (r niche lower - p niche lower) - if positive, assign 0 and if negative, multiply by -1
  ## cold fundamental over: (r niche lower - ctmin) - if positive, assign 0 and if positive, leave
  ## cold fundamental under: (r niche lower - ctmin) - if negative, assign 0 and if positive, multiply by -1
  niche_filling <- niche_lims %>%
    mutate(warm_under = ifelse((p_niche_upper - r_niche_upper) < 0, 0, 
                               -(p_niche_upper - r_niche_upper))) %>%
    mutate(warm_over = ifelse((r_niche_upper - p_niche_upper) < 0, 0, (r_niche_upper - p_niche_upper))) %>%
    mutate(f_warm_over = ifelse((r_niche_upper - ctmax) < 0, 0, (r_niche_upper - ctmax))) %>%
    mutate(f_warm_under = ifelse((r_niche_upper - ctmax) > 0, 0, -(r_niche_upper - ctmax))) %>%
    mutate(cold_under = ifelse((p_niche_lower - r_niche_lower) > 0, 0, 
                               (p_niche_lower - r_niche_lower))) %>%
    mutate(cold_over = ifelse((r_niche_lower - p_niche_lower) > 0, 0, 
                              -(r_niche_lower - p_niche_lower))) %>%
    mutate(f_cold_over = ifelse((r_niche_lower - ctmin) > 0, 0, (r_niche_upper - ctmin))) %>%
    mutate(f_cold_under = ifelse((r_niche_lower - ctmin) < 0, 0, -(r_niche_lower - ctmin))) 
  
  return(niche_filling)
}

