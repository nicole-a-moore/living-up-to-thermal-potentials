## calculating potential range underfilling and overfilling 
library(tidyverse)
library(raster)
select <- dplyr::select

## read in rasterized realized ranges:
rrs <- readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds")
names(rrs) <- str_replace_all(names(rrs), '\\.', '_')

## read in potential ranges:
potential_ranges_te <- readRDS("data-processed/potential-ranges/all-realms/potential_ranges_filling_script.rds") 
potential_ranges_acc <- readRDS("data-processed/potential-ranges/all-realms/potential_ranges_filling_script_acc.rds")

## read in traits:
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = "_"))

calculate_range_filling <- function(potential_ranges, type) {
  prs <- potential_ranges[[1]]
  
  # ## get underlying mean temperatures:
  # mean_Te <- stack("large-files/operative-temperatures/Te_mean.grd")
  # mean_marine <- stack("data-processed/temperature-data/marine/raster_marine_mean.grd")
  # mean_intertidal <- stack("data-processed/temperature-data/intertidal/raster_intertidal_mean.grd")

  ## get underlying hot/cold extreme temperatures:
  hot_temps_Te <- stack("large-files/operative-temperatures/Te_shade_max.grd")
  cold_temps_Te <- stack("large-files/operative-temperatures/Te_shade_min.grd")
  hot_temps_marine <- stack('data-processed/temperature-data/marine/raster_marine_high.grd')
  cold_temps_marine <- stack('data-processed/temperature-data/marine/raster_marine_low.grd')
  hot_temps_intertidal <- stack('data-processed/temperature-data/intertidal/raster_intertidal_high.grd')
  cold_temps_intertidal <- stack('data-processed/temperature-data/intertidal/raster_intertidal_low.grd')
  
  i = 1
  while (i < nlayers(prs) + 1) {
    ## for each pair of realized and potential thermal ranges:
    pr <- prs[[i]]
    
    ## get info about ranges:
    split <- str_split_fixed(names(prs)[[i]], '_', n = 3) 
    species <- paste(split[1,1], split[1,2], sep = '_')
    range_id <- paste(species, split[1,3], sep = '_')
    realm <- traits$Realm[first(which(species == traits$genus_species))]
    
    rr <- rrs[[which(names(rrs) == range_id)]]
    
    ## make sure potential range is not empty:
    pr_empty <- length(which(is.na(values(pr)) == FALSE)) <= 0
    
    ## if potential range is not empty: 
    if (!pr_empty) {
      
      ## calculate total # of cells in potential + realized range: 
      pr_count = length(which(!is.na(values(pr))))
      rr_count = length(which(!is.na(values(rr))))
      union_of_cells <- pr_count + rr_count
      
      ## divide cells in potential + realized range by total # of cells
      pr = pr/union_of_cells
      rr = rr/union_of_cells
      
      ## assign cells in rr not in pr a value of 0 and vice versa 
      union_pr <- pr
      union_rr <- rr
      union_pr[is.na(union_pr) & !is.na(union_rr)] <- 0
      union_rr[is.na(union_rr) & !is.na(union_pr)] <- 0
      
      ## subtract values in rr from pr and take absolute value:
      # this gets rid of areas of 'filling'
      diff <- abs(union_pr - union_rr)
      
      ## 1. calculate range overfilling
      ## mask by realized range to get rid of areas of overfilling:
      D_o_grid <- mask(diff, rr)
      #plot(D_o_grid)
      
      ## count cells
      cells_o = length(which(!is.na(values(D_o_grid)) & values(D_o_grid) != 0))
      
      ## 2. calculate range underfilling
      ## mask by potential range to get rid of areas of overfilling:
      D_u_grid <- mask(diff, pr)
      #plot(D_u_grid)
      
      ## count cells
      cells_u = length(which(!is.na(values(D_u_grid)) & values(D_u_grid) != 0))
      
      ## 3. calculate range filling in equatorward and poleward range halves
      ## split cells in half at midpoint of filled potential range
      diff[diff != 0] <- NA
      if (length(which(!is.na(values(diff)))) != 0) {
        mid <- data.frame(rasterToPoints(diff))
        mid_lat <- (max(mid$y) + min(mid$y))/2
        
        fill_rast <- D_u_grid
        potential_rast <- fill_rast
        potential_rast[!is.na(potential_rast)] <- 1
        fill_rast[fill_rast != 0] <- NA # make filled areas have a 1, everything else NA
        f_pts <- data.frame(rasterToPoints(fill_rast))
        p_pts <- data.frame(rasterToPoints(potential_rast))
        
        ## if above equator:
        if (mid_lat > 0) {
          ## poleward filling = # cells filled in poleward / number cells total in poleward
          pol <- length(which(f_pts$y >= mid_lat))/length(which(p_pts$y >= mid_lat))
          ## equatorward:
          equ <- length(which(f_pts$y <= mid_lat))/length(which(p_pts$y <= mid_lat))
        }
        else {
          ## poleward filling = # cells filled in poleward / number cells total in poleward
          pol <- length(which(f_pts$y <= mid_lat))/length(which(p_pts$y <= mid_lat))
          ## equatorward:
          equ <- length(which(f_pts$y >= mid_lat))/length(which(p_pts$y >= mid_lat))
        }
      }
      else {
        pol = equ = NA 
      }
      
      ## 4. extract temps in areas of range underfilling and filling
      ## combine temps in rr and pr:
      if (realm == 'Terrestrial') {
        # temps = mean_Te[[which(names(mean_Te) == species)]]
        hot_temps = hot_temps_Te[[which(names(hot_temps_Te) == species)]]
        cold_temps = cold_temps_Te[[which(names(cold_temps_Te) == species)]]
      }
      else if (realm == 'Marine') {
        # temps = mean_marine
        hot_temps = hot_temps_marine
        cold_temps = cold_temps_marine
      }
      else if (realm == 'Intertidal') {
        # temps = mean_intertidal
        hot_temps = hot_temps_intertidal
        cold_temps = cold_temps_intertidal
      }
      
      ## get rid of areas of filling in underfilling grid:
      D_u_grid[D_u_grid == 0] <- NA
      ## get rid of areas of overfilling in overfilling grid (to get only filled areas)
      D_f_grid <- rr
      D_f_grid[D_f_grid != pr] <- 0
      D_f <- sum(values(D_f_grid),na.rm = TRUE)
      D_u <- sum(values(D_u_grid),na.rm = TRUE)
      
      # temps_f <- values(temps)[which(!is.na(values(mask(temps, D_f_grid))))]
      # temps_u <- values(temps)[which(!is.na(values(mask(temps, D_u_grid))))]
      
      hot_temps_f <- values(hot_temps)[which(!is.na(values(mask(hot_temps, D_f_grid))))]
      hot_temps_u <- values(hot_temps)[which(!is.na(values(mask(hot_temps, D_u_grid))))]
      
      cold_temps_f <- values(cold_temps)[which(!is.na(values(mask(cold_temps, D_f_grid))))]
      cold_temps_u <- values(cold_temps)[which(!is.na(values(mask(cold_temps, D_u_grid))))]
      
      ## deal with ranges that have no over/underfilling or no known temps in the areas:
      if (D_f == 0) {
        # temps_f = 
        hot_temps_f = cold_temps_f = NA
      }
      if (D_u == 0) {
        # temps_u = 
        hot_temps_u = cold_temps_u = NA
      }

    }
    else {
      #temps_f = temps_u = 
      D_f = D_u = hot_temps_f = cold_temps_f = hot_temps_u = cold_temps_u = 
        cells_o = cells_u = rr_count = pr_count = pol = equ = NA 
    }
    if (i == 1) {
      rangefilling <- data.frame(range = range_id, 
                                 species = species, 
                                 type = type,
                                 o_cells = cells_o,
                                 u_cells = cells_u,
                                 rr_cells = rr_count,
                                 pr_cells = pr_count, 
                                 pol_fill = pol, 
                                 equ_fill = equ)
      
      temps_under <- data.frame(range = range_id, 
                                species = species, 
                                type = type,
                                u_or_o = "under",
                               # temps = temps_u, 
                                hot_temps = hot_temps_u, 
                                cold_temps = cold_temps_u, 
                                o_cells = cells_o,
                                u_cells = cells_u,
                                rr_cells = rr_count,
                                pr_cells = pr_count)
      
      temps_filled <- data.frame(range = range_id, 
                               species = species,
                               type = type,
                               u_or_o = "over",
                              # temps = temps_f, 
                               hot_temps = hot_temps_f, 
                               cold_temps = cold_temps_f, 
                               o_cells = cells_o,
                               u_cells = cells_u,
                               rr_cells = rr_count,
                               pr_cells = pr_count)
    }
    else {
      rangefilling <- rbind(rangefilling, data.frame(range = range_id, 
                                                     species = species, 
                                                     type = type,
                                                     o_cells = cells_o,
                                                     u_cells = cells_u,
                                                     rr_cells = rr_count,
                                                     pr_cells = pr_count,
                                                     pol_fill = pol, 
                                                     equ_fill = equ))
      
      temps_under <- rbind(temps_under, data.frame(range = range_id, 
                                                   species = species, 
                                                   type = type,
                                                   u_or_o = "under",
                                                  # temps = temps_u, 
                                                   hot_temps = hot_temps_u, 
                                                   cold_temps = cold_temps_u, 
                                                   o_cells = cells_o,
                                                   u_cells = cells_u,
                                                   rr_cells = rr_count,
                                                   pr_cells = pr_count))
      
      temps_filled <- rbind(temps_filled, data.frame(range = range_id,
                                                 species = species, 
                                                 type = type,
                                                 u_or_o = "over",
                                                # temps = temps_f, 
                                                 hot_temps = hot_temps_f, 
                                                 cold_temps = cold_temps_f, 
                                                 o_cells = cells_o,
                                                 u_cells = cells_u,
                                                 rr_cells = rr_count,
                                                 pr_cells = pr_count))
    }
    print(paste("Finished range number:", i))
    i = i + 1 
  }
  
  if (type == "acclimatized") {
    acc_data <- read.csv("data-processed/traits/acclimation-data.csv")
    rangefilling <- filter(rangefilling, species %in% acc_data$genus_species)
    temps_under <- filter(temps_under, species %in% acc_data$genus_species)
    temps_filled <- filter(temps_filled, species %in% acc_data$genus_species)
  }
  
  return(list(rangefilling, temps_under, temps_filled))
}

rf_te <- calculate_range_filling(potential_ranges = potential_ranges_te, type = "Te")
rf_acc <- calculate_range_filling(potential_ranges = potential_ranges_acc,
                                  type = "acclimatized")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Exporting data for analysis          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## combine data from all types 
rangefilling <- rbind(rf_te[[1]], rf_acc[[1]])
temps_under <- rbind(rf_te[[2]], rf_acc[[2]])
temps_filled <- rbind(rf_te[[3]], rf_acc[[3]])

niche_data <- read.csv("data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics.csv") %>%
  select(range, realm, species, source, range_area_km2, lat_mp, type, ctmax, ctmin) %>%
  filter(!duplicated(.)) %>%
  mutate(species = str_replace_all(species, " ", '_'))

## add other information about ranges:
rangefilling <- left_join(rangefilling, niche_data)
temps_under <- left_join(temps_under, niche_data)
temps_filled <- left_join(temps_filled, niche_data)

## save:
write.csv(rangefilling, "data-processed/potential-ranges/range-filling/rangefilling.csv", row.names = FALSE)
write.csv(temps_under, "data-processed/potential-ranges/range-filling/temps_under.csv", row.names = FALSE)
write.csv(temps_filled, "data-processed/potential-ranges/range-filling/temps_filled.csv", row.names = FALSE)
