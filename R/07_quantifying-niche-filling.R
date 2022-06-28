### calculating potential thermal niche filling
library(tidyverse)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Read in data           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = " "))


## calculate thermal niche filling for all sensitivity sets separately 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Operative temperatures           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
### assumes species are always in the shade at the hottest hour of the hottest day
### assumes one thermal limit applies across entire range of species 
niche_Te <- read.csv("data-processed/thermal-niches/all-realms/niche_Te.csv")

## split potential and realized
r_niches <- niche_Te %>%
  filter(type == 'realized_Te')
p_niches <- niche_Te %>%
  filter(type == 'potential_Te')

niche_lims_Te <- extract_niche_limits(r_niches, p_niches, thermal_limits, type = "Te")

filling_Te <- calculate_filling(niche_lims_Te)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Operative temperatures with behavioural thermoregulation           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
### tests assumption that species are always in the shade at the hottest hour of the hottest day
### lets species thermoregulate between the sun and shade to reach their Tpref 
## read in Tpref data:
tpref <- read.csv("data-processed/traits/Tpref_clean.csv")

niche_Te_in_sun <- read.csv("data-processed/thermal-niches/all-realms/niche_Te_in_sun.csv")

## split potential and realized
r_niches <- niche_Te_in_sun %>%
  filter(type == 'realized_Te_in_sun')
p_niches <- niche_Te_in_sun %>%
  filter(type == 'potential_Te_in_sun')

niche_lims_Te_in_sun <- extract_niche_limits(r_niches, p_niches, thermal_limits, type = "Te_in_sun")

filling_Te_in_sun <- calculate_filling(niche_lims_Te_in_sun)

## subset niches to only species we have Tpref for
niche_lims_Te_in_sun_sub <- niche_lims_Te_in_sun %>%
  mutate(genus_species = paste(str_split_fixed(range, "\\_", n=3)[,1], 
                               str_split_fixed(range, "\\_", n=3)[,2], sep = "_")) %>%
  filter(genus_species %in% tpref$genus_species)

niche_lims_Te_sub <- niche_lims_Te %>%
  mutate(genus_species = paste(str_split_fixed(range, "\\_", n=3)[,1], 
                               str_split_fixed(range, "\\_", n=3)[,2], sep = "_")) %>%
  filter(genus_species %in% tpref$genus_species)

## merge with tpref data:
niche_lims_Te_sub <- tpref %>%
  mutate(limit_type = ifelse(limit_type == 'max, min', 'tmax, tmin', 
                             ifelse(limit_type == 'min', 'tmin', 'tmax'))) %>%
  left_join(niche_lims_Te_sub, .) %>%
  filter(!is.infinite(r_niche_upper))

## get max Te sun
maxTesun <- niche_lims_Te_in_sun_sub %>%
  select(range, r_niche_upper) %>%
  rename('max_Te_sun' = r_niche_upper) %>%
  unique()
  
niche_lims_t <- left_join(niche_lims_Te_sub, maxTesun)

## we are going to adjust r niche upper and p niche upper based on Tpref
niche_lims_tpref <- niche_lims_t %>%
  # If Tpref > max Te shade, assume animal’s body temperature is Tpref 
  mutate(r_niche_upper = ifelse(Tpref > r_niche_upper, Tpref, r_niche_upper)) %>%
  mutate(p_niche_upper = ifelse(Tpref > p_niche_upper, Tpref, p_niche_upper)) %>%
  # If Tpref > max Te sun, assume animal’s body temperature is max Te sun 
  mutate(r_niche_upper = ifelse(Tpref > max_Te_sun, max_Te_sun, r_niche_upper)) %>%
  mutate(p_niche_upper = ifelse(Tpref > max_Te_sun, max_Te_sun, p_niche_upper)) %>%
  filter(!is.na(Tpref))


filling_Te_tpref <- calculate_filling(niche_lims_tpref)
filling_Te_tpref <- filling_Te_tpref[,-c(12:20)]
filling_Te_tpref$type = "Te_tpref"


## we are going to adjust r niche upper and p niche upper based on Tpref/Tb
niche_lims_tpreftb <- niche_lims_t %>%
  # If Tpref/Tb > max Te shade, assume animal’s body temperature is Tpref 
  mutate(r_niche_upper = ifelse(!is.na(Tpref) & (Tpref > r_niche_upper), Tpref,
                                ifelse(!is.na(Tb) & (Tb > r_niche_upper), Tb, 
                                       r_niche_upper))) %>%
  mutate(p_niche_upper = ifelse(!is.na(Tpref) & (Tpref > p_niche_upper), Tpref,
                                ifelse(!is.na(Tb) & (Tb > p_niche_upper), Tb,
                                       p_niche_upper))) %>%
  # If Tpref/Tb > max Te sun, assume animal’s body temperature is max Te sun 
  mutate(r_niche_upper = ifelse(!is.na(Tpref) & (Tpref > max_Te_sun), max_Te_sun,
                                ifelse(!is.na(Tb) & (Tb > max_Te_sun), Tb, 
                                       r_niche_upper))) %>%
  mutate(p_niche_upper = ifelse(!is.na(Tpref) & (Tpref > max_Te_sun), Tpref,
                                ifelse(!is.na(Tb) & (Tb > max_Te_sun), Tb,
                                       p_niche_upper))) 

filling_Te_tpreftb <- calculate_filling(niche_lims_tpreftb)
filling_Te_tpreftb <- filling_Te_tpreftb[,-c(12:20)]
filling_Te_tpreftb$type = "Te_tpreftb"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Acclimatized thermal limits           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
### tests assumption that one thermal limit applies across entire range of species 
### potential ranges in this version of the niche were created by allowing tmax and tmin to vary spatially accoridng to environmental acclimatization temperatures  
niche_acc <- read.csv("data-processed/thermal-niches/all-realms/niche_acclimatized.csv")

## split potential and realized
r_niches <- niche_acc %>%
  filter(type == 'realized_acclimatized')
p_niches <- niche_acc %>%
  filter(type == 'potential_acclimatized')
thermal_limits <- niche_acc %>%
  filter(type == 'tlims')

niche_lims_acc <- extract_niche_limits(r_niches, p_niches, thermal_limits, type = "acclimatized")

filling_acc <- calculate_filling(niche_lims_acc)

## filter to only the species that could had at least one thermal limit acclimatized:
acc_data <- read.csv("data-processed/traits/acclimation-data.csv")

filling_acc <- filling_acc %>%
  mutate(genus_species = paste(str_split_fixed(range, '_', n = 3)[,1], 
                               str_split_fixed(range, '_', n = 3)[,2], sep = "_")) %>%
  filter(genus_species %in% acc_data$genus_species) %>%
  select(-genus_species)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Exporting data for analysis          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
filling_all <- rbind(filling_Te, filling_Te_tpref, 
                     filling_Te_tpreftb, filling_acc)

## add realized range latitudinal midpoint and geographic area:
realized_ranges <- st_read("large-files/realized-ranges/realized-ranges.shp")
realized_ranges$range <- paste(str_replace_all(realized_ranges$species, " ", "_"), 
                               realized_ranges$source, sep = "_") 

colnames(realized_ranges)[4] <- "range_area_km2"

filling_all <- left_join(filling_all, realized_ranges, by = c("range", "realm")) %>%
  select(-geometry)

write.csv(filling_all, "data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics.csv", row.names = FALSE)



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
