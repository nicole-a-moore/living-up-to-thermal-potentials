## creating minimum reproducible dataset

## read in niche filling data 
cold = readRDS("data-processed/thermal-niches/niche-filling/cold-complete-cases_unscaled.rds")
warm = readRDS("data-processed/thermal-niches/niche-filling/warm-complete-cases_unscaled.rds")

## read in range filling data 
range = read.csv("data-processed/potential-ranges/range-filling/range-complete-cases_unscaled.csv")

colnames(warm)
colnames(cold)
colnames(range)

## rename some columns 
cold <- cold %>%
  rename("NicheFilling_Cool" = filling_value,
         "Metric_Tmin" = metric)
warm <- warm %>%
  rename("NicheFilling_Warm" = filling_value,
         "Metric_Tmax" = metric)
range <- range %>%
  rename("LogRangeFilling" = log_prop_occupied,
         "EquatorwardBiasInUnderfilling" = bias_in_uf) %>%
  select(-log_pol_equ_diff, -metric, range) 

## join them by species 
niche = full_join(cold, warm)

full = full_join(range, niche)

length(unique(paste(full$Genus, full$Species)))
## note: some species have multiple rows since different realized range was used to calculate rangefilling 

## now add: 
## warm niche filling given behaviour 
behav <- readRDS("data-processed/potential-ranges/range-filling/warm-complete-cases_unscaled_behaviour.rds")

colnames(behav)

behav <- behav %>%
  rename("NicheFilling_Warm_Behaviour" = filling_value,
         "Metric_Tmax" = metric) %>%
  select(range, Genus, Species, Order, Family, Class, NicheFilling_Warm_Behaviour,
         Metric_Tmax)

full <- left_join(full, behav)

## niche filling given acclimation 
acc_warm <- readRDS("data-processed/potential-ranges/range-filling/warm-complete-cases_unscaled_acc.rds")
acc_cold <- readRDS("data-processed/potential-ranges/range-filling/cold-complete-cases_unscaled_acc.rds")

acc_warm <- acc_warm %>%
  rename("NicheFilling_Warm_Acclimation" = filling_value,
         "Metric_Tmax" = metric) %>%
  select(range, Genus, Species, Order, Family, Class, NicheFilling_Warm_Acclimation,
         Metric_Tmax)

acc_cold <- acc_cold %>%
  rename("NicheFilling_Cool_Acclimation" = filling_value,
         "Metric_Tmin" = metric) %>%
  select(range, Genus, Species, Order, Family, Class, NicheFilling_Cool_Acclimation,
         Metric_Tmin)

full <- left_join(full, acc_warm) %>%
  left_join(., acc_cold)

## range filling given acc
range_acc <- readRDS("data-processed/potential-ranges/range-filling/range-complete-cases_unscaled_acc.rds")

range_acc <- range_acc %>%
  rename("LogRangeFilling_Acclimation" = log_prop_occupied,
         "EquatorwardBiasInUnderfilling_Acclimation" = bias_in_uf) %>%
  select(-metric, range) 

full <- left_join(full, range_acc)

## check sample sizes - full data
length(which(!is.na(full$NicheFilling_Cool))) ## 227 cold niche filling
length(which(!is.na(full$NicheFilling_Warm))) ## 382 cold niche filling
length(which(!is.na(full$LogRangeFilling))) ## 156 range filling
length(which(!is.na(full$EquatorwardBiasInUnderfilling))) ## 156 range filling

## check sample sizes - subsets
length(which(!is.na(full$NicheFilling_Warm_Behaviour))) ## 219 warm niche filling behaviour subset
length(which(!is.na(full$NicheFilling_Warm_Acclimation))) ## 212 warm niche filling acc subset
length(which(!is.na(full$NicheFilling_Cool_Acclimation))) ## 134 warm niche filling acc subset
length(which(!is.na(full$EquatorwardBiasInUnderfilling_Acclimation))) ## 90 warm niche filling acc subset

## clean up and reorder the column names 
colnames(full)

full <- full %>%
  select(-log_range_area) %>%
  rename("Realm" = realm,
         "AbsLatitude" = abs_lat_mp,
         "DispersalDistance" = dispersal_distance_continuous,
         "RangeSize" = rr_cells,
         "ColdSeasonDormancy" = cold_season_dormancy_,
         "HotSeasonDormancy" = hot_season_dormancy_,
         "RealizedRangeSource" = range, 
         "LogBodySize" = log_maximum_body_size) %>%
  mutate(RealizedRangeSource = str_split_fixed(.$RealizedRangeSource, "\\_", 3)[,3]) 

## add tmin and tmax 
thermal <- thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv")%>%
  select(genus_species, thermal_limit, type, metric) %>%
  mutate(Genus = str_split_fixed(.$genus_species, "\\_", 2)[,1],
         Species = str_split_fixed(.$genus_species, "\\_", 2)[,2],
         Metric_Tmax = ifelse(type == "max", metric, NA),
         Metric_Tmin = ifelse(type == "min", metric, NA)) %>%
  spread(type, thermal_limit) %>%
  select(-genus_species, -metric) %>%
  unique(.) %>%
  rename("Tmax" = max, "Tmin" = min) %>%
  group_by(Genus, Species) %>%
  fill(c(Tmin, Tmax, Metric_Tmax, Metric_Tmin), .direction = "updown")

full <- left_join(full, thermal)

## reorder columns
full <- full %>%
  select(Realm, Class, Order, Family, Genus, Species, Tmin, Tmax, Metric_Tmin, Metric_Tmax, RealizedRangeSource, AbsLatitude, 
         RangeSize, LogBodySize, DispersalDistance, HotSeasonDormancy, ColdSeasonDormancy, 
         NicheFilling_Warm, NicheFilling_Cool, LogRangeFilling, EquatorwardBiasInUnderfilling,
         NicheFilling_Warm_Behaviour,
         NicheFilling_Warm_Acclimation, NicheFilling_Cool_Acclimation, LogRangeFilling_Acclimation,
         EquatorwardBiasInUnderfilling_Acclimation) %>%
  arrange(Realm, Class, Order, Family, Genus, Species, RealizedRangeSource)
 
# na <- full %>%
#   filter(is.na(Metric_Tmax) & is.na(Metric_Tmin))

## fill NAs
full = full %>%
  group_by(Genus, Species, Class, Order, Family) %>%
  fill(Metric_Tmax, Metric_Tmin,
       DispersalDistance, HotSeasonDormancy, ColdSeasonDormancy,
       .direction = "updown")


write.csv(full, "data-processed/minimum-dataset.csv", row.names = FALSE)
