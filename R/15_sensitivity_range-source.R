## range source sensitivity
library(tidyverse)
select <- dplyr::select

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####      Prepping data       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
uofill <- read.csv("data-processed/niche-filling/thermal-niche-filling-metrics.csv") 
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv") %>%
  rename("realm" = Realm) %>%
  select(-limit_type)


thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv")%>%
  select(genus_species, type, metric) %>%
  unique(.) %>%
  rename("limit_type" = type) %>%
  mutate(limit_type = ifelse(limit_type == 'max, min', 
                             'tmax, tmin',
                             ifelse(limit_type ==  'max', 
                                    'tmax',
                                    ifelse(limit_type ==  'min', 
                                           'tmin',NA))))

## make one column for both cold and warm filling, with edge_type specifying whether cold or warm
## make species and source columns 
data <- uofill %>%
  mutate(genus_species = paste(str_split_fixed(.$range, '_', 3)[,1], 
                               str_split_fixed(.$range, '_', 3)[,2], sep = '_')) %>%
  mutate(source = str_split_fixed(.$range, '_', 3)[,3]) %>%
  select(range, genus_species, source, everything()) %>%
  mutate(source = factor(.$source, levels = c("IUCN", "GBIF", "GARD"), ordered = TRUE)) %>%
  arrange(genus_species, type, source) %>%
  rename("sensitivity_type" = type) %>%
  filter(sensitivity_type == "Te")

data <- data %>%
  mutate(perfect_warm = ifelse(warm_over == 0 & warm_under == 0, 0, NA)) %>%
  mutate(perfect_cold = ifelse(cold_over == 0 & cold_under == 0, 0, NA)) %>%
  gather(key = "filling_type", value = "filling_value", c(warm_under, warm_over, cold_under, 
                                                          cold_over, perfect_warm, perfect_cold)) %>%
  filter(filling_value != 0 | filling_value == 0 & (filling_type %in% 
                                                      c("perfect_cold","perfect_warm"))) %>%
  select(range, sensitivity_type, limit_type,
         genus_species, source, realm, dormancy, filling_type, filling_value, 
         lat_mp, range_area_km2) %>%
  mutate(edge_type = ifelse(str_detect(filling_type, "warm"), "warm", "cold")) %>%
  mutate(type = ifelse(str_detect(filling_type, "warm"), "max", "min")) %>%
  arrange(genus_species, filling_type, sensitivity_type, source) %>%
  filter(!is.infinite(filling_value)) ## get rid of infinite filling values:



## plot correlation 
## cold
cold <- data %>%
  filter(edge_type == "cold") %>%
  spread(key = source, value = "filling_value") %>%
  mutate(source = ifelse(is.na(IUCN), "GARD", "IUCN")) %>%
  group_by(genus_species) %>%
  mutate(expert_informed = ifelse(!is.na(IUCN), max(IUCN, na.rm=T), 
                                  max(GARD, na.rm=T)),
         GBIF = max(GBIF, na.rm=T)) %>%
  filter(!is.infinite(expert_informed), !is.infinite(GBIF))
  
  
cold_plot <- cold %>%
  ggplot(., aes(x = expert_informed, y = GBIF, shape = source)) + theme_bw() + 
  geom_point(colour = "steelblue") +
  geom_abline(intercept = 0, slope = 1) + 
  theme(panel.grid = element_blank()) + 
  labs(x = "Cool niche filling (expert-informed range)",
       y = "Cool niche filling (occurrence-inferred range)") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "none")

warm <- data %>%
  filter(edge_type == "warm") %>%
  spread(key = source, value = "filling_value") %>%
  mutate(source = ifelse(is.na(IUCN), "GARD", "IUCN")) %>%
  group_by(genus_species) %>%
  mutate(expert_informed = ifelse(!is.na(IUCN), max(IUCN, na.rm=T), 
                                  max(GARD, na.rm=T)),
         GBIF = max(GBIF, na.rm=T)) %>%
  filter(!is.infinite(expert_informed), !is.infinite(GBIF)) 


warm_plot <- warm %>%
  ggplot(., aes(x = expert_informed, y = GBIF, shape = source)) + theme_bw() + 
  geom_point(colour = '#b45346') +
  geom_abline(intercept = 0, slope = 1) + 
  theme(panel.grid = element_blank()) + 
  labs(x = "Warm niche filling (expert-informed range)",
       y = "Warm niche filling (occurrence-inferred range)") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "none")
  

## now do the same for range filling values 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####      Prepping data       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in under/over/filling data:
rf <- read.csv("data-processed/range-filling/rangefilling.csv")
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv")
## get rid of species in traits that are not in uofill
traits <- filter(traits, genus_species %in% rf$species)

## merge with traits and thermal limit info
traits = traits %>%
  rename("realm" = Realm) %>%
  select(-limit_type) 

thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv") %>%
  select(genus_species, type, metric) %>%
  unique(.) %>%
  rename("limit_type" = type) %>%
  mutate(limit_type = ifelse(limit_type == 'max, min', 
                             'tmax, tmin',
                             ifelse(limit_type ==  'max', 
                                    'tmax',
                                    ifelse(limit_type ==  'min', 
                                           'tmin',NA)))) %>%
  group_by(genus_species) %>%
  mutate(limit_type = paste(limit_type, collapse = ", "))%>%
  mutate(n_distinct = n_distinct(metric)) %>%
  ungroup(.) %>%
  mutate(metric = ifelse(n_distinct == 1, as.character(metric), 'mixed')) %>%
  unique(.) %>%
  select(-n_distinct) 

rf$filling <- log((rf$pr_cells - rf$u_cells) / rf$pr_cells)

rf <- rf %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  left_join(., thermal_limits, by = c('species' = 'genus_species')) %>%
  left_join(., traits, by = c('species' = 'genus_species', 'realm')) 


range <- rf %>%
  spread(key = source, value = filling) %>%
  mutate(source = ifelse(is.na(IUCN), "GARD", "IUCN")) %>%
  group_by(species) %>%
  mutate(expert_informed = ifelse(!is.na(IUCN), max(IUCN, na.rm=T), 
                                  max(GARD, na.rm=T)),
         GBIF = max(GBIF, na.rm=T)) %>%
  filter(!is.infinite(expert_informed), !is.infinite(GBIF)) 

range_plot <- range %>%
  ggplot(., aes(x = expert_informed, y = GBIF, shape = source)) + theme_bw() + 
  geom_point(colour = 'darkgrey') +
  geom_abline(intercept = 0, slope = 1) + 
  theme(panel.grid = element_blank()) + 
  labs(x = "Log(range filling) (expert-informed range)",
       y =  "Log(range filling) (occurrence-inferred range)", 
       shape = "Expert range source") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none")


## arrange them 
rsource <- ggdraw() +
  draw_plot(cold_plot, x = 0, y = 0, height = 1, width = 0.33) +
  draw_plot(warm_plot, x = 0.33, y = 0, height = 1, width = 0.33) +
  draw_plot(range_plot, x = 0.66, y = 0, height = 1, width = 0.33) +
  draw_plot_label(label = c("a)", "b)", "c)"),
                  x = c(0, 0.33, 0.66),
                  y = c(1, 1, 1), size = 10, 
                  color = "grey30")

ggsave(rsource, path = "figures", 
       filename = "range-source-plot.png", 
       device = "png", width = 12, height = 4)


## legend
range %>%
  ggplot(., aes(x = expert_informed, y = GBIF, shape = source)) + theme_bw() + 
  geom_point(colour = 'black') +
  geom_abline(intercept = 0, slope = 1) + 
  theme(panel.grid = element_blank()) + 
  labs(x = "Log(range filling) (expert-informed range)",
       y =  "Log(range filling) (occurrence-inferred range)", 
       shape = "Expert range source") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") 

