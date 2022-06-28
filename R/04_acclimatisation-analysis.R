## this script wrangles the data needed to estimate 'acclimatised' thermal limits for each species 
library(tidyverse)
library(broom)
rename = dplyr::rename
summarise = dplyr::summarise
count = dplyr::count

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######     Calculate global ARR as weighted mean          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######     Data with pre-calculated ARRs           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## collate all ARR data together, making sure none are duplicated:
##  GUNDERSON AND STILLMAN 
gunderson <- read.csv("data-raw/arr-data/g&s_data.csv")[1:287,] %>%
  mutate(ref_number = str_split_fixed(str_split_fixed(ref, "\\[", n=2)[,2], "\\]", n=2)[,1])
gunderson_refs <- data.frame(ref = read.csv("data-raw/arr-data/g&s_refs.csv")[,1])%>%
  mutate(ref_number = str_split_fixed(ref, "\\.", n=2)[,1],
         ref_full = str_split_fixed(ref, "\\.", n=2)[,2]) %>%
  select(-ref)
  
gunderson <- left_join(gunderson, gunderson_refs, by = c("ref_number"))

gunderson$ref_number[100] <- "18, 19"
gunderson$ref_full[100] <- paste(gunderson$ref_full[18], gunderson$ref_full[18], sep = ", ")

gunderson <- gunderson %>%
  mutate(compilation = "Gunderson", genus_species = paste(genus, species, sep = "_"), 
         ref = ref_full, lifestage = NA) %>%
  dplyr::rename("tmax" = CTmax.ARR, "tmin" = CTmin.ARR, "Genus" = genus, "Species" = species) %>%
  gather(key = "parameter_tmax_or_tmin", value = "ARR", c(tmax, tmin)) %>%
  filter(!is.na(ARR)) %>%
  select(genus_species, compilation, ARR, Genus, Species, ref, lifestage,
         parameter_tmax_or_tmin)

##  MORLEY ET AL. 
morley <- read_csv("data-raw/arr-data/morley.csv")[8:380,]
colnames(morley) <- morley[1,]
morley <- morley[-1,1:46]

morley <- morley %>%
  rename(ARR = `CTMax ARR`) %>%
  mutate(genus_species = paste(Genus, Species, sep = "_"), compilation = 'Morley', Order = NA, 
         lifestage = NA, parameter_tmax_or_tmin = "tmax", "ref" = Reference,  ARR = as.numeric(ARR)) %>%
  select(genus_species, compilation, ARR, Genus, Species, ref, lifestage,
         parameter_tmax_or_tmin) %>%
  filter(!is.na(ARR))

##  ROHR ET AL.  
rohr <- read_csv("data-raw/arr-data/rohr.csv") 
colnames(rohr) <- rohr[2,]
rohr <- rohr[-c(1:2, 254:259),1:20]

rohr <- read_csv("data-raw/arr-data/amphib-rohr.csv") %>%
  select(Genus3, Species3, Stage, Reference) %>%
  mutate(genus_species = paste(Genus3, Species3, sep = "_")) %>%
  group_by(genus_species) %>%
  mutate(Reference = paste(Reference, collapse = " "), Stage = paste(Stage, collapse = " ")) %>%
  ungroup(.) %>%
  filter(!duplicated(.)) %>%
  left_join(rohr, .) 

rohr <- rohr %>%
  rename("ARR" = `Average of ARR`, "Genus" = Genus3, "Species" = Species3, "ref" = Reference,
         "lifestage" = Stage) %>%
  mutate(genus_species = paste(Genus, Species, sep = "_"), compilation = "Rohr",
         parameter_tmax_or_tmin = "tmax", ARR = as.numeric(ARR)) %>%
  select(genus_species, compilation, ARR, Genus, Species,lifestage, 
         ref, parameter_tmax_or_tmin) %>%
  filter(!is.na(ARR))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######      Data with raw thermal limits           #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##   INTRATHERM
intratherm <- read.csv("data-raw/arr-data/intratherm-with-elev.csv") %>%
  filter(original_compilation != "Comte") %>% ## get rid of Comte duplicated data
  mutate(genus_species = paste(genus, species, sep = '_'), compilation = "Intratherm") %>%
  select(genus_species, genus, species, acclim_temp, ref, parameter_value, 
         parameter_tmax_or_tmin, metric_type, life_stage.x, compilation) %>%
  dplyr::rename("lifestage" = life_stage.x, "Genus" = genus, "Species" = species)

##   COMTE AND OLDEN 
comte <- read.csv("data-raw/arr-data/comte-all.csv") %>%
  mutate(genus_species = str_replace_all(Species, ' ', "_"),  compilation = "Comte") %>%
  select(genus_species, Species, Temperature.of.acclimation...C., Thermal.limit...C., Source,
         Methodology, Life.stage, compilation) %>%
  mutate(parameter_tmax_or_tmin = "tmax", "Genus" = str_split_fixed(genus_species, "\\_", n=2)[,1],
         "Species" = str_split_fixed(genus_species, "\\_", n=2)[,2]) %>%
  mutate(metric_type = ifelse(Methodology == "dynamic", "critical", "lethal")) %>%  # dynamic = CTmax, static = LT
  dplyr::rename("acclim_temp" = Temperature.of.acclimation...C., "parameter_value" = Thermal.limit...C., 
                "ref" = Source, "lifestage" = Life.stage) %>%
  select(-Methodology) 

## combine all datasets:
arr <- rbind(intratherm, comte) 

## filter to species with known acclimation temperatures, parameter values and at least two estimates
arr <- filter(arr, !is.na(acclim_temp), 
              !is.na(parameter_value), 
              genus_species %in% genus_species[which(duplicated(genus_species))])

## split into upper and lower thermal limits:
arr_upper <- filter(arr, parameter_tmax_or_tmin == "tmax") %>% 
  group_by(genus_species) %>%
  mutate(ref = paste(ref, collapse = ", "))
arr_lower <- filter(arr, parameter_tmax_or_tmin == "tmin")  %>% 
  group_by(genus_species) %>%
  mutate(ref = paste(ref, collapse = ", "))

l <- str_split(arr_upper$ref, "\\, ")
for (i in 1:length(l)) {
  l[[i]] <- paste(unique(l[[i]]), collapse = ", ")
}
arr_upper$ref <- unlist(l)
l <- str_split(arr_lower$ref, "\\, ")
for (i in 1:length(l)) {
  l[[i]] <- paste(unique(l[[i]]), collapse = ", ")
}
arr_lower$ref <- unlist(l)


## fit lm to thermal limits of each spp, extract ARR (slope of lm)
arr_upper_slopes <- arr_upper %>%
  do(tidy(lm(parameter_value ~ acclim_temp, data = .), conf.int = TRUE)) %>% #fit lm to each group
  filter(term == "acclim_temp") %>% #extract just the slopes - get rid of intercept 
  rename(ARR = estimate) %>% #rename these ARR
  select(genus_species, ARR) %>%
  left_join(select(arr_upper, -c(acclim_temp, parameter_value, metric_type)),
            ., by = c("genus_species")) %>%
  filter(!duplicated(genus_species, ARR), !is.na(ARR)) 

arr_lower_slopes <- arr_lower %>%
  do(tidy(lm(parameter_value ~ acclim_temp, data = .), conf.int = TRUE)) %>% #fit lm to each group
  filter(term == "acclim_temp") %>% #extract just the slopes - get rid of intercept 
  rename(ARR = estimate) %>% #rename these ARR
  select(genus_species, ARR)  %>%
  left_join(select(arr_lower, -c(acclim_temp, parameter_value, metric_type)),
            ., by = c("genus_species")) %>%
  filter(!duplicated(genus_species, ARR), !is.na(ARR)) 

arr_all <- bind_rows(rohr, morley, gunderson, arr_upper_slopes, arr_lower_slopes) %>%
  mutate(id = 1:nrow(.)) 



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######          Remove duplicates                #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## find which are duplicated: 
dup <- arr_all %>%
  mutate(genus_species = paste(genus_species, parameter_tmax_or_tmin, sep = " ")) %>%
  filter(genus_species %in% .$genus_species[duplicated(.$genus_species)]) %>%
  group_by(genus_species) %>%
  filter(length(unique(compilation)) >= 2)%>%
  ungroup() 

## search for first word of Morley/Rohr/Gunderson reference in reference list for duplicates in Intratherm/Comte
## if duplicate reference is found, put a 1  
m_r_g <- dup %>%
  filter(compilation %in% c("Rohr", "Gunderson", "Morley")) %>%
  mutate(first_word = ifelse(compilation %in% c("Morley", "Gunderson"),
                str_split_fixed(ref, " ", n= 2)[,1],
                ifelse(compilation == "Rohr", 
                       str_split_fixed(ref, "\\,", n= 2)[,1], 
                       NA))) %>%
  select(genus_species, first_word, id, compilation, ref) 

## fix hard ones
m_r_g$first_word[which(m_r_g$id %in% c(577, 578))] <- 'Davis'
m_r_g$first_word[which(m_r_g$id %in% c(254))] <- "Jumbam"
m_r_g$first_word[which(m_r_g$id %in% c(229, 230))] <- "Mitchell"
m_r_g$first_word[which(m_r_g$id %in% c(218))] <-"Terblanche"
m_r_g$first_word[which(m_r_g$id %in% c(371, 372))] <- "Young"
m_r_g$first_word[which(m_r_g$id %in% c(193))] <- "Rajaguru"
m_r_g$first_word[which(m_r_g$id %in% c(206,219,220))] <- "Simons"
m_r_g$first_word[which(m_r_g$id %in% c(140))] <- "Davreau"
m_r_g$first_word[which(m_r_g$id %in% c(97))] <- "Howard"
m_r_g$first_word[which(m_r_g$id %in% c(124))] <- "Chen"

i_c <- dup %>% 
  filter(compilation %in% c("Intratherm", "Comte")) %>%
  select(genus_species, id, compilation, ref) %>%
  rename("ref_ic" = ref) 
  
## check for duplicates, if duplicated keep intratherm/comte data 
merge <- left_join(i_c, m_r_g, by = c("genus_species")) %>%
  mutate(i_c = ifelse(str_detect(ref_ic, first_word),1,0))

## search for duplicates within leftovers:
leftovers <- m_r_g %>%
  filter(!id %in% merge$id.y)

merge <- merge %>%
  filter(i_c == 0) 

leftovers <- leftovers %>%
  rbind(., filter(m_r_g, id %in% merge$id.y)) 

count <- leftovers %>%
# prioritize: rohr, then morley, then gunderson since that is order of publication 
  group_by(genus_species) %>%
  summarize(n = length(unique(compilation))) 

deathrow = leftovers %>% 
  filter(genus_species %in% count$genus_species[which(count$n > 1)])

leftovers = leftovers %>% 
  filter(genus_species %in% count$genus_species[which(count$n == 1)])

keepers <- deathrow %>%
  group_by(genus_species) %>%
  mutate(r = any(compilation == "Rohr"), g = any(compilation == "Gunderson"),
         m = any(compilation == "Morley")) %>%
  mutate(keep_from = ifelse(r, "Rohr", 
                            ifelse(g, "Gunderson",
                                   ifelse(m, "Morley", NA)))) %>%
  filter(keep_from == compilation) %>% 
  select(-r, -m, -g, -first_word, -keep_from) 
  
## combine all data from dup to keep:
keepsies <- i_c %>%
  rename("ref" = ref_ic) %>%
  bind_rows(., keepers, leftovers) 

keepsies <- filter(dup, id %in% keepsies$id) %>%
  mutate(parameter_tmax_or_tmin = str_split_fixed(genus_species, " ", n=2)[,2], 
         genus_species = str_split_fixed(genus_species, " ", n=2)[,1])
  
new_arr <- arr_all %>%
  filter(!id %in% dup$id) %>% ## get rid of all duplicates
  bind_rows(., keepsies) ## add back properly filtered duplicates 
  

## for species with multiple ARRs from one/more sources, calculate mean so each species contributes 1 estimate 
new_arr <- new_arr %>%
  group_by(genus_species, parameter_tmax_or_tmin) %>%
  summarise(ARR = mean(ARR), 
         id = paste(id, collapse = ", "),
         ref = paste(ref, collapse = ", "), 
         compilation = paste(compilation, collapse = ", "),
         genus_species = first(genus_species)) %>%
  mutate(genus_species = str_replace_all(genus_species, " ", "_"))

t <- str_split(new_arr$compilation, "\\, ")
for (i in 1:length(t)) {
  t[[i]] <- paste(unique(t[[i]]), collapse = ", ")
}
v <- unlist(t)
new_arr$compilation <- v


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######             Group species                 #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(taxize)

## get species' higher taxonomy
taxonomy <- tax_name(query = unique(str_replace_all(new_arr$genus_species, "\\_", " "))
                      [c(1:156, 158, 161:303, 305:374, 376:690)],
                     get = c("kingdom", "phylum", "class","order","family"), 
                     db = "both")

#saveRDS(taxonomy, "data-processed/intermediate-files/arr_taxized.rds")
taxonomy <- readRDS("data-processed/intermediate-files/arr_taxized.rds")

unique = unique(str_replace_all(new_arr$genus_species, "_", " "))
missing <- data.frame(genus_species = 
                        str_replace(unique[!unique %in% taxonomy$query], " ", "_"),
                      phylum = NA, family = NA, order = NA, class = NA, kingdom = NA)


taxonomy <- taxonomy %>%
  group_by(query) %>%
  mutate(kingdom = first(na.omit(kingdom)),
         class = first(na.omit(class)),
         order = first(na.omit(order)),
         family = first(na.omit(family)),
         phylum = first(na.omit(phylum)),
         genus_species = str_replace_all(query, " ", "_")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-db, -query) 

taxonomy <- rbind(taxonomy, missing)

## fill in missing:
missing <- taxonomy %>%
  filter(is.na(kingdom))

write.csv(missing, "data-processed/intermediate-files/arr_taxized_missing.csv", row.names = FALSE)

missing <- read.csv("data-raw/traits/arr_taxized_missing_filled.csv")

taxonomy <- taxonomy %>%
  filter(!genus_species %in% missing$genus_species) %>%
  rbind(., missing)

## fix missing phylums, make kindgom consistent, get rid of non-animals 
taxonomy <- taxonomy %>% 
  mutate(kingdom = ifelse(kingdom == "Metazoa", "Animalia", as.character(kingdom))) %>%
  filter(kingdom == "Animalia") 

## join and assess what we can do:
new_arr <- left_join(taxonomy, new_arr)

tally_arr <- new_arr %>%
  filter(parameter_tmax_or_tmin == "tmax") %>%
  count(class)

thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = "_"))

tally_tlims <- thermal_limits %>%
  count(Class)%>%
  rename("n_tlims" = n, "class" = Class)

upp_count = left_join(tally_tlims, tally_arr)

tally_arr <- new_arr %>%
  filter(parameter_tmax_or_tmin == "tmin") %>%
  count(class)

tally_tlims <- thermal_limits %>%
  count(Class) %>%
  rename("n_tlims" = n, "class" = Class)

low_count = left_join(tally_tlims, tally_arr)

## calculate aggregated ARR at the class level
mean_ARRs <- new_arr %>%
  filter(class %in% thermal_limits$Class) %>%
  filter(ARR > -0.15 & ARR < 2) %>%
  group_by(class, parameter_tmax_or_tmin) %>%
  summarise(mean_ARR = mean(ARR), sd_ARR = sd(ARR), n_ARR = n()) %>%
  filter(n_ARR > 1)

new_arr_class <- new_arr %>% 
  count(class) %>%
  left_join(., new_arr) %>%
  mutate(class_label = paste(class, " (n=", n, ")", sep = "")) 

library(PNWColors)
pal = pnw_palette(name="Sunset2", n=8, type="continuous")

class_arr <- new_arr_class %>%
  filter(n != 1) %>%
  mutate(label = ifelse(parameter_tmax_or_tmin == "tmax", "Upper thermal limit", 
                        "Lower thermal limit")) %>%
  filter(class %in% thermal_limits$Class) %>%
  filter(ARR > -0.15 & ARR < 2, class != "Arachnida") %>%
  ggplot(., aes(y = class_label, x = ARR, col = class)) + geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 1) +
  theme_light() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  labs(y = 'Class', x = "") +
  facet_wrap(~label) +
  scale_color_manual(values = pal)
    
realms <- thermal_limits %>%
  select(realm, Class) %>%
  mutate(realm = ifelse(Class == "Arachnida", 'Terrestrial', as.character(realm))) %>%
  mutate(realm = ifelse(realm %in% c("Intertidal", "Marine"), 'Aquatic', as.character(realm))) %>%
  filter(!duplicated(.)) 

realms = rbind(realms, data.frame(realm = c("Aquatic"), Class = c("Teleostei")))

new_arr_realms <- new_arr %>%
  left_join(., realms, by = c("class" = "Class")) 

mean_ARRs_realm <- new_arr_realms %>%
  filter(ARR > -0.15 & ARR < 2) %>%
  group_by(realm, parameter_tmax_or_tmin) %>%
  summarise(mean_ARR = mean(ARR), sd_ARR = sd(ARR), n_ARR = n()) %>%
  mutate(realm = ifelse(realm == 'Aquatic', "Marine", as.character(realm)))

mean_ARRs_realm[5:6,] <- mean_ARRs_realm[1:2,]
mean_ARRs_realm$realm[5:6] <- "Intertidal"

new_arr_realms <- new_arr_realms %>% 
  count(realm) %>%
  left_join(., new_arr_realms) %>%
  mutate(realm_label = paste(realm, " (n=", n, ")", sep = ""))

pal = pnw_palette(name="Sunset", n=2, type="continuous") 

realm_arr <- new_arr_realms %>%
  mutate(label = ifelse(parameter_tmax_or_tmin == "tmax", "Upper thermal limit", 
                        "Lower thermal limit")) %>%
  filter(!is.na(realm)) %>%
  filter(ARR > -0.15 & ARR < 2) %>%
  ggplot(., aes(y = realm_label, x = ARR, col = realm)) + geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.1), size = 1) +
  theme_light() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  labs(y = 'Realm', x = "Acclimation response ratio (ARR)") +
  facet_wrap(~label) +
  scale_color_manual(values = pal) 

## combine into figure:
library(cowplot)
library(gridExtra)
combined <- grid.arrange(class_arr, realm_arr, nrow = 2)

ggsave(combined, path = "figures/additional-figures/acclimatisation-analysis", filename = "ARR-by-class-and-realm.png", 
       device = "png", height = 4.8, width = 4.8)


## solve for species-specific intercepts using class- or realm-specific mean ARR:
mean_ARRs$class_specific = "TRUE"
acc_data <- thermal_limits %>%
  mutate(parameter_tmax_or_tmin = ifelse(type == "max", "tmax", "tmin")) %>%
  left_join(., mean_ARRs_realm, by = c("realm", "parameter_tmax_or_tmin")) %>%
  rename("realm_mean_ARR" = mean_ARR, "realm_sd_ARR" = sd_ARR, "realm_n_ARR" = n_ARR) %>%
  left_join(., mean_ARRs, by = c("Class" = "class", "parameter_tmax_or_tmin")) %>%
  rename("class_mean_ARR" = mean_ARR, "class_sd_ARR" = sd_ARR, "class_n_ARR" = n_ARR) %>%
  mutate("ARR_equ_slope" = ifelse(!is.na(class_mean_ARR), class_mean_ARR, realm_mean_ARR))%>%
  mutate(ARR_equ_int = thermal_limit - ARR_equ_slope*acclimation_temperature) %>%
  group_by(genus_species) %>%
  filter(any(!is.na(ARR_equ_int))) %>%
  ungroup(.) %>%
  filter(!is.na(acclimation_temperature)) 

## plot global arr and all of our data:
acc_data <- acc_data %>%
  group_by(Class, type) %>%
  mutate(mean_int_class = mean(ARR_equ_int)) %>%
  group_by(realm, type) %>%
  mutate(mean_int_realm = mean(ARR_equ_int)) %>%
  ungroup() %>%
  mutate(int = ifelse(!is.na(mean_int_class), mean_int_class, mean_int_realm)) 

upper <- filter(acc_data, type == "max")
lower <- filter(acc_data, type == "min")

upper_plot <- acc_data %>%
  filter(type ==  "max") %>%
  ggplot(., aes(x = acclimation_temperature, y = thermal_limit)) + 
  geom_abline(slope = upper$ARR_equ_slope, intercept = upper$ARR_equ_int,
              colour = "grey", size = 0.1) + 
  geom_point(col = "darkred") +
  labs(x = "Acclimation temperature (°C)", y = "Upper thermal limit (°C)") +
  theme_light() +
  theme(panel.grid = element_blank())

lower_plot <- acc_data %>%
  filter(type ==  "min") %>%
  ggplot(., aes(x = acclimation_temperature, y = thermal_limit)) + 
  geom_abline(slope = lower$ARR_equ_slope, intercept = lower$ARR_equ_int,
              col = "grey", size = 0.1) + 
  geom_point(col = "darkblue") +
  labs(x = "Acclimation temperature (°C)", y = "Lower thermal limit (°C)") +
  theme_light() +
  theme(panel.grid = element_blank())

slopey <- grid.arrange(upper_plot, lower_plot, nrow = 2)

ggsave(slopey, path = "figures/additional-figures/acclimatisation-analysis", filename = "ARRs-all-species.png", 
       device = "png", height = 4.8, width = 4.8)

acc_methods <- ggdraw() + 
  draw_plot(combined, x = 0, y = 0, width = 0.6, height = 1) +
  draw_plot(slopey, x = 0.6, y = 0, width = 0.4, height = 1) +  
  draw_plot_label(label = c("a)", "b)"),
                  x = c(0, 0.6),
                  y = c(1, 1), size = 10, 
                  color = "grey30") 

ggsave(acc_methods, path = "figures/additional-figures/acclimatisation-analysis", filename = "ARR-methods.png", 
       device = "png", height = 4.8, width = 8)

write.csv(acc_data, "data-processed/traits/acclimation-data.csv", row.names = FALSE)
