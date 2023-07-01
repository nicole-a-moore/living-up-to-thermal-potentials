## this script brings together IUCN, GARD, and GBIF realized range maps for ectothermic animal species in Globtherm
library(tidyverse)
library(rgdal)
library(sf)
library(rnaturalearth)
library(latticeExtra)
library(taxize)
library(plyr)

thermal_limits <- read_csv("data-raw/globtherm_full_dataset_2019.csv") %>%
  filter(thermy == "ectotherm")

## do a bit of name cleaning
thermal_limits$Genus[which(str_detect(thermal_limits$Genus, "<ca>") == TRUE)] <- "Trachemys"
thermal_limits$Species[which(str_detect(thermal_limits$Species, "<ca>") == TRUE)] <- c("quadracus", "atra","spiloptera", "spiloptera")
thermal_limits$Species[which(str_detect(thermal_limits$Species, "parvulus_parvulus") == TRUE)] <- "parvulus"
thermal_limits$genus_species <- paste(thermal_limits$Genus, thermal_limits$Species, sep = "_") 

## use taxize to get a list of species synonyms and the most taxonomically correct names:
thermal_species <- unique(paste(thermal_limits$Genus, thermal_limits$Species, sep = " "))
taxa_tt <- data.frame(binomial = thermal_species)
tsn_search_tt <- get_tsn(as.character(taxa_tt$binomial), accepted = FALSE)
#saveRDS(tsn_search_tt, "data-processed/intermediate-files/tsn_search_tt.rds")
tsn_search_tt <- readRDS("data-processed/intermediate-files/tsn_search_tt.rds")

tsns_tt <- data.frame(tsn_search_tt)
tsns_tt$binomial <- taxa_tt$binomial

found <- tsns_tt %>%
  subset(match == "found")  ## get only found found spp

## get synonyms
syns <- synonyms(tsn_search_tt)
#saveRDS(syns, "data-processed/intermediate-files/syns_tt.rds")
syns <- readRDS("data-processed/intermediate-files/syns_tt.rds")

syns_df <- ldply(syns, data.frame) %>%
  select(-c(1:2))
syns_df <- left_join(syns_df, found, by = c("sub_tsn" = "ids")) %>%
  filter(!is.na(sub_tsn)) %>%
  select(binomial, everything())

## add back species that were not found: 
not_found <- tsns_tt %>%
  subset(match == "not found") 

no_syn <- tsns_tt %>%
  subset(match == "found") %>%
  filter(!ids %in% syns_df$sub_tsn)
  
syns_df <- bind_rows(syns_df, no_syn) %>% 
  select(-class, -ids) %>%
  bind_rows(., not_found) %>%
  mutate(acc_name = ifelse(is.na(acc_name), as.character(binomial), as.character(acc_name)))

names <- syns_df %>%
  select(binomial, acc_name) %>%
  filter(!duplicated(.))

## save a key to reference correct names later:
write.csv(names, "data-processed/intermediate-files/globtherm_taxize-key.csv", row.names = FALSE)

## add column for corrected name:
thermal_limits$temp <- paste(thermal_limits$Genus, thermal_limits$Species, sep = " ")
thermal_limits <- left_join(thermal_limits, names, 
                            by = c("temp" = "binomial"))

colnames(names)[1] = "thermal_species"

thermal_limits <- thermal_limits %>%
  select(-temp) 

thermal_limits$metric[which(thermal_limits$metric %in% c('DLT50', 'LT50'))] = 'lt'

write.csv(thermal_limits, "data-processed/traits/globtherm_taxized.csv", row.names = FALSE)

length(unique(thermal_limits$genus_species)) #988

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######      IUCN RANGES        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## figure out which fish polygon groups to download from IUCN
families <- c("Pomacanthidae", "Blenniidae"	, "Albulidae", "Elopidae", "Megalopidae", "Chaetodontidae", "Epinephelidae", "Tetraodontidae", "Sparidae", "Centracanthidae", "Acanthuridae", "Syngnathidae", "Aulostomidae", "Centriscidae", "Fistulariidae", "Solenostomidae", "Istiophoridae", "Scombridae" , "Xiphiidae", "Labridae", "Scaridae")
classes <- c("Chondrichthyes", "Myxini")
orders <- c("Clupeiformes")

families_overlap <- families[which(families %in% thermal_limits$Family)]
classes_overlap <- classes[which(classes %in% thermal_limits$Class)]
orders_overlap <- orders[which(orders %in% thermal_limits$Order)]
## download sets that include: "Blenniidae", "Tetraodontidae", "Sparidae", "Labridae","Chondrichthyes", "Clupeiformes"

## sort through downloaded range maps, filter out species we do not have in thermal tolerance data
files <- c("large-files/IUCN/AMPHIBIANS/AMPHIBIANS.shp", 
           "large-files/IUCN/BLENNIES/BLENNIES.shp",
           "large-files/IUCN/CLUPEIFORMES/CLUPEIFORMES.shp",
           "large-files/IUCN/PUFFERFISH/PUFFERFISH.shp",
           "large-files/IUCN/REPTILES/REPTILES.shp",
           "large-files/IUCN/SEABREAMS_PORGIES_PICARELS/SEABREAMS_PORGIES_PICARELS.shp",
           "large-files/IUCN/SHARKS_RAYS_CHIMAERAS/SHARKS_RAYS_CHIMAERAS.shp",
           "large-files/IUCN/WRASSES_PARROTFISHES/WRASSES_PARROTFISHES.shp")
i = 1
while (i < length(files) + 1) {
  all_spp <- st_read(files[i])
  
  spp <- thermal_species[which(thermal_species %in% all_spp$binomial)]
  overlap <- all_spp %>%
    filter(binomial %in% spp)
  
  # look for iucn range names listed under synonyms
  spp_syns <- syns_df[which(syns_df$syn_name %in% all_spp$binomial),] %>%
    mutate(species = binomial) %>%
    select(syn_name, species) %>%
    filter(!duplicated(.))
  
  # change iucn range name from the synonym to its name in the thermal tolerance data 
  overlap_syns <- all_spp %>%
    filter(binomial %in% spp_syns$syn_name) %>%
    left_join(., spp_syns, by = c("binomial" = "syn_name")) %>%
    mutate(binomial = species) %>%
    select(-species)
  
  ##combine:
  if (i == 1) {
    IUCN <- rbind(overlap, overlap_syns)
  }
  else {
    IUCN <- rbind(IUCN, overlap, overlap_syns)
  } 
  
  i = i + 1
}

## collect same ID number (species) into one MULTIPOLYGON:
IUCN <- aggregate(IUCN, list(IUCN$id_no), function(x) x[1]) ## 322 ranges

## write out to file:
st_write(IUCN, "data-raw/polygons/IUCN-ectotherms.shp", driver = "ESRI Shapefile", 
         append = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######      GARD RANGES        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## checking the Meiri data for more reptile distributions 
## data from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.83s7k

## read in ranges:
GARD <- st_read("large-files/realized-ranges/GARD1.1_dissolved_ranges/modeled_reptiles.shp")

## remove species with IUCN ranges since GARD took ranges from IUCN when available
IUCN <- st_read("large-files/realized-ranges/IUCN-ectotherms.shp")
GARD <- GARD[which(!GARD$Binomial %in% IUCN$binomial),]

## look for GARD species in synonym database: 
syns_gard <- syns_df[which(syns_df$syn_name %in% GARD$Binomial),] %>%
  select(syn_name, binomial)
perfmatch_gard <- syns_df[which(syns_df$binomial %in% GARD$Binomial),] %>%
  select(binomial) %>%
  filter(!duplicated(.))

GARD_syns <- left_join(syns_gard, GARD, by = c("syn_name" = "Binomial")) %>%
  mutate(binomial = ifelse(!is.na(binomial), as.character(binomial), as.character(Binomial))) %>%
  select(-syn_name)

GARD_perfmatch <- left_join(perfmatch_gard, GARD, by = c("binomial" = "Binomial")) 

GARD_ol <- rbind(GARD_syns, GARD_perfmatch) %>% 
  select(binomial, geometry) ## 57 ranges!

## write: 
st_write(GARD_ol, "data-raw/polygons/GARD_ol.shp", driver = "ESRI Shapefile", 
         append = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######      GBIF RANGES        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
GBIF_og <- st_read("large-files/realized-ranges/GBIF polygons/Filtered occurences ectotherm animals_020817.shp")

# get species synonyms from files:
path = "large-files/filtered-occurrences/dropbox filtered occurrences for sWEEP/Filtered occurences ectotherm animals 16 12 31 copy/"
files <- list.files(path, pattern="*.csv")
filenames <- paste(path, files, sep="/")

syn_match <- data.frame(filenames = str_split_fixed(files, "\\.", n=2)[,1])

i = 1 
og_names <- c()
while (i < length(filenames)+1) {
  cur <- read.csv(filenames[i]) 
  
  og_names <- append(og_names, as.character(cur$species[1]))
  
  i = i + 1
}

syn_match$og_names <- og_names

## fix species name mistakes in range file:
GBIF <- GBIF_og 
GBIF$species <- as.character(GBIF$species)
GBIF$species[which(GBIF$species == "Clubiona triviali")] <- "Clubiona trivialis"

## replace species names in range file with the names originally used for their filenames
rearranged <- data.frame(species = GBIF$species)
rearranged <- left_join(rearranged, syn_match, by = c("species" = "og_names")) %>%
  mutate(filenames = ifelse(is.na(filenames), as.character(species), as.character(filenames)))
GBIF$species <- rearranged$filenames

## look for GBIF species in synonym database: 
syns_gbif <- syns_df[which(syns_df$syn_name %in% GBIF$species),] %>%
  select(syn_name, binomial)
perfmatch_gbif <- syns_df[which(syns_df$binomial %in% GBIF$species),] %>%
  select(binomial) %>%
  filter(!duplicated(.))

GBIF_syns <- left_join(syns_gbif, GBIF, by = c("syn_name" = "species")) %>%
  mutate(binomial = ifelse(!is.na(binomial), as.character(binomial), as.character(species))) %>%
  select(-syn_name)
GBIF_perfmatch <- left_join(perfmatch_gbif, GBIF, by = c("binomial" = "species")) 

which(GBIF_syns$binomial %in% GBIF_perfmatch$binomial) ## oh -- all syns were also in gbif as correct name, so do not add twice!

## merge to add column with corrected name:
GBIF_ol <- GBIF_perfmatch %>% 
  select(binomial, geometry) ## 234 ranges!


## Note: some GBIF species go missing here because they do not have thermal limits
## should not happen - look into why they are not in globtherm data
## Norops should be Anolis 
missing <- GBIF$species[which(!GBIF$species %in% GBIF_ol$binomial)]
GBIF_missing <- GBIF[which(!GBIF$species %in% GBIF_ol$binomial), ]

GBIF_missing$single_i <- str_replace_all(missing, "ii", "i")

i_gbif <- syns_df[which(syns_df$binomial %in% GBIF_missing$single_i),] %>%
  select(binomial) %>%
  filter(!duplicated(.))

GBIF_i <- GBIF_missing[which(GBIF_missing$single_i %in% i_gbif$binomial),] %>%
  mutate(binomial = single_i) %>%
  select(binomial, geometry) %>%
  as.data.frame(.)

## merge
GBIF_ol <- rbind(GBIF_ol, GBIF_i)  ## 240 ranges!
GBIF_missing <- GBIF_missing[which(!GBIF_missing$single_i %in% GBIF_i$binomial),] 
GBIF_missing$species ## 17 still missing

st_write(GBIF_ol, "large-files/realized-ranges/gbif_error-checked.shp", append = FALSE)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######      Combine GBIF, GARD and IUCN ranges        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
IUCN <- st_read("data-raw/polygons/IUCN-ectotherms.shp") %>%
  select(binomial, geometry)
GBIF <- st_read("large-files/gbif_error-checked.shp") 
GARD <- st_read("data-raw/polygons/GARD_ol.shp") 

## add source column to identify which ranges are from GBIF vs which are from IUCN:
GBIF$source <- "GBIF"
IUCN$source <- "IUCN"
GARD$source <- "GARD"

## combine:
realized_ranges <- rbind(IUCN, GBIF, GARD) %>%
  dplyr::rename("species" = binomial) ## okay, we have 623 ectotherm ranges


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### Take inventory of unique and overlapping species ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
length(unique(realized_ranges$species)) ## 491 unique species 
length(which(duplicated(realized_ranges$species))) ## 128 are duplicated

length(unique(IUCN$binomial)) ## 319 spp
length(unique(GBIF$binomial)) ## 240 spp
length(unique(GARD$binomial)) ## 56 spp 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######  Add extra information to realized ranges      #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## add extra information to realized ranges about realm, which hemisphere species is in, whether realized range crosses the equator, latitudinal midpoint, geographic area 
realms <- thermal_limits %>%
  select(genus_species, realm) %>%
  mutate(species = str_replace_all(genus_species, "_", " ")) %>%
  select(-genus_species) %>%
  filter(!duplicated(.))

realized_ranges <-  merge(realized_ranges, realms, by = "species")

## calculate latitudinal midpoint and geographic area:
realized_ranges <- realized_ranges %>%
  mutate(range_area_km2 = round(units::set_units(st_area(.), km^2), digits = 0)) 

## get bbox extents to approximate lat midpoint:
bbox <- as.data.frame(do.call("rbind", lapply(st_geometry(realized_ranges), st_bbox)))
lat_mp <- (bbox$ymin + bbox$ymax)/2
realized_ranges$lat_mp <- lat_mp

## investigate ranges that are duplicated:
rr <- realized_ranges
rr$range_id <- paste(rr$species, rr$source, sep = "_")
dups <- rr[rr$range_id %in% rr$range_id[which(duplicated(rr$range_id))], ] # 4 spp have two ranges from the same source

## Anolis lemurinus - keep larger (sec)
## Crotalus atrox - keep larger (first)
## Liolaemus fitzingerii - keep larger (sec)
## Liolaemus kriegi - keep larger (sec) 

dups <- dups[c(2,3,5,8),]
dups <- select(dups, -range_id)

realized_ranges <- realized_ranges %>%
  filter(!species %in% dups$species) %>%
  rbind(., dups)

## get rid of Freshwater spp:
realized_ranges <- filter(realized_ranges, realm != "Freshwater")

st_write(realized_ranges, "large-files/realized-ranges/realized-ranges.shp", append = FALSE,
         overwrite = TRUE)
##NOTE: some range areas too large to write 

## write new thermal limits database with only species that we have ranges for 
thermal_limits_new <- thermal_limits[paste(thermal_limits$Genus, thermal_limits$Species, sep = " ") %in% realized_ranges$species,]
length(unique(thermal_limits_new$genus_species)) ## 474

write.csv(thermal_limits_new, "data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv", row.names = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######     Investigating overlap        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## look at species with ranges in IUCN and made by Greta
## all 85 are squamata 
countries <- map_data("world")

thermal_limits_new <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv")
duplicated <- IUCN$binomial[which(IUCN$binomial %in% GBIF$binomial)]

duplicate <- filter(realized_ranges, species == as.character(duplicated[1]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[1]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
legend(st_bbox(duplicate)$xmin-5, st_bbox(duplicate)$ymin+2, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)

duplicate <- filter(realized_ranges, species == as.character(duplicated[2]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[2]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
legend(st_bbox(duplicate)$xmin-5, st_bbox(duplicate)$ymin+10, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)

duplicate <- filter(realized_ranges, species == as.character(duplicated[3]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[3]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
legend(st_bbox(duplicate)$xmin+8, st_bbox(duplicate)$ymin+9, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)

duplicate <- filter(realized_ranges, species == as.character(duplicated[4]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[4]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
legend(st_bbox(duplicate)$xmin+28, st_bbox(duplicate)$ymin+22, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)