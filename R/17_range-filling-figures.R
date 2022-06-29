## range filling figures 
library(tidyverse)
library(raster)
library(cowplot)
library(sf)
library(PNWColors)
select = dplyr::select
rename = dplyr::rename

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       Prepping data for figure making     ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## niche:
uofill <- read.csv("data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics.csv") 
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv")
thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv")%>%
  select(genus_species, type, metric) %>%
  unique(.)

## split range id into species and source:
uofill <- uofill %>%
  mutate(species = paste(str_split_fixed(.$range, '_', 3)[,1], 
                         str_split_fixed(.$range, '_', 3)[,2], sep = ' ')) %>%
  mutate(source = str_split_fixed(.$range, '_', 3)[,3]) %>%
  select(range, species, source, everything()) %>%
  mutate(source = factor(.$source, levels = c("IUCN", "GARD", "GBIF"), ordered = TRUE)) %>%
  arrange(species, type, source)

## select data using IUCN realized range if species has realized range from multiple sources 
uofill_iucn <- uofill %>%
  mutate(temp = paste(species, type, sep = '')) %>%
  filter(!duplicated(temp)) %>%
  select(-temp)

## split by sensitivity analysis groups
types <- group_split(uofill_iucn, type)

Te_acclimatized <- types[[1]]
Te <- types[[2]]
Te_tpref <- types[[3]]
Te_tpreftb <- types[[4]]
Te_subset_tpreftb <- Te[which(Te$species %in% Te_tpreftb$species),]
Te_subset_tpref <- filter(Te, Te$species %in% Te_tpref$species)
Te_subset_acclimatized <- filter(Te, Te$species %in% Te_acclimatized$species)

types = list(Te_acclimatized, Te, Te_tpref, Te_tpreftb,
             Te_subset_acclimatized, Te_subset_tpreftb, Te_subset_tpref)
names(types) = c("Te_acclimatized", "Te", "Te_tpref", "Te_tpreftb",
                 "Te_subset_acclimatized", "Te_subset_tpreftb",
                 "Te_subset_tpref")
## range: 
## read in results:
rf <- read.csv("data-processed/potential-ranges/range-filling/rangefilling-metrics_model-ready.csv")
nf <- read.csv("data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics_model-ready.csv")

## calculate the proportion of occupied cells in the potential range
rf$prop_occupied <- (rf$pr_cells - rf$u_cells) / rf$pr_cells
rf$log_prop_occupied <- log(rf$prop_occupied)
rf <- filter(rf, !is.infinite(log_prop_occupied))

## re-order factors to give desired contrasts
rf$realm <- relevel(factor(rf$realm), ref = "Terrestrial")
rf$Trophic_position <- relevel(factor(rf$Trophic_position), ref = "insectivore")
rf$Trophic_position <- relevel(factor(rf$Trophic_position), ref = "omnivore")
rf$Trophic_position <- relevel(factor(rf$Trophic_position), ref = "herbivore")

## split by type:
types <- group_split(rf, type)

acc <- types[[1]]
te <- types[[2]]
te_subset_acc <- filter(te, te$species %in% acc$species)
acc <- filter(acc, acc$species %in% te_subset_acc$species)

## give priority to IUCN ranges, then GARD, then GBIF
te_subset_acc <- te_subset_acc %>%
  select(range, species, source, everything()) %>%
  mutate(source = factor(.$source, levels = c("IUCN", "GARD", "GBIF"), ordered = TRUE)) %>%
  arrange(species, type, source) %>% 
  filter(!duplicated(species)) 

acc <- acc %>%
  select(range, species, source, everything()) %>%
  mutate(source = factor(.$source, levels = c("IUCN", "GARD", "GBIF"), ordered = TRUE)) %>%
  arrange(species, type, source) %>% 
  filter(!duplicated(species)) 


##~~~~~~~~~~~~~~~~~~~~~~~##
#####       Plots    ######
##~~~~~~~~~~~~~~~~~~~~~~~##
## plot collection point on map 
## colour by difference between 95% percentile of cold body temperature across range and species' coldest body temperature at collection point 
library(rgeos)
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

op_temps <- stack("large-files/operative-temperatures/Te_shade_min.grd")
mar_temps <- stack("data-processed/temperature-data/marine/raster_marine_low.grd")
int_temps <- stack("data-processed/temperature-data/intertidal/raster_intertidal_low.grd") 

## read in niche limits to get species list 
filling <- read.csv("data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics_model-ready.csv") %>%
  mutate(species = str_replace_all(genus_species, "_", ".")) %>%
  mutate(range = paste(species, source, sep = "_")) %>%
  filter(sensitivity_type == "Te") %>%
  select(range, genus_species, realm) %>%
  unique(.)

## read in realized ranges and subset to species in analysis 
ranges <- readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds")

ranges <- ranges[[which(names(ranges) %in% filling$range)]]

## read in collection points 
thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv")

## subset to ones on our model
locs <- select(thermal_limits, genus_species, collection_latitude, collection_longitude) %>%
  unique(.) %>% 
  left_join(filling, .) %>%
  filter(!is.na(collection_longitude), !is.na(collection_latitude)) %>%
  unique(.)


# loop through list of species and:
sp <- 1
data <- c()
while(sp <= nrow(locs)) {
  species <- locs[sp,]
    
  # 1. pick out realized range
  range <- ranges[[which(names(ranges) == species$range)]]
  
  # 2. calculate 95% percentile temp across range
  # 3. extract mean temp at collection location
  if(species$realm == "Intertidal") {
    min <- min(values(mask(terr_temps, range)), na.rm=T)
    
    xy <- data.frame(x = species$collection_longitude, y = species$collection_latitude)
    pt_temp <- data.frame(extract(int_temps, SpatialPoints(xy), sp = T))[1,1]
    
    if(is.na(pt_temp)) {
      ## buffer collection points and find nearest
      xy <- gBuffer(SpatialPoints(xy), width = 15)
      pt_temp <- max(extract(int_temps, xy, sp = T)[[1]], na.rm = T)
    }
  }
  else if(species$realm == "Marine") {
    min <- min(values(mask(mar_temps, range)), na.rm=T)
    
    xy <- data.frame(x = species$collection_longitude, y = species$collection_latitude)
    pt_temp <- data.frame(extract(mar_temps, SpatialPoints(xy), sp = T))[1,1]
   
    if(is.na(pt_temp)) {
      ## buffer collection points and find nearest
      xy <- gBuffer(SpatialPoints(xy), width = 15)
      pt_temp <- max(extract(mar_temps, xy, sp = T)[[1]], na.rm = T)
    }
  }
  else {
    ## get species' mean operative temps
    terr_temps <- op_temps[[which(names(op_temps) == species$genus_species)]]
    
    min <- min(values(mask(terr_temps, range)), na.rm=T)
    
    xy <- data.frame(x = species$collection_longitude, y = species$collection_latitude)
    pt_temp <- data.frame(extract(terr_temps, SpatialPoints(xy), sp = T))[1,1]
    
    if(is.na(pt_temp)) {
      ## buffer collection points and find nearest
      xy <- gBuffer(SpatialPoints(xy), width = 15)
      pt_temp <- max(extract(terr_temps, xy, sp = T)[[1]], na.rm = T)
    }
  }

  ## save:
  data <- rbind(data, data.frame(range = species$range, min_temp = min, pt_temp = pt_temp))
  
  sp = sp + 1
}

data <- left_join(locs, data)
ggplot(data, aes(x = min_temp, y = pt_temp)) + geom_point() + geom_abline(intercept = 0, slope = 1)

data$diff <- data$pt_temp - data$min_temp
hist(data$diff)

data <- filter(data, !is.na(pt_temp), !is.infinite(pt_temp))

## plot a map
countries <- map_data("world")
col_map <- ggplot(countries, aes(x=long, y=lat, group = group)) + theme_minimal() + 
  geom_polygon(fill = "grey") + 
  coord_fixed() + 
  geom_point(data = data, aes(y = collection_latitude, x = collection_longitude, 
                              colour = diff), 
             inherit.aes = F) + 
  labs(x =  "Latitude", y = "Longitude", 
       colour = "Difference between\ncoldest body temperature\nat collection point\nand coldest realized\nbody temperature (°C)") +
  scale_colour_gradient2(mid = "white", high = "#b45346", low = 'steelblue', midpoint = 0) +
  scale_x_continuous(breaks = c(-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150)) + 
  scale_y_continuous(breaks = c(-90, -60, -30, 0, 30, 60, 90)) +
  theme(legend.title = element_text(size = 9))
  
ggsave(col_map, path = "figures/extended-data/", filename = "collection-point-map.png", width = 9, height = 7)



## map for ink scape:
map <- as.data.frame(raster('data-processed/masks/raster_terr_mask_nichemapr.grd'), xy = TRUE)
colnames(map)[1:3] <- c("longitude", "latitude", "mask")

rrs <- as.data.frame(readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds"), xy=TRUE)
colnames(rrs)[1:2] <- c("longitude", "latitude")

prs <- as.data.frame(stack(readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_acclimatized.rds")[[1]]), 
                     xy=TRUE)
colnames(prs)[1:2] <- c("longitude", "latitude")

range = 'Elgaria_multicarinata_IUCN'

## get realized range:
species <- range %>%
  str_replace_all("_", ".")

split <- str_split_fixed(range, "\\_", n = 3)
source <- split[1,3]

## get species name and realm
species <- paste(split[1,1], split[1,2], sep = ".")

rr_index <- which(str_detect(colnames(rrs), species) & str_detect(colnames(rrs), source))

r <- rrs[,c(1:2, rr_index)] 
colnames(r)[3] <- "rr"
r <- left_join(r, prs[,c(1:2, which(colnames(prs) == range))]) 
colnames(r)[4] <- "pr"
r$rr <- as.factor(r$rr)

rangemap = ggplot(r, aes(x = longitude, y = latitude)) +
  xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = rr)) + 
  scale_fill_manual(values = "#b0b0b0", na.value = "transparent") +
  annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.12,
           fill = r$pr) +
  annotate(geom="raster", x=map$longitude, y=map$latitude, alpha = 0.1,
           fill = map$mask) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

ggsave(rangemap, path = 'figures/main/didactic/', filename = 'range-map-ggplot.svg', 
       height = 6, width = 10, units = "in", device = "svg")

map <- stack('data-processed/temperature-data/terrestrial/raster_terr_high.grd')[[1]]
mask <- raster('data-processed/masks/raster_terr_mask_nichemapr.grd')
map <- mask(map, mask)

map <- as.data.frame(map, xy = T)
colnames(map)[1:2] <- c("longitude", "latitude")

base <- as.data.frame(mask, xy = T)
colnames(base)[1:2] <- c("longitude", "latitude")

rrs <- as.data.frame(readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds"), xy=TRUE)
colnames(rrs)[1:2] <- c("longitude", "latitude")

prs <- as.data.frame(stack(readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_acclimatized.rds")[[1]]), xy=TRUE)
colnames(prs)[1:2] <- c("longitude", "latitude")

range = 'Elgaria_multicarinata_IUCN'

## get realized range:
species <- range %>%
  str_replace_all("_", ".")

split <- str_split_fixed(range, "\\_", n = 3)
source <- split[1,3]

## get species name and realm
species <- paste(split[1,1], split[1,2], sep = ".")

rr_index <- which(str_detect(colnames(rrs), species) & str_detect(colnames(rrs), source))

r <- rrs[,c(1:2, rr_index)] 
colnames(r)[3] <- "rr"
r <- left_join(r, prs[,c(1:2, which(colnames(prs) == range))]) 
colnames(r)[4] <- "pr"
r$rr <- as.factor(r$rr)
r <- left_join(map[,c(1:3)], r) 
colnames(r)[3] <- "map"

r$normalized = (r$map-min(r$map, na.rm = TRUE))/(max(r$map, na.rm = TRUE)-min(r$map, na.rm = TRUE))

r <- r %>%
  mutate(rr = as.numeric(rr)) %>%
  mutate(normalized = ifelse(is.na(pr) & is.na(rr), normalized, NA)) %>%
  mutate(rr = as.factor(rr), pr = as.factor(pr))

saved = r

rangemap_coloured = ggplot(r, aes(x = longitude, y = latitude)) +
  xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = rr)) + 
  scale_fill_manual(values = "yellow", na.value = "transparent") +
  annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.3,
           fill = r$pr) +
  annotate(geom="raster", x=r$longitude, y=r$latitude,
           fill = scales::colour_ramp(c("#0083ccff","#d40000ff"))(r$normalized), alpha = 0.15) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

ggsave(rangemap_coloured, path = 'figures/main/didactic', filename = 'range-map-coloured-ggplot.svg', 
       height = 6, width = 10, units = "in", device = "svg")

r <- rrs[,c(1:2, rr_index)] 
colnames(r)[3] <- "rr"
r <- left_join(r, prs[,c(1:2, which(colnames(prs) == range))]) 
colnames(r)[4] <- "pr"
r$rr <- as.factor(r$rr)
r <- left_join(map[,c(1:3)], r) 
colnames(r)[3] <- "map"

r$normalized = (r$map-min(r$map, na.rm = TRUE))/(max(r$map, na.rm = TRUE)-min(r$map, na.rm = TRUE))

r <- r %>%
  mutate(rr = as.numeric(rr)) %>%
  mutate(normalized = ifelse(is.na(rr), normalized, NA)) %>%
  mutate(rr = as.factor(rr), pr = as.factor(pr))

rangemap_no_pr = ggplot(r, aes(x = longitude, y = latitude)) +
  xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = rr)) + 
  scale_fill_manual(values = "#b0b0b0", na.value = "transparent") +
  annotate(geom="raster", x=r$longitude, y=r$latitude,
           fill = scales::colour_ramp(c("#0083ccff","#d40000ff"))(r$normalized), alpha = 0.2) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

ggsave(rangemap_no_pr, path = 'figures/main/didactic', filename = 'range-map-only-rr.svg', 
       height = 6, width = 10, units = "in", device = "svg")

# get rid of areas of range overfilling
r <- rrs[,c(1:2, rr_index)] 
colnames(r)[3] <- "rr"
r <- left_join(r, prs[,c(1:2, which(colnames(prs) == range))]) 
colnames(r)[4] <- "pr"
r$rr <- as.factor(r$rr)
r <- left_join(map[,c(1:3)], r) 
colnames(r)[3] <- "map"

r$normalized = (r$map-min(r$map, na.rm = TRUE))/(max(r$map, na.rm = TRUE)-min(r$map, na.rm = TRUE))

r <- r %>%
  mutate(rr = as.numeric(rr)) %>%
  mutate(rr = ifelse(is.na(pr), NA, rr)) %>%
  mutate(normalized = ifelse(is.na(rr), normalized, NA)) %>%
  mutate(rr = as.factor(rr), pr = as.factor(pr))

rangemap_no_pr = ggplot(r, aes(x = longitude, y = latitude)) +
  xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = rr)) + 
  scale_fill_manual(values = "#b0b0b0", na.value = "transparent") +
  annotate(geom="raster", x=r$longitude, y=r$latitude,
           fill = scales::colour_ramp(c("#0083ccff","#d40000ff"))(r$normalized), alpha = 0.2) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

ggsave(rangemap_no_pr, path = 'figures/main/didactic', filename = 'range-map-only-rr-no-overfill.svg', 
       height = 6, width = 10, units = "in", device = "svg")


## make rr polygon with stripey pattern
library(ggpattern)

pr_poly <- rasterFromXYZ(r[,c(1,2,5)]) %>%
  rasterToPolygons(.) %>%
  aggregate(.,dissolve = TRUE) %>%
  broom::tidy()

holes <- filter(pr_poly, hole == TRUE) 
notholes <- filter(pr_poly, hole == FALSE)

stripeyboi <- ggplot(notholes, aes(x = long, y = lat, group = group)) +
  theme_void() +
  geom_polygon_pattern(
    pattern = 'stripe',
    fill = 'transparent',
    pattern_fill = 'black',
    color = 'black',
    pattern_spacing = 0.02,
    pattern_density = 0.005) +
  coord_fixed() +
  geom_polygon(data = holes, 
               fill = 'white',
               colour = 'black') 

ggsave(stripeyboi, path = 'figures/main/didactic', 
       filename = 'range-map-only-pr-stripey.png', 
       height = 6, width = 10, units = "in", device = "png", bg = "transparent")

pr_poly <- rasterFromXYZ() %>%
  rasterToPolygons(.) %>%
  aggregate(.,dissolve = TRUE)m%>%
  broom::tidy()

box <- cbind(x=c(360,0,0,360), y=c(0, 0, 180, 180))
box <- SpatialPolygons(list(Polygons(list(Polygon(box)),"1")))
plot(box)

box <-broom::tidy(box)

stripeybox <- ggplot(box, aes(x = long, y = lat, group = group)) +
  theme_void() +
  geom_polygon_pattern(
    pattern = 'stripe',
    fill = 'transparent',
    pattern_fill = 'black',
    color = 'black',
    pattern_spacing = 0.01,
    pattern_density = 0.005) +
  coord_fixed() 

ggsave(stripeybox, path = 'figures/main/didactic', 
       filename = 'range-map-stripey-box.png', 
       height = 6, width = 10, units = "in", device = "png", bg = "transparent")

### try one that shows the underlying distribution of temperatutes
rrs <- as.data.frame(stack(readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_Te.rds")[[4]]), xy=TRUE)
colnames(rrs)[1:2] <- c("longitude", "latitude")

prs <- as.data.frame(stack(readRDS("data-processed/potential-ranges/terrestrial/prs_terrestrial_acclimatized.rds")[[2]]), xy=TRUE)
colnames(prs)[1:2] <- c("longitude", "latitude")

r <- rrs[,c(1:2, rr_index)] 
colnames(r)[3] <- "rr"
r <- left_join(r, prs[,c(1:2, which(colnames(prs) == range))]) 
colnames(r)[4] <- "pr"
r$rr <- as.factor(r$rr)

rangemap = ggplot(r, aes(x = longitude, y = latitude)) +
  xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = rr)) + 
  #scale_fill_manual(values = "yellow", na.value = "transparent") +
  annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.4,
           fill = r$pr) +
  annotate(geom="raster", x=map$longitude, y=map$latitude, alpha=.1,
           fill = map$seasonal_high_temp) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 