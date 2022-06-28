## range filling figures 
library(tidyverse)
library(raster)
select <- dplyr::select

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       Prepping data for figure making     ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## niche:
uofill <- read.csv("data-processed/niche-filling/thermal-niche-filling-metrics.csv") 
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
rf <- read.csv("data-processed/range-filling/rangefilling-metrics_model-ready.csv")
nf <- read.csv("data-processed/niche-filling/thermal-niche-filling-metrics_model-ready.csv")

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


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       Plot range filling versus niche filling     ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## get rid of niche filling data for species without range filling values 
rf <- te
rf$filling <- log(rf$prop_occupied <- (rf$pr_cells - rf$u_cells) / rf$pr_cells)
rf <- rf[!is.na(rf$filling),]

nf <- Te
nf <- nf[nf$range %in% rf$range,]
nf$species <- str_replace_all(nf$species, " ", "_")

## combine data 
data <- left_join(nf, rf)
## make columns for warm and cold niche filling 
data$warm_filling <- ifelse(data$warm_over == 0, data$warm_under, data$warm_over) 
data$cold_filling <- ifelse(data$cold_over == 0, data$cold_under, data$cold_over) 

## get rid of niche overfilling, since we got rid of areas of overfilling:
data$warm_filling <- ifelse(data$warm_filling > 0, 0, data$warm_filling) 
data$cold_filling <- ifelse(data$cold_filling > 0, 0, data$cold_filling) 

## make positive so axes can say underfilling
data <- data %>%
  mutate(warm_underfilling = abs(warm_filling),
         cold_underfilling = abs(cold_filling))

## plot range filling versus warm niche filling:
ggplot(data, aes(x = filling, y = warm_underfilling)) + geom_point() +
  theme_minimal() +
  labs(y = "Warm niche underfilling (°C)", x = "Range filling")

ggplot(data, aes(x = filling, y = cold_underfilling)) + geom_point() +
  theme_minimal() +
  labs(y = "Cold niche underfilling (°C)", x = "Range filling")

ggplot(data, aes(x = cold_underfilling, y = warm_underfilling, size = filling)) + geom_point() +
  theme_minimal() +
  labs(y = "Warm niche underfilling (°C)", x = "Cold niche underfilling (°C)", size = "Range filling")

ggplot(data, aes(x = log(filling), y = warm_underfilling)) + geom_point() +
  theme_minimal() +
  labs(y = "Warm niche underfilling (°C)", x = "Ln range filling")

ggplot(data, aes(x = log(filling), y = cold_underfilling)) + geom_point() +
  theme_minimal() +
  labs(y = "Cold niche underfilling (°C)", x = "Ln range filling")

ggplot(data, aes(x = cold_underfilling, y = warm_underfilling, size = log(filling))) + 
  geom_point(position = position_jitter()) +
  theme_minimal() +
  labs(y = "Warm niche underfilling (°C)", x = "Cold niche underfilling (°C)", 
       size = "Ln range filling") +
  scale_size(range = c(0.5, 2)) 





## where is the realized niche position within the potential niche?
data <- data %>%
  mutate(r_niche_upper = ifelse(r_niche_upper > p_niche_upper, p_niche_upper, r_niche_upper)) %>%
  mutate(r_niche_lower = ifelse(r_niche_lower < p_niche_lower, p_niche_lower, r_niche_lower))

data$realized_position = (data$r_niche_upper - data$r_niche_lower)/2
data$potential_position = (data$p_niche_upper - data$p_niche_lower)/2

data$relative_niche_position = data$potential_position - data$realized_position


ggplot(data, aes(x = filling, y = relative_niche_position)) + geom_point() +
  labs(x = "Occupied proportion of potential range", 
       y = "Relative niche position\n(potential - realized niche midpoint)") +
  theme_minimal()

ggplot(data, aes(x = cold_underfilling, y = warm_underfilling, size = log(filling),
                 colour = relative_niche_position)) + 
  geom_point(position = position_jitter()) +
  theme_minimal() +
  labs(y = "Warm niche underfilling (°C)", x = "Cold niche underfilling (°C)", 
       size = "Ln range filling") +
  scale_size(range = c(0.5, 2)) 



## make nice plots:
data <- gather(data, key = "facet", value = "underfilling_value", 
                     c(warm_underfilling, cold_underfilling))
  

warm_range <- data %>%
  filter(facet == "warm_underfilling") %>%
  ggplot(., aes(x = log(filling), y = underfilling_value)) + 
  geom_point(aes(colour = underfilling_value)) +
  labs(y = "Shortfall of temperatures occupied (°C)", x = "Ln proportion of occupied potential range") +
  scale_colour_gradient(high = "#b45346", low = "white") +
  theme_bw() +
  theme(axis.line.x.top = element_line(),
        axis.line.y.right = element_line(), 
        panel.border = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank())  +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 1) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 1) +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(-7.5, 0)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  geom_point(shape = 1, aes(x = log(filling), y = underfilling_value), 
             inherit.aes = F, stroke = 0.5) 

cold_range <- data %>%
  filter(facet == "cold_underfilling") %>%
  ggplot(., aes(x = log(filling), y = underfilling_value)) + 
  geom_point(aes(colour = underfilling_value)) +
  labs(y = "Shortfall of temperatures occupied (°C)", x = "Ln proportion of occupied potential range") +
  scale_colour_gradient(high = "steelblue", low = "white") +
  theme_bw() +
  theme(axis.line.x.top = element_line(),
        axis.line.y.right = element_line(), 
        panel.border = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank())  +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 1) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 1) +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(-7.5, 0)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  geom_point(shape = 1, aes(x = log(filling), y = underfilling_value), 
             inherit.aes = F, stroke = 0.5) 


ggplot(data, aes(x = log(filling), y = underfilling_value, colour = facet)) + 
  geom_point() +
  labs(y = "Shortfall of temperatures occupied (°C)", 
       x = "Ln proportion of occupied potential range") +
  facet_wrap(~facet) +
  scale_colour_manual(values = c("steelblue", "#b45346")) +
  theme_bw() +
  theme(axis.line.x.top = element_line(),
        axis.line.y.right = element_line(), 
        panel.border = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank())  +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 1) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 1) +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(-7.5, 0)) +
  scale_y_continuous(expand = c(0.01,0.01)) 
  


rel_niche_pos <- ggplot(data, aes(x = log(filling), y = relative_niche_position, colour = relative_niche_position)) + 
  geom_point() +
  labs(y = "Potential - realized niche\nmidpoint (°C)", x = "Ln proportion of occupied potential range") +
  theme_bw() +
  theme(axis.line.x.top = element_line(),
        axis.line.y.right = element_blank(), 
        panel.border = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank())  +
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  scale_x_continuous(expand = c(0.01,0.01), limits = c(-7.5, 0)) +
  scale_y_continuous(expand = c(0.01,0.01), limits = c(-11, 11)) +
  scale_colour_gradient(high = "#b45346", low = "white") +
  geom_point(shape = 1, aes(x = log(filling), y = relative_niche_position), 
             inherit.aes = F, stroke = 0.5)


ggsave(rel_niche_pos, path ="figures/range-filling", filename = "rel-niche-position.png",
       width = 4, height = 3, device = "png")
ggsave(cold_range, path ="figures/range-filling", filename = "cold-range-filling.png",
       width = 2, height = 2, device = "png")
ggsave(warm_range, path ="figures/range-filling", filename = "warm-range-filling.png",
       width = 2, height = 2, device = "png")

cold_range = cold_range + labs(y = "", x = "")

warm_range = warm_range + labs(y = "", x = "")

grob <- grid.arrange(cold_range, warm_range, nrow = 2)


grid <- grid.arrange(arrangeGrob(grob, 
                                  left = textGrob("Shortfall of temperatures occupied (°C)", 
                                                  gp=gpar(size = 6), rot=90)))

ggsave(grid, path ="figures/range-filling", filename = "warm-cold-range-filling.png",
       width = 4, height = 3, device = "png")



## plot collection point on map 
## colour by difference between 95% percentile of cold body temperature across range and species' coldest body temperature at collection point 
library(rgeos)
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

op_temps <- stack("data-processed/temperature-data/operative-temperatures/Te_shade_min.grd")
mar_temps <- stack("data-processed/temperature-data/raster_marine_low.grd")
int_temps <- stack("data-processed/temperature-data/raster_intertidal_low.grd") 

## read in niche limits to get species list 
filling <- read.csv("data-processed/niche-filling/thermal-niche-filling-metrics_model-ready.csv") %>%
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
  
ggsave(col_map, path = "figures/extended-data/", filename = "collection_map.png", width = 9, height = 7)



## map for ink scape:
map <- as.data.frame(raster('data-processed/masks/raster_terr_mask_nichemapr.grd'), xy = T)
colnames(map)[1:3] <- c("longitude", "latitude", "mask")

rrs <- as.data.frame(stack("data-processed/realized-ranges/rasterized_rrs_nichemapr.grd"), xy=TRUE)
colnames(rrs)[1:2] <- c("longitude", "latitude")

prs <- as.data.frame(stack(readRDS("data-processed/potential-ranges/prs_terrestrial_acclimatized.rds")[[1]]), xy=TRUE)
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


ggsave(rangemap, path = 'figures/didactic/', filename = 'range-map-ggplot.svg', 
       height = 6, width = 10, units = "in", device = "svg")

map <- stack('data-processed/temperature-data/raster_terr_high.grd')[[1]]
mask <- raster('data-processed/masks/raster_terr_mask_nichemapr.grd')
map <- mask(map, mask)

map <- as.data.frame(map, xy = T)
colnames(map)[1:2] <- c("longitude", "latitude")

base <- as.data.frame(mask, xy = T)
colnames(base)[1:2] <- c("longitude", "latitude")

rrs <- as.data.frame(stack("data-processed/realized-ranges/rasterized_rrs_nichemapr.grd"), xy=TRUE)
colnames(rrs)[1:2] <- c("longitude", "latitude")

prs <- as.data.frame(stack(readRDS("data-processed/potential-ranges/prs_terrestrial_acclimatized.rds")[[1]]), xy=TRUE)
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

ggsave(rangemap_coloured, path = 'figures/didactic/', filename = 'range-map-coloured-ggplot.svg', 
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

ggsave(rangemap_no_pr, path = 'figures/didactic/', filename = 'range-map-only-rr.svg', 
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

ggsave(rangemap_no_pr, path = 'figures/didactic/', filename = 'range-map-only-rr-no-overfill.svg', 
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
  

ggsave(stripeyboi, path = 'figures/didactic/', 
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

ggsave(stripeybox, path = 'figures/didactic/', 
       filename = 'range-map-stripey-box.png', 
       height = 6, width = 10, units = "in", device = "png", bg = "transparent")




### try one that shows the underlying distribution of temperatutes
rrs <- as.data.frame(stack(readRDS("data-processed/potential-ranges/prs_terrestrial_Te.rds")[[4]]), xy=TRUE)
colnames(rrs)[1:2] <- c("longitude", "latitude")

prs <- as.data.frame(stack(readRDS("data-processed/potential-ranges/prs_terrestrial_acclimatized_Te.rds")[[2]]), xy=TRUE)
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





### for new figure showing restrictions
## read in saved data for Elgaria_multicarinata_IUCN:
pr <- readRDS("data-processed/potential_raster_Em.rds")
elev <- readRDS("data-processed/elev_Em.rds")
realm <- readRDS("data-processed/realm_Em.rds")
pr_restricted <- readRDS("data-processed/pr-restricted_Em.rds")
rniche_res <- readRDS("data-processed/p-restricted_rniche.rds")

pr <- as.data.frame(pr, xy = T)
colnames(pr)[1:3] <- c("longitude", "latitude", "pr")

elev <- as.data.frame(elev, xy = T)
colnames(elev)[1:3] <- c("longitude", "latitude", "elev")

realm <- as.data.frame(realm, xy=TRUE)
colnames(realm)[1:3] <- c("longitude", "latitude", "realm")

pr_restricted <- as.data.frame(pr_restricted, xy=TRUE)
colnames(pr_restricted)[1:3] <- c("longitude", "latitude", "pr_res")

rniche_res <- as.data.frame(rniche_res, xy=TRUE)
colnames(rniche_res)[1:3] <- c("longitude", "latitude", "rniche_res")

map <- stack('data-processed/temperature-data/raster_terr_high.grd')[[1]] 
map <- as.data.frame(map, xy = T)
colnames(map)[1:3] <- c("longitude", "latitude", "map")

map$normalized = (map$map-min(map$map, na.rm = TRUE))/(max(map$map, na.rm = TRUE)-min(map$map, na.rm = TRUE))

mask <- map[1:3]
colnames(mask)[1:3] <- c("longitude", "latitude", "mask")
mask = mask %>%
  mutate(mask = ifelse(!is.na(mask), 1, mask))

rr_index <- which(str_detect(colnames(rrs), species) & str_detect(colnames(rrs), source))

r <- rrs[,c(1:2, rr_index)] 
colnames(r)[3] <- "rr"

test <- map

map <- raster("data-processed/temperature-data/raster_marine_mean.grd")
map <- as.data.frame(map[[1]], xy = T)
colnames(map)[1:3] <- c("longitude", "latitude", "map")
map$normalized = (map$map-min(map$map, na.rm = TRUE))/(max(map$map, na.rm = TRUE)-min(map$map, na.rm = TRUE))
test <- map

one <- ggplot(test, aes(x = longitude, y = latitude, alpha = 0.4)) +
  xlim(-180, 180) + ylim(-90, 90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = normalized)) + 
  scale_fill_gradient2(low = "#0083ccff", high = "#d40000ff", mid = "white", midpoint = 0.75,
                       na.value = "transparent") +
  # annotate(geom="raster", x=mask$longitude, y=mask$latitude, alpha = 0.1,
  #          fill = mask$mask) +
  # annotate(geom="raster", x=r$longitude, y=r$latitude, alpha = 0.5,
  #          fill = r$rr) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

test <- left_join(pr, map[,c(1:4)]) 
test <- filter(test, !is.na(pr))

two <- ggplot(test, aes(x = longitude, y = latitude, alpha = 0.4)) +
  xlim(-180, 180) + ylim(-90, 90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = normalized)) + 
  scale_fill_gradient2(low = "#0083ccff", high = "#d40000ff", mid = "white", midpoint = 0.75,
                                                            na.value = "transparent") +
  annotate(geom="raster", x=mask$longitude, y=mask$latitude, alpha = 0.1,
           fill = mask$mask) +
  annotate(geom="raster", x=r$longitude, y=r$latitude, alpha = 0.5,
           fill = r$rr) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

test <- left_join(elev, map[,c(1:4)]) 
test <- filter(test, !is.na(elev))

three <- ggplot(test, aes(x = longitude, y = latitude, alpha = 0.3)) +
  xlim(-180, 180) + ylim(-90, 90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = normalized)) + 
  scale_fill_gradient2(low = "#0083ccff", high = "#d40000ff", mid = "white", midpoint = 0.75,
                       na.value = "transparent") +
  annotate(geom="raster", x=mask$longitude, y=mask$latitude, alpha = 0.1,
           fill = mask$mask) +
  annotate(geom="raster", x=r$longitude, y=r$latitude, alpha = 0.5,
           fill = r$rr) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

test <- left_join(pr_restricted, map[,c(1:4)]) 
test <- filter(test, !is.na(pr_res))

four <- ggplot(test, aes(x = longitude, y = latitude, alpha = 0.5)) +
  xlim(-180, 180) + ylim(-90, 90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = normalized)) + 
  scale_fill_gradient2(low = "#0083ccff", high = "#d40000ff", mid = "white", midpoint = 0.75,
                       na.value = "transparent") +
  annotate(geom="raster", x=mask$longitude, y=mask$latitude, alpha = 0.03,
           fill = mask$mask) +
  annotate(geom="raster", x=r$longitude, y=r$latitude, alpha = 0.5,
           fill = r$rr) +
  annotate(geom="raster", x=realm$longitude, y=realm$latitude, alpha = 0.07,
           fill = realm$realm) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

ggsave(one, path = 'figures/didactic/', 
       filename = 'range-making_1.png', 
       height = 6, width = 10, units = "in", device = "png", bg = "transparent")

ggsave(two, path = 'figures/didactic/', 
       filename = 'range-making_2.png', 
       height = 6, width = 10, units = "in", device = "png", bg = "transparent")

ggsave(three, path = 'figures/didactic/', 
       filename = 'range-making_3.png', 
       height = 6, width = 10, units = "in", device = "png", bg = "transparent")

ggsave(four, path = 'figures/didactic/', 
       filename = 'range-making_4.png', 
       height = 6, width = 10, units = "in", device = "png", bg = "transparent")

pr_poly <- rasterToPolygons(rasterFromXYZ(rniche_res)) %>%
  aggregate(., dissolve = TRUE) %>%
  broom::tidy()

holes <- filter(pr_poly, hole == TRUE) 
notholes <- filter(pr_poly, hole == FALSE)

cropped_test <- filter(test, latitude > 0, longitude < -50)
cropped_mask <- filter(mask, latitude > 0, longitude < -50)
cropped_r <- filter(r, latitude > 0, longitude < -50)
cropped_realm <- filter(realm, latitude > 0, longitude < -50)

lat_mp <- cropped_r %>%
  filter(!is.na(rr)) 

split_range <- cropped_test %>%
  ggplot(., aes(x = longitude, y = latitude, alpha = 0.5)) +
  xlim(-180, 180) + ylim(-90, 90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = normalized)) + 
  scale_fill_gradient2(low = "white", high = "white", mid = "white", midpoint = 0.75,
                       na.value = "transparent") +
  annotate(geom="raster", x=cropped_mask$longitude, y=cropped_mask$latitude, alpha = 0.1,
           fill = cropped_mask$mask) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") + 
  geom_polygon_pattern(data = notholes, aes(x = long, y = lat, fill = hole,
                                         group = piece),
    pattern = 'stripe',
    fill = 'grey90',
    pattern_fill = 'black',
    color = 'black',
    pattern_spacing = 0.03,
    pattern_density = 0.005,
    size = 0.5) +
  annotate(geom="raster", x=cropped_r$longitude, y=cropped_r$latitude, alpha = 0.70,
           fill = cropped_r$rr) +
  geom_segment(aes(y = (max(lat_mp$latitude) + min(lat_mp$latitude))/2,
               yend = (max(lat_mp$latitude) + min(lat_mp$latitude))/2,
               x = -130, xend = -68))

ggsave(split_range, path = 'figures/didactic/', 
       filename = 'split-range.png', 
       height = 4, width = 4, units = "in", device = "png", bg = "transparent")


cropped_test <- filter(test, latitude  > 15,latitude  < 80, longitude < -75, longitude > -140)
cropped_mask <- filter(mask, latitude > 15, latitude  < 80, longitude < -75, longitude > -140)
cropped_r <- filter(r, latitude > 15, latitude  < 80, longitude < -75, longitude > -140)
buffer <- buffer(rasterFromXYZ(r), width = 250000)
plot(buffer)

buffer <- as.data.frame(buffer, xy=TRUE)
colnames(buffer)[1:3] <- c("longitude", "latitude", "buffer")
buffer <- filter(buffer, latitude > 15, latitude  < 80, longitude < -75, longitude > -140)

zoomy <- cropped_test %>%
  ggplot(., aes(x = longitude, y = latitude, alpha = 0.5)) +
  xlim(-180, 180) + ylim(-90, 90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = normalized)) + 
  scale_fill_gradient2(low = "white", high = "white", mid = "white", midpoint = 0.75,
                       na.value = "transparent") +
  annotate(geom="raster", x=cropped_mask$longitude, y=cropped_mask$latitude, alpha = 0.1,
           fill = cropped_mask$mask) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") +
  annotate(geom="raster", x=cropped_r$longitude, y=cropped_r$latitude, alpha = 0.5,
           fill = cropped_r$rr) +
  annotate(geom="raster", x=buffer$longitude, y=buffer$latitude, alpha = 0.5,
           fill = buffer$buffer) +
  geom_segment(aes(y = (max(lat_mp$latitude) + min(lat_mp$latitude))/2,
                   yend = (max(lat_mp$latitude) + min(lat_mp$latitude))/2,
                   x = -130, xend = -110))

ggsave(zoomy, path = 'figures/didactic/', 
       filename = 'zoomed-range.png', 
       height = 4, width = 4, units = "in", device = "png", bg = "transparent")

test <- map
cropped_test <- filter(test, latitude  > 15,latitude  < 80, longitude < -75, longitude > -140)

range_potential <- test %>%
  ggplot(., aes(x = longitude, y = latitude)) +
  xlim(-180, 0) + ylim(0, 90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = normalized)) +
  scale_fill_gradient2(low = "#0083ccff", high = "#d40000ff", mid = "white", midpoint = 0.6,
                       na.value = "transparent") +
  # annotate(geom="raster", x=cropped_mask$longitude, y=cropped_mask$latitude, alpha = 0.1,
  #          fill = cropped_mask$mask)  +
  theme_void() +
  theme(legend.position = "none") 
+ 
  geom_polygon(data = notholes, aes(x = long, y = lat, group = piece), fill = "transparent", colour = "black")



ggplot(test, aes(x = longitude, y = latitude, alpha = 0.4)) +
  xlim(-180, 180) + ylim(-90, 90) + coord_fixed(ratio = 1) +
  geom_raster(aes(fill = normalized)) + 
  scale_fill_gradient2(low = "#0083ccff", high = "#d40000ff", mid = "white", midpoint = 0.75,
                       na.value = "transparent") +
  annotate(geom="raster", x=mask$longitude, y=mask$latitude, alpha = 0.1,
           fill = mask$mask) +
  annotate(geom="raster", x=r$longitude, y=r$latitude, alpha = 0.5,
           fill = r$rr) +
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) + theme_void() +
  theme(legend.position = "none") 

ggsave(range_potential, path = 'figures/didactic/', 
       filename = 'temp.png', 
       height = 4, width = 4, units = "in", device = "png", bg = "transparent")
       