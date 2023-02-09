### model selection and averaging for niche filling
library(tidyverse)
library(nlme)
library(car)
library(MuMIn)
library(knitr)
library(cowplot)
library(raster)
library(grid)
library(gridExtra)
library(ggridges)
select <- dplyr::select
rename <-dplyr::rename

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Prepping for model fitting          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
uofill <- read.csv("data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics.csv") 
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
  arrange(genus_species, type, source) 

## give priority to IUCN ranges (but test sensitivity later)
data <- data %>%
  mutate(temp = paste(genus_species, type, sep = '')) %>%
  filter(!duplicated(temp)) %>%
  select(-temp) %>%
  rename("sensitivity_type" = type) 

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


## merge with species traits and geographic predictors:
data <- left_join(data, traits) %>%
  left_join(., thermal_limits) %>%
  unique(.)

colnames(data)[30] <- "dispersal_distance_category"
colnames(data)[32] <- "dispersal_ability_category"

data$sensitivity_type <- factor(data$sensitivity_type, 
                                levels = c("Te", "Te_tpreftb", "acclimatized", "Te_tpref"),
                                ordered = TRUE)

data <- data %>%
  mutate(dispersal_distance_continuous = ifelse(dispersal_distance_category == "0-1", 1, 
                                                ifelse(dispersal_distance_category == "1-10", 10,
                                                       ifelse(dispersal_distance_category == "10-100", 100,
                                                              ifelse(dispersal_distance_category == "100+", 1000, 
                                                                     NA)))))
data$abs_lat_mp <- abs(data$lat_mp)
data$log_maximum_body_size <- log(data$maximum_body_size_SVL_HBL_cm_)

data$dispersal_ability_category_general <- ifelse(str_detect(data$dispersal_ability_category, 
                                                             "non-pelagic"), "non-pelagic", 
                                                  as.character(data$dispersal_ability_category)) 

data$dispersal_ability_category_general <- ifelse(str_detect(data$dispersal_ability_category_general,
                                                             "pelagic") & 
                                                    !str_detect(data$dispersal_ability_category_general,
                                                                "non-pelagic"),
                                                  'pelagic', 
                                                  as.character(data$dispersal_ability_category_general)) 

## get realized range area as number of cells
rrs <- readRDS("data-processed/realized-ranges/rasterized_rrs_nichemapr.rds")
names(rrs) <- str_replace_all(names(rrs), '\\.', '_')
for (i in 1:nlayers(rrs)) {
  rr <- rrs[[i]]

  if (i == 1) {
    range_area <- data.frame(range = names(rrs)[i],
                             rr_cells = length(which(!is.na(values(rr)))))
  }
  else {
    range_area <- rbind(range_area, data.frame(range = names(rrs)[i],
                                               rr_cells = length(which(!is.na(values(rr))))))
  }
}

data <- left_join(data, range_area) %>%
  mutate(log_range_area = log(rr_cells)) %>%
  filter(!is.infinite(log_range_area))

## write out:
write.csv(data, "data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics_model-ready.csv", row.names = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####            Niche filling models             #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in the data
data <- read.csv("data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics_model-ready.csv")

thermal_limits <- read.csv("data-processed/traits/thermal-limits_ectotherms-with-ranges_taxized.csv")%>%
  select(genus_species, type, metric) %>%
  unique(.)

## re-order factors to give desired contrasts
data$realm <- relevel(factor(data$realm), ref = "Terrestrial")

## split data by sensitivity type
sens_types <- group_split(data, sensitivity_type)

Te <- sens_types[[2]]
Te_tpreftb <- sens_types[[4]]
acc <- sens_types[[1]]
Te_subset_tpreftb <- Te[which(Te$genus_species %in% Te_tpreftb$genus_species),]
Te_subset_acc <- filter(Te, Te$genus_species %in% acc$genus_species)

## fitting to cold and warm filling separately
cold <- Te %>%
  filter(edge_type == "cold") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "min")  
## assign metric type used to measure cold limit 

warm <- Te %>%
  filter(edge_type == "warm") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "max") 
## assign metric type used to measure warm limit 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####           The ~trait~ models                #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#First, look for multicolinearity among our variables
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#check for multicolinearity among variables (in a linear model)
modvif_cold <- lm(filling_value ~ abs_lat_mp + log_range_area + realm  + dispersal_ability_category +
                    dispersal_distance_continuous  +
                    log_maximum_body_size + hot_season_dormancy_ + 
                    cold_season_dormancy_ + metric,
                  data = cold)
vif(modvif_cold)
## removed:  dispersal_ability_category

modvif_warm <- lm(filling_value ~ abs_lat_mp + log_range_area + realm + dispersal_ability_category +
                    dispersal_distance_continuous +
                    log_maximum_body_size + hot_season_dormancy_ + 
                    cold_season_dormancy_ + metric,
                  data = warm)
vif(modvif_warm)
## removed: dispersal_ability_category 

#Second, rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
cold <- select(cold, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size, cold_season_dormancy_,metric,
                       Class, Order, Family, Genus, Species))
warm <- select(warm, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size, hot_season_dormancy_,metric, 
                       Class, Order, Family, Genus, Species))

## write list of species
write.csv(unique(paste(cold$Genus, cold$Species, sep=" ")), 
          "data-processed/thermal-niches/splist_cold-niche-filling.csv",
          row.names = FALSE)
write.csv(unique(paste(warm$Genus, warm$Species, sep=" ")), 
          "data-processed/thermal-niches/splist_warm-niche-filling.csv",
          row.names = FALSE)

# get complete cases
cold <- subset(cold, complete.cases(cold))
# dim(cold)
warm <- subset(warm, complete.cases(warm))
# dim(warm)

## write list of species
write.csv(unique(paste(cold$Genus, cold$Species, sep=" ")), 
          "data-processed/thermal-niches/splist_cold-niche-filling_model.csv",
          row.names = FALSE)
write.csv(unique(paste(warm$Genus, warm$Species, sep=" ")), 
          "data-processed/thermal-niches/splist_warm-niche-filling_model.csv",
          row.names = FALSE)

#scale the continuous variables by subtracting mean and dividing by sd
means <- sds <- c()
cold = as.data.frame(cold)
coldnorm <- cold
x = 1
for (i in c(2,3,5,6)) {
  coldnorm[,i] <- (cold[,i]-mean(c(cold[,i])))/sd(cold[,i])
  means[x] = mean(c(cold[,i]))
  sds[x] = sd(c(cold[,i]))
  x = x+1
}
## save scaling factors:
scalers = data.frame(type = "cold", var = colnames(cold)[c(2,3,5,6)], means = means, sds = sds)

warm = as.data.frame(warm)
warmnorm <- warm
x = 1
for (i in c(2,3,5,6)) {
  warmnorm[,i] <- (warm[,i]-mean(c(warm[,i])))/sd(warm[,i])
  means[x] = mean(c(warm[,i]))
  sds[x] = sd(c(warm[,i]))
  x = x+1
}
scalers = rbind(scalers,
                data.frame(type = "warm", var = colnames(cold)[c(2,3,5,6)], means = means, sds = sds))


warm = warmnorm
cold = coldnorm


#Now, run the full models and check out the residuals:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#run full lme model
model_cold <- lme(filling_value ~ abs_lat_mp*realm + log_range_area +
                    dispersal_distance_continuous  +
                    log_maximum_body_size + metric,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = cold
                  
)
#check residuals
hist(resid(model_cold))
plot(model_cold)

#against predictors:
E1 <- resid(model_cold)
plot(E1 ~ abs_lat_mp, data = cold)
plot(E1 ~ log_range_area, data = cold)
plot(E1 ~ realm, data = cold)
plot(E1 ~ dispersal_distance_continuous, data = cold)

model_warm <- lme(filling_value ~ abs_lat_mp*realm + log_range_area +
                    dispersal_distance_continuous  +
                    log_maximum_body_size + metric,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = warm
)

#check residuals
hist(resid(model_warm))
plot(model_warm)

E1 <- resid(model_warm)
plot(E1 ~ abs_lat_mp, data = warm)
plot(E1 ~ log_range_area, data = warm)
plot(E1 ~ realm, data = warm)
plot(E1 ~ dispersal_distance_continuous, data = warm)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_cold <- dredge(model_cold, extra="R^2"))
confset.95p_cold <- get.models(allmodels_cold, subset = cumsum(weight)<=.95) #get confidence set, 6 models
avgm_cold <- model.avg(confset.95p_cold) #do averaging
summary(avgm_cold)

sum <- summary(avgm_cold)
df <- as.data.frame(sum$coefmat.full) 
CI <- as.data.frame(confint(avgm_cold, full=T, level = 0.95)) # get confidence intervals for full model
df$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df$CI.max <-CI$`97.5 %`
data.table::setDT(df, keep.rownames = "coefficient") #put rownames into column
df$how_conf = ifelse((df$CI.min < 0 & df$CI.max < 0) | (df$CI.min > 0 & df$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0

whisker_cold <- ggplot(data=df, aes(x=fct_rev(coefficient), y=Estimate)) + 
  geom_hline(yintercept = 0) + 
  theme_classic() +
  geom_pointrange(aes(ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = '#4682b4', size = 1, shape = 21) + 
  coord_flip() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("white", "#4682b4")) +
  labs(y = "Effect of variable on cool niche filling", x = "") + 
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal",
                              "Metric: lethal",
                              "Realized range size (log no. cells)", 
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Abs. realized range latitudinal midpoint ", 
                              "Reference"), drop = FALSE) +
  scale_y_continuous(limits = c(-18,18))


# plotting estimates (fixed effects) 
suppressWarnings(allmodels_warm <- dredge(model_warm, extra="R^2")) 
confset.95p_warm <- get.models(allmodels_warm, subset = cumsum(weight)<=.95) #get confidence set,6 models
avgm_warm <- model.avg(confset.95p_warm) #do averaging
summary(avgm_warm)

# plotting estimates (fixed effects) 
sum <- summary(avgm_warm)
df <- as.data.frame(sum$coefmat.full) 
CI <- as.data.frame(confint(avgm_warm, full=T, level = 0.95)) # get confidence intervals for full model
df$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df$CI.max <-CI$`97.5 %`
data.table::setDT(df, keep.rownames = "coefficient") #put rownames into column
df$how_conf = ifelse((df$CI.min < 0 & df$CI.max < 0) | (df$CI.min > 0 & df$CI.max > 0), "so_conf", "not")# add column to show which predictors have CIs that don't cross 0

whisker_warm <- ggplot(data=df, aes(x=fct_rev(coefficient), y=Estimate)) + 
  geom_hline(yintercept = 0) + 
  theme_classic() +
  geom_pointrange(aes(ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = '#b45346', size = 1, shape = 21) + 
  coord_flip() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("white", "#b45346")) +
  labs(y = "Effect of variable on warm niche filling", x = "") + 
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal",
                              "Metric: lethal",
                              "Realized range size (log no. cells)", 
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "Abs. realized range latitudinal midpoint  x realm: subtidal",
                              "Abs. realized range latitudinal midpoint  x realm: intertidal",
                              "Abs. realized range latitudinal midpoint ", 
                              "Reference")) +
  scale_y_continuous(limits = c(-18,18))

## save figures:
saveRDS(whisker_warm, "data-processed/intermediate-files/whisker_warm.rds")
saveRDS(whisker_cold, "data-processed/intermediate-files/whisker_cold.rds")

ggsave(whisker_warm, path ="figures/extended-data", 
       filename = "whisker-plot_warm-niche-filling.png",
       width = 6, height = 3, device = "png")
ggsave(whisker_cold, path ="figures/extended-data", 
       filename = "whisker-plot_cold-niche-filling.png",
       width = 6, height = 3, device = "png")

## save complete cases
saveRDS(coldnorm, "data-processed/thermal-niches/niche-filling/cold-complete-cases.rds")
saveRDS(warmnorm, "data-processed/thermal-niches/niche-filling/warm-complete-cases.rds")


### assess model averaging by looking at top model set 
library(dotwhisker)
library(RColorBrewer)
## add average model to lists so it can be visualized 
cold_list <- append(confset.95p_cold, list(avgm_cold))
names(cold_list)[7] <- "Average"

colours <- colorRampPalette(brewer.pal(8, "Set2"))(7)
colours[1] <- "black"

cold_dwplot <- dwplot(cold_list,  vline = geom_vline(xintercept = 0, colour = "grey50"),
                      show_intercept = TRUE) + 
  scale_y_discrete(labels = c("Dispersal distance (km)",
                              "Maximum body size (log cm)",
                              "SD (observations)",
                              "SD (intercept)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Realm: subtidal",
                              "Realm: intertidal",
                              "Metric: lethal",
                              "Realized range size (log no. cells)", 
                              "Abs. realized range latitudinal midpoint ", 
                              "Reference")) +
  labs(colour = '', x = "Effect of variable on cool niche filling") +
  theme_light() +
  scale_color_manual(values = colours)  +
  theme(panel.grid = element_blank()) + 
  scale_x_continuous(limits = c(-18,18))

saveRDS(cold_dwplot, "data-processed/intermediate-files/whisker_cold_model-avg.rds")

warm_list <- append(confset.95p_warm, list(avgm_warm))
names(warm_list)[7] <- "Average"

colours <- colorRampPalette(brewer.pal(8, "Set2"))(11)
colours[1] <- "black"

warm_dwplot <- dwplot(warm_list, 
                      vline = geom_vline(xintercept = 0, colour = "grey50"),
                      show_intercept = TRUE) + 
  scale_y_discrete(labels = c("Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "SD (observations)",
                              "SD (intercept)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Realm: subtidal",
                              "Realm: intertidal",
                              "Metric: lethal",
                              "Realized range size (log no. cells)",
                              "Abs. realized range latitudinal midpoint ",
                              "Reference")) +
  labs(colour = '', x = "Effect of variable on warm niche filling") +
  theme_light() +
  scale_color_manual(values = colours) +
  theme(panel.grid = element_blank())  + 
  scale_x_continuous(limits = c(-18,18))

saveRDS(warm_dwplot, "data-processed/intermediate-files/whisker_warm_model-avg.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####                Exporting tables           ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_warm)
sum_warm <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", "slope", "slope", "slope", 
              "intercept","intercept",
              "intercept")
coefs <- sum_warm[c(1,2,7,8,3,9,10,5,6,4),1] # reorder
f_effects <-  c("(Intercept)",
                           "Abs. realized range latitudinal midpoint ",
                           "Abs. realized range latitudinal midpoint x realm: intertidal)",
                           "Abs. realized range latitudinal midpoint x realm: subtidal)",
                           "Realized range size (log no.cells)",
                           "Dispersal distance (km)",
                           "Maximum body size (log cm)",
                           "Realm: intertidal",
                           "Realm: subtidal",
                           "Metric: lethal")
std_err <- sum_warm[c(1,2,7,8,3,9,10,5,6,4),2]  # reorder fixed effects
z_val <- sum_warm[c(1,2,7,8,3,9,10,5,6,4),4]
p_val <- sum_warm[c(1,2,7,8,3,9,10,5,6,4),5]
names(coefs) <- names(std_err) <- names(z_val) <- names(p_val) <- NULL

## put all into a table:
results_w <- data.frame("fixed effects" = f_effects, 
                        "effect type" = eff_type,
                        "estimate" = coefs,
                        "s.e." = std_err,
                        "z-value" = z_val,
                        "p-value" = p_val)
colnames(results_w) = c("fixed effects",	"effect type",	"estimate", "s.e.",'z-value',	"p-value")

sum <- summary(avgm_cold)
sum_cold <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", "slope", "slope", "slope", 
              "intercept","intercept",
              "intercept")
coefs <- sum_cold[c(1,2,7,8,3,10,9,5,6,4),1] # reorder
f_effects = c("(Intercept)",
              "Abs. realized range latitudinal midpoint ",
              "Abs. realized range latitudinal midpoint x realm: intertidal)",
              "Abs. realized range latitudinal midpoint x realm: subtidal)",
              "Realized range size (log no. cells)",
              "Dispersal distance (km)",
              "Maximum body size (log cm)",
              "Realm: intertidal",
              "Realm: subtidal",
              "Metric: lethal")
std_err <- sum_cold[c(1,2,7,8,3,10,9,5,6,4),2] # reorder fixed effects
z_val <- sum_cold[c(1,2,7,8,3,10,9,5,6,4),4]
p_val <- sum_cold[c(1,2,7,8,3,10,9,5,6,4),5]
names(coefs) <- names(std_err) <- names(z_val) <- names(p_val) <- NULL

results_c <- data.frame("fixed effects" = f_effects, 
                        "effect type" = eff_type,
                        "estimate" = coefs,
                        "s.e." = std_err,
                        "z-value" = z_val,
                        "p-value" = p_val)
colnames(results_c) = c("fixed effects",	"effect type",	"estimate", "s.e.",'z-value',	"p-value")

results <- rbind(results_w, results_c)
results[c(c(3,4,5))] <- round(results[,c(3,4,5)], digits = 2) # round estimates to 2 decimal places
results[6] <- round(results[,6], digits = 3) # round p values to 3 decimal places
## add significance indicators
results$`p-value` <- ifelse(results$`p-value` == 0, "<0.001 **", 
       ifelse(results$`p-value` <= 0.01, paste(results$`p-value`, " **", sep= ""), 
              ifelse(results$`p-value` <= 0.05, paste(results$`p-value`, " *", sep= ""),
                     ifelse(results$`p-value` <= 0.1, paste(results$`p-value`, " ", sep= ""),
                            as.character(results$`p-value`)))))
table <- results %>%
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
  htmlTable(rnames = c("Warm niche edge", rep("", 10), 
                                       "Cold niche edge", rep("", 10)))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####          Plotting main results           ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# function to find mode of each predictor
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####             Across range area            ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
cold = coldnorm

int_cd <- cold %>% filter(realm == "Intertidal")
mar_cd <- cold %>% filter(realm == "Marine")
terr_cd <- cold %>% filter(realm == "Terrestrial")

new_data <- data.frame(expand_grid(abs_lat_mp = c(median(terr_cd$abs_lat_mp, na.rm=T),
                                                  median(mar_cd$abs_lat_mp, na.rm=T),
                                                  median(int_cd$abs_lat_mp, na.rm=T)),
                                   dispersal_distance_continuous = 
                                     unique(cold$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(terr_cd$log_maximum_body_size,na.rm=T),
                                       median(int_cd$log_maximum_body_size,na.rm=T),
                                       median(mar_cd$log_maximum_body_size,na.rm=T)) ,
                                   log_range_area = seq(min(cold$log_range_area, 
                                                                na.rm = T), 
                                                            max(cold$log_range_area, na.rm = T), 
                                                            length.out = 1000),
                                   realm = c("Terrestrial", 'Marine', 'Intertidal'),
                                   metric = c("ct","lt")))

new_data$realm <- factor(new_data$realm, levels = c("Terrestrial", 'Intertidal', 'Marine'))

pred_cold <- predict(avgm_cold, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_cold <- new_data %>%
  mutate(filling_value = pred_cold$fit) %>%
  mutate(filling_value_SE = pred_cold$se.fit)

## Convert scaled prediction to original data scale of log_range_area
fitted_pred_cold$log_range_area <- fitted_pred_cold$log_range_area * scalers$sd[2] + scalers$mean[2]
cold$log_range_area <- cold$log_range_area * scalers$sd[2] + scalers$mean[2]
terr_cd$log_range_area <- terr_cd$log_range_area * scalers$sd[2] + scalers$mean[2]
int_cd$log_range_area <-int_cd$log_range_area * scalers$sd[2] + scalers$mean[2]
mar_cd$log_range_area <- mar_cd$log_range_area * scalers$sd[2] + scalers$mean[2]

realm_cold_terr <- fitted_pred_cold %>%
  filter(realm == "Terrestrial") %>%
  filter(metric == getmode(terr_cd$metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(terr_cd$dispersal_distance_continuous)) %>%
  filter(log_range_area >= min(terr_cd$log_range_area, na.rm=T) & 
           log_range_area <= max(terr_cd$log_range_area, na.rm=T))%>%
  filter(abs_lat_mp == median(terr_cd$abs_lat_mp,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(terr_cd$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(niche_edge = "Cold")

realm_cold_int <- fitted_pred_cold %>%
  filter(realm == "Intertidal") %>%
  filter(metric == getmode(int_cd$metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(int_cd$dispersal_distance_continuous)) %>%
  filter(log_range_area >= min(int_cd$log_range_area, na.rm=T) & 
           log_range_area <= max(int_cd$log_range_area, na.rm=T)) %>%
  filter(abs_lat_mp == median(int_cd$abs_lat_mp,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(int_cd$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(niche_edge = "Cold") 

realm_cold_mar <- fitted_pred_cold %>%
  filter(realm == "Marine") %>%
  filter(metric == getmode(mar_cd$metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(mar_cd$dispersal_distance_continuous)) %>%
  filter(log_range_area >= min(mar_cd$log_range_area, na.rm=T) & 
           log_range_area <= max(mar_cd$log_range_area, na.rm=T)) %>%
  filter(abs_lat_mp == median(mar_cd$abs_lat_mp,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(mar_cd$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(niche_edge = "Cold")

warm = warmnorm

int_cd <- warm %>% filter(realm == "Intertidal")
mar_cd <- warm %>% filter(realm == "Marine")
terr_cd <- warm %>% filter(realm == "Terrestrial")

new_data <- data.frame(expand_grid(abs_lat_mp = c(median(terr_cd$abs_lat_mp, na.rm=T),
                                                  median(mar_cd$abs_lat_mp, na.rm=T),
                                                  median(int_cd$abs_lat_mp, na.rm=T)),
                                   dispersal_distance_continuous =
                                     unique(warm$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(c(median(terr_cd$log_maximum_body_size,na.rm=T),
                                         median(int_cd$log_maximum_body_size,na.rm=T),
                                         median(mar_cd$log_maximum_body_size,na.rm=T))), 
                                   log_range_area = seq(min(warm$log_range_area, 
                                                                na.rm = T), 
                                                            max(warm$log_range_area, na.rm = T), 
                                                            length.out = 1000),
                                   realm = c("Terrestrial", 'Marine', 'Intertidal'),
                                   metric = c("ct","lt")))

new_data$realm <- factor(new_data$realm, levels = c("Terrestrial", 'Intertidal', 'Marine'))

pred_warm <- predict(avgm_warm, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

## Convert scaled prediction to original data scale of log_range_area
fitted_pred_warm$log_range_area <- fitted_pred_warm$log_range_area * scalers$sd[6] + scalers$mean[6]
warm$log_range_area <- warm$log_range_area * scalers$sd[6] + scalers$mean[6]
terr_cd$log_range_area <- terr_cd$log_range_area * scalers$sd[6] + scalers$mean[6]
int_cd$log_range_area <-int_cd$log_range_area * scalers$sd[6] + scalers$mean[6]
mar_cd$log_range_area <- mar_cd$log_range_area * scalers$sd[6] + scalers$mean[6]

realm_warm_terr <- fitted_pred_warm %>%
  filter(realm == "Terrestrial") %>%
  filter(metric == getmode(terr_cd$metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(terr_cd$dispersal_distance_continuous)) %>%
  filter(log_range_area >= min(terr_cd$log_range_area, na.rm=T) & 
           log_range_area <= max(terr_cd$log_range_area, na.rm=T)) %>%
  filter(abs_lat_mp == median(terr_cd$abs_lat_mp,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(terr_cd$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(niche_edge = "Warm")

realm_warm_mar <- fitted_pred_warm %>%
  filter(realm == "Marine") %>%
  filter(metric == getmode(mar_cd$metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(mar_cd$dispersal_distance_continuous)) %>%
  filter(log_range_area >= min(mar_cd$log_range_area, na.rm=T) & 
           log_range_area <= max(mar_cd$log_range_area, na.rm=T)) %>%
  filter(abs_lat_mp == median(mar_cd$abs_lat_mp,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(mar_cd$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(niche_edge = "Warm")

realm_warm_int <- fitted_pred_warm %>%
  filter(realm == "Intertidal") %>%
  filter(metric == getmode(int_cd$metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(int_cd$dispersal_distance_continuous)) %>%
  filter(log_range_area >= min(int_cd$log_range_area, na.rm=T) & 
           log_range_area <= max(int_cd$log_range_area, na.rm=T)) %>%
  filter(abs_lat_mp == median(int_cd$abs_lat_mp,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(int_cd$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(niche_edge = "Warm")

## plot predictions
warm$cold_season_dormancy_ = NA
cold$hot_season_dormancy_ = NA
warm$niche_edge = "Warm"
cold$niche_edge = "Cold"

both <- rbind(warm, cold)
both$niche_edge <- factor(both$niche_edge, levels = c("Warm", "Cold"), ordered = T)

##~~~~ RANGE AREA ~~~~~~##
realm_cold_no_dormancy <- rbind(realm_cold_mar, realm_cold_int) %>%
  rbind(.,  realm_cold_terr) %>%
  mutate(facet = ifelse(.$realm %in% c("Terrestrial", "Intertidal"), "Air-interacting", "Subtidal"))

realm_warm_no_dormancy <- rbind(realm_warm_mar, realm_warm_int) %>%
  rbind(.,  realm_warm_terr) %>%
  mutate(facet = ifelse(.$realm %in% c("Terrestrial", "Intertidal"), "Air-interacting", "Subtidal"))

## try new colour palette
#c( "#74add1", "#45b0b4","#313695", "#d16280", "#f46d43", "#a50026")
palette <- c( "#74add1", "#45b0b4","#313695", "#f57396", "#f58d42", "#a50026")
linetype <- c("dotdash", "longdash", "solid", "dotdash", "longdash", "solid")
shapes <- c(17,15,19,17,15,19)
labels <- c("Intertidal marine", "Subtidal marine", "Terrestrial", "Intertidal marine", "Subtidal marine", "Terrestrial")

realm_warm_no_dormancy$colour <- factor(paste(realm_warm_no_dormancy$niche_edge,
                                              realm_warm_no_dormancy$realm), 
                                        levels = c("Cold Terrestrial", "Cold Intertidal", "Cold Marine",
                                                   "Warm Terrestrial", "Warm Intertidal", "Warm Marine"), 
                                        ordered = T)
realm_cold_no_dormancy$colour <- factor(paste(realm_cold_no_dormancy$niche_edge, realm_cold_no_dormancy$realm), 
                                        levels = c("Cold Terrestrial", "Cold Intertidal", "Cold Marine",
                                                   "Warm Terrestrial", "Warm Intertidal", "Warm Marine"), 
                                        ordered = T)
both$colour <- factor(paste(both$niche_edge, both$realm), 
                      levels = c("Cold Intertidal", "Cold Marine", "Cold Terrestrial",
                                 "Warm Intertidal", "Warm Marine", "Warm Terrestrial"), 
                      ordered = T)

range_area <- both %>%
  ggplot(., aes(x = log_range_area, y = filling_value,
                shape = colour, colour = colour)) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  geom_point(alpha = 0.6) +
  labs(col = "", shape = "", x =  "Realized range size (log no. cells)", 
       fill = "", linetype = "",
       y =  "Niche filling (°C)") +
  geom_line(data = realm_warm_no_dormancy, inherit.aes = F, size = 2,
            aes(x = log_range_area, y = filling_value, colour = colour,
                group = realm)) +
  geom_line(data = realm_cold_no_dormancy, inherit.aes = F,  size = 2,
            aes(x = log_range_area, y = filling_value, colour = colour,
                group = realm)) +
  scale_color_manual(values = palette, labels = labels) +
  scale_fill_manual(values = palette, labels = labels) +
  geom_ribbon(data = realm_warm_no_dormancy, inherit.aes = F,
              aes(x = log_range_area,
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = realm, fill = colour),
              alpha=0.3) +
  geom_ribbon(data = realm_cold_no_dormancy, inherit.aes = F, 
              aes(x = log_range_area, 
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = realm, fill = colour),
              alpha=0.15) +
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(t = 1, b = 0, r = 1, l = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(), 
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(-23, 36), breaks = c(-20, 0, 20, 40)) +
  scale_x_continuous(limits = c(-0.35, 8)) +
  scale_linetype_manual(values = linetype,
                        labels = labels) +
  scale_shape_manual(values = shapes,
                     labels = labels) +
  facet_wrap(~niche_edge, nrow = 2) 

saveRDS(range_area, "data-processed/intermediate-files/predictions_niche_range-area.rds")

ggsave(range_area, path ="figures/additional-figures/predictions", filename = "predictions_niche_range-area.png",
       width = 5.5, height = 6, device = "png")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####             Across latitude              ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
cold = coldnorm

int_cd <- cold %>% filter(realm == "Intertidal")
mar_cd <- cold %>% filter(realm == "Marine")
terr_cd <- cold %>% filter(realm == "Terrestrial")

new_data <- data.frame(expand_grid(abs_lat_mp = seq(min(cold$abs_lat_mp), max(cold$abs_lat_mp), 
                                                    length.out = 180),
                                   dispersal_distance_continuous = 
                                     unique(cold$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(terr_cd$log_maximum_body_size,
                                              na.rm=T),
                                       median(int_cd$log_maximum_body_size,
                                              na.rm=T),
                                       median(mar_cd$log_maximum_body_size,
                                              na.rm=T)), 
                                   log_range_area = c(median(terr_cd$log_range_area,
                                                             na.rm=T),
                                                      median(int_cd$log_range_area,
                                                             na.rm=T),
                                                      median(mar_cd$log_range_area,
                                                             na.rm=T)), 
                                   realm = c("Terrestrial", 'Marine', 'Intertidal'),
                                   metric = c("ct","lt")))

new_data$realm <- factor(new_data$realm, levels = c("Terrestrial", 'Intertidal', 'Marine'))

pred_cold <- predict(avgm_cold, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_cold <- new_data %>%
  mutate(filling_value = pred_cold$fit) %>%
  mutate(filling_value_SE = pred_cold$se.fit)

## Convert scaled prediction to original data scale of log_range_area
fitted_pred_cold$abs_lat_mp <- fitted_pred_cold$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
cold$abs_lat_mp <- cold$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
terr_cd$abs_lat_mp <- terr_cd$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
int_cd$abs_lat_mp <-int_cd$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
mar_cd$abs_lat_mp <- mar_cd$abs_lat_mp * scalers$sd[1] + scalers$mean[1]

realm_cold_terr <- fitted_pred_cold %>%
  filter(realm == "Terrestrial") %>%
  filter(metric == getmode(metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(terr_cd$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(terr_cd$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(terr_cd$abs_lat_mp, na.rm=T) ) %>%
  filter(log_range_area == median(terr_cd$log_range_area,
                                  na.rm=T)) %>%
  filter(log_maximum_body_size == median(terr_cd$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(niche_edge = "Cold")

realm_cold_mar <- fitted_pred_cold %>%
  filter(realm == "Marine") %>%
  filter(metric == getmode(metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(mar_cd$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(mar_cd$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(mar_cd$abs_lat_mp, na.rm=T) ) %>%
  filter(log_range_area == median(mar_cd$log_range_area,
                                      na.rm=T)) %>%
  filter(log_maximum_body_size == median(mar_cd$log_maximum_body_size,
                                         na.rm=T))  %>%
  mutate(niche_edge = "Cold")

realm_cold_int <- fitted_pred_cold %>%
  filter(realm == "Intertidal") %>%
  filter(metric == getmode(metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(int_cd$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(int_cd$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(int_cd$abs_lat_mp, na.rm=T) ) %>%
  filter(log_range_area == median(int_cd$log_range_area,
                                      na.rm=T)) %>%
  filter(log_maximum_body_size == median(int_cd$log_maximum_body_size,
                                         na.rm=T))  %>%
  mutate(niche_edge = "Cold")

warm = warmnorm

int_cd <- warm %>% filter(realm == "Intertidal")
mar_cd <- warm %>% filter(realm == "Marine")
terr_cd <- warm %>% filter(realm == "Terrestrial")

new_data <- data.frame(expand_grid(abs_lat_mp = seq(min(warm$abs_lat_mp), max(warm$abs_lat_mp), 
                                                    length.out = 180),
                                   dispersal_distance_continuous = 
                                     unique(warm$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(terr_cd$log_maximum_body_size,
                                                                    na.rm=T),
                                                             median(int_cd$log_maximum_body_size,
                                                                    na.rm=T),
                                                             median(mar_cd$log_maximum_body_size,
                                                                    na.rm=T)), 
                                   log_range_area = c(median(terr_cd$log_range_area,
                                                             na.rm=T),
                                                      median(int_cd$log_range_area,
                                                             na.rm=T),
                                                      median(mar_cd$log_range_area,
                                                             na.rm=T)), 
                                   realm = c("Terrestrial", 'Marine', 'Intertidal'),
                                   metric = c("ct","lt")))

new_data$realm <- factor(new_data$realm, levels = c("Terrestrial", 'Intertidal', 'Marine'))

pred_warm <- predict(avgm_warm, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

## Convert scaled prediction to original data scale of log_range_area
fitted_pred_warm$abs_lat_mp <- fitted_pred_warm$abs_lat_mp * scalers$sd[5] + scalers$mean[5]
warm$abs_lat_mp <- warm$abs_lat_mp * scalers$sd[5] + scalers$mean[5]
terr_cd$abs_lat_mp <- terr_cd$abs_lat_mp * scalers$sd[5] + scalers$mean[5]
int_cd$abs_lat_mp <-int_cd$abs_lat_mp * scalers$sd[5] + scalers$mean[5]
mar_cd$abs_lat_mp <- mar_cd$abs_lat_mp * scalers$sd[5] + scalers$mean[5]

realm_warm_terr <- fitted_pred_warm %>%
  filter(realm == "Terrestrial") %>%
  filter(metric == getmode(metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(terr_cd$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(terr_cd$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(terr_cd$abs_lat_mp, na.rm=T)) %>%
  filter(log_range_area == median(terr_cd$log_range_area,
                                      na.rm=T)) %>%
  filter(log_maximum_body_size == median(terr_cd$log_maximum_body_size,
                                      na.rm=T))  %>%
  mutate(niche_edge = "Warm")

realm_warm_mar <- fitted_pred_warm %>%
  filter(realm == "Marine") %>%
  filter(metric == getmode(metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(mar_cd$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(mar_cd$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(mar_cd$abs_lat_mp, na.rm=T) ) %>%
  filter(log_range_area == median(mar_cd$log_range_area,
                                      na.rm=T)) %>%
  filter(log_maximum_body_size == median(mar_cd$log_maximum_body_size,
                                         na.rm=T))  %>%
  mutate(niche_edge = "Warm")

realm_warm_int <- fitted_pred_warm %>%
  filter(realm == "Intertidal") %>%
  filter(metric == getmode(metric)) %>%
  filter(dispersal_distance_continuous ==
           getmode(int_cd$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(int_cd$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(int_cd$abs_lat_mp, na.rm=T) ) %>%
  filter(log_range_area == median(int_cd$log_range_area,
                                      na.rm=T)) %>%
  filter(log_maximum_body_size == median(int_cd$log_maximum_body_size,
                                         na.rm=T))  %>%
  mutate(niche_edge = "Warm")

## plot predictions
warm$cold_season_dormancy_ = NA
cold$hot_season_dormancy_ = NA
warm$niche_edge = "Warm"
cold$niche_edge = "Cold"

both <- rbind(warm, cold)
both$niche_edge <- factor(both$niche_edge, levels = c("Warm", "Cold"), ordered = T)

##~~~~~~ LATITUDE ~~~~~~##
realm_cold_no_dormancy <- rbind(realm_cold_int, realm_cold_mar) %>%
  rbind(., realm_cold_terr) %>%
  mutate(facet = ifelse(.$realm %in% c("Terrestrial", "Intertidal"), "Air-interacting", "Subtidal"))

realm_warm_no_dormancy <- rbind(realm_warm_int, realm_warm_mar) %>%
  rbind(., realm_warm_terr) %>%
  mutate(facet = ifelse(.$realm %in% c("Terrestrial", "Intertidal"), "Air-interacting", "Subtidal"))

## try new colour palette
palette <- c( "#74add1", "#45b0b4","#313695", "#d16280", "#f46d43", "#a50026")
linetype <- c("dotdash", "longdash", "solid", "dotdash", "longdash", "solid")
shapes <- c(17,15,19,17,15,19)
labels <- c("Intertidal marine", "Subtidal marine", "Terrestrial", "Intertidal marine", "Subtidal marine", "Terrestrial")

realm_warm_no_dormancy$colour <- factor(paste(realm_warm_no_dormancy$niche_edge,
                                              realm_warm_no_dormancy$realm), 
                                        levels = c("Cold Terrestrial", "Cold Intertidal", "Cold Marine",
                                                   "Warm Terrestrial", "Warm Intertidal", "Warm Marine"), 
                                        ordered = T)
realm_cold_no_dormancy$colour <- factor(paste(realm_cold_no_dormancy$niche_edge, realm_cold_no_dormancy$realm), 
                                        levels = c("Cold Terrestrial", "Cold Intertidal", "Cold Marine",
                                                   "Warm Terrestrial", "Warm Intertidal", "Warm Marine"), 
                                        ordered = T)
both$colour <- factor(paste(both$niche_edge, both$realm), 
                      levels = c("Cold Intertidal", "Cold Marine", "Cold Terrestrial",
                                 "Warm Intertidal", "Warm Marine", "Warm Terrestrial"), 
                      ordered = T)

latitude <- both %>%
  ggplot(., aes(x = abs_lat_mp, y = filling_value,
                 colour = colour,
                #shape = colour
                )) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  geom_point(alpha = 0.6) +
  labs(col = "", shape = "",
       x =  "Absolute realized range latitudinal midpoint (°N/S)", 
       fill = "", linetype = "",
       y =  "Niche filling (°C)") +
  geom_line(data = realm_warm_no_dormancy, inherit.aes = F, size = 1,
            aes(x = abs_lat_mp, y = filling_value, colour = colour,
                group = realm)) +
  geom_line(data = realm_cold_no_dormancy, inherit.aes = F, size = 1,
            aes(x = abs_lat_mp, y = filling_value, colour = colour,
                group = realm)) +
  scale_color_manual(values = palette,
                     labels = labels)+
  geom_ribbon(data = realm_warm_no_dormancy, inherit.aes = F, 
              aes(x = abs_lat_mp, 
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = realm, fill = colour),
              alpha=0.3) +
  geom_ribbon(data = realm_cold_no_dormancy, inherit.aes = F, 
              aes(x = abs_lat_mp, 
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = realm, fill = colour),
              alpha=0.3) +
  scale_fill_manual(values = palette, 
                    labels = labels) +
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(t = 1, b = 0, r = 1, l = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(-23, 36), breaks = c(-20, 0, 20, 40)) +  
  scale_x_continuous(limits = c(-2.5, 66)) +
  scale_linetype_manual(values = linetype, 
                        labels = labels) +
  scale_shape_manual(values = shapes,
                     labels = labels) + facet_wrap(niche_edge~realm) 

latitude <- latitude + theme(legend.position = "none")

saveRDS(latitude, "data-processed/intermediate-files/predictions_niche_latitude.rds")

ggsave(latitude, path ="figures/main", filename = "predictions_niche_latitude.png",
       width = 6.5, height = 4, device = "png")

latitude_talk <- both %>%
  ggplot(., aes(x = abs_lat_mp, y = filling_value,
                shape = colour, colour = colour)) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  geom_point(alpha = 0.6) +
  labs(col = "", shape = "",
       x =  "Absolute realized range latitudinal midpoint (°N/S)", 
       fill = "", linetype = "",
       y =  "Niche filling (°C)") +
  geom_line(data = realm_warm_no_dormancy, inherit.aes = F, size = 1,
            aes(x = abs_lat_mp, y = filling_value, colour = colour,
                group = realm)) +
  geom_line(data = realm_cold_no_dormancy, inherit.aes = F, size = 1,
            aes(x = abs_lat_mp, y = filling_value, colour = colour,
                group = realm)) +
  scale_color_manual(values = palette,
                     labels = labels)+
  geom_ribbon(data = realm_warm_no_dormancy, inherit.aes = F, 
              aes(x = abs_lat_mp, 
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = realm, fill = colour),
              alpha=0.3) +
  geom_ribbon(data = realm_cold_no_dormancy, inherit.aes = F, 
              aes(x = abs_lat_mp, 
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = realm, fill = colour),
              alpha=0.3) +
  scale_fill_manual(values = palette, 
                    labels = labels) +
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(t = 1, b = 0, r = 1, l = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(-26, 36), breaks = c(-20, 0, 20, 40)) +  
  scale_x_continuous(limits = c(-2.5, 66)) +
  scale_linetype_manual(values = linetype, 
                        labels = labels) +
  scale_shape_manual(values = shapes,
                     labels = labels) + facet_wrap(niche_edge~realm) 
latitude_talk <- latitude_talk + theme(legend.position = "none")

ggsave(latitude_talk, path ="figures/additional-figures/for-talk", filename = "predictions_niche_latitude_talk.png",
       width = 6, height = 4, device = "png")


latitude_blank <- both %>%
  ggplot(., aes(x = abs_lat_mp, y = filling_value,
                shape = colour, colour = colour)) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  labs(col = "", shape = "",
       x =  "Absolute realized range latitudinal midpoint (°N/S)", 
       fill = "", linetype = "",
       y =  "Niche filling (°C)") +
  scale_color_manual(values = palette,
                     labels = labels)+
  scale_fill_manual(values = palette, 
                    labels = labels) +
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(t = 1, b = 0, r = 1, l = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(-26, 36), breaks = c(-20, 0, 20, 40)) +  
  scale_x_continuous(limits = c(-2.5, 66)) +
  scale_linetype_manual(values = linetype, 
                        labels = labels) +
  scale_shape_manual(values = shapes,
                     labels = labels) + facet_wrap(niche_edge~realm) 
latitude_blank <- latitude_blank + theme(legend.position = "none")

ggsave(latitude_blank, path ="figures/additional-figures/for-talk", filename = "predictions_niche_latitude_talk_blank.png",
       width = 6, height = 4, device = "png")


## save legend:
latitude <- readRDS("data-processed/intermediate-files/predictions_niche_latitude.rds")
latitude <- latitude + theme(legend.position = "right")

legend = get_legend(latitude)

ggsave(legend, path ="figures/main", filename = "predictions_niche_legend.png",
       width = 2, height = 2, device = "png")



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####                Realm distribution plots           ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
warm = warmnorm

new_data <- data.frame(expand_grid(abs_lat_mp = 0,
                                   dispersal_distance_continuous = 0,
                                   log_maximum_body_size = 0, 
                                   log_range_area = 0,
                                   realm = c("Terrestrial", 'Marine', 'Intertidal'),
                                   metric = c("ct","lt")))

new_data$realm <- factor(new_data$realm, levels = c("Terrestrial", 'Intertidal', 'Marine'))

pred_warm <- predict(avgm_warm, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

## get mean and 95 CI of realms, across all levels of discrete variables and holding lat/body size/range area constant
means_cold <- fitted_pred_cold %>%
  mutate(CI = 1.96*filling_value_SE) %>%
  group_by(realm) %>%
  mutate(mean = mean(filling_value, na.rm = T),
         mean_CI = mean(CI)) %>%
  ungroup() %>%
  select(realm, mean, mean_CI) %>%
  unique(.) %>%
  mutate(edge_type = "Cold niche edge")

means_warm <- fitted_pred_warm %>%
  mutate(CI = 1.96*filling_value_SE) %>%
  group_by(realm) %>%
  mutate(mean = mean(filling_value, na.rm = T),
         mean_CI = mean(CI)) %>%
  ungroup() %>%
  select(realm, mean, mean_CI) %>%
  unique(.) %>%
  mutate(edge_type = "Warm niche edge")

cold$edge_type = "Cold niche edge"
warm$edge_type = "Warm niche edge"
warm$cold_season_dormancy_ = NA
cold$hot_season_dormancy_ = NA
warm$niche_edge = "Warm"
both <- rbind(cold, warm)


palette <- c( "#74add1","#f57396","#45b0b4", "#f58d42","#313695","#a50026")
both$colour <- paste(both$realm, both$edge_type)

realm_distributions <- ggplot(both, aes(x = filling_value, y = realm, 
                            fill = colour, colour = colour)) +
  geom_density_ridges(aes(height = ..density..), stat = "density", scale = 0.9,
                      alpha = 0.6, trim = TRUE)  +
  geom_point(shape = 95, size = 4) +
  theme_ridges() + 
  theme(legend.position = "none") +
  coord_flip() +
  scale_fill_manual(values = palette) +
  scale_colour_manual(values = palette)+
  labs(y = "", 
       x = "", #x = "Niche filling (°C)"
  ) +
  scale_y_discrete(labels = c("Terrestrial", "Intertidal\nmarine",
                              "Subtidal\nmarine")) +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10, hjust = 0.5),
        axis.text.y = element_text(size = 8),
        panel.border = element_rect(colour = "black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())   +
  annotate(geom = "text", label = "overfilling", y = 0.5, x = 12, angle = 90) +
  annotate(geom = "text", label = "underfilling", y = 0.5, x = -14, angle = 90) +
  geom_vline(xintercept = 0, size = 0.4) +
  scale_x_continuous(breaks = c(-20, -10, 0, 10, 20, 30), 
                     labels = c("-20°C", "-10°C", 
                                "0°C", "10°C",
                                "20°C", "30°C"),
                     limits = c(-25, 40)) +
  facet_wrap(~edge_type) 

ggsave(realm_distributions, "figures/main", 
       filename = "distributions_niche_realm.png",
       width = 5, height = 3, device = "png")

## save dataset 
full_withdormancy <- saveRDS(both, "data-processed/intermediate-files/full_withdormancy.rds") 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####             Plot distributions              #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## join acclimation complete case data subset to whole complete case data
acc_cc_warm <- readRDS("data-processed/intermediate-files/acc-warm-complete-cases.rds")
acc_cc_cold <- readRDS("data-processed/intermediate-files/acc-cold-complete-cases.rds")

acc_cc_warm$edge_type = "Warm niche edge"
acc_cc_cold$edge_type = "Cold niche edge"
acc_cc_warm$cold_season_dormancy_ = NA
acc_cc_cold$hot_season_dormancy_ = NA
acc_cc_warm$colour <- "acclimatized"
acc_cc_cold$colour <- "acclimatized"

both$colour = both$edge_type
acc_cc_warm$niche_edge = "Warm"
acc_cc_cold$niche_edge = "Cold"

acc_both <- rbind(both, acc_cc_warm) %>%
  rbind(., acc_cc_cold)

## plot em all together
ghostly <- ggplot(acc_both, aes(x = filling_value, y = realm, 
                                fill = colour, colour = colour)) +
  geom_density_ridges(aes(height = ..density..), stat = "density", scale = 0.9,
                      alpha = 0.6, trim = TRUE)  +
  geom_point(data = both, shape = 95, size = 4) +
  theme_ridges() +
  theme(legend.position = "none") +
  coord_flip() +
  scale_fill_manual(values = c("darkgrey", "steelblue", "#b45346")) +
  scale_colour_manual(values =c("darkgrey", "steelblue", "#b45346")) +
  labs(y = "", 
       x = "", #x = "Niche filling (°C)"
  ) +
  scale_y_discrete(labels = c("Terrestrial", "Marine\n(intertidal)",
                              "Marine\n(subtidal)")) +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10, hjust = 0.5),
        axis.text.y = element_text(size = 8),
        panel.border = element_rect(colour = "black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())   +
  annotate(geom = "text", label = "overfilling", y = 0.5, x = 12, angle = 90) +
  annotate(geom = "text", label = "underfilling", y = 0.5, x = -14, angle = 90) +
  geom_vline(xintercept = 0, size = 0.4) +
  scale_x_continuous(breaks = c(-20, -10, 0, 10, 20, 30), 
                     labels = c("-20°C", "-10°C", 
                                "0°C", "10°C",
                                "20°C", "30°C"),
                     limits = c(-25, 40)) +
  facet_wrap(~edge_type) 

saveRDS(ghostly, "data-processed/intermediate-files/niche_realm-with-ghostly-acc.rds")

ggsave(ghostly, "figures/main", filename = "distributions_niche_realm-with-acc.png",
       width = 5, height = 3, device = "png")



## plot only ghosts and overlay in inkscape
acc_no_int <- acc_both %>%
  mutate(filling_value = ifelse(niche_edge == "Cold" & realm == "Intertidal" & colour == "acclimatized",
                                NA,
                                filling_value))

acc_int <- acc_both %>%
  mutate(filling_value = ifelse(!(niche_edge == "Cold" & realm == "Intertidal" & colour == "acclimatized"),
                                NA,
                                filling_value))
  
ghostly_only <- acc_both %>%
  filter(colour == "acclimatized") %>%
  ggplot(., aes(x = filling_value, y = realm, 
                fill = colour, colour = colour)) +
  geom_density_ridges(aes(height = ..density..), stat = "density", scale = 0.9,
                      alpha = 0.6, trim = TRUE) + 
  theme_ridges() +
  theme(legend.position = "none") +
  coord_flip() +
  scale_fill_manual(values = c("darkgrey"), labels = c('Acclimatized')) +
  scale_colour_manual(values =c("darkgrey"), labels = c('Acclimatized')) +
  labs(y = "", 
       x = "", #x = "Niche filling (°C)"
  ) +
  scale_y_discrete(labels = c("Terrestrial", "Marine\n(intertidal)",
                              "Marine\n(subtidal)")) +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10, hjust = 0.5),
        axis.text.y = element_text(size = 8),
        panel.border = element_rect(colour = "black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())   +
  annotate(geom = "text", label = "overfilling", y = 0.5, x = 12, angle = 90,
           colour = "transparent") +
  annotate(geom = "text", label = "underfilling", y = 0.5, x = -14, angle = 90,
           colour = "transparent") +
  scale_x_continuous(breaks = c(-20, -10, 0, 10, 20, 30), 
                     labels = c("-20°C", "-10°C", 
                                "0°C", "10°C",
                                "20°C", "30°C"),
                     limits = c(-25, 40)) +
  facet_wrap(~edge_type)

ggsave(ghostly_only, "figures/main", filename = "distributions_niche_realm-acc-only_both.png",
       width = 5, height = 3, device = "png")

## save legend of ghost:
ghost_leg <- ghostly_only + theme(legend.position = 'left') + labs(colour = "",
                                                                   fill = "") 
ghost_leg <- get_legend(ghost_leg)

ggsave(ghost_leg, "figures/main", filename = "distributions_niche_realm-acc-legend.png",
       width = 2, height = 1, device = "png")
