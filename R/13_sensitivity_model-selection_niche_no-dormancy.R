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
rename <- dplyr::rename

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
  filter(type == "min") %>%  
  filter(cold_season_dormancy_ == "No",hot_season_dormancy_ == "No") ## GET RID OF DORMANT SPECIES
## assign metric type used to measure cold limit 

warm <- Te %>%
  filter(edge_type == "warm") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "max") %>%
  filter(hot_season_dormancy_ == "No", cold_season_dormancy_ == "No") ## GET RID OF DORMANT SPECIES
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
                    log_maximum_body_size  + metric,
                  data = cold)
vif(modvif_cold)
## removed:  dispersal_ability_category

modvif_warm <- lm(filling_value ~ abs_lat_mp + log_range_area + realm + dispersal_ability_category +
                    dispersal_distance_continuous +
                    log_maximum_body_size + metric,
                  data = warm)
vif(modvif_warm)
## removed: dispersal_ability_category 

#Second, rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
cold <- select(cold, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size,metric,
                       Class, Order, Family, Genus, Species))
warm <- select(warm, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size,metric, 
                       Class, Order, Family, Genus, Species))

# get complete cases
cold <- subset(cold, complete.cases(cold))
# dim(cold)
warm <- subset(warm, complete.cases(warm))
# dim(warm)

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
                    log_maximum_body_size  + metric,
                  
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
confset.95p_cold <- get.models(allmodels_cold, subset = cumsum(weight)<=.95) #get confidence set, 12 models
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
  labs(y = "Effect of variable on cold niche filling", x = "") + 
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal",
                              "Metric: lethal",
                              "Realized range size (log no. cells)", 
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Abs. realized range latitudinal midpoint ", 
                              "Reference"), drop = FALSE)


# plotting estimates (fixed effects) 
suppressWarnings(allmodels_warm <- dredge(model_warm, extra="R^2")) 
confset.95p_warm <- get.models(allmodels_warm, subset = cumsum(weight)<=.95) #get confidence set
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
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Abs. realized range latitudinal midpoint ", 
                              "Reference"))

## combine and save figure:
whisker <- ggdraw() + 
  draw_plot(whisker_warm, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(whisker_cold, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(label = c("a)", "b)"),
                  x = c(0, 0.5),
                  y = c(1, 1), size = 10, 
                  color = "grey30") 

ggsave(whisker, width = 10, height = 3.5, path = "figures/additional-figures/dormancy/", 
       filename = "whisker-plot_niche_no-dormants.png", 
       device = "png")

saveRDS(whisker_warm, "data-processed/intermediate-files/whisker_warm_no-dormants.rds")
saveRDS(whisker_cold, "data-processed/intermediate-files/whisker_cold_no-dormants.rds")

ggsave(whisker_warm, path ="figures/additional-figures/dormancy/", 
       filename = "niche_model-results-warm_no-dormants.png",
       width = 5.5, height = 3, device = "png")
ggsave(whisker_cold, path ="figures//additional-figures/dormancy/", 
       filename = "niche_model-results-cold_no-dormants.png",
       width = 5.5, height = 3, device = "png")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####          Plotting main results           ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# function to find mode of each predictor
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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
warm$niche_edge = "Warm"
cold$niche_edge = "Cold"

both <- rbind(warm, cold)
both$niche_edge <- factor(both$niche_edge, levels = c("Warm", "Cold"), ordered = T)


##~~~~~~ LATITUDE ~~~~~~##
realm_cold_no_dormancy <- rbind(realm_cold_int, realm_cold_mar) %>%
  rbind(., realm_cold_terr) 

realm_warm_no_dormancy <- rbind(realm_warm_int, realm_warm_mar) %>%
  rbind(., realm_warm_terr)  

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
                shape = colour, colour = colour)) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  geom_point(alpha = 0.6) +
  labs(col = "", shape = "",
       x =  "Realized range latitudinal midpoint (°N/S)", 
       fill = "", linetype = "",
       y =  "Niche filling (°C)") +
  geom_line(data = realm_warm_no_dormancy, inherit.aes = F, size = 2,
            aes(x = abs_lat_mp, y = filling_value, colour = colour,
                group = realm)) +
  geom_line(data = realm_cold_no_dormancy, inherit.aes = F, size = 2,
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
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(-29, 40), breaks = c(-20, 0, 20, 40)) +  
  scale_x_continuous(limits = c(-2.5, 60)) +
  scale_linetype_manual(values = linetype, 
                        labels = labels) +
  scale_shape_manual(values = shapes,
                     labels = labels) +
  facet_wrap(~niche_edge)

saveRDS(latitude, "data-processed/intermediate-files/predictions_latitude_no-dormants.rds")

ggsave(latitude, path ="figures/additional-figures/dormancy/", filename = "predictions_niche_latitude_no-dormants.png",
       width = 8, height = 3.5, device = "png")

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

cold$edge_type = "Cold niche edge"
warm$edge_type = "Warm niche edge"
warm$niche_edge = "Warm"
cold$niche_edge = "Cold"
both <- rbind(cold, warm)

palette <- c( "#74add1","#d16280", "#45b0b4", "#f46d43","#313695","#a50026")
both$colour <- paste(both$realm, both$edge_type)

## add ghosts of full data other distributions 
full_withdormancy <- readRDS("data-processed/intermediate-files/full_withdormancy.rds") # read data

## join data
full_withdormancy$colour = paste(full_withdormancy$realm, full_withdormancy$edge_type)
both <- select(both, -niche_edge)
full_withdormancy <- select(full_withdormancy, -cold_season_dormancy_, -hot_season_dormancy_, -niche_edge)
both$colour = paste("no dormancy")

dorm_both <- rbind(both, full_withdormancy) 

realm_distributions <- ggplot(full_withdormancy, aes(x = filling_value, y = realm, 
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

ggsave(realm_distributions, "figures/extended-data", 
       filename = "distributions_niche_no-dormants.png",
       width = 5, height = 3, device = "png")

full_withdormancy$colour = paste("full data", full_withdormancy$edge_type)

## plot only ghosts and overlay in inkscape
ghostly_only <- both %>%
  filter(colour == "no dormancy") %>%
  ggplot(., aes(x = filling_value, y = realm, 
                                fill = colour, colour = colour)) +
  geom_density_ridges(aes(height = ..density..), stat = "density", scale = 0.9,
                      alpha = 0.6, trim = TRUE) +
  geom_point(data = both, shape = 95, size = 4) +
  theme_ridges() +
  theme(legend.position = "none") +
  coord_flip() +
  scale_fill_manual(values = c("darkgrey"), labels = c('No dormancy')) +
  scale_colour_manual(values =c("darkgrey"), labels = c('No dormancy')) +
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

ggsave(ghostly_only, "figures/extended-data", filename = "distributions_niche_no-dormants_grey.png",
       width = 5, height = 3, device = "png")

## save legend of ghost:
ghost_leg <- ghostly_only + theme(legend.position = 'left') + labs(colour = "",
                                                                   fill = "") 
ghost_leg <- get_legend(ghost_leg)

ggsave(ghost_leg, "figures/extended-data", filename = "distributions_niche_no-dormants_legend.png",
       width = 2, height = 1, device = "png")