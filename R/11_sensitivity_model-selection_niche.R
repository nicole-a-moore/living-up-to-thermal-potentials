## model selection and averaging for sensitivity analysis subsets
library(tidyverse)
library(nlme)
library(car)
library(MuMIn)
library(knitr)
library(gridExtra)
library(cowplot)
library(ggridges)
select <- dplyr::select
rename <- dplyr::rename

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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Tpref/Tb subset                   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Subset WITH adjustment ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
warm <- Te_tpreftb %>%
  filter(edge_type == "warm") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "max") 
## assign metric type used to measure warm limit 

#Rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
warm <- select(warm, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size, hot_season_dormancy_,metric, 
                       Class, Order, Family, Genus, Species))

# get complete cases
warm <- subset(warm, complete.cases(warm))
# dim(warm)

#scale the continuous variables by subtracting mean and dividing by sd
means <- sds <- c()
warm = as.data.frame(warm)
warmnorm <- warm
x = 1
for (i in c(2,3,5,6)) {
  warmnorm[,i] <- (warm[,i]-mean(c(warm[,i])))/sd(warm[,i])
  means[x] = mean(c(warm[,i]))
  sds[x] = sd(c(warm[,i]))
  x = x+1
}
scalers = data.frame(type = "warm", var = colnames(warm)[c(2,3,5,6)], means = means, sds = sds)


warm = warmnorm

#Now, run the full models
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
model_warm <- lme(filling_value ~ abs_lat_mp + log_range_area +
                    dispersal_distance_continuous  +
                    log_maximum_body_size + metric,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = warm
)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_warm <- dredge(model_warm, extra="R^2")) 
confset.95p_warm <- get.models(allmodels_warm, subset = cumsum(weight)<=.95) #get confidence set
avgm_warm <- model.avg(confset.95p_warm) #do averaging
summary(avgm_warm)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Exporting tables                  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_warm)
sum_warm <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", "slope", 
              "intercept")
coefs <- sum_warm[c(1,2,4,6,3,5),1] # reorder
f_effects <- c("(Intercept)",
               "Abs. realized range latitudinal midpoint",
               "Realized range size (log no. cells)",
               "Dispersal distance (km)",
               "Maximum body size (log cm)",
               "Metric: lethal")
std_err <- sum_warm[c(1,2,4,6,3,5),2]  # reorder fixed effects
z_val <- sum_warm[c(1,2,4,6,3,5),4]
p_val <- sum_warm[c(1,2,4,6,3,5),5]
names(coefs) <- names(std_err) <- names(z_val) <- names(p_val) <- NULL

## put all into a table:
results_w <- data.frame("fixed effects" = f_effects, 
                        "effect type" = eff_type,
                        "estimate" = coefs,
                        "s.e." = std_err,
                        "z-value" = z_val,
                        "p-value" = p_val)
colnames(results_w) = c("fixed effects",	"effect type",	"estimate", "s.e.",'z-value',	"p-value")

results_w[c(c(3,4,5))] <- round(results_w[,c(3,4,5)], digits = 2) # round estimates to 2 decimal places
results_w[6] <- round(results_w[,6], digits = 3) # round p values to 3 decimal places
## add significance indicators
results_w$`p-value` <- ifelse(results_w$`p-value` == 0, "<0.001 **", 
                            ifelse(results_w$`p-value` <= 0.01, paste(results_w$`p-value`, " **", sep= ""), 
                                   ifelse(results_w$`p-value` <= 0.05, paste(results_w$`p-value`, " *", sep= ""),
                                          ifelse(results_w$`p-value` <= 0.1, paste(results_w$`p-value`, " ", sep= ""),
                                                 as.character(results_w$`p-value`)))))


table <- results_w %>%
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
  htmlTable(rnames = c("Warm niche edge", rep("", 10), 
                       "Cold niche edge", rep("", 10)))

## save model results for plotting
avgm_warm_behav <- avgm_warm
warm_behav <- warm

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Subset WITHOUT adjustment ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## fitting to cold and warm filling separately
warm <- Te_subset_tpreftb %>%
  filter(edge_type == "warm") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "max") 
## assign metric type used to measure warm limit 


#Rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
warm <- select(warm, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size, hot_season_dormancy_,metric, 
                       Class, Order, Family, Genus, Species))

# get complete cases
warm <- subset(warm, complete.cases(warm))
# dim(warm)

#scale the continuous variables by subtracting mean and dividing by sd
means <- sds <- c()
warm = as.data.frame(warm)
warmnorm <- warm
x = 1
for (i in c(2,3,5,6)) {
  warmnorm[,i] <- (warm[,i]-mean(c(warm[,i])))/sd(warm[,i])
  means[x] = mean(c(warm[,i]))
  sds[x] = sd(c(warm[,i]))
  x = x+1
}
scalers = data.frame(type = "warm", var = colnames(warm)[c(2,3,5,6)], means = means, sds = sds)


warm = warmnorm

#Now, run the full models
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
model_warm <- lme(filling_value ~ abs_lat_mp + log_range_area +
                    dispersal_distance_continuous  +
                    log_maximum_body_size + metric,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = warm
)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_warm <- dredge(model_warm, extra="R^2")) 
confset.95p_warm <- get.models(allmodels_warm, subset = cumsum(weight)<=.95) #get confidence set
avgm_warm <- model.avg(confset.95p_warm) #do averaging
summary(avgm_warm)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Exporting tables                  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_warm)
sum_warm <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", "slope", 
              "intercept")
coefs <- sum_warm[c(1,2,3,6,5,4),1] # reorder
f_effects <- c("(Intercept)",
               "Abs. realized range latitudinal midpoint",
               "Realized range size (log no. cells)",
               "Dispersal distance (km)",
               "Maximum body size (log cm)",
               "Metric: lethal")
std_err <- sum_warm[c(1,2,3,6,5,4),2]  # reorder fixed effects
z_val <- sum_warm[c(1,2,3,6,5,4),4]
p_val <- sum_warm[c(1,2,3,6,5,4),5]
names(coefs) <- names(std_err) <- names(z_val) <- names(p_val) <- NULL

## put all into a table:
results_w <- data.frame("fixed effects" = f_effects, 
                        "effect type" = eff_type,
                        "estimate" = coefs,
                        "s.e." = std_err,
                        "z-value" = z_val,
                        "p-value" = p_val)
colnames(results_w) = c("fixed effects",	"effect type",	"estimate", "s.e.",'z-value',	"p-value")

results_w[c(c(3,4,5))] <- round(results_w[,c(3,4,5)], digits = 2) # round estimates to 2 decimal places
results_w[6] <- round(results_w[,6], digits = 3) # round p values to 3 decimal places
## add significance indicators
results_w$`p-value` <- ifelse(results_w$`p-value` == 0, "<0.001 **", 
                              ifelse(results_w$`p-value` <= 0.01, paste(results_w$`p-value`, " **", sep= ""), 
                                     ifelse(results_w$`p-value` <= 0.05, paste(results_w$`p-value`, " *", sep= ""),
                                            ifelse(results_w$`p-value` <= 0.1, paste(results_w$`p-value`, " ", sep= ""),
                                                   as.character(results_w$`p-value`)))))

table <- results_w %>%
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
  htmlTable(rnames = c("Warm niche edge", rep("", 10), 
                       "Cold niche edge", rep("", 10)))

## save model results for plotting
avgm_warm_no_behav <- avgm_warm
warm_no_behav <- warm

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####   Plotting predictions with vs. without behaviour     ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# function to find mode of each predictor
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
## across range area
warm = warm_behav
new_data <- data.frame(expand_grid(abs_lat_mp = c(median(warm$abs_lat_mp, na.rm=T)),
                                   dispersal_distance_continuous = 
                                     median(warm$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(warm$log_maximum_body_size,na.rm=T)) ,
                                   log_range_area = seq(min(warm$log_range_area, 
                                                            na.rm = T), 
                                                        max(warm$log_range_area, na.rm = T), 
                                                        length.out = 1000),
                                   metric = c("ct","lt")))

pred_warm <- predict(avgm_warm_behav, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm_behav <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

## Convert scaled prediction to original data scale 
fitted_pred_warm_behav$log_range_area <- fitted_pred_warm_behav$log_range_area * scalers$sd[2] + scalers$mean[2]
warm_behav$log_range_area <- warm_behav$log_range_area * scalers$sd[2] + scalers$mean[2]

realm_warm_behav <- fitted_pred_warm_behav %>%
  filter(metric == getmode(warm$metric)) %>%
  filter(dispersal_distance_continuous ==
           median(warm$dispersal_distance_continuous)) %>%
  filter(log_range_area >= min(warm_behav$log_range_area, na.rm=T) & 
           log_range_area <= max(warm_behav$log_range_area, na.rm=T))%>%
  filter(abs_lat_mp == median(warm$abs_lat_mp,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(warm$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(behaviour = "no") 

warm = warm_no_behav
new_data <- data.frame(expand_grid(abs_lat_mp = c(median(warm$abs_lat_mp, na.rm=T)),
                                   dispersal_distance_continuous = 
                                     median(warm$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(warm$log_maximum_body_size,na.rm=T)) ,
                                   log_range_area = seq(min(warm$log_range_area, 
                                                            na.rm = T), 
                                                        max(warm$log_range_area, na.rm = T), 
                                                        length.out = 1000),
                                   metric = c("ct","lt")))

pred_warm <- predict(avgm_warm_no_behav, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm_no_behav <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

## Convert scaled prediction to original data scale 
fitted_pred_warm_no_behav$log_range_area <- fitted_pred_warm_no_behav$log_range_area * scalers$sd[2] + scalers$mean[2]
warm_no_behav$log_range_area <- warm_no_behav$log_range_area * scalers$sd[2] + scalers$mean[2]

realm_warm_no_behav <- fitted_pred_warm_no_behav %>%
  filter(metric == getmode(warm$metric)) %>%
  filter(dispersal_distance_continuous ==
           median(warm$dispersal_distance_continuous)) %>%
  filter(log_range_area >= min(warm_no_behav$log_range_area, na.rm=T) & 
           log_range_area <= max(warm_no_behav$log_range_area, na.rm=T))%>%
  filter(abs_lat_mp == median(warm$abs_lat_mp,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(warm$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(behaviour = "no") 

### plot 
warm_behav$behaviour = "yes"
warm_no_behav$behaviour = "no"
warm_both <- rbind(warm_no_behav, warm_behav)

realm_warm_behav$behaviour = "yes"
realm_warm_no_behav$behaviour = "no"
realm_warm_both <- rbind(realm_warm_no_behav, realm_warm_behav)

range_area <- warm_both %>%
  ggplot(., aes(x = log_range_area, 
                y = filling_value, colour = behaviour, 
                #shape = behaviour, linetype = behaviour
                )) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  geom_point() +
  labs(col = "", shape = "", fill = "", x =  "Realized range size (log no. cells)", 
       y =  "Warm niche filling (°C)",
       linetype = "") +
  geom_line(data = realm_warm_both, inherit.aes = F, size = 2,
            aes(x = log_range_area, y = filling_value, colour = behaviour,
                group = behaviour, 
                #linetype = behaviour
                )) +
  scale_color_manual(values = c("black", "darkgrey"), 
                     labels = c("No behavioural\nthermoregulation", 
                                "Behavioural\nthermoregulation")) +
  geom_ribbon(data = realm_warm_both, inherit.aes = F, 
              aes(x = log_range_area, 
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = behaviour, fill = behaviour),                                                                                                                                                                                                                                                
              alpha=0.3) +
  scale_fill_manual(values = c("black", "darkgrey"),
                    labels =c("No behavioural\nthermoregulation", 
                              "Behavioural\nthermoregulation")) +
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(t = 1, b = 0, r = 1, l = 1),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  # scale_linetype_discrete(labels =c("No behavioural\nthermoregulation", 
  #                                   "Behavioural\nthermoregulation")) +
  # scale_shape_discrete(labels = c("No behavioural\nthermoregulation", 
  #                                 "Behavioural\nthermoregulation"))  +
  scale_y_continuous(limits = c(-29, 40), breaks = c(-20, 0, 20, 40)) +
  scale_x_continuous(limits = c(-0.35, 8))

saveRDS(range_area, "data-processed/intermediate-files/niche_behav-thermo_range-area.rds")

## across latitude
warm = warm_behav

new_data <- data.frame(expand_grid(abs_lat_mp = seq(min(warm$abs_lat_mp, 
                                                        na.rm = T), 
                                                    max(warm$abs_lat_mp, na.rm = T), 
                                                    length.out = 1000),
                                   dispersal_distance_continuous = 
                                     median(warm$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(warm$log_maximum_body_size,na.rm=T)) ,
                                   log_range_area = c(median(warmnorm$log_range_area, na.rm=T)),
                                   metric = c("ct","lt")))

pred_warm <- predict(avgm_warm_behav, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm_behav <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

## Convert scaled prediction to original data scale 
fitted_pred_warm_behav$abs_lat_mp <- fitted_pred_warm_behav$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
warm_behav$abs_lat_mp <- warm_behav$abs_lat_mp * scalers$sd[1] + scalers$mean[1]

realm_warm_behav <- fitted_pred_warm_behav %>%
  filter(metric == getmode(warm$metric)) %>%
  filter(dispersal_distance_continuous ==
           median(warm$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(warm_behav$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(warm_behav$abs_lat_mp, na.rm=T))%>%
  filter(log_range_area == median(warmnorm$log_range_area,
                              na.rm=T)) %>%
  filter(log_maximum_body_size == median(warm$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(behaviour = "no") 

warm = warm_no_behav
new_data <- data.frame(expand_grid(abs_lat_mp = seq(min(warm$abs_lat_mp, 
                                                        na.rm = T), 
                                                    max(warm$abs_lat_mp, na.rm = T), 
                                                    length.out = 1000),
                                   dispersal_distance_continuous = 
                                     median(warm$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(warm$log_maximum_body_size,na.rm=T)) ,
                                   log_range_area = c(median(warmnorm$log_range_area, na.rm=T)),
                                   metric = c("ct","lt")))

pred_warm <- predict(avgm_warm_no_behav, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm_no_behav <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

## Convert scaled prediction to original data scale 
fitted_pred_warm_no_behav$abs_lat_mp <- fitted_pred_warm_no_behav$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
warm_no_behav$abs_lat_mp <- warm_no_behav$abs_lat_mp * scalers$sd[1] + scalers$mean[1]

realm_warm_no_behav <- fitted_pred_warm_no_behav %>%
  filter(metric == getmode(warm$metric)) %>%
  filter(dispersal_distance_continuous ==
           median(warm$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(warm_no_behav$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(warm_no_behav$abs_lat_mp, na.rm=T))%>%
  filter(log_range_area == median(warmnorm$log_range_area,
                                  na.rm=T)) %>%
  filter(log_maximum_body_size == median(warm$log_maximum_body_size,
                                         na.rm=T)) %>%
  mutate(behaviour = "no")

### plot 
warm_behav$behaviour = "yes"
warm_no_behav$behaviour = "no"
warm_both <- rbind(warm_no_behav, warm_behav)

realm_warm_behav$behaviour = "yes"
realm_warm_no_behav$behaviour = "no"
realm_warm_both <- rbind(realm_warm_no_behav, realm_warm_behav)


latitude = warm_both %>%
  ggplot(., aes(x = abs_lat_mp,
                y = filling_value, colour = behaviour, 
               # shape = behaviour, linetype = behaviour
                )) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  geom_point() +
  labs(col = "", fill = "", shape = "", x =  "Absolute realized range latitudinal midpoint (°N/S)", 
       y =  "Warm niche filling (°C)",
       linetype = "") +
  geom_line(data = realm_warm_both, inherit.aes = F, size = 2,
            aes(x = abs_lat_mp, y = filling_value, colour = behaviour,
                group = behaviour, 
                #linetype = behaviour
                )) +
  scale_color_manual(values = c("black", "darkgrey"), 
                     labels = c("No behavioural\nthermoregulation", 
                                "Behavioural\nthermoregulation")) +
  geom_ribbon(data = realm_warm_both, inherit.aes = F, 
              aes(x = abs_lat_mp, 
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = behaviour, fill = behaviour),
              alpha=0.3) +
  scale_fill_manual(values = c("black", "darkgrey"),
                    labels = c("No behavioural\nthermoregulation", 
                               "Behavioural\nthermoregulation")) +
  theme_bw() +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(t = 1, b = 0, r = 1, l = 1),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(-29, 40), breaks = c(-20, 0, 20, 40)) +
  scale_x_continuous(limits = c(-2.5, 66))

saveRDS(latitude, "data-processed/intermediate-files/niche_behav-thermo_latitude.rds")

## break apart and save 
legend <- get_legend(latitude)

latitude <- latitude + guides(shape = "none", colour = "none", linetype = "none", 
                              fill = 'none') 

ggsave(legend, path ="figures/extended-data", 
       filename = "predictions_warm-niche_behav-across-lat_legend.png",
       width = 2, height = 3, device = "png")

ggsave(latitude, path ="figures/extended-data", 
       filename = "predictions_warm-niche_behav-across-lat.png",
       width = 4, height = 3, device = "png")

## show how distribution of warm niche filling values changes when behaviour is considered
warm = warmnorm
new_data <- data.frame(expand_grid(abs_lat_mp = median(warm$abs_lat_mp, na.rm=T),
                                   dispersal_distance_continuous = 
                                     median(warm$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(warm$log_maximum_body_size,na.rm=T)) ,
                                   log_range_area = c(median(warmnorm$log_range_area, na.rm=T)),
                                   metric = c("ct","lt")))

pred_warm <- predict(avgm_warm_behav, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm_behav <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

new_data <- data.frame(expand_grid(abs_lat_mp = median(warm$abs_lat_mp, na.rm=T),
                                   dispersal_distance_continuous = 
                                     median(warm$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(warm$log_maximum_body_size,na.rm=T)) ,
                                   log_range_area = c(median(warmnorm$log_range_area, na.rm=T)),
                                   metric = c("ct","lt")))

pred_warm <- predict(avgm_warm_no_behav, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_warm_no_behav <- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

means_behav <- fitted_pred_warm_behav %>%
  mutate(CI = 1.96*filling_value_SE) %>%
  mutate(mean = mean(filling_value, na.rm = T),
         mean_CI = mean(CI)) %>%
  select(mean, mean_CI) %>%
  unique(.) %>%
  mutate(behaviour = "Subset (behaviour)")

means_no_behav <- fitted_pred_warm_no_behav %>%
  mutate(CI = 1.96*filling_value_SE) %>%
  mutate(mean = mean(filling_value, na.rm = T),
         mean_CI = mean(CI)) %>%
  select(mean, mean_CI) %>%
  unique(.) %>%
  mutate(behaviour = "Subset (no behaviour)")

warm <- readRDS("data-processed/thermal-niches/niche-filling/warm-complete-cases.rds")
warm$behaviour = "Full dataset\n(no behaviour)"

## add full dataset:
warm_both <- rbind(warm_both, warm)

distrib <- warm_both %>%
  ggplot(., aes(x = filling_value, y = behaviour, fill = behaviour, colour = behaviour)) +
  geom_density_ridges(aes(height = ..density..), stat = "density", scale = 0.9,
                      alpha = 0.6, trim = TRUE) +
  geom_point(shape = 95, size = 4) +
  theme_ridges() + 
  theme(legend.position = "none") +
  coord_flip() +
  scale_fill_manual(values = c("#b45346", "black", "darkgrey")) +
  scale_colour_manual(values = c("#b45346", "black", "darkgrey")) +
  labs(y = "", 
       x = "", #x = "Warm niche filling (°C)"
  )  +
  scale_y_discrete(labels = c("Full dataset\n(no behaviour)", "Subset\n(no behaviour)",
                              "Subset\n(behaviour)")) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10, hjust = 0.5),
        axis.text.y = element_text(size = 8),
        panel.border = element_rect(colour = "black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank())  +
  annotate(geom = "text", label = "overfilling", y = 0.5, x = 7, angle = 90) +
  annotate(geom = "text", label = "underfilling", y = 0.5, x = -8, angle = 90) +
  geom_vline(xintercept = 0, size = 0.4) +
  scale_x_continuous(breaks = c(-20, -10, 0, 10), labels = c("-20°C", "-10°C", 
                                                             "0°C", "10°C")) 

ggsave(distrib, path ="figures/extended-data", filename = "distributions_warm-niche_behav.png",
       width = 5.5, height = 4, device = "png")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####          Acclimatization subset            ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Subset WITH adjustment ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## fitting to cold and warm filling separately
cold <- acc %>%
  filter(edge_type == "cold") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "min")  
## assign metric type used to measure cold limit 

warm <- acc %>%
  filter(edge_type == "warm") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "max") 
## assign metric type used to measure warm limit 


#Rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
cold <- select(cold, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size, cold_season_dormancy_,metric,
                       Class, Order, Family, Genus, Species))
warm <- select(warm, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size, hot_season_dormancy_,metric, 
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

#Now, run the full models
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
model_cold <- lme(filling_value ~ abs_lat_mp*realm + log_range_area +
                     dispersal_distance_continuous  +
                     log_maximum_body_size  + metric,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = cold
                  
)
model_warm <- lme(filling_value ~ abs_lat_mp*realm + log_range_area +
                    dispersal_distance_continuous  +
                    log_maximum_body_size + metric,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = warm
)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_cold <- dredge(model_cold, extra="R^2"))
confset.95p_cold <- get.models(allmodels_cold, subset = cumsum(weight)<=.95) #get confidence set
avgm_cold <- model.avg(confset.95p_cold) #do averaging
summary(avgm_cold)

suppressWarnings(allmodels_warm <- dredge(model_warm, extra="R^2")) 
confset.95p_warm <- get.models(allmodels_warm, subset = cumsum(weight)<=.95) #get confidence set
avgm_warm <- model.avg(confset.95p_warm) #do averaging
summary(avgm_warm)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Exporting tables                  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_warm)
sum_warm <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", "slope", "slope", "slope", 
              "intercept","intercept",
              "intercept")
coefs <- sum_warm[c(1,2,7,8,3,10,9,5,6,4),1] # reorder
f_effects <-  c("(Intercept)",
                "Abs. realized range latitudinal midpoint",
                "Abs. realized range latitudinal midpoint x realm: intertidal)",
                "Abs. realized range latitudinal midpoint x realm: subtidal)",
                "Realized range size (log no.cells)",
                "Dispersal distance (km)",
                "Maximum body size (log cm)",
                "Realm: intertidal",
                "Realm: subtidal",
                "Metric: lethal")
std_err <- sum_warm[c(1,2,7,8,3,10,9,5,6,4),2]  # reorder fixed effects
z_val <- sum_warm[c(1,2,7,8,3,10,9,5,6,4),4]
p_val <- sum_warm[c(1,2,7,8,3,10,9,5,6,4),5]
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
coefs <- sum_cold[c(1,2,8,9,4,10,3,6,7,5),1] # reorder

f_effects = c("(Intercept)",
              "Abs. realized range latitudinal midpoint",
              "Abs. realized range latitudinal midpoint x realm: intertidal)",
              "Abs. realized range latitudinal midpoint x realm: subtidal)",
              "Realized range size (log km2)",
              "Dispersal distance (km)",
              "Maximum body size (log cm)",
              "Realm: intertidal",
              "Realm: subtidal",
              "Metric: lethal")
std_err <- sum_cold[c(1,2,8,9,4,10,3,6,7,5),2] # reorder fixed effects
z_val <- sum_cold[c(1,2,8,9,4,10,3,6,7,5),4]
p_val <- sum_cold[c(1,2,8,9,4,10,3,6,7,5),5]
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


## save model results for plotting
avgm_warm_acc <- avgm_warm
avgm_cold_acc <- avgm_cold

## save complete cases
saveRDS(warm, "data-processed/intermediate-files/acc-warm-complete-cases.rds")
saveRDS(cold, "data-processed/intermediate-files/acc-cold-complete-cases.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Subset WITHOUT adjustment ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## fitting to cold and warm filling separately
cold <- Te_subset_acc %>%
  filter(edge_type == "cold") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "min")  
## assign metric type used to measure cold limit 

warm <- Te_subset_acc %>%
  filter(edge_type == "warm") %>%
  select(-metric) %>%
  left_join(., thermal_limits, by = c("genus_species", "type")) %>%
  filter(type == "max") 
## assign metric type used to measure warm limit 


#Rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
cold <- select(cold, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size, cold_season_dormancy_,metric,
                       Class, Order, Family, Genus, Species))
warm <- select(warm, c(filling_value, abs_lat_mp, log_range_area, realm,
                       dispersal_distance_continuous,
                       log_maximum_body_size, hot_season_dormancy_,metric, 
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

#Now, run the full models
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
model_cold <- lme(filling_value ~ abs_lat_mp*realm + log_range_area +
                    dispersal_distance_continuous  +
                    log_maximum_body_size + metric,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = cold
                  
)
model_warm <- lme(filling_value ~ abs_lat_mp*realm + log_range_area +
                    dispersal_distance_continuous  +
                    log_maximum_body_size + metric,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = warm
)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_cold <- dredge(model_cold, extra="R^2"))
confset.95p_cold <- get.models(allmodels_cold, subset = cumsum(weight)<=.95) #get confidence set
avgm_cold <- model.avg(confset.95p_cold) #do averaging
summary(avgm_cold)

suppressWarnings(allmodels_warm <- dredge(model_warm, extra="R^2")) 
confset.95p_warm <- get.models(allmodels_warm, subset = cumsum(weight)<=.95) #get confidence set
avgm_warm <- model.avg(confset.95p_warm) #do averaging
summary(avgm_warm)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Exporting tables                  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_warm)
sum_warm <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", "slope", "slope", "slope", 
              "intercept","intercept",
              "intercept")
coefs <- sum_warm[c(1,2,7,8,3,10,9,5,6,4),1] # reorder
f_effects <-  c("(Intercept)",
                "Abs. realized range latitudinal midpoint",
                "Abs. realized range latitudinal midpoint x realm: intertidal",
                "Abs. realized range latitudinal midpoint x realm: subtidal",
                "Realized range size (log no.cells)",
                "Dispersal distance (km)",
                "Maximum body size (log cm)",
                "Realm: intertidal",
                "Realm: subtidal",
                "Metric: lethal")
std_err <- sum_warm[c(1,2,7,8,3,10,9,5,6,4),2]  # reorder fixed effects
z_val <- sum_warm[c(1,2,7,8,3,10,9,5,6,4),4]
p_val <- sum_warm[c(1,2,7,8,3,10,9,5,6,4),5]
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
coefs <- sum_cold[c(1,2,8,9,4,10,3,6,7,5),1] # reorder

f_effects = c("(Intercept)",
              "Abs. realized range latitudinal midpoint",
              "Abs. realized range latitudinal midpoint x realm: intertidal)",
              "Abs. realized range latitudinal midpoint x realm: subtidal)",
              "Realized range size (log km2)",
              "Dispersal distance (km)",
              "Maximum body size (log cm)",
              "Realm: intertidal",
              "Realm: subtidal",
              "Metric: lethal")
std_err <- sum_cold[c(1,2,8,9,4,10,3,6,7,5),2] # reorder fixed effects
z_val <- sum_cold[c(1,2,8,9,4,10,3,6,7,5),4]
p_val <- sum_cold[c(1,2,8,9,4,10,3,6,7,5),5]
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


## save model results for plotting
avgm_warm_no_acc <- avgm_warm
avgm_cold_no_acc <- avgm_cold

saveRDS(warm, "data-processed/intermediate-files/no-acc-warm-complete-cases.rds")
saveRDS(cold, "data-processed/intermediate-files/no-acc-cold-complete-cases.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####                Making distribution comparison plots                 #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
warm_no_acc <- readRDS("data-processed/intermediate-files/no-acc-warm-complete-cases.rds")
warm_acc <- readRDS("data-processed/intermediate-files/acc-warm-complete-cases.rds")
cold_no_acc <- readRDS("data-processed/intermediate-files/no-acc-cold-complete-cases.rds")
cold_acc <- readRDS("data-processed/intermediate-files/acc-cold-complete-cases.rds")
warm <- readRDS("data-processed/thermal-niches/niche-filling/warm-complete-cases.rds")
cold <- readRDS("data-processed/thermal-niches/niche-filling/cold-complete-cases.rds")

warm_acc <- warm_acc[,which(colnames(warm_acc) %in% colnames(warm_no_acc))]

warm_acc$type = cold_acc$type = "acclimatized"
warm_no_acc$type = cold_no_acc$type = "not acclimatized"
warm$type = cold$type = "full data"
cold_no_acc$edge_type = cold_acc$edge_type = cold$edge_type = "Cold niche filling"
warm_no_acc$edge_type = warm_acc$edge_type = warm$edge_type = "Warm niche filling"
cold_no_acc$hot_season_dormancy_ = cold_acc$hot_season_dormancy_ = cold$hot_season_dormancy_ = NA
warm_no_acc$cold_season_dormancy_ = warm_acc$cold_season_dormancy_ = warm$cold_season_dormancy_ =NA

## combine them all
both <- rbind(warm_acc, warm_no_acc) %>%
  rbind(., warm) %>%
  rbind(., cold_acc) %>%
  rbind(., cold_no_acc) %>%
  rbind(., cold) %>%
  mutate(colour = ifelse(type == "not acclimatized", "not acclimatized", 
                         ifelse(type == "acclimatized", "acclimatized", edge_type))) %>%
  mutate(y = paste(type, realm)) %>%
  mutate(y = factor(.$y, levels = c( "full data Terrestrial", "not acclimatized Terrestrial", "acclimatized Terrestrial",
                                     "full data Intertidal", "not acclimatized Intertidal", "acclimatized Intertidal",
                                     "full data Marine", "not acclimatized Marine", "acclimatized Marine"), 
                    ordered = T))

no_int <- both %>%
  mutate(filling_value = ifelse(y %in% c("not acclimatized Intertidal", "acclimatized Intertidal",
                                         "full data Intertidal") &
                                  edge_type == "Cold niche filling",
                                NA,
                                filling_value))
int_only <- both  %>%
  mutate(filling_value = ifelse(!(y %in% c("not acclimatized Intertidal", "acclimatized Intertidal",
                                         "full data Intertidal" ) &
                                  edge_type == "Cold niche filling"),
                                NA,
                                filling_value))

## plot:
acclimatized_dist <- no_int %>%
  ggplot(., aes(x = filling_value, y = y, fill = colour, colour = colour)) +
  geom_density_ridges(aes(height = ..density..), stat = "density", scale = 0.9, 
                      alpha = 0.6, trim = TRUE)  +
  geom_density_ridges(data = int_only, aes(height = ..density..), stat = "density", scale = 0.9, 
                      alpha = 0.6, trim = TRUE) +
  geom_point(data = both, shape = 95, size = 4) +
  theme_ridges() + 
  theme(legend.position = "none") +
  coord_flip() +
  facet_wrap(~edge_type) +
  scale_fill_manual(values = c("darkgrey", "steelblue", "black", "#b45346")) +
  scale_colour_manual(values = c("darkgrey", "steelblue", "black", "#b45346")) +
  labs(y = "", 
       #x = "Niche filling (°C)"
       x = "") +
  scale_y_discrete(labels = c("Terrestrial", "", "", "Intertidal\n marine", "", "",
                              "Subtidal\n marine", "","")) +
  theme(strip.text = element_text(size = 10, vjust = 1),
        axis.text.x = element_text(size = 10, hjust = 0),
        axis.title.y = element_text(size = 10, hjust = 0.5),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())  +
  annotate(geom = "text", label = "overfilling", y = 0.55, x = 14, angle = 90) +
  annotate(geom = "text", label = "underfilling", y = 0.55, x = -15, angle = 90) +
  geom_vline(xintercept = 0, size = 0.4) +
  scale_x_continuous(breaks = c(-20, -10, 0, 10, 20), labels = c("-20°C", "-10°C", 
                                                                 "0°C", "10°C",
                                                                 "20°C"))

ggsave(acclimatized_dist, path ="figures/extended-data", filename = "distributions_niche_acc.png",
       width = 8, height = 4, device = "png")

## make legend:
acclimatized_legend <- acclimatized_dist + theme(legend.position = "right") +
  guides(colour = "none") +
  scale_fill_manual(values = c("black", "darkgrey", "steelblue", "#b45346"), 
                    labels = c("Subset (no acclimatization)", "Subset (acclimatization)", 
                    "Full dataset (no acclimatization)", "")) +
  labs(fill = "")
legend <- get_legend(acclimatized_legend)


ggsave(legend, path ="figures/extended-data", filename = "distributions_niche_acc_legend.png",
       width = 4, height = 2, device = "png")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Making whisker plots               #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(dotwhisker)

## behavioural thermoreg vs. no behavioural thermoreg:
## make model list
behav_warm <- list(avgm_warm_behav, avgm_warm_no_behav)
names(behav_warm) <- c("Behavioural thermoregulation", "No behavioural thermoregulation")

behav_dw <- dwplot(behav_warm, 
                   vline = geom_vline(xintercept = 0, colour = "grey50"),
                   show_intercept = TRUE) +
  scale_y_discrete(labels = c("Dispersal distance (km)",
                              "Maximum body size (log cm)",
                              "Metric: lethal",
                              "Realized range size (log no. cells)", 
                              "Abs. realized range latitudinal midpoint", 
                              "Reference")) +
  labs(colour = '', x = "Effect of variable on warm niche filling") +
  theme_light() +
  scale_colour_manual(values = c("#b45346", "black")) + 
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.75, 'lines'),
        axis.text.y = element_text(size = 8), 
        legend.position = "none", panel.grid = element_blank())

saveRDS(behav_dw, "data-processed/intermediate-files/whisker_behav.rds") 

ggsave(behav_dw, width = 4.5, height = 2.5, path = "figures/extended-data", 
       filename = "whisker-plot_behaviour.png", 
       device = "png")


## acclimatization vs. no acclimatization:
acc_warm <- list(avgm_warm_acc, avgm_warm_no_acc)
acc_cold <- list(avgm_cold_acc, avgm_cold_no_acc)
names(acc_warm) <- c("Acclimatiztion", "No acclimatization")
names(acc_cold) <- c("Acclimatiztion", "No acclimatization")

acc_warm_dw <- dwplot(acc_warm, 
                      vline = geom_vline(xintercept = 0, colour = "grey50"),
                      show_intercept = TRUE) +
  scale_y_discrete(labels = c("Dispersal distance (km)",
                              "Maximum body size (log cm)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Realm: subtidal",
                              "Realm: intertidal",
                              "Metric: lethal",
                              "Realized range size (log no. cells)", 
                              "Abs. realized range latitudinal midpoint", 
                              "Reference")) +
  labs(colour = '', x = "Effect of variable on warm niche filling") +
  theme_light() +
  scale_colour_manual(values = c("#b45346", "black")) + 
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.75, 'lines'),
        axis.text.y = element_text(size = 8),  legend.position = "bottom",
        panel.grid = element_blank())  + 
  scale_x_continuous(limits = c(-24, 30))
  
acc_cold_dw <- dwplot(acc_cold,
                      vline = geom_vline(xintercept = 0, colour = "grey50"),
                      show_intercept = TRUE) +
  scale_y_discrete(labels = c("Dispersal distance (km)",
                              "Maximum body size (log cm)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Realm: subtidal",
                              "Realm: intertidal",
                              "Metric: lethal",
                              "Realized range size (log no. cells)", 
                              "Abs. realized range latitudinal midpoint", 
                              "Reference")) +
  labs(colour = '', x = "Effect of variable on cool niche filling") +
  theme_light() +
  scale_colour_manual(values = c("steelblue", "black")) + 
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.75, 'lines'),
          axis.text.y = element_text(size = 8),
        legend.position = "bottom", panel.grid = element_blank()) + 
  scale_x_continuous(limits = c(-24, 30)) 

saveRDS(acc_cold_dw, "data-processed/intermediate-files/whisker_acc_cold.rds")
saveRDS(acc_warm_dw, "data-processed/intermediate-files/whisker_acc_warm.rds")

acc_dw <- ggdraw() +
  draw_plot(acc_cold_dw, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(acc_warm_dw, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(label = c("a)", "b)"),
                  x = c(0, 0.5),
                  y = c(1, 1), size = 10, 
                  color = "grey30") +
  theme(plot.background = element_rect(fill="white", color = "white"))
            
ggsave(acc_dw, width = 12, height = 5, path = "figures/extended-data", 
       filename = "whisker-plot_acclimatization_niche.png", 
       device = "png")

