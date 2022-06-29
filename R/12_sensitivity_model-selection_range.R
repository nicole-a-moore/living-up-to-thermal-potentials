## creating null model for niche filling in geographic space
library(tidyverse)
library(sf)
library(raster)
library(nlme)
library(MuMIn)
library(car)
select <- dplyr::select


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Calculate proportion of potential range that is occupied  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in results:
rf <- read.csv("data-processed/potential-ranges/range-filling/rangefilling-metrics_model-ready.csv")
nf <- read.csv("data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics_model-ready.csv")

## calculate the proportion of occupied cells in the potential range
rf$prop_occupied <- (rf$pr_cells - rf$u_cells) / rf$pr_cells
rf$log_prop_occupied <- log(rf$prop_occupied)
rf <- filter(rf, !is.infinite(log_prop_occupied))

## re-order factors to give desired contrasts
rf$realm <- relevel(factor(rf$realm), ref = "Terrestrial")

## calculate difference between proportion of cells underfilled in poleward versus equatorward range half
rf$bias_in_uf <- (1- rf$equ_fill) - (1-rf$pol_fill)

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

## add traits
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv") %>%
  rename("realm" = Realm) %>%
  select(-limit_type) %>%
  mutate(species = paste(Genus, Species, sep = "_"))
te_subset_acc = left_join(te_subset_acc, traits)
acc = left_join(acc, traits)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Subset WITH adjustment   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#check for multicolinearity among variables (in a linear model)
modvif <- lm(log_prop_occupied ~ abs_lat_mp  + realm  +
               dispersal_distance_continuous  +
               log_maximum_body_size,
             data = acc)
vif(modvif)
## removed:  dispersal_ability_category, metric, aliased

# get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
acc <- select(acc, c(log_prop_occupied, abs_lat_mp, realm, 
                   dispersal_distance_continuous,
                   log_maximum_body_size, rr_cells, metric,
                   Class, Order, Family, Genus, Species, range, bias_in_uf))


# get complete cases
acc <- subset(acc, complete.cases(acc))
dim(acc)

## re-order factors to give desired contrasts
acc$realm <- relevel(factor(acc$realm), ref = "Terrestrial")


#Second, rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
uf <- acc

# rescale the continuous variables
uf = as.data.frame(uf)
ufnorm <- uf
for (i in c(2,4,5)) {
  ufnorm[,i] <- (uf[,i]-mean(c(uf[,i]), na.rm = T))/sd(uf[,i], na.rm = T)
}
dim(ufnorm)

uf = ufnorm

#Now, run the full models and check out the residuals:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#run full lme model
model_uf <- lme(log_prop_occupied ~ abs_lat_mp + realm + dispersal_distance_continuous +
                  log_maximum_body_size,
                
                random = ~1|Class/Order/Family/Genus,
                
                data = uf
)
#check residuals
hist(resid(model_uf))
plot(model_uf)

#against predictors:
E1 <- resid(model_uf)
plot(E1 ~ abs_lat_mp, data = uf)
plot(E1 ~ log_range_area, data = uf)
plot(E1 ~ realm, data = uf)
plot(E1 ~ dispersal_distance_continuous, data = uf)
plot(E1 ~ log_maximum_body_size, data = uf)

## check for heteroskedasticity
library(lmtest)
bptest(model_uf)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_uf <- dredge(model_uf, extra="R^2"))
confset.95p_uf <- get.models(allmodels_uf, subset = cumsum(weight)<=.95) #get confidence set
avgm_uf <- model.avg(confset.95p_uf) #do averaging
summary(avgm_uf)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####                Exporting tables           ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_uf)
sum_warm <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", 
              "intercept","intercept")
coefs <- sum_warm[c(1,6,4,5,2,3),1] # reorder
f_effects <- c("(Intercept)",
               "Abs. realized range latitudinal midpoint",
               "Dispersal distance (km)",
               "Maximum body size (log cm)",
               "Realm: intertidal",
               "Realm: subtidal")
std_err <- sum_warm[c(1,6,4,5,2,3),2]  # reorder fixed effects
z_val <- sum_warm[c(1,6,4,5,2,3),4]
p_val <- sum_warm[c(1,6,4,5,2,3),5]
names(coefs) <- names(std_err) <- names(z_val) <- names(p_val) <- NULL

## put all into a table:
results <- data.frame("fixed effects" = f_effects, 
                        "effect type" = eff_type,
                        "estimate" = coefs,
                        "s.e." = std_err,
                        "z-value" = z_val,
                        "p-value" = p_val)
colnames(results) = c("fixed effects",	"effect type",	"estimate", "s.e.",'z-value',	"p-value")

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
  htmlTable(., rnames = rep("", 6))

## save model for whisker plot:
avgm_uf_acc <- avgm_uf

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Subset WITHOUT adjustment ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#check for multicolinearity among variables (in a linear model)
modvif <- lm(log_prop_occupied ~ abs_lat_mp  + realm  +
               dispersal_distance_continuous + 
               log_maximum_body_size,
             data = te_subset_acc)
vif(modvif)
## removed:  dispersal_ability_category, metric, aliased

# get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
te_subset_acc <- select(te_subset_acc, c(log_prop_occupied, abs_lat_mp, realm, 
                     dispersal_distance_continuous,
                     log_maximum_body_size, rr_cells, metric,
                     Class, Order, Family, Genus, Species, range, bias_in_uf))


# get complete cases
te_subset_acc <- subset(te_subset_acc, complete.cases(te_subset_acc))
dim(te_subset_acc)

## re-order factors to give desired contrasts
te_subset_acc$realm <- relevel(factor(te_subset_acc$realm), ref = "Terrestrial")

#Second, rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
uf <- te_subset_acc

# rescale the continuous variables
uf = as.data.frame(uf)
ufnorm <- uf
for (i in c(2,4,5)) {
  ufnorm[,i] <- (uf[,i]-mean(c(uf[,i]), na.rm = T))/sd(uf[,i], na.rm = T)
}
dim(ufnorm)

uf = ufnorm

#Now, run the full models and check out the residuals:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#run full lme model
model_uf <- lme(log_prop_occupied ~ abs_lat_mp + realm + dispersal_distance_continuous +
                  log_maximum_body_size,
                
                random = ~1|Class/Order/Family/Genus,
                
                data = uf
)
#check residuals
hist(resid(model_uf))
plot(model_uf)

#against predictors:
E1 <- resid(model_uf)
plot(E1 ~ abs_lat_mp, data = uf)
plot(E1 ~ log_range_area, data = uf)
plot(E1 ~ realm, data = uf)
plot(E1 ~ dispersal_distance_continuous, data = uf)
plot(E1 ~ log_maximum_body_size, data = uf)

## check for heteroskedasticity
library(lmtest)
bptest(model_uf)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_uf <- dredge(model_uf, extra="R^2"))
confset.95p_uf <- get.models(allmodels_uf, subset = cumsum(weight)<=.95) #get confidence set
avgm_uf <- model.avg(confset.95p_uf) #do averaging
summary(avgm_uf)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####                Exporting tables           ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_uf)
sum_warm <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", 
              "intercept","intercept")
coefs <- sum_warm[c(1,4,5,6,2,3),1] # reorder
f_effects <- c("(Intercept)",
               "Abs. realized range latitudinal midpoint",
               "Realized range size (log no. cells)",
               "Maximum body size (log cm)",
               "Realm: intertidal",
               "Realm: subtidal")
std_err <- sum_warm[c(1,4,5,6,2,3),2]  # reorder fixed effects
z_val <- sum_warm[c(1,4,5,6,2,3),4]
p_val <- sum_warm[c(1,4,5,6,2,3),5]
names(coefs) <- names(std_err) <- names(z_val) <- names(p_val) <- NULL

## put all into a table:
results <- data.frame("fixed effects" = f_effects, 
                      "effect type" = eff_type,
                      "estimate" = coefs,
                      "s.e." = std_err,
                      "z-value" = z_val,
                      "p-value" = p_val)
colnames(results) = c("fixed effects",	"effect type",	"estimate", "s.e.",'z-value',	"p-value")

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
  htmlTable(., rname = rep("", 6))

## save model for whisker plot:
avgm_uf_no_acc <- avgm_uf


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Making whisker plots               #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(dotwhisker)

## acclimatization vs. no acclimatization:
acc_uf <- list(avgm_uf_acc, avgm_uf_no_acc)
names(acc_uf) <- c("Acclimatization", "No acclimatization")

acc_uf_dw <- dwplot(acc_uf, 
                    vline = geom_vline(xintercept = 0, colour = "grey50"),
                    show_intercept = TRUE) + 
  scale_y_discrete(labels = c("x realm: intertidal Abs. realized range latitudinal midpoint", 
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "Realm: subtidal",
                              "Realm: intertidal",
                              "Reference")) +
  labs(colour = "", x = "Effect of variable on range filling") +
  theme_light() +
  scale_colour_manual(values = c("grey", "black")) + 
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.75, 'lines'),
        axis.text.y = element_text(size = 8),  legend.position = "bottom",
        panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-5, 6.25))

ggsave(acc_uf_dw, width = 6, height = 3, path = "figures/extended-data", 
       filename = "whisker-plot_acclimatization_range.png", 
       device = "png")

saveRDS(acc_uf_dw, "data-processed/intermediate-files/whisker_acc_range.rds")

## combine with other acc plots:
acc_cold_dw <- readRDS("data-processed/intermediate-files/whisker_acc_cold.rds") 
acc_warm_dw <- readRDS("data-processed/intermediate-files/whisker_acc_warm.rds")

acc_cold_dw <- acc_cold_dw + theme(legend.position = "none")
acc_warm_dw <- acc_warm_dw + theme(legend.position = "none")
acc_uf_dw <- acc_uf_dw + theme(legend.position = "none")

acc_dw <-  ggdraw() + 
  draw_plot(acc_warm_dw, x = 0, y = 0.63, width = 1, height = 0.37) +
  draw_plot(acc_cold_dw, x = 0, y = 0.26, width = 1, height = 0.37) + 
  draw_plot(acc_uf_dw, x = 0, y = 0, width = 1, height = 0.26) +
  draw_plot_label(label = c("b)", "c)", "d)"),
                  x = c(0, 0,0),
                  y = c(1, 0.63, 0.26), size = 10, 
                  color = "grey30") 

ggsave(acc_dw, width = 6, height = 4.5, path = "figures/extended-data", 
       filename = "whisker-plot_acclimatization_all.png", 
       device = "png")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Subset WITH adjustment   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## bias model:
uf = ufnorm

#Now, run the full models and check out the residuals:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#run full lme model
model_asym <- lme(bias_in_uf ~ abs_lat_mp*realm + dispersal_distance_continuous +
                  log_maximum_body_size,
                
                random = ~1|Class/Order/Family/Genus,
                
                data = uf
)
#check residuals
hist(resid(model_asym))
plot(model_asym)

#against predictors:
E1 <- resid(model_asym)
plot(E1 ~ abs_lat_mp, data = uf)
plot(E1 ~ realm, data = uf)
plot(E1 ~ dispersal_distance_continuous, data = uf)
plot(E1 ~ log_maximum_body_size, data = uf)

## check for heteroskedasticity
library(lmtest)
bptest(model_asym)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_asym <- dredge(model_asym, extra="R^2", subset = (`abs_lat_mp:realm`)))
confset.95p_asym <- get.models(allmodels_asym, subset = cumsum(weight)<=.95) #get confidence set
avgm_asym <- model.avg(confset.95p_asym) #do averaging
summary(confset.95p_asym[[1]])
avgm_asym <- confset.95p_asym[[1]]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####                Exporting tables           ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_asym)
sum_warm <- sum$coefficients$fixed
eff_type <- c("intercept", "slope", "slope", "slope", 
              "intercept","intercept")
coefs <- sum_warm[c(1,2,5,6,3,4)] # reorder
f_effects <- c("(Intercept)",
               "Abs. realized range latitudinal midpoint",
               "Abs. realized range latitudinal midpoint x realm: intertidal",
               "Abs. realized range latitudinal midpoint x realm: subtidal",
               "Realm: intertidal",
               "Realm: subtidal")
ttable <- as.data.frame(sum$tTable)
std_err <- ttable[c(1,2,5,6,3,4),2]  # reorder fixed effects
z_val <- ttable[c(1,2,5,6,3,4),4]
p_val <- ttable[c(1,2,5,6,3,4),5]
names(coefs) <- names(std_err) <- names(z_val) <- names(p_val) <- NULL

## put all into a table:
results <- data.frame("fixed effects" = f_effects, 
                      "effect type" = eff_type,
                      "estimate" = coefs,
                      "s.e." = std_err,
                      "z-value" = z_val,
                      "p-value" = p_val)
colnames(results) = c("fixed effects",	"effect type",	"estimate", "s.e.",'z-value',	"p-value")

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
  htmlTable(., rname = rep("", 6))

## save model for whisker plot:
avgm_asym_acc <- avgm_asym

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Subset WITHOUT adjustment   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#Now, run the full models and check out the residuals:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#run full lme model
model_asym <- lme(bias_in_uf ~ abs_lat_mp*realm,
                  
                  random = ~1|Class/Order/Family/Genus,
                  
                  data = te_subset_acc
)
#check residuals
hist(resid(model_asym))
plot(model_asym)

#against predictors:
E1 <- resid(model_asym)
plot(E1 ~ abs_lat_mp, data = te_subset_acc)
plot(E1 ~ realm, data = te_subset_acc)
plot(E1 ~ dispersal_distance_continuous, data = te_subset_acc)
plot(E1 ~ log_maximum_body_size, data = te_subset_acc)

## check for heteroskedasticity
library(lmtest)
bptest(model_asym)

avgm_asym <- model_asym

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####                Exporting tables           ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_asym)
sum_warm <- sum$coefficients$fixed
eff_type <- c("intercept", "slope", "slope", "slope", 
              "intercept","intercept")
coefs <- sum_warm[c(1,2,5,6,3,4)] # reorder
f_effects <- c("(Intercept)",
               "Abs. realized range latitudinal midpoint",
               "Abs. realized range latitudinal midpoint x realm: intertidal",
               "Abs. realized range latitudinal midpoint x realm: subtidal",
               "Realm: intertidal",
               "Realm: subtidal")
ttable <- as.data.frame(sum$tTable)
std_err <- ttable[c(1,2,5,6,3,4),2]  # reorder fixed effects
z_val <- ttable[c(1,2,5,6,3,4),4]
p_val <- ttable[c(1,2,5,6,3,4),5]
names(coefs) <- names(std_err) <- names(z_val) <- names(p_val) <- NULL

## put all into a table:
results <- data.frame("fixed effects" = f_effects, 
                      "effect type" = eff_type,
                      "estimate" = coefs,
                      "s.e." = std_err,
                      "z-value" = z_val,
                      "p-value" = p_val)
colnames(results) = c("fixed effects",	"effect type",	"estimate", "s.e.",'z-value',	"p-value")

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
  htmlTable(., rname = rep("", 6))

## save model for whisker plot:
avgm_asym_no_acc <- avgm_asym

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Making whisker plots               #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(dotwhisker)

## acclimatization vs. no acclimatization:
acc_uf <- list(avgm_uf_acc, avgm_uf_no_acc)
names(acc_uf) <- c("Acclimatization", "No acclimatization")

acc_uf_dw <- dwplot(acc_uf, 
                    vline = geom_vline(xintercept = 0, colour = "grey50"),
                    show_intercept = TRUE) + 
  scale_y_discrete(labels = c("intertidal Abs. realized range latitudinal midpoint", 
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "Realm: subtidal",
                              "Realm: intertidal",
                              "Reference")) +
  labs(colour = "", x = "Effect of variable on range filling") +
  theme_light() +
  scale_colour_manual(values = c("grey", "black")) + 
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.75, 'lines'),
        axis.text.y = element_text(size = 8),  legend.position = "bottom",
        panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-5, 6.25))

ggsave(acc_uf_dw, width = 6, height = 3, path = "figures/extended-data", 
       filename = "whisker-plot_acclimatization_range.png", 
       device = "png")

saveRDS(acc_uf_dw, "data-processed/intermediate-files/whisker_acc_range.rds")

## bias 
## acclimatization vs. no acclimatization:
acc_asym <- list(avgm_asym_acc, avgm_asym_no_acc)
names(acc_asym) <- c("Acclimatization", "No acclimatization")

acc_asym_dw <- dwplot(acc_asym, 
                    vline = geom_vline(xintercept = 0, colour = "grey50"),
                    show_intercept = TRUE) + 
  scale_y_discrete(labels = c("SD obs.",
                              "SD int", "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Realm: subtidal",
                              "Realm: intertidal",
                              "Abs. realized range latitudinal midpoint", 
                              "Reference")) +
  labs(colour = "", x = "Effect of variable on asymmetry in underfilling") +
  theme_light() +
  scale_colour_manual(values = c("grey", "black")) + 
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.75, 'lines'),
        axis.text.y = element_text(size = 8),  legend.position = "bottom",
        panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-0.65, 0.8))

ggsave(acc_asym_dw, width = 6, height = 3, path = "figures/extended-data", 
       filename = "whisker-plot_acclimatization_asym.png", 
       device = "png")

saveRDS(acc_asym_dw, "data-processed/intermediate-files/whisker_acc_asym.rds")


## combine whisker plots:
acc_cold_dw <- readRDS("data-processed/intermediate-files/whisker_acc_cold.rds") 
acc_warm_dw <- readRDS("data-processed/intermediate-files/whisker_acc_warm.rds")

acc_cold_dw <- acc_cold_dw + theme(legend.position = "none")
acc_warm_dw <- acc_warm_dw + theme(legend.position = "none")
acc_uf_dw <- acc_uf_dw + theme(legend.position = "none")
acc_asym_dw <- acc_asym_dw + theme(legend.position = "none")

acc_dw <- ggdraw() + 
  draw_plot(acc_warm_dw, x = 0, y = 0.73, width = 1, height = 0.27) +
  draw_plot(acc_cold_dw, x = 0, y = 0.46, width = 1, height = 0.27) + 
  draw_plot(acc_uf_dw, x = 0, y = 0.25, width = 1, height = 0.21) +
  draw_plot(acc_asym_dw, x = 0, y = 0, width = 1, height = 0.25) +
  draw_plot_label(label = c("b)", "c)", "d)", "e)"),
                  x = c(0, 0, 0, 0),
                  y = c(1, 0.73, 0.46, 0.25), size = 10, 
                  color = "grey30") 

ggsave(acc_dw, width = 6.2, height = 6.5, path = "figures/extended-data", 
       filename = "whisker-plot_acclimatization_all.png", 
       device = "png")

