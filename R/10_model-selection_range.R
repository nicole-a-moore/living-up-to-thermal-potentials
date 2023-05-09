## model selection and averaging for range filling
library(tidyverse)
library(sf)
library(raster)
library(nlme)
library(MuMIn)
library(car)
library(cowplot)
library(gridExtra)
library(grid)
library(PNWColors)
select <- dplyr::select
rename <- dplyr::rename

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       Prep for model fitting              #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in under/over/filling data:
rf <- read.csv("data-processed/potential-ranges/range-filling/rangefilling.csv")
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

rf <- rf %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  left_join(., thermal_limits, by = c('species' = 'genus_species')) %>%
  left_join(., traits, by = c('species' = 'genus_species', 'realm')) 

## give priority to IUCN ranges, then GARD, then GBIF
rf <- rf %>%
  select(range, species, source, everything()) %>%
  mutate(source = factor(.$source, levels = c("IUCN", "GARD", "GBIF"), ordered = TRUE)) %>%
  arrange(species, type, source) %>% 
  filter(!duplicated(paste(species, type)))

## make dispersal distance continuous, take log of range size and absolute value of latitude
rf <- rf %>%
  mutate(log_range_area = log(rr_cells)) %>%
  mutate(log_maximum_body_size = log(maximum_body_size_SVL_HBL_cm_)) %>%
  mutate(abs_lat_mp = abs(lat_mp)) %>%
  mutate(dispersal_distance_continuous = ifelse(dispersal_distance_category == "0-1", 1, 
                                                ifelse(dispersal_distance_category == "1-10", 10,
                                                       ifelse(dispersal_distance_category == "10-100", 100,
                                                              ifelse(dispersal_distance_category == "100+", 1000, 
                                                                     NA)))))


## calculate the proportion of unoccupied cells in the potential range
rf$prop_unoccupied <- rf$u_cells / rf$pr_cells

## write out:
write.csv(rf, "data-processed/potential-ranges/range-filling/rangefilling-metrics_model-ready.csv", row.names = FALSE)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Calculate proportion of potential range that is occupied  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in results:
rf <- read.csv("data-processed/potential-ranges/range-filling/rangefilling-metrics_model-ready.csv")
nf <- read.csv("data-processed/thermal-niches/niche-filling/thermal-niche-filling-metrics_model-ready.csv")

## calculate the proportion of occupied cells in the potential range
rf$prop_occupied <- (rf$pr_cells - rf$u_cells) / rf$pr_cells
rf$log_prop_occupied <- log(rf$prop_occupied)

## calculate difference between proportion of cells occupied in poleward versus equatorward range half
rf$log_pol_equ_diff <- log(rf$pol_fill) - log(rf$equ_fill)
rf$log_prop_occupied_equ <- log(rf$equ_fill)

## calculate difference between proportion of cells underfilled in poleward versus equatorward range half
rf$bias_in_uf <- (1- rf$equ_fill) - (1-rf$pol_fill)

## get rid of spp with no filling:
rf <- filter(rf, !is.infinite(log_prop_occupied))

## re-order factors to give desired contrasts
rf$realm <- relevel(factor(rf$realm), ref = "Terrestrial")

## split by type:
types <- group_split(rf, type)

acc <- types[[1]]
te <- types[[2]]

## give priority to IUCN ranges, then GARD, then GBIF
te <- te %>%
  select(range, species, source, everything()) %>%
  mutate(source = factor(.$source, levels = c("IUCN", "GARD", "GBIF"), ordered = TRUE)) %>%
  arrange(species, type, source) %>% 
  filter(!duplicated(species)) 

## add traits
traits <- read.csv("data-processed/traits/rangetherm-traits_all-spp.csv") %>%
  rename("realm" = Realm) %>%
  select(-limit_type) %>%
  mutate(species = paste(Genus, Species, sep = "_"))
te = left_join(te, traits)

## plot values 
pal = pnw_palette("Sunset", 3, type = "discrete")

hist <- te %>%
  filter(prop_occupied > 0) %>%
  ggplot(., aes(x = prop_occupied, fill = realm)) + geom_histogram() + theme_bw() +
  labs(x = "Range filling", y = "Frequency") +
  scale_x_continuous(limits = c(0,1)) + 
  scale_y_continuous(limits = c(0,33)) +
  theme(panel.grid = element_blank(),
        legend.position = "none") + 
  scale_fill_manual(values = pal) 

hist_log <- te %>%
  ggplot(., aes(x = log_prop_occupied, fill = realm)) + geom_histogram() + theme_bw() +
  labs(x = "Log(range filling)", y = "", fill = "") +
  theme(panel.grid = element_blank()) + 
  scale_fill_manual(values = pal, labels = c("Marine (intertidal)", "Marine (subtidal)", 
                                             "Terrestrial")) 

histos <- ggdraw() + 
  draw_plot(hist, x = 0, y = 0, width = 0.4, height = 1) +
  draw_plot(hist_log, x = 0.4, y = 0, width = 0.6, height = 1) +
  draw_plot_label(label = c("a)", "b)"),
                  x = c(0, 0.4),
                  y = c(1, 1), size = 10, 
                  color = "grey30") 

ggsave(histos, path = "figures/additional-figures", 
       filename = "distributions_range_no-log-vs-log.png", 
       width = 9, height = 3)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####           The ~trait~ models                #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#check for multicolinearity among variables (in a linear model)
modvif <- lm(log_prop_occupied ~ abs_lat_mp  + realm  +
               dispersal_distance_continuous + 
               log_maximum_body_size,
             data = te)
vif(modvif)
## removed:  dispersal_ability_category, metric, aliased

# get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## select only columns we care about:
te <- select(te, c(log_prop_occupied, log_pol_equ_diff, abs_lat_mp, realm,
                               dispersal_distance_continuous, 
                               log_maximum_body_size, rr_cells, metric,
                               Class, Order, Family, Genus, Species, range, bias_in_uf))

## write out list of species
splist_range <- select(te, Genus, Species, Class, Order, Family, Genus) %>%
  mutate(genus_species = paste(Genus, Species, sep=" "))

write.csv(splist_range, 
          "data-processed/thermal-niches/splist_range-filling.csv",
          row.names = FALSE)

# get complete cases
te <- subset(te, complete.cases(te))
dim(te)

## write out list of species
write.csv(unique(paste(te$Genus, te$Species, sep=" ")), 
          "data-processed/thermal-niches/splist_range-filling_model.csv",
          row.names = FALSE)

## re-order factors to give desired contrasts
te$realm <- relevel(factor(te$realm), ref = "Terrestrial")

#Second, rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
uf <- te

# rescale the continuous variables
uf = as.data.frame(uf)
means <- sds <- c()
ufnorm <- uf
x = 1
for (i in c(3,5,6)) {
  ufnorm[,i] <- (uf[,i]-mean(c(uf[,i]), na.rm = T))/sd(uf[,i], na.rm = T)
  means[x] = mean(c(uf[,i]))
  sds[x] = sd(c(uf[,i]))
  x = x+1
}
dim(ufnorm)
## save scaling factors:
scalers = data.frame(var = colnames(uf)[c(3,5,6)], means = means, sds = sds)

uf = ufnorm

#Now, run the full models and check out the residuals:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#run full lme model
model_uf <- lme(log_prop_occupied ~ abs_lat_mp*realm + dispersal_distance_continuous +
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
confset.95p_uf <- get.models(allmodels_uf, subset = cumsum(weight)<=.95) #get confidence set, 8 models
avgm_uf <- model.avg(confset.95p_uf) #do averaging
summary(avgm_uf)

sum <- summary(avgm_uf)
df <- as.data.frame(sum$coefmat.full) 
CI <- as.data.frame(confint(avgm_uf, full=T, level = 0.95)) # get confidence intervals for full model
df$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df$CI.max <-CI$`97.5 %`
data.table::setDT(df, keep.rownames = "coefficient") #put rownames into column
df$how_conf = ifelse((df$CI.min < 0 & df$CI.max < 0) | (df$CI.min > 0 & df$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0

whisker_uf <- ggplot(data=df, aes(x=fct_rev(coefficient), y=Estimate)) + 
  geom_hline(yintercept = 0) + 
  theme_classic() +
  geom_pointrange(aes(ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = '#b0b0b0', size = 1, shape = 21) + 
  coord_flip() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("white", "#b0b0b0")) +
  labs(y = "Effect of variable on range filling", x = "") +
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal",
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Abs. realized range latitudinal midpoint ",
                              "Reference")) +
  scale_y_continuous(limits = c(-10,10))

saveRDS(whisker_uf, "data-processed/intermediate-files/whisker_range.rds")

ggsave(whisker_uf, width = 6, height = 2.5, path = "figures/extended-data", 
       filename = "whisker-plot_range-filling.png", 
       device = "png")

### assess model averaging by looking at top model set 
library(dotwhisker)
library(RColorBrewer)
## add average model to lists so it can be visualized 
uf_list <- append(confset.95p_uf, list(avgm_uf))
names(uf_list)[9] <- "Average"

colours <- colorRampPalette(brewer.pal(8, "Set2"))(9)
colours[1] <- "black"

range_dwplot <- dwplot(uf_list, 
                       vline = geom_vline(xintercept = 0, colour = "grey50"),
                       show_intercept = TRUE) + 
  scale_y_discrete(labels = c("Dispersal distance (km)",
                              "Maximum body size (log cm)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "SD (observations)",
                              "SD (intercept)",
                              "Realm: subtidal",
                              "Realm: intertidal",
                              "Abs. realized range latitudinal midpoint ",
                              "Reference")) +
  labs(colour = '', x = "Effect of variable on range filling") +
  theme_light() +
  scale_color_manual(values = colours) +
  theme(panel.grid = element_blank())  + 
  scale_x_continuous(limits = c(-10,10))

saveRDS(range_dwplot, "data-processed/intermediate-files/whisker_range_model-avg.rds")

## combine with other plots
warm_dwplot <- readRDS("data-processed/intermediate-files/whisker_warm_model-avg.rds")
cold_dwplot <- readRDS("data-processed/intermediate-files/whisker_cold_model-avg.rds")
range_dwplot <- readRDS("data-processed/intermediate-files/whisker_range_model-avg.rds")

warm_dwplot <- warm_dwplot + theme(legend.text = element_text(size = 6), legend.key.size = unit(0.75, 'lines'),
                                   axis.text.y = element_text(size = 8))
cold_dwplot <- cold_dwplot + theme(legend.text = element_text(size = 6), legend.key.size = unit(0.75, 'lines'),
                                   axis.text.y = element_text(size = 8))
range_dwplot <- range_dwplot + theme(legend.text = element_text(size = 6), legend.key.size = unit(0.75, 'lines'),
                                     axis.text.y = element_text(size = 8))

model_check <- ggdraw() + 
  draw_plot(warm_dwplot, x = 0, y = 0.63, width = 1, height = 0.37) +
  draw_plot(cold_dwplot, x = 0, y = 0.26, width = 1, height = 0.37) + 
  draw_plot(range_dwplot, x = 0, y = 0, width = 1, height = 0.26) +
  draw_plot_label(label = c("a)", "b)", "c)"),
                  x = c(0, 0,0),
                  y = c(1, 0.63, 0.26), size = 10, 
                  color = "grey30") 

ggsave(model_check, height = 6, width = 6.5, path = "figures/extended-data", 
       filename = "whisker-plot_model-averaging-check.png", 
       device = "png")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Exporting tables                  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(htmlTable)

## make data frame of details I want to include in our table:
sum <- summary(avgm_uf)
sum_uf <- sum$coefmat.full
eff_type <- c("intercept", "slope", "slope", "slope", "slope","slope","intercept","intercept")
coefs <- sum_uf[c(1,2,5,6,8,7,3,4),1] # reorder
f_effects <- c("(Intercept)",
               "Abs. realized range latitudinal midpoint",
               "Abs. realized range latitudinal midpoint x realm: subtidal",
               "Abs. realized range latitudinal midpoint x realm: intertidal",
               "Dispersal distance (km)",
               "Maximum body size (log cm)",
               "Realm: intertidal",
               "Realm: marine")
std_err <- sum_uf[c(1,2,5,6,8,7,3,4),2]  # reorder fixed effects
z_val <- sum_uf[c(1,2,5,6,8,7,3,4),4]
p_val <- sum_uf[c(1,2,5,6,8,7,3,4),5]
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
  htmlTable(., rnames = rep("", 8))

## write out species list 
splist_asymm <- select(uf, Genus, Species, Class, Order, Family, Genus) %>%
  mutate(genus_species = paste(Genus, Species, sep=" "))

write.csv(splist_asymm, 
          "data-processed/thermal-niches/splist_asymm-range-filling.csv",
          row.names = FALSE)

## fit model to bias in range underfilling 
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


#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_asym <- dredge(model_asym, extra="R^2", subset = (`abs_lat_mp:realm`)))
confset.95p_asym <- get.models(allmodels_asym, subset = cumsum(weight)<=.95) #get confidence set, 5 models
avgm_asym <- model.avg(confset.95p_asym) # only 1 model so cannot average 
summary(confset.95p_asym[[1]])
avgm_asym <- confset.95p_asym[[1]]
  
sum <- summary(confset.95p_asym[[1]])
df <- as.data.frame(sum$coefficients$fixed) 
intervals <- as.data.frame(intervals(confset.95p_asym[[1]], which = "fixed")$fixed)
df$CI.min <- intervals$lower
df$CI.max <- intervals$upper  # get confidence intervals for model
data.table::setDT(df, keep.rownames = "coefficient") #put rownames into column
df$how_conf = ifelse((df$CI.min < 0 & df$CI.max < 0) | (df$CI.min > 0 & df$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0

whisker_asym <- ggplot(data=df, aes(x=fct_rev(coefficient), y=`sum$coefficients$fixed`)) + 
  geom_hline(yintercept = 0) + 
  theme_classic() +
  geom_pointrange(aes(ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = '#b0b0b0', size = 1, shape = 21) + 
  coord_flip() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("white", "#b0b0b0")) +
  labs(y = "Effect of variable on equatorward bias in underfilling", x = "") +
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Abs. realized range latitudinal midpoint ",
                              "Reference")) +
  scale_y_continuous(limits = c(-0.8,0.8))

saveRDS(whisker_asym, "data-processed/intermediate-files/whisker_asym-underfilling.rds")

ggsave(whisker_asym, width = 7, height = 2.5, path = "figures/extended-data", 
       filename = "whisker-plot_asym-underfilling.png", 
       device = "png")
 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Exporting tables                  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## make data frame of details I want to include in our table:
sum <- summary(confset.95p_asym[[1]])
sum_uf <- sum$coefficients$fixed
eff_type <- c("intercept", "slope", "slope", "slope","intercept","intercept")
coefs <- sum_uf[c(1,2,5,6,3,4)] # reorder
f_effects <- c("(Intercept)",
               "Abs. realized range latitudinal midpoint",
               "Abs. realized range latitudinal midpoint x realm: subtidal",
               "Abs. realized range latitudinal midpoint x realm: intertidal",
               "Realm: intertidal",
               "Realm: marine")
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
  htmlTable(., rnames = rep("", 8))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         Plotting model predictions         ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# function to find mode of each predictor
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####                Asymmetry                 ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## plot predictions
uf <- ufnorm

int <- uf %>% filter(realm == "Intertidal")
mar <- uf %>% filter(realm == "Marine")
terr <- uf %>% filter(realm == "Terrestrial")

new_data <- data.frame(expand_grid(abs_lat_mp = seq(min(uf$abs_lat_mp), max(uf$abs_lat_mp), 
                                                    length.out = 1000),
                                   realm = c("Terrestrial", 'Marine', 'Intertidal')))

new_data$realm <- factor(new_data$realm, levels = c("Terrestrial", 'Intertidal', 'Marine'))

pred <- predict(avgm_asym, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred <- new_data %>%
  mutate(filling_value = pred$fit) %>%
  mutate(filling_value_SE = pred$se.fit)

## Convert scaled prediction to original data scale of log_range_area
fitted_pred$abs_lat_mp <- fitted_pred$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
uf$abs_lat_mp <- uf$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
terr$abs_lat_mp <- terr$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
int$abs_lat_mp <-int$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
mar$abs_lat_mp <- mar$abs_lat_mp * scalers$sd[1] + scalers$mean[1]

realm_terr <- fitted_pred %>%
  filter(realm == "Terrestrial") %>%
  filter(abs_lat_mp >= min(terr$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(terr$abs_lat_mp, na.rm=T) ) 

realm_mar <- fitted_pred %>%
  filter(realm == "Marine") %>%
  filter(abs_lat_mp >= min(mar$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(mar$abs_lat_mp, na.rm=T) ) 

realm_int <- fitted_pred %>%
  filter(realm == "Intertidal") %>%
  filter(abs_lat_mp >= min(int$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(int$abs_lat_mp, na.rm=T) )

## plot predictions
realms_all <- rbind(realm_int, realm_mar) %>%
  rbind(., realm_terr)

asym_lat <- uf %>%
  ggplot(., aes(x = abs_lat_mp, y = bias_in_uf,
                shape = realm, colour = bias_in_uf, fill = bias_in_uf, group = realm, linetype = realm)) +
  labs(col = "", shape = "", x =  "Absolute realized range latitudinal midpoint (°N/S)", 
       fill = "", linetype = "",
       y =  "Equatorward bias in range underfilling") +
  scale_fill_gradient2(midpoint = 0, high = "#b45346", low = "steelblue", mid = "white",
                       limits = c(-1, 1)) +
  geom_ribbon(data = realms_all, inherit.aes = F, 
              aes(x = abs_lat_mp, 
                  ymin = (filling_value-1.96*filling_value_SE),
                  ymax = filling_value+1.96*filling_value_SE,
                  group = realm), 
              fill = "darkgrey",
              alpha=0.15) +
  geom_point(aes(fill = bias_in_uf), size = 2.5, stroke = 0.5, colour = "black") +
  geom_line(data = realms_all, inherit.aes = F,
            aes(x = abs_lat_mp, y = filling_value, 
                group = realm, linetype = realm)) +
  theme_bw() +
  theme(plot.margin = margin(t = 1, b = 0, r = 1, l = 1),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  scale_linetype_discrete(labels = c("Terrestrial", "Intertidal marine", "Subtidal marine")) +
  scale_shape_manual(labels = c("Terrestrial", "Intertidal marine", "Subtidal marine"),
                     values = c(21,22,24)) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  scale_y_continuous(limits = c(-0.5,1)) 

saveRDS(asym_lat, "data-processed/intermediate-files/asym-underfilling_latitude.rds")

ggsave(asym_lat, path ="figures/main", filename = "predictions_asym_latitude.png",
       width = 6.5, height = 3.5, device = "png")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####             Range filling                ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
uf <- ufnorm

int <- uf %>% filter(realm == "Intertidal")
mar <- uf %>% filter(realm == "Marine")
terr <- uf %>% filter(realm == "Terrestrial")

new_data <- data.frame(expand_grid(abs_lat_mp = seq(min(uf$abs_lat_mp), max(uf$abs_lat_mp), 
                                                    length.out = 1000),
                                   dispersal_distance_continuous = 
                                     unique(uf$dispersal_distance_continuous),
                                   log_maximum_body_size = 
                                     c(median(terr$log_maximum_body_size,
                                              na.rm=T),
                                       median(int$log_maximum_body_size,
                                              na.rm=T),
                                       median(mar$log_maximum_body_size,
                                              na.rm=T)), 
                                   realm = c("Terrestrial", 'Marine', 'Intertidal')))

new_data$realm <- factor(new_data$realm, levels = c("Terrestrial", 'Intertidal', 'Marine'))

pred <- predict(avgm_uf, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred <- new_data %>%
  mutate(filling_value = pred$fit) %>%
  mutate(filling_value_SE = pred$se.fit)

## Convert scaled prediction to original data scale of log_range_area
fitted_pred$abs_lat_mp <- fitted_pred$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
uf$abs_lat_mp <- uf$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
terr$abs_lat_mp <- terr$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
int$abs_lat_mp <-int$abs_lat_mp * scalers$sd[1] + scalers$mean[1]
mar$abs_lat_mp <- mar$abs_lat_mp * scalers$sd[1] + scalers$mean[1]

realm_terr <- fitted_pred %>%
  filter(realm == "Terrestrial") %>%
  filter(dispersal_distance_continuous ==
           getmode(terr$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(terr$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(terr$abs_lat_mp, na.rm=T) ) %>%
  filter(log_maximum_body_size == median(terr$log_maximum_body_size,
                                         na.rm=T)) 

realm_mar <- fitted_pred %>%
  filter(realm == "Marine") %>%
  filter(dispersal_distance_continuous ==
           getmode(mar$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(mar$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(mar$abs_lat_mp, na.rm=T) ) %>%
  filter(log_maximum_body_size == median(mar$log_maximum_body_size,
                                         na.rm=T))  

realm_int <- fitted_pred %>%
  filter(realm == "Intertidal") %>%
  filter(dispersal_distance_continuous ==
           getmode(int$dispersal_distance_continuous)) %>%
  filter(abs_lat_mp >= min(int$abs_lat_mp, na.rm=T) & 
           abs_lat_mp <= max(int$abs_lat_mp, na.rm=T) ) %>%
  filter(log_maximum_body_size == median(int$log_maximum_body_size,
                                         na.rm=T))  

## plot predictions
realms_all <- rbind(realm_int, realm_mar) %>%
  rbind(., realm_terr)

latitude <- uf %>%
  mutate(realm = ifelse(realm == "Intertidal", "Intertidal marine", 
                        ifelse(realm == "Marine", "Subtidal marine",
                               "Terrestrial"))) %>%
  mutate(realm = factor(.$realm, levels = c("Terrestrial", "Intertidal marine", "Subtidal marine"), 
                        ordered = TRUE)) %>%
  ggplot(., aes(x = abs_lat_mp, y = log_prop_occupied,
                shape = realm, colour = bias_in_uf, fill = bias_in_uf, group = realm)) +
  labs(col = "", shape = "", x =  "Absolute realized range latitudinal midpoint (°N/S)", 
       fill = "", linetype = "",
       y =  "Log proportion\nof range filled") +
  theme_bw() +
  scale_fill_gradient2(midpoint = 0, high = "#b45346", low = "steelblue", mid = "white",
                       limits = c(-1, 1)) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  theme(plot.margin = margin(t = 1, b = 0, r = 1, l = 1),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  geom_point(aes(fill = bias_in_uf), size = 2.5, stroke = 0.5, colour = "black") +
  scale_linetype_discrete(labels = c("Terrestrial", "Intertidal marine", "Subtidal marine")) +
  scale_shape_manual(labels = c("Terrestrial", "Intertidal marine", "Subtidal marine"),
                     values = c(21,22,24)) +
  scale_y_continuous(limits = c(-7.2, 0.1)) +
  facet_wrap(~realm) + 
  theme(legend.position = "none")

saveRDS(latitude, "data-processed/intermediate-files/range_latitude.rds")

ggsave(latitude, path ="figures/main", filename = "predictions_range_latitude.png",
       width = 6.5, height = 2.5, device = "png")

## combine all whisker plots and save figure:
whisker_uf <- readRDS("data-processed/intermediate-files/whisker_range.rds")
whisker_warm <- readRDS("data-processed/intermediate-files/whisker_warm.rds")
whisker_cold <- readRDS("data-processed/intermediate-files/whisker_cold.rds")
whisker_asym <- readRDS("data-processed/intermediate-files/whisker_asym-underfilling.rds")

whisker <- ggdraw() + 
  draw_plot(whisker_warm, x = 0, y = 0.73, width = 1, height = 0.27) +
  draw_plot(whisker_cold, x = 0, y = 0.46, width = 1, height = 0.27) + 
  draw_plot(whisker_uf, x = 0, y = 0.23, width = 1, height = 0.23) +
  draw_plot(whisker_asym, x = 0, y = 0.03, width = 1, height = 0.20) +
  draw_plot_label(label = c("a)", "b)", "c)", "d)"),
                  x = c(0, 0, 0, 0),
                  y = c(1, 0.73, 0.46, 0.23), size = 10, 
                  color = "grey30") 

ggsave(whisker, width = 7, height = 8, path = "figures/extended-data", 
       filename = "whisker-plot_warm-cold-range.png", 
       device = "png")


### make table of species names for appendix:
c = read.csv("data-processed/thermal-niches/splist_cold-niche-filling.csv")
w = read.csv("data-processed/thermal-niches/splist_warm-niche-filling.csv") 
r = read.csv("data-processed/thermal-niches/splist_range-filling.csv") 
ur = read.csv("data-processed/thermal-niches/splist_asymm-range-filling.csv") 

c_sp = select(c, genus_species) %>% rename("x" = genus_species)
w_sp = select(w, genus_species) %>% rename("x" = genus_species)
r_sp = select(r, genus_species) %>% rename("x" = genus_species)

c_m = read.csv("data-processed/thermal-niches/splist_cold-niche-filling_model.csv") %>%
  mutate("Cool niche filling" = "x")
w_m = read.csv("data-processed/thermal-niches/splist_warm-niche-filling_model.csv") %>%
  mutate("Warm niche filling" = "x")
r_m = read.csv("data-processed/thermal-niches/splist_range-filling_model.csv") %>%
  mutate("Range filling" = "x")
ur = read.csv("data-processed/thermal-niches/splist_asymm-range-filling.csv") %>%
  select(genus_species) %>% 
  rename("x" = genus_species) %>%
  mutate("Bias in underfilling" = "x") 
  
appendix <- full_join(c_sp, w_sp) %>%
  full_join(., r_sp) %>%
  left_join(., c_m)%>%
  left_join(., w_m)%>%
  left_join(., r_m)%>%
  left_join(., ur) %>%
  filter(x != 'Tarentola boettgeri')

write.csv(appendix, "data-processed/AppendixA.csv")

### make table showing taxonomic breakdown for appendix:
appendix2 <- full_join(c, w) %>%
  full_join(., r) %>%
  left_join(., c) %>%
  filter(genus_species != 'Tarentola boettgeri')

## plot by class/clade 
appendix2 %>%
  ggplot(aes(x = Class)) + geom_bar() +
  coord_flip() 

appendix2 <- appendix2 %>% 
  group_by(Class) %>%
  tally() %>%
  arrange(-n) %>%
  rename("Number of species" = "n", 
         "Class or clade" = "Class") 

write.csv(appendix2, "data-processed/AppendixA2.csv", row.names = F)

