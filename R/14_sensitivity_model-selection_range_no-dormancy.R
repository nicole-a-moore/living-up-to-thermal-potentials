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
select <- dplyr::select
rename <- dplyr::rename

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

## calculate difference between proportion of cells underfilled in poleward versus equatorward range half
rf$bias_in_uf <- (1- rf$equ_fill) - (1-rf$pol_fill)

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

## get rid of species with dormancy 
te <- filter(te, cold_season_dormancy_ != "Yes", hot_season_dormancy_ != "Yes")

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
te <- select(te, c(log_prop_occupied, abs_lat_mp, realm, 
                               dispersal_distance_continuous, 
                               log_maximum_body_size, rr_cells, metric,
                               Class, Order, Family, Genus, Species, range, bias_in_uf))


# get complete cases
te <- subset(te, complete.cases(te))
dim(te)

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
for (i in c(2,4,5)) {
  ufnorm[,i] <- (uf[,i]-mean(c(uf[,i]), na.rm = T))/sd(uf[,i], na.rm = T)
  means[x] = mean(c(uf[,i]))
  sds[x] = sd(c(uf[,i]))
  x = x+1
}
dim(ufnorm)
## save scaling factors:
scalers = data.frame(var = colnames(uf)[c(2,4,5)], means = means, sds = sds)

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
                              "Abs. realized range latitudinal midpoint", 
                              "Reference"))

saveRDS(whisker_uf, "data-processed/intermediate-files/whisker_range_no-dormants.rds")

ggsave(whisker_uf, width = 5, height = 3.5, path = "figures/extended-data/", 
       filename = "whisker-plot_range_no-dormants.png", 
       device = "png")


## bias in underfilling:
#Second, rescale variables and get complete cases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
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
bptest(model_uf)

#Now, fit all the models using dredge, select the confidence set and average the models:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#fit all models
suppressWarnings(allmodels_asym <- dredge(model_asym, extra="R^2", subset = (`abs_lat_mp:realm`)))
confset.95p_asym <- get.models(allmodels_asym, subset = cumsum(weight)<=.95) #get confidence set
avgm_asym <- model.avg(confset.95p_asym) #do averaging
summary(confset.95p_asym[[1]])
avgm_asym <- confset.95p_asym[[1]]

sum <- summary(avgm_asym)
df <- as.data.frame(sum$coefficients$fixed) 
intervals <- as.data.frame(intervals(confset.95p_asym[[1]], which = "fixed")$fixed)
df$CI.min <- intervals$lower
df$CI.max <- intervals$upper  # get confidence intervals for model
data.table::setDT(df, keep.rownames = "coefficient") #put rownames into column
df$how_conf = ifelse((df$CI.min < 0 & df$CI.max < 0) | (df$CI.min > 0 & df$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0

whisker_asym <- ggplot(data=df, aes(x=fct_rev(coefficient), y=sum$coefficients$fixed)) + 
  geom_hline(yintercept = 0) + 
  theme_classic() +
  geom_pointrange(aes(ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = '#b0b0b0', size = 1, shape = 21) + 
  coord_flip() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("white", "#b0b0b0")) +
  labs(y = "Effect of variable on asymmetry in underfilling", x = "") +
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal", 
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "Abs. realized range latitudinal midpoint", 
                              "Reference"))

saveRDS(whisker_asym, "data-processed/intermediate-files/whisker_asym_no-dormants.rds")

ggsave(whisker_asym, width = 5, height = 3.5, path = "figures/extended-data/", 
       filename = "whisker-plot_asym_no-dormants.png", 
       device = "png")

## combine all whisker plots and save figure:
whisker_uf <- readRDS("data-processed/intermediate-files/whisker_range_no-dormants.rds")
whisker_warm <- readRDS("data-processed/intermediate-files/whisker_warm_no-dormants.rds")
whisker_cold <- readRDS("data-processed/intermediate-files/whisker_cold_no-dormants.rds")

whisker_cold <- whisker_cold  +
  scale_y_continuous(limits = c(-19,19)) +
  labs(y = "Effect of variable on cool niche filling")

whisker_warm <- whisker_warm + 
  scale_y_continuous(limits = c(-18,18)) 

whisker_uf <- whisker_uf + scale_x_discrete(labels = c("Realm: subtidal",
                                                       "Realm: intertidal",
                                                       "Realized range size (log no. cells)", 
                                                       "Maximum body size (log cm)",
                                                       "Dispersal distance (km)",
                                                       "x realm: intertidal Abs. realized range latitudinal midpoint", 
                                                       "Reference")) +
  scale_y_continuous(limits = c(-10,10))

whisker_asym <- whisker_asym + 
  scale_y_continuous(limits = c(-1,1)) + 
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal", 
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "x realm: intertidal Abs. realized range latitudinal midpoint", 
                              "Reference"))

whisker <- ggdraw() + 
  draw_plot(whisker_warm, x = 0, y = 0.72, width = 1, height = 0.28) +
  draw_plot(whisker_cold, x = 0, y = 0.46, width = 1, height = 0.28) + 
  draw_plot(whisker_uf, x = 0, y = 0.22, width = 1, height = 0.22) + 
  draw_plot(whisker_asym, x = 0, y = 0, width = 1, height = 0.22)

ggsave(whisker, width = 6, height = 5.5, path = "figures/extended-data", 
       filename = "whisker-plot_no-dormants.png", 
       device = "png")

ggsave(whisker_warm, width = 6, height = 3.5, path = "figures/extended-data", 
       filename = "whisker-plot_warm-niche_no-dormants.png", 
       device = "png")

ggsave(whisker_cold, width = 6, height = 3.5, path = "figures/extended-data", 
       filename = "whisker-plot_cold-niche_no-dormants.png", 
       device = "png")

ggsave(whisker_uf, width = 5, height = 4, path = "figures/extended-data", 
       filename = "whisker-plot_range-filling_no-dormants.png", 
       device = "png")