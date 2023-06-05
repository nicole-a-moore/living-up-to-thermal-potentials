## running analysis to see if results are sensitive to use of PGLS versus LME with nested taxonomy 
library(tidyverse)
library(ape)
library(MCMCglmm)
library(RRphylo)
library(nlme)
library(geiger)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####     prep phylogenetic tree          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in phylogeny from TimeTree
tree <- read.tree("data-raw/phylo/TimeTreePhylo.nwk")
plot(tree)

is.ultrametric(tree)
is.binary(tree) ## binary - so only 2 branches per node 
sum(tree$edge.length < 0) ## no negative branch lengths
sum(tree$edge.length == 0) ## 18 zero branch lengths

## add small, small quantity to all branches so that none are zero
tree$edge.length <- tree$edge.length + 1e-10


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## fit phylogenetic generalized least squares ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use gls() from nlme package
## argument "correlation" allows you to specifiy a variance-covariance matrix derived from a phylogeny

## since not all species have TimeTree information, we will have to run analyses on data subsets 
## so, compare pgls output to lme output analyzed on the data subset to ensure changes in results are not due to data subsetting 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##    warm niche filling     ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in data to model
warm <- readRDS("data-processed/thermal-niches/niche-filling/warm-complete-cases.rds")
warm$genus_species <- paste(warm$Genus, warm$Species, sep = "_")

## subset to species in phylogeny
length(unique(warm$genus_species)) ## 382 spp
warm <- filter(warm, genus_species %in% tree$tip.label)
length(unique(warm$genus_species)) ## 311 spp

## remove species in tree that are not in data
row.names(warm) <- warm$genus_species
nc <- name.check(tree, warm)
nc ## 66 spp

tree_warm <- drop.tip(tree, nc$tree_not_data)
#plot(tree_warm)

nc <- name.check(tree_warm, warm)
nc ## OK

## create corClass (correlation structure class) object 
## use corPagel 
corr = corPagel(value = 1, 
            phy = tree_warm, 
            form = ~genus_species)

corr = Initialize(corr, data = warm) # initialize object

## extract the correlation matrix and inspect 
matrix <- matrix(corMatrix(corr),
                           nrow = length(unique(warm$genus_species)))
#View(matrix)
## looks good!
## make sure diagonal = 1
diag(matrix) == 1

## now run gls with correlation matrix
## same fixed effect structure as lme 
gls_warm = gls(filling_value ~ abs_lat_mp*realm + log_range_area +
      dispersal_distance_continuous  +
      log_maximum_body_size + metric,
    
      correlation = corr, 
      
    data = warm)

## get summary
summary(gls_warm)

#check residuals
hist(resid(gls_warm))
plot(gls_warm)

## plot em!
sum <- summary(gls_warm)
coefs <- sum$coefficients
df_warm <- as.data.frame(sum$coefficients) 
df_warm <- rename(df_warm, "coefs" = "sum$coefficients")
CI <- as.data.frame(confint(gls_warm, full=T, level = 0.95)) # get confidence intervals for full model
df_warm$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df_warm$CI.max <-CI$`97.5 %`
data.table::setDT(df_warm, keep.rownames = "coefficient") #put rownames into column
df_warm$how_conf = ifelse((df_warm$CI.min < 0 & df_warm$CI.max < 0) | (df_warm$CI.min > 0 & df_warm$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0

## now lme warm:
lme_warm = lme(filling_value ~ abs_lat_mp*realm + log_range_area +
                 dispersal_distance_continuous  +
                 log_maximum_body_size + metric,
               
               random = ~1|Class/Order/Family/Genus,
               
               data = warm)

## get summary
summary(lme_warm)

#check residuals
hist(resid(lme_warm))
plot(lme_warm)

## plot em!
sum <- summary(lme_warm)
coefs <- sum$coefficients$fixed
df2_warm <- as.data.frame(sum$coefficients$fixed) 
df2_warm <- rename(df2_warm, "coefs" = "sum$coefficients$fixed")
ints <- intervals(lme_warm)[[1]]
df2_warm$CI.min <- ints[,1]
df2_warm$CI.max <- ints[,3]
data.table::setDT(df2_warm, keep.rownames = "coefficient") #put rownames into column
df2_warm$how_conf = ifelse((df2_warm$CI.min < 0 & df2_warm$CI.max < 0) | (df2_warm$CI.min > 0 & df2_warm$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0

warm_plot <- ggplot(data=df2_warm, aes(x=fct_rev(coefficient), y=coefs)) + 
  geom_hline(yintercept = 0) + 
  geom_pointrange(data=df_warm, inherit.aes = FALSE,
                  aes(x=fct_rev(coefficient), y=coefs, ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = 'black', size = 1, shape = 21,
                  position = position_nudge(x = 0.15))+
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
                              "Reference"), drop = FALSE) +
  scale_y_continuous(limits = c(-20.5,24))

saveRDS(warm_plot, "data-processed/intermediate-files/PGLS_warm-plot.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##    cold niche filling     ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in data to model
cold <- readRDS("data-processed/thermal-niches/niche-filling/cold-complete-cases.rds")
cold$genus_species <- paste(cold$Genus, cold$Species, sep = "_")

## subset to species in phylogeny
length(unique(cold$genus_species)) ## 227 spp
cold <- filter(cold, genus_species %in% tree$tip.label)
length(unique(cold$genus_species)) ## 201 spp

## remove species in tree that are not in data
row.names(cold) <- cold$genus_species
nc <- name.check(tree, cold)
nc ## 176 spp

tree_cold <- drop.tip(tree, nc$tree_not_data)
#plot(tree_cold)

nc <- name.check(tree_cold, cold)
nc ## OK

## create corClass (correlation structure class) object 
## use corPagel
corr = corPagel(value = 1, 
                   phy = tree_cold, 
                   form = ~genus_species)

corr = Initialize(corr, data = cold) # initialize object

## extract the correlation matrix and inspect 
matrix <- matrix(corMatrix(corr),
                 nrow = length(unique(cold$genus_species)))
#View(matrix)
## looks good!
## make sure diagonal = 1
diag(matrix) == 1

## now run gls with correlation matrix
## same fixed effect structure as lme 
gls_cold = gls(filling_value ~ abs_lat_mp*realm + log_range_area +
                 dispersal_distance_continuous  +
                 log_maximum_body_size + metric,
               
               correlation = corr, 
               
               data = cold)

## get summary
summary(gls_cold)

#check residuals
hist(resid(gls_cold))
plot(gls_cold)

## plot em!
sum <- summary(gls_cold)
coefs <- sum$coefficients
df_cold <- as.data.frame(sum$coefficients) 
df_cold <- rename(df_cold, "coefs" = "sum$coefficients")
CI <- as.data.frame(confint(gls_cold, full=T, level = 0.95)) # get confidence intervals for full model
df_cold$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df_cold$CI.max <-CI$`97.5 %`
data.table::setDT(df_cold, keep.rownames = "coefficient") #put rownames into column
df_cold$how_conf = ifelse((df_cold$CI.min < 0 & df_cold$CI.max < 0) | (df_cold$CI.min > 0 & df_cold$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0

## now lme cold:
lme_cold = lme(filling_value ~ abs_lat_mp*realm + log_range_area +
                 dispersal_distance_continuous  +
                 log_maximum_body_size + metric,
               
               random = ~1|Class/Order/Family/Genus,
               
               data = cold)

## get summary
summary(lme_cold)

#check residuals
hist(resid(lme_cold))
plot(lme_cold)

## plot em!
sum <- summary(lme_cold)
coefs <- sum$coefficients$fixed
df2_cold <- as.data.frame(sum$coefficients$fixed) 
df2_cold <- rename(df2_cold, "coefs" = "sum$coefficients$fixed")
ints <- intervals(lme_cold, which = "fixed")[[1]]
df2_cold$CI.min <- ints[,1]
df2_cold$CI.max <- ints[,3]
data.table::setDT(df2_cold, keep.rownames = "coefficient") #put rownames into column
df2_cold$how_conf = ifelse((df2_cold$CI.min < 0 & df2_cold$CI.max < 0) | (df2_cold$CI.min > 0 & df2_cold$CI.max > 0), 
                      "so_conf", "not")# add column showing which predictor CIs don't cross 0

cold_plot <- ggplot(data=df2_cold, aes(x=fct_rev(coefficient), y=coefs)) + 
  geom_hline(yintercept = 0) + 
  geom_pointrange(data=df_cold, inherit.aes = FALSE,
                  aes(x=fct_rev(coefficient), y=coefs, ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = 'black', size = 1, shape = 21,
                  position = position_nudge(x = 0.15)) +
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
  scale_y_continuous(limits = c(-20.5,30.1))

saveRDS(cold_plot, "data-processed/intermediate-files/PGLS_cold-plot.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         range  filling          #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in data to model
rf <- read.csv("data-processed/potential-ranges/range-filling/range-complete-cases.csv")
rf$genus_species <- paste(rf$Genus, rf$Species, sep = "_")

## subset to species in phylogeny
length(unique(rf$genus_species)) ## 156 spp
rf <- filter(rf, genus_species %in% tree$tip.label)
length(unique(rf$genus_species)) ## 137 spp

## remove species in tree that are not in data
row.names(rf) <- rf$genus_species
nc <- name.check(tree, rf)
nc ## 240 spp

tree_rf <- drop.tip(tree, nc$tree_not_data)
#plot(tree_rf)

nc <- name.check(tree_rf, rf)
nc ## OK

## create corClass (correlation structure class) object 
## use corPagel
corr = corPagel(value = 1, 
                phy = tree_rf, 
                form = ~genus_species)

corr = Initialize(corr, data = rf) # initialize object

## extract the correlation matrix and inspect 
matrix <- matrix(corMatrix(corr),
                 nrow = length(unique(rf$genus_species)))
#View(matrix)
## looks good!
## make sure diagonal = 1
diag(matrix) == 1

## now run gls with correlation matrix
## same fixed effect structure as lme 
gls_rf = gls(log_prop_occupied ~ abs_lat_mp*realm + dispersal_distance_continuous +
               log_maximum_body_size,
               
               correlation = corr, 
               
               data = rf)

## get summary
summary(gls_rf)

#check residuals
hist(resid(gls_rf))
plot(gls_rf)

## plot em!
sum <- summary(gls_rf)
coefs = sum$coefficients
df_rf <- as.data.frame(sum$coefficients) 
df_rf <- rename(df_rf, "coefs" = "sum$coefficients")
CI <- as.data.frame(confint(gls_rf, full=T, level = 0.95)) # get confidence intervals for full model
df_rf$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df_rf$CI.max <-CI$`97.5 %`
data.table::setDT(df_rf, keep.rownames = "coefficient") #put rownames into column
df_rf$how_conf = ifelse((df_rf$CI.min < 0 & df_rf$CI.max < 0) | (df_rf$CI.min > 0 & df_rf$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0
## now lme rf:
lme_rf = lme(log_prop_occupied ~ abs_lat_mp*realm + dispersal_distance_continuous +
               log_maximum_body_size,
               
               random = ~1|Class/Order/Family/Genus,
               
               data = rf)

## get summary
summary(lme_rf)

#check residuals
hist(resid(lme_rf))
plot(lme_rf)

## plot em!
sum <- summary(lme_rf)
coefs <- sum$coefficients$fixed
df2_rf <- as.data.frame(sum$coefficients$fixed) 
df2_rf <- rename(df2_rf, "coefs" = "sum$coefficients$fixed")
ints <- intervals(lme_rf, which = "fixed")[[1]]
df2_rf$CI.min <- ints[,1]
df2_rf$CI.max <- ints[,3]
data.table::setDT(df2_rf, keep.rownames = "coefficient") #put rownames into column
df2_rf$how_conf = ifelse((df2_rf$CI.min < 0 & df2_rf$CI.max < 0) | (df2_rf$CI.min > 0 & df2_rf$CI.max > 0), 
                      "so_conf", "not")# add column showing which predictor CIs don't cross 0

rf_plot <- ggplot(data=df2_rf, aes(x=fct_rev(coefficient), y=coefs)) + 
  geom_hline(yintercept = 0) + 
  geom_pointrange(data=df_rf, 
                  aes(x=fct_rev(coefficient), y=coefs, ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = 'black', size = 1, shape = 21,
                  position = position_nudge(x = 0.15)) +
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
  scale_y_continuous(limits = c(-13,13))

saveRDS(rf_plot, "data-processed/intermediate-files/PGLS_rf-plot.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####    bias in range underfilling   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## now run gls with correlation matrix
## same fixed effect structure as lme 
gls_uf = gls(bias_in_uf ~ abs_lat_mp*realm + dispersal_distance_continuous +
               log_maximum_body_size,
             
             correlation = corr, 
             
             data = rf)

## get summary
summary(gls_uf)

#check residuals
hist(resid(gls_uf))
plot(gls_uf)

## plot em!
sum <- summary(gls_uf)
coefs <- sum$coefficients
df_uf <- as.data.frame(sum$coefficients) 
CI <- as.data.frame(confint(gls_uf, full=T, level = 0.95)) # get confidence intervals for full model
df_uf$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df_uf$CI.max <-CI$`97.5 %`
data.table::setDT(df_uf, keep.rownames = "coefficient") #put rownames into column
df_uf$how_conf = ifelse((df_uf$CI.min < 0 & df_uf$CI.max < 0) | (df_uf$CI.min > 0 &  df_uf$CI.max > 0), 
                     "so_conf", "not")# add column showing which predictor CIs don't cross 0

## now lme for bias:
lme_uf = lme(bias_in_uf ~ abs_lat_mp*realm + dispersal_distance_continuous +
               log_maximum_body_size,
             
             random = ~1|Class/Order/Family/Genus,
             
             data = rf)

## get summary
summary(lme_uf)

#check residuals
hist(resid(lme_uf))
plot(lme_uf)

## plot em!
sum <- summary(lme_uf)
coefs <- sum$coefficients$fixed
df2_uf <- as.data.frame(sum$coefficients$fixed) 
ints <- intervals(lme_uf, which = "fixed")[[1]]
df2_uf$CI.min <- ints[,1]
df2_uf$CI.max <- ints[,3]
data.table::setDT(df2_uf, keep.rownames = "coefficient") #put rownames into column
df2_uf$how_conf = ifelse((df2_uf$CI.min < 0 & df2_uf$CI.max < 0) | (df2_uf$CI.min > 0 & df2_uf$CI.max > 0), 
                      "so_conf", "not")# add column showing which predictor CIs don't cross 0

uf_plot <- ggplot(data=df2_uf, aes(x=fct_rev(coefficient), y=coefs)) + 
  geom_hline(yintercept = 0) + 
  geom_pointrange(data=df_uf, 
                  aes(x=fct_rev(coefficient), y=coefs, ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = 'black', size = 1, shape = 21,
                  position = position_nudge(x = 0.15)) +
  theme_classic() +
  geom_pointrange(aes(ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = '#b0b0b0', size = 1, shape = 21) + 
  coord_flip() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("white", "#b0b0b0")) +
  labs(y = "Effect of variable on equatorward bias in underfilling", x = "") +
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal",
                              "Maximum body size (log cm)",
                              "Dispersal distance (km)",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Abs. realized range latitudinal midpoint ",
                              "Reference")) +
  scale_y_continuous(limits = c(-2.1,2.1)) 

saveRDS(uf_plot, "data-processed/intermediate-files/PGLS_uf-plot.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####   combine the figures   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

warm_plot <- readRDS("data-processed/intermediate-files/PGLS_warm-plot.rds")

composite <- ggdraw() + 
  draw_plot(warm_plot, x = 0, y = 0.75, width = 1, height = 0.25) +
  draw_plot(cold_plot, x = 0, y = 0.5, width = 1, height = 0.25) + 
  draw_plot(rf_plot, x = 0, y = 0.25, width = 1, height = 0.25) +
  draw_plot(uf_plot, x = 0, y = 0, width = 1, height = 0.25) +
  draw_plot_label(label = c("a)", "b)", "c)", "d)"),
                  x = c(0, 0, 0, 0),
                  y = c(1, 0.75, 0.5, 0.25),
                  size = 10, 
                  color = "grey30") 


ggsave(composite, height = 7, width = 7, path = "figures/extended-data", 
       filename = "pgls-vs-lme.png", 
       device = "png")

