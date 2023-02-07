## this script makes a bunch of figures to visualize the niche filling data
library(tidyverse)
library(cowplot)
library(lme4)
library(nlme)

thermal =thermal %>%
  select(genus_species, Genus, Species, realm, Phylum, Class, Order, Family, type, 
         thermal_limit, metric, 
         collection_latitude, collection_longitude, 
         reference)

write.csv(thermal, "files for gabriel/species metadata/species_information.csv", row.names = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       Prepping data for figure making     ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
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

Te <- types[[2]]
Te_tpreftb <- types[[4]]
Te_subset_tpreftb <- Te[which(Te$species %in% Te_tpreftb$species),]
Te_acclimatized <- types[[1]]
Te_subset_acclimatized <- filter(Te, Te$species %in% Te_acclimatized$species)

types = list(Te, Te_tpreftb,Te_subset_tpreftb,Te_acclimatized,
             Te_subset_acclimatized)
names(types) = c("Te", "Te_tpreftb","Te_subset_tpreftb", "Te_acclimatized", 
                 "Te_subset_acclimatized")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####              Making figures               ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## make loop to plot each type:
type = 1 
while (type < length(types) + 1) {
  folder = names(types)[type]
  path = paste("figures/additional-figures/niche-filling/", folder, sep = "")
  data <- types[[type]]
  
  ## prep dataset for figure making:
  fig <- data %>%
    filter(!is.infinite(r_niche_upper),
           !is.infinite(r_niche_lower),
           !is.infinite(p_niche_upper),
           !is.infinite(r_niche_lower)) %>%
    mutate(upper_inner = ifelse(warm_under < 0, "r_niche_upper", "p_niche_upper")) %>%
    mutate(lower_inner = ifelse(cold_under < 0, "r_niche_lower", "p_niche_lower")) %>%
    gather(key = "niche_limit_type", value = "limit_value", c(r_niche_upper, r_niche_lower, p_niche_upper, 
                                                              p_niche_lower, ctmax, ctmin))%>%
    filter(!is.na(limit_value)) %>%
    mutate(r_p_c = ifelse(str_detect(niche_limit_type, "r_niche"), "r_niche", ifelse(str_detect
                                                                                     (niche_limit_type, 
                                                                                      "p_niche"),"p_niche",
                                                                                     "ctlim"))) %>%
    mutate(w_c = ifelse(str_detect(niche_limit_type, "upper") | str_detect(niche_limit_type, "max"), 
                        "warm", "cold")) %>%
    mutate(line_group1 = ifelse((cold_under < 0) & (w_c == "cold") & (r_p_c != 'ctlim'), 'cold_under', 
                                ifelse((cold_over > 0) & (w_c == "cold") & (r_p_c != 'ctlim'), 'cold_over',
                                       ifelse((warm_under < 0) & (w_c == "warm") & (r_p_c != 'ctlim'),
                                              'warm_under',
                                              ifelse((warm_over > 0) & (w_c == "warm") & 
                                                       (r_p_c != 'ctlim'),
                                                     "warm_over", NA))))) %>%
    mutate(line_type = ifelse(str_detect(line_group1, "over") & (r_p_c != "ctlim"), 
                              "Overfilling", ifelse(str_detect(line_group1, "under") & (r_p_c != "ctlim"), 
                                                    "Underfilling", NA))) %>%
    mutate(line_group2 = paste(range, line_group1)) %>%
    mutate(thin_line1 = ifelse((upper_inner == "r_niche_upper") & (niche_limit_type == "r_niche_upper"),"inner",
                               ifelse((upper_inner == "p_niche_upper") & (niche_limit_type == "p_niche_upper"),
                                      "inner", ifelse((lower_inner == "p_niche_lower") & 
                                                        (niche_limit_type == "p_niche_lower"), "inner", 
                                                      ifelse((lower_inner == "r_niche_lower") & 
                                                               (niche_limit_type == "r_niche_lower"), 
                                                             "inner", NA))))) %>%
    mutate(thin_line2 = paste(range, thin_line1)) 
  
  fig_filling <- fig %>%
    filter(!is.na(thin_line1))
  
  fig_overfilling <- fig %>%
    filter(!is.na(line_group1)) %>%
    filter(line_type == "Overfilling") 
  
  fig_underfilling <- fig %>%
    filter(!is.na(line_group1)) %>%
    filter(line_type == "Underfilling") 
  
  ## plot filling across realms:
  beauty <- ggplot(data = fig, aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
    scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
    guides(color = FALSE) +
    geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5) +
    geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
    geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
    scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                      "Realized")) + 
    labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
         x = "Realized range latitudinal midpoint (°N)") +
    geom_point(fill = "white") + 
    theme_bw() +
    scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) + 
    theme(panel.grid = element_blank())
  
  ggsave(beauty, path = path, filename = "niche-filling-across-lat.png", 
         device = "png", width = 11, height = 6)
  
  ## plot filling across realms:
 beauty_talk <- fig %>%
    filter(!r_p_c %in% c("ctlim")) %>%
    filter(niche_limit_type %in% c("p_niche_upper", "p_niche_lower")) %>%
    ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
    scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
    guides(color = FALSE) +
    #geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5, colour = "transparent") +
    #geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
    #geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
    scale_shape_manual(values=c(21, 19), labels = c("Potential",
                                                       "Realized")) + 
    labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
         x = "Realized range latitudinal midpoint (°N)") +
    geom_point(fill = "white") + 
    theme_bw() +
    scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) + 
    theme(panel.grid = element_blank()) +
    scale_y_continuous(limits = c(-45, 50))
  
  ggsave(beauty_talk, path = path, filename = "niche-filling-across-lat_potential.png", 
         device = "png", width = 8, height = 4)
  
  beauty_abs <- ggplot(data = fig, aes(x = abs(lat_mp), y = limit_value, col = w_c, shape = r_p_c)) +
    scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
    guides(color = FALSE) +
    geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5) +
    geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
    geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
    scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                      "Realized")) + 
    labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
         x = "Realized range latitudinal midpoint (°N/S)") +
    geom_point(fill = "white") + 
    theme_bw() +
    scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(0, 66)) + 
    theme(panel.grid = element_blank())
  
  ggsave(beauty_abs, path = path, filename = "niche-filling-across-abs-lat.png", 
         device = "png", width = 11, height = 6)

  if(names(types)[type] %in% c("Te_acclimatized", "Te")) {
    no_leg <- ggplot(data = fig, aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5) +
      geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
      geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                        "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
           x = "Realized range latitudinal midpoint (°N)") +
      geom_point(fill = "white") + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) + 
      theme(panel.grid = element_blank())
    
    fig_filling_t <- fig_filling %>%
      filter(realm == "Terrestrial") 
    
    fig_overfilling_t <- fig_overfilling %>%
      filter(realm == "Terrestrial")
    
    fig_underfilling_t <- fig_underfilling %>%
      filter(realm == "Terrestrial")
    
    terr <- fig %>%
      filter(realm == "Terrestrial") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_t, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.4) +
      geom_line(data = fig_overfilling_t, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.4) +
      geom_line(data = fig_filling_t, aes(group = thin_line2), col = "grey80", size = 0.5*0.4) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                        "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.75, fill = "white", stroke = 0.25) + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57)) + 
      theme(panel.grid = element_blank())
    
    fig_filling_m <- fig_filling %>%
      filter(realm == "Marine") 
    
    fig_overfilling_m <- fig_overfilling %>%
      filter(realm == "Marine")
    
    fig_underfilling_m <- fig_underfilling %>%
      filter(realm == "Marine")
    
    marine <- fig %>%
      filter(realm == "Marine") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_m, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.4) +
      geom_line(data = fig_overfilling_m, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.4) +
      geom_line(data = fig_filling_m, aes(group = thin_line2), col = "grey80", size = 0.5*0.4) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                        "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.75, fill = "white", stroke = 0.25) + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57))+ 
      theme(panel.grid = element_blank())
    
    fig_filling_i <- fig_filling %>%
      filter(realm == "Intertidal") 
    
    fig_overfilling_i <- fig_overfilling %>%
      filter(realm == "Intertidal")
    
    fig_underfilling_i <- fig_underfilling %>%
      filter(realm == "Intertidal")
    
    int <- fig %>%
      filter(realm == "Intertidal") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_i, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.4) +
      geom_line(data = fig_overfilling_i, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.4) +
      geom_line(data = fig_filling_i, aes(group = thin_line2), col = "grey80", size = 0.5*0.4) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                        "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.75, fill = "white", stroke = 0.25) + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57)) + 
      theme(panel.grid = element_blank())
    
    ## ## ## ## ## ## ## ## ## ## ## 
    ## make a didactic legend:
    legend <- get_legend(beauty) 
    
    which(unique(fig$range) =='Liolaemus_kingii_IUCN')
    
    didac_u <- filter(fig_underfilling, range == unique(fig$range)[211]) 
    didac_o <- filter(fig_overfilling, range == unique(fig$range)[211])
    didac_f <- filter(fig_filling, range == unique(fig$range)[211])
    
    label <- fig %>%
      filter(range == unique(fig$range)[211]) 
    label = label[1:4,] %>%
      mutate(label_y = c(41.54995, 33.4209, 7.010000, 0.005046)) %>%
      mutate(label = c("unavailable warm niche space", 'warm underfilling', 
                       'no unavailable cold niche space', 'cold overfilling')) 
    
    didactic <- fig %>%
      filter(range == unique(fig$range)[211]) %>% # select a good example of warm underfilling, cold overfilling 
      ggplot(., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) +
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = didac_u, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.812) +
      geom_line(data = didac_o, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1.2) +
      geom_line(data = didac_f, aes(group = thin_line2), col = "grey80", size = 0.5) + 
      scale_shape_manual(values=c(3, 21, 19)) + 
      geom_point(size = 4, stroke = 1, fill = "white") + 
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) + theme_void() +
      geom_text(data = label, inherit.aes = F, aes(label = label, x = lat_mp +8, y = label_y), 
                hjust = 0, vjust=0.5, size = 3) 
    
    beauty_with_didactic <- ggdraw() +
      draw_plot(no_leg, 0, 0, 0.6, 1) +
      draw_plot(marine, 0.6, 2/3, 0.2, 1/3) +
      draw_plot(int, 0.6, 1/3, 0.2, 1/3) +
      draw_plot(terr, 0.6, 0, 0.2, 1/3) +
      draw_plot(legend, 0.78, 1/2, 0.2, 1/2) +
      draw_plot(didactic, 0.79, 0.1, 0.2, 0.4) +
      draw_plot_label(label = c("a)", "b)", "c)", "d)"),
                      x = c(0, 0.6, 0.6, 0.6),
                      y = c(0.99, 0.99, (2/3-0.01), (1/3-0.01)), size = 10, 
                      color = "grey30") +
      theme(plot.background = element_rect(fill="white", color = "white"))
    
    ggsave(beauty_with_didactic, path = path, 
           filename = "niche-filling-across-lat-and-realm-with-didactic.png",
           device = "png", width = 10, height = 5)
      
    terr <- fig %>%
      filter(realm == "Terrestrial") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_t, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.3) +
      geom_line(data = fig_overfilling_t, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.2) +
      geom_line(data = fig_filling_t, aes(group = thin_line2), col = "grey80", size = 0.5*0.2) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                         "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.45, fill = "white", stroke = 0.25) + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57)) + 
      theme(panel.grid = element_blank())
    
    int <- fig %>%
      filter(realm == "Intertidal") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_i, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.3) +
      geom_line(data = fig_overfilling_i, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.2) +
      geom_line(data = fig_filling_i, aes(group = thin_line2), col = "grey80", size = 0.5*0.2) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                         "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.45, fill = "white", stroke = 0.25) + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57)) + 
      theme(panel.grid = element_blank())
    
    marine <- fig %>%
      filter(realm == "Marine") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_m, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.3) +
      geom_line(data = fig_overfilling_m, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.2) +
      geom_line(data = fig_filling_m, aes(group = thin_line2), col = "grey80", size = 0.5*0.2) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                         "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.45, fill = "white", stroke = 0.25) + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57))+ 
      theme(panel.grid = element_blank())
    
    ggsave(int, path = path, 
           filename = "niche-filling-int.png",
           device = "png", width = 2, height = 1.7)
    ggsave(terr, path = path, 
           filename = "niche-filling-ter.png",
           device = "png", width = 2, height = 1.7)
    ggsave(marine, path = path, 
           filename = "niche-filling-mar.png",
           device = "png", width = 2, height = 1.7)
    
    ## ## ## ## ## ## ## ## ## ## ## 
    
    beauty_comb <- ggdraw() +
      draw_plot(no_leg, 0, 0, 0.75, 1) +
      draw_plot(marine, 0.75, 2/3, 0.25, 1/3) +
      draw_plot(int, 0.75, 1/3, 0.25, 1/3) +
      draw_plot(terr, 0.75, 0, 0.25, 1/3) +
      draw_plot_label(label = c("a)", "b)", "c)", "d)"),
                       x = c(0, 0.75, 0.75, 0.75),
                      y = c(0.99, 0.99, (2/3-0.01), (1/3-0.01)), size = 10, 
                      color = "grey30")
    
    if(type == 1) {
      ggsave(beauty_comb, path = "figures/main/",
             filename = "niche-filling-across-lat-and-realm.png", 
             device = "png", width = 8, height = 4.5)
    }
    else {
      ggsave(beauty_comb, path = path, 
           filename = "niche-filling-across-lat-and-realm.png", 
           device = "png", width = 8, height = 4.5)
    }
    
  }
  
  ## distribution of filling vals:
  filling <- data %>%
    mutate(perfect_warm = ifelse(warm_over == 0 & warm_under == 0, 0, NA)) %>%
    mutate(perfect_cold = ifelse(cold_over == 0 & cold_under == 0, 0, NA)) %>%
    gather(key = "filling_type", value = "filling_value", c(warm_under, warm_over, cold_under, 
                                                            cold_over, perfect_warm, perfect_cold)) %>%
    filter(filling_value != 0 | filling_value == 0 & (filling_type %in% c("perfect_cold","perfect_warm"))) %>%
    select(range, species, source, realm, dormancy, filling_type, filling_value, lat_mp) %>%
    mutate(warm_cold= ifelse(str_detect(filling_type, "warm"), "warm", "cold")) 
  
  perf_fill <- filling %>%
    filter(filling_value == 0)
  
  warm_fig <- filling %>%
    filter(warm_cold == "warm") %>%
    filter(!is.infinite(filling_value)) %>%
    ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + 
    geom_segment(aes(x = abs(lat_mp),  xend = abs(lat_mp), y = 0,yend = filling_value)) +
    scale_color_manual(values = c('#b45346')) + theme_bw() +
    guides(colour = FALSE) + 
    labs(y = "", x = "Absolute realized range latitudinal midpoint (°N)") + 
    scale_y_continuous(limits = c(-24, 45)) + 
    theme(panel.grid = element_blank())
  #+geom_point(data = perf_fill, shape = 5)
  
  cold_fig <- filling %>%
    filter(warm_cold == "cold") %>%
    filter(!is.infinite(filling_value)) %>%
    ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + 
    geom_segment(aes(x = abs(lat_mp), xend = abs(lat_mp), y = 0, yend = filling_value)) +
    scale_color_manual(values = c("steelblue")) + theme_bw()  +
    guides(colour = FALSE) +
    labs(y = "Shortfall or excess of temperatures occupied (°C)", 
         x = "Absolute realized range latitudinal midpoint (°N)") + 
    scale_y_continuous(limits = c(-24, 45)) + 
    theme(panel.grid = element_blank())
  
  fig_fill <- ggdraw() +
    draw_plot(cold_fig, 0, 0, 0.5, 1) +
    draw_plot(warm_fig, 0.5, 0, 0.5, 1) +
    draw_plot_label(label = c("a)", "b)"),
                    x = c(0, 0.5),
                    y = c(0.99, 0.99), size = 10, 
                    color = "grey30")
  
  ggsave(fig_fill, path = path, 
         filename = "niche-filling-distribution-cold-warm.png", 
         device = "png", width = 9, height = 5)
  
  
  warm_dots <- filling %>%
    filter(warm_cold == "warm") %>%
    filter(!is.infinite(filling_value)) %>%
    ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
    scale_color_manual(values = c('#b45346')) + theme_bw() +
    guides(colour = FALSE) + 
    labs(y = "", x = "Absolute realized range latitudinal midpoint (°N/S)") + 
    scale_y_continuous(limits = c(-30, 30)) + 
    theme(panel.grid = element_blank())
  
  cold_dots <- filling %>%
    filter(warm_cold == "cold") %>%
    filter(!is.infinite(filling_value)) %>%
    ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
    scale_color_manual(values = c('steelblue')) + theme_bw() +
    guides(colour = FALSE) + 
    labs(x = "Absolute realized range latitudinal midpoint (°N/S)", 
         y = "Shortfall or excess of temperatures occupied (°C)") + 
    scale_y_continuous(limits = c(-40, 40)) + 
    theme(panel.grid = element_blank())
  
  fig_dots <- ggdraw() +
    draw_plot(cold_dots, 0, 0, 0.5, 1) +
    draw_plot(warm_dots, 0.5, 0, 0.5, 1) +
    draw_plot_label(label = c("a)", "b)"),
                    x = c(0, 0.5),
                    y = c(0.99, 0.99), size = 10, 
                    color = "grey30")
  
  ggsave(fig_dots, path = path, 
         filename = "niche-filling-distribution-cold-warm-dots.png", 
         device = "png", width = 8, height = 4.5)
  
  
  ## underfilling fundamental thermal niche:
  fund_warm = data %>%
    mutate(f_filling = ifelse(is.infinite(f_warm_under) | is.na(f_warm_under),
                              NA,
                              ifelse(f_warm_under == 0,
                                     f_warm_over, 
                                     -f_warm_under))) %>%
    ggplot(., aes(y = f_filling, x = abs(lat_mp))) + geom_point(colour = c('#b45346')) + 
    theme_bw() + geom_smooth(method = "lm", colour = c('#b45346')) + 
    labs(x = "Absolute realized range latitudinal midpoint (°N/S)",
         y = "Shortfall/excess of temperatures occupied beyond warm niche limit (°C)")
  
  fund_cold = data %>%
    mutate(f_filling = ifelse(is.infinite(f_cold_under) | is.na(f_cold_under),
                              NA,
                              ifelse(f_cold_under == 0,
                                     f_cold_over, 
                                     -f_cold_under))) %>%
    ggplot(., aes(y = f_filling, x = abs(lat_mp))) + geom_point(colour = c('steelblue')) + 
    theme_bw() + geom_smooth(method = "lm", colour = c('steelblue')) + 
    labs(x = "Absolute realized range latitudinal midpoint (°N/S)",
         y = "Shortfall/excess of temperatures occupied beyond cold niche limit (°C)")
  
  
  ## mean and distribution 
  filling %>%
    ggplot(., aes(y = filling_value, x = realm, colour = warm_cold)) + geom_violin(outlier.shape = NA) + geom_jitter(width = 0.2) + theme_bw() +
    facet_wrap(~warm_cold) +
    scale_color_manual(values = c("steelblue", '#b45346')) +
    guides(colour = F) +
    labs(x = "", y = "Shortfall/excess of temperatures occupied (°C)") +
    stat_summary(fun.y=mean, geom="point", size=2, colour = "black")+ 
    theme(panel.grid = element_blank())

  ## just warm:
  if (names(types)[type] %in% c("Te_acclimatized", "Te_subset_acclimatized")) {
    warm_dots <- filling %>%
      filter(warm_cold == "warm") %>%
      filter(!is.infinite(filling_value)) %>%
      ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
      scale_color_manual(values = c('#b45346')) + theme_bw() +
      guides(colour = FALSE) + 
      labs(y = "Shortfall or excess of temperatures occupied (°C)",
           x = "Absolute realized range latitudinal midpoint (°N/S)") + 
      scale_y_continuous(limits = c(-25, 12)) + 
      theme(panel.grid = element_blank())
    
    cold_dots <- filling %>%
      filter(warm_cold == "cold") %>%
      filter(!is.infinite(filling_value)) %>%
      ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
      scale_color_manual(values = c('steelblue')) + theme_bw() +
      guides(colour = FALSE) + 
      labs(y = "Shortfall or excess of temperatures occupied (°C)",
           x = "Absolute realized range latitudinal midpoint (°N/S)") + 
      scale_y_continuous(limits = c(-30, 30)) + 
      theme(panel.grid = element_blank())
  }
  else {
    warm_dots <- filling %>%
      filter(warm_cold == "warm") %>%
      filter(realm == "Terrestrial") %>%
      filter(!is.infinite(filling_value)) %>%
      ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
      scale_color_manual(values = c('#b45346')) + theme_bw() +
      guides(colour = FALSE) + 
      labs(y = "Shortfall or excess of temperatures occupied (°C)",
           x = "Absolute realized range latitudinal midpoint (°N/S)") + 
      scale_y_continuous(limits = c(-25, 12)) + 
      theme(panel.grid = element_blank())
    
    cold_dots <- filling %>%
      filter(warm_cold == "cold") %>%
      filter(!is.infinite(filling_value)) %>%
      ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
      scale_color_manual(values = c('steelblue')) + theme_bw() +
      guides(colour = FALSE) + 
      labs(y = "Shortfall or excess of temperatures occupied (°C)",
           x = "Absolute realized range latitudinal midpoint (°N/S)") + 
      scale_y_continuous(limits = c(-30, 30)) + 
      theme(panel.grid = element_blank())
    
  }
  ggsave(warm_dots, path = path, 
         filename = "niche-filling-distribution-warm-dots.png", 
         device = "png", width = 4.5, height = 4)
  ggsave(cold_dots, path = path, 
         filename = "niche-filling-distribution-cold-dots.png", 
         device = "png", width = 4.5, height = 4)
  
  ## plot only one thermal limit ones
  fig_one_limit <- fig %>%
    filter(is.na(lower_inner) | is.na(upper_inner))
  
  fig_filling <- fig_one_limit %>%
    filter(!is.na(thin_line1))
  
  fig_overfilling <- fig_one_limit %>%
    filter(!is.na(line_group1)) %>%
    filter(line_type == "Overfilling") 
  
  fig_underfilling <- fig_one_limit %>%
    filter(!is.na(line_group1)) %>%
    filter(line_type == "Underfilling") 
  
  ## plot filling across realms:
  beauty_one_limit <- ggplot(data = fig_one_limit, aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
    scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
    guides(color = FALSE) +
    geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5) +
    geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
    geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
    scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                      "Realized")) + 
    labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
         x = "Realized range latitudinal midpoint (°N)") +
    geom_point(fill = "white") + 
    theme_bw() +
    scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) + 
    theme(panel.grid = element_blank())
  
  ggsave(beauty_one_limit, path = path, filename = "niche-filling-across-lat-one-limit.png", 
         device = "png", width = 11, height = 6)
  
  beauty_abs_one_limit <- ggplot(data = fig_one_limit, aes(x = abs(lat_mp), y = limit_value, col = w_c, shape = r_p_c)) +
    scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
    guides(color = FALSE) +
    geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5) +
    geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
    geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
    scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                      "Realized")) + 
    labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
         x = "Realized range latitudinal midpoint (°N/S)") +
    geom_point(fill = "white") + 
    theme_bw() +
    scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(0, 66)) + 
    theme(panel.grid = element_blank())
  
  ggsave(beauty_abs_one_limit, path = path, filename = "niche-filling-across-abs-lat-one-limit.png", 
         device = "png", width = 11, height = 6)
  
  if(names(types)[type] %in% c("Te_acclimatized", "Te")) {
    no_leg <- ggplot(data = fig_one_limit, aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5) +
      geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
      geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                        "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
           x = "Realized range latitudinal midpoint (°N)") +
      geom_point(fill = "white") + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) + 
      theme(panel.grid = element_blank())
    
    fig_filling_t <- fig_filling %>%
      filter(realm == "Terrestrial") 
    
    fig_overfilling_t <- fig_overfilling %>%
      filter(realm == "Terrestrial")
    
    fig_underfilling_t <- fig_underfilling %>%
      filter(realm == "Terrestrial")
    
    terr <- fig_one_limit %>%
      filter(realm == "Terrestrial") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_t, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.5) +
      geom_line(data = fig_overfilling_t, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.5) +
      geom_line(data = fig_filling_t, aes(group = thin_line2), col = "grey80", size = 0.5*0.5) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                        "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.5, fill = "white") + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57)) + 
      theme(panel.grid = element_blank())
    
    fig_filling_m <- fig_filling %>%
      filter(realm == "Marine") 
    
    fig_overfilling_m <- fig_overfilling %>%
      filter(realm == "Marine")
    
    fig_underfilling_m <- fig_underfilling %>%
      filter(realm == "Marine")
    
    marine <- fig_one_limit %>%
      filter(realm == "Marine") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_m, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.5) +
      geom_line(data = fig_overfilling_m, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.5) +
      geom_line(data = fig_filling_m, aes(group = thin_line2), col = "grey80", size = 0.5*0.5) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                        "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.5, fill = "white") + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57)) + 
      theme(panel.grid = element_blank())
    
    fig_filling_i <- fig_filling %>%
      filter(realm == "Intertidal") 
    
    fig_overfilling_i <- fig_overfilling %>%
      filter(realm == "Intertidal")
    
    fig_underfilling_i <- fig_underfilling %>%
      filter(realm == "Intertidal")
    
    int <- fig_one_limit %>%
      filter(realm == "Intertidal") %>%
      ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
      scale_colour_manual(values = c("steelblue", '#b45346'), labels = c("Cold limit", "Warm limit")) +
      guides(color = FALSE, shape = FALSE, linetype = FALSE) +
      geom_line(data = fig_underfilling_i, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.5) +
      geom_line(data = fig_overfilling_i, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.5) +
      geom_line(data = fig_filling_i, aes(group = thin_line2), col = "grey80", size = 0.5*0.5) + 
      scale_shape_manual(values=c(3, 21, 19), labels = c("Fundamental", "Potential",
                                                        "Realized")) + 
      labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
           x = "") +
      geom_point(size = 0.5, fill = "white") + 
      theme_bw() +
      scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
      scale_y_continuous(limits = c(-45, 57)) + 
      theme(panel.grid = element_blank())
    
    if (type == 1) {
      beauty_comb <- ggdraw() +
        draw_plot(no_leg, 0, 0, 0.75, 1) +
        draw_plot(marine, 0.75, 1/2, 0.25, 0.5) +
        draw_plot(terr, 0.75, 0, 0.25, 0.5) +
        draw_plot_label(label = c("a)", "b)", "c)"),
                        x = c(0, 0.75, 0.75),
                        y = c(0.99, 0.99, (0.5-0.01)), size = 10, 
                        color = "grey30")
    }
    else if(type == 4) {
        beauty_comb <- ggdraw() +
          draw_plot(no_leg, 0, 0, 0.75, 1) +
          draw_plot(marine, 0.75, 2/3, 0.25, 1/3) +
          draw_plot(terr, 0.75, 0, 0.25, 1/3) +
          draw_plot_label(label = c("a)", "b)", "c)", "d)"),
                          x = c(0, 0.75, 0.75, 0.75),
                          y = c(0.99, 0.99, (2/3-0.01), (1/3-0.01)), size = 10, 
                          color = "grey30")
    }
    else {
      beauty_comb <- ggdraw() +
        draw_plot(no_leg, 0, 0, 0.75, 1) +
        draw_plot(marine, 0.75, 2/3, 0.25, 1/3) +
        draw_plot(int, 0.75, 1/3, 0.25, 1/3) +
        draw_plot(terr, 0.75, 0, 0.25, 1/3) +
        draw_plot_label(label = c("a)", "b)", "c)", "d)"),
                        x = c(0, 0.75, 0.75, 0.75),
                        y = c(0.99, 0.99, (2/3-0.01), (1/3-0.01)), size = 10, 
                        color = "grey30")
    }
    
    ggsave(beauty_comb, path = path, 
           filename = "niche-filling-across-lat-and-realm-one-limit.png", 
           device = "png", width = 11, height = 6)
    
  }
  
  ## distribution of filling vals:
  filling <- data %>%
    filter(limit_type %in% c('tmax', 'tmin')) %>%
    mutate(perfect_warm = ifelse(warm_over == 0 & warm_under == 0, 0, NA)) %>%
    mutate(perfect_cold = ifelse(cold_over == 0 & cold_under == 0, 0, NA)) %>%
    gather(key = "filling_type", value = "filling_value", c(warm_under, warm_over, cold_under, 
                                                            cold_over, perfect_warm, perfect_cold)) %>%
    filter(filling_value != 0 | filling_value == 0 & (filling_type %in% c("perfect_cold","perfect_warm"))) %>%
    select(range, species, source, realm, dormancy, filling_type, filling_value, lat_mp) %>%
    mutate(warm_cold= ifelse(str_detect(filling_type, "warm"), "warm", "cold")) 
  
  perf_fill <- filling %>%
    filter(filling_value == 0)
  
  warm_fig <- filling %>%
    filter(warm_cold == "warm") %>%
    filter(!is.infinite(filling_value)) %>%
    ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + 
    geom_segment(aes(x = abs(lat_mp),  xend = abs(lat_mp), y = 0,yend = filling_value)) +
    scale_color_manual(values = c('#b45346')) + theme_bw() +
    guides(colour = FALSE) + 
    labs(y = "", x = "Absolute realized range latitudinal midpoint (°N)") + 
    scale_y_continuous(limits = c(-24, 45)) + 
    theme(panel.grid = element_blank())
  #+geom_point(data = perf_fill, shape = 5)
  
  cold_fig <- filling %>%
    filter(warm_cold == "cold") %>%
    filter(!is.infinite(filling_value)) %>%
    ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + 
    geom_segment(aes(x = abs(lat_mp), xend = abs(lat_mp), y = 0, yend = filling_value)) +
    scale_color_manual(values = c("steelblue")) + theme_bw()  +
    guides(colour = FALSE) +
    labs(y = "Shortfall or excess of temperatures occupied (°C)", 
         x = "Absolute realized range latitudinal midpoint (°N)") + 
    scale_y_continuous(limits = c(-24, 45)) + 
    theme(panel.grid = element_blank())
  
  fig_fill <- ggdraw() +
    draw_plot(cold_fig, 0, 0, 0.5, 1) +
    draw_plot(warm_fig, 0.5, 0, 0.5, 1) +
    draw_plot_label(label = c("a)", "b)"),
                    x = c(0, 0.5),
                    y = c(0.99, 0.99), size = 10, 
                    color = "grey30")
  
  ggsave(fig_fill, path = path, 
         filename = "niche-filling-distribution-cold-warm-one-limit.png", 
         device = "png", width = 9, height = 5)
  
  
  warm_dots <- filling %>%
    filter(warm_cold == "warm") %>%
    filter(!is.infinite(filling_value)) %>%
    ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
    scale_color_manual(values = c('#b45346')) + theme_bw() +
    guides(colour = FALSE) + 
    labs(y = "", x = "Absolute realized range latitudinal midpoint (°N/S)") + 
    scale_y_continuous(limits = c(-30, 30)) + 
    theme(panel.grid = element_blank())
  
  cold_dots <- filling %>%
    filter(warm_cold == "cold") %>%
    filter(!is.infinite(filling_value)) %>%
    ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
    scale_color_manual(values = c('steelblue')) + theme_bw() +
    guides(colour = FALSE) + 
    labs(x = "Absolute realized range latitudinal midpoint (°N/S)", 
         y = "Shortfall or excess of temperatures occupied (°C)") + 
    scale_y_continuous(limits = c(-40, 40)) + 
    theme(panel.grid = element_blank())
  
  fig_dots <- ggdraw() +
    draw_plot(cold_dots, 0, 0, 0.5, 1) +
    draw_plot(warm_dots, 0.5, 0, 0.5, 1) +
    draw_plot_label(label = c("a)", "b)"),
                    x = c(0, 0.5),
                    y = c(0.99, 0.99), size = 10, 
                    color = "grey30")
  
  ggsave(fig_dots, path = path, 
         filename = "niche-filling-distribution-cold-warm-dots-one-limit.png", 
         device = "png", width = 8, height = 4.5)
  
  ## just warm:
  if (names(types)[type] %in% c("Te_acclimatized", "Te_subset_acclimatized")) {
    warm_dots <- filling %>%
      filter(warm_cold == "warm") %>%
      filter(!is.infinite(filling_value)) %>%
      ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
      scale_color_manual(values = c('#b45346')) + theme_bw() +
      guides(colour = FALSE) + 
      labs(y = "Shortfall or excess of temperatures occupied (°C)",
           x = "Absolute realized range latitudinal midpoint (°N/S)") + 
      scale_y_continuous(limits = c(-25, 12)) + 
      theme(panel.grid = element_blank())
    
    cold_dots <- filling %>%
      filter(warm_cold == "cold") %>%
      filter(!is.infinite(filling_value)) %>%
      ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
      scale_color_manual(values = c('steelblue')) + theme_bw() +
      guides(colour = FALSE) + 
      labs(y = "Shortfall or excess of temperatures occupied (°C)",
           x = "Absolute realized range latitudinal midpoint (°N/S)") + 
      scale_y_continuous(limits = c(-30, 30)) + 
      theme(panel.grid = element_blank())
  }
  else {
    warm_dots <- filling %>%
      filter(warm_cold == "warm") %>%
      filter(realm == "Terrestrial") %>%
      filter(!is.infinite(filling_value)) %>%
      ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
      scale_color_manual(values = c('#b45346')) + theme_bw() +
      guides(colour = FALSE) + 
      labs(y = "Shortfall or excess of temperatures occupied (°C)",
           x = "Absolute realized range latitudinal midpoint (°N/S)") + 
      scale_y_continuous(limits = c(-25, 12)) + 
      theme(panel.grid = element_blank())
    
    cold_dots <- filling %>%
      filter(warm_cold == "cold") %>%
      filter(!is.infinite(filling_value)) %>%
      ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
      scale_color_manual(values = c('steelblue')) + theme_bw() +
      guides(colour = FALSE) + 
      labs(y = "Shortfall or excess of temperatures occupied (°C)",
           x = "Absolute realized range latitudinal midpoint (°N/S)") + 
      scale_y_continuous(limits = c(-30, 30)) + 
      theme(panel.grid = element_blank())
    
  }
  ggsave(warm_dots, path = path, 
         filename = "niche-filling-distribution-warm-dots-one-limit.png", 
         device = "png", width = 4.5, height = 4)
  ggsave(cold_dots, path = path, 
         filename = "niche-filling-distribution-cold-dots-one-limit.png", 
         device = "png", width = 4.5, height = 4)

  type = type + 1
}