## this script makes a bunch of figures to visualize the niche filling data
library(tidyverse)
library(cowplot)
library(lme4)
library(nlme)
select = dplyr::select

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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####   Look at how thermal breadth changes with latitude    ######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## filter to species with only both thermal tolerance limits 
both_lims <- Te %>%
  filter(!is.na(ctmax) & !is.na(ctmin)) %>%
  mutate(abs_lat_mp = abs(lat_mp)) %>%
  mutate(realm = factor(.$realm, levels = c("Terrestrial", "Intertidal", "Marine")),
                                            ordered = TRUE) %>% ## reorder realm
  mutate(perfect_warm = ifelse(warm_over == 0 & warm_under == 0, 0, NA)) %>%
  gather(key = "filling_type", value = "filling_value", c(warm_under, warm_over, perfect_warm)) %>%
  filter(filling_value != 0 | filling_value == 0 & (filling_type %in% c("perfect_warm"))) %>%
  filter(!is.infinite(filling_value), !is.infinite(filling_value)) %>%
  filter(!is.na(filling_value), !is.infinite(filling_value))

## calculate thermal breadth 
both_lims$tolerance_breadth = as.numeric(as.character(both_lims$ctmax - both_lims$ctmin)) 

## plot 
both_lims %>%
  ggplot(aes(x = abs_lat_mp, y = tolerance_breadth, colour = realm)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_smooth(method = "lm") +
  labs(x = "Realized range latitudinal midpoint (°N/S)",
       y = "Thermal tolerance breadth (°C)", 
       colour = "") +
  scale_colour_manual(values = c("#a50026", "#f57396", "#f58d42"),
                      labels = c("Terrestrial", "Intertidal marine", "Subtidal marine"))

mod = lm(tolerance_breadth ~ realm*abs_lat_mp,
  data = both_lims)

summary(mod)

ter = filter(both_lims, realm == "Terrestrial")
newdata_ter = expand.grid(realm = unique(ter$realm),
                      abs_lat_mp = seq(from = min(ter$abs_lat_mp), 
                                       to = max(ter$abs_lat_mp),
                                       by = 0.1))
int = filter(both_lims, realm == "Intertidal")
newdata_int = expand.grid(realm = unique(int$realm),
                          abs_lat_mp = seq(from = min(int$abs_lat_mp), 
                                           to = max(int$abs_lat_mp),
                                           by = 0.1))
mar = filter(both_lims, realm == "Marine")
newdata_mar = expand.grid(realm = unique(mar$realm),
                          abs_lat_mp = seq(from = min(mar$abs_lat_mp), 
                                           to = max(mar$abs_lat_mp),
                                           by = 0.1))
newdata = rbind(newdata_int, newdata_ter, newdata_mar)

preds = predict(mod, newdata = newdata, se.fit = T, level = 0)

preds <- newdata %>%
  mutate(pred_breadth = preds$fit) %>%
  mutate(pred_breadth_SE = preds$se.fit) %>%
  mutate(realm = factor(.$realm, levels = c("Terrestrial", "Intertidal", "Marine")),
         ordered = TRUE) 

## colour by warm range edge underfilling 
ttb_realms <- both_lims %>%
  mutate(realm = factor(.$realm, levels = c("Terrestrial", "Intertidal", "Marine")),
         ordered = TRUE) %>%
  ggplot(aes(x = abs_lat_mp, y = tolerance_breadth, fill = filling_value, shape = realm)) +
  geom_point(size = 2.5, stroke = 0.5, colour = "black") +
  scale_linetype_discrete(labels = c("Terrestrial", "Intertidal marine", "Subtidal marine")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Absolute realized range latitudinal midpoint (°N/S)",
       y = "Thermal tolerance breadth (°C)", 
       colour = "Warm niche\nunderfilling (°C)",
       fill = "Warm niche\nunderfilling (°C)",
       linetype = "",
       shape = "") +
  scale_fill_gradient(low = "#b45346", high = "white") +
  geom_line(data = preds, aes(x = abs_lat_mp, y = pred_breadth, group = realm,
                              linetype = realm),
            inherit.aes = F) +
  scale_shape_manual(labels = c("Terrestrial", "Intertidal marine", "Subtidal marine"),
                     values = c(21,22,24)) +
  geom_ribbon(data = preds, aes(x = abs_lat_mp,
                                ymin = (pred_breadth-1.96*pred_breadth_SE),
                                ymax = pred_breadth+1.96*pred_breadth_SE,
                                group = realm),
              fill = "darkgrey",
              alpha = 0.15,
              colour = NA,
              inherit.aes = F) 

ggsave(ttb_realms, path ="figures/extended-data", filename = "predictions_tolerance-breadth_latitude.png",
       width = 5, height = 3.25, device = "png")


## see if warm underfilling increases with thermal tolerance breadth
mod_warm_breadth = lm(data = both_lims, 
   filling_value ~ tolerance_breadth*realm)

summary(mod_warm_breadth)

ter = filter(asym_lat, realm == "Terrestrial")
newdata_ter = expand.grid(realm = unique(ter$realm),
                          tolerance_breadth = seq(from = min(ter$tolerance_breadth), 
                                           to = max(ter$tolerance_breadth),
                                           by = 0.1))
int = filter(asym_lat, realm == "Intertidal")
newdata_int = expand.grid(realm = unique(int$realm),
                          tolerance_breadth = seq(from = min(int$tolerance_breadth), 
                                           to = max(int$tolerance_breadth),
                                           by = 0.1))
mar = filter(both_lims, realm == "Marine")
newdata_mar = expand.grid(realm = unique(mar$realm),
                          tolerance_breadth = seq(from = min(mar$tolerance_breadth), 
                                           to = max(mar$tolerance_breadth),
                                           by = 0.1))
newdata = rbind(newdata_int, newdata_ter, newdata_mar)

preds = predict(mod_warm_breadth, newdata = newdata, se.fit = T, level = 0)

preds <- newdata %>%
  mutate(pred_warm = preds$fit) %>%
  mutate(pred_warm_SE = preds$se.fit) %>%
  mutate(realm = factor(.$realm, levels = c("Terrestrial", "Intertidal", "Marine")),
         ordered = TRUE) 

ttb_wnf <- both_lims %>%
  mutate(realm = factor(.$realm, levels = c("Terrestrial", "Intertidal", "Marine")),
         ordered = TRUE) %>%
  ggplot(aes(x = tolerance_breadth, y = filling_value, fill = filling_value, shape = realm)) +
  scale_fill_gradient(low = "#b45346", high = "white") +
  geom_point(size = 2.5, stroke = 0.5, colour = "black") +
  scale_linetype_discrete(labels = c("Terrestrial", "Intertidal marine", "Subtidal marine")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Thermal tolerance breadth (°C)",
       y = "Warm niche filling (°C)", 
       colour = "",
       fill = "Warm niche\nunderfilling (°C)", 
       shape = "",
       linetype =  "") +
  geom_line(data = preds, aes(x = tolerance_breadth, y = pred_warm, group = realm,
                              linetype = realm),
            inherit.aes = F) +
  scale_shape_manual(labels = c("Terrestrial", "Intertidal marine", "Subtidal marine"),
                     values = c(21,22,24)) +
  geom_ribbon(data = preds, aes(x = tolerance_breadth,
                                ymin = (pred_warm-1.96*pred_warm_SE),
                                ymax = pred_warm+1.96*pred_warm_SE,
                                group = realm),
              fill = "darkgrey",
              alpha = 0.15,
              colour = NA,
              inherit.aes = F) +
  geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(15, 55))

ggsave(ttb_wnf, path ="figures/extended-data", filename = "predictions_tolerance-breadth_warm-underfilling.png",
       width = 5.05, height = 3.25, device = "png")


## make dot whisker plots

# plotting estimates (fixed effects) 
sum <- summary(mod)
df <- as.data.frame(sum$coefficients) 
CI <- as.data.frame(confint(mod, full=T, level = 0.95)) # get confidence intervals for full model
df$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df$CI.max <-CI$`97.5 %`
data.table::setDT(df, keep.rownames = "coefficient") #put rownames into column
df$how_conf = ifelse((df$CI.min < 0 & df$CI.max < 0) | (df$CI.min > 0 & df$CI.max > 0), "so_conf", "not")# add column to show which predictors have CIs that don't cross 0

## reorder rows 
df$coefficient <- factor(df$coefficient, levels = c("(Intercept)", "abs_lat_mp", "realmIntertidal:abs_lat_mp",
                                                    "realmMarine:abs_lat_mp", "realmIntertidal", "realmMarine"), 
                         ordered = TRUE)

whisker_lat <- ggplot(data=df, aes(x=fct_rev(coefficient), y=Estimate)) + 
  geom_hline(yintercept = 0) + 
  theme_classic() +
  geom_pointrange(aes(ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = '#b45346', size = 1, shape = 21) + 
  coord_flip() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("white", "#b45346")) +
  labs(y = "Effect of variable on thermal tolerance breadth", x = "") + 
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal",
                              "Abs. realized range latitudinal midpoint x realm: subtidal",
                              "Abs. realized range latitudinal midpoint x realm: intertidal",
                              "Abs. realized range latitudinal midpoint", 
                              "Reference"))


sum <- summary(mod_warm_breadth)
df2 <- as.data.frame(sum$coefficients) 
CI <- as.data.frame(confint(mod_warm_breadth, full=T, level = 0.95)) # get confidence intervals for full model
df2$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df2$CI.max <-CI$`97.5 %`
data.table::setDT(df2, keep.rownames = "coefficient") #put rownames into column
df2$how_conf = ifelse((df2$CI.min < 0 & df2$CI.max < 0) | (df2$CI.min > 0 & df2$CI.max > 0), "so_conf", "not")# add column to show which predictors have CIs that don't cross 0

## reorder rows 
df2$coefficient <- factor(df2$coefficient, levels = c("(Intercept)", "tolerance_breadth", "rtolerance_breadth:realmIntertidal",
                                                    "tolerance_breadth:realmMarine", "realmIntertidal", "realmMarine"), 
                         ordered = TRUE)

whisker_wniche <- ggplot(data=df2, aes(x=fct_rev(coefficient), y=Estimate)) + 
  geom_hline(yintercept = 0) + 
  theme_classic() +
  geom_pointrange(aes(ymin=CI.min, ymax=CI.max,  
                      fill = how_conf), colour = '#b45346', size = 1, shape = 21) + 
  coord_flip() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("white", "#b45346")) +
  labs(y = "Effect of variable on warm niche underfilling", x = "") + 
  scale_x_discrete(labels = c("Realm: subtidal",
                              "Realm: intertidal",
                              "Thermal tolerance breadth x realm: subtidal",
                              "Thermal tolerance breadth x realm: intertidal",
                              "Thermal tolerance breadth", 
                              "Reference"))

## write out to file 
ggsave(whisker_wniche, path ="figures/extended-data", filename = "whisker-plot_tolerance-breadth_warm-underfilling.png",
       width = 6.05, height = 1.5, device = "png")
ggsave(whisker_lat, path ="figures/extended-data", filename = "whisker-plot_tolerance-breadth_latitude.png",
       width = 6.75, height = 1.5, device = "png")


  