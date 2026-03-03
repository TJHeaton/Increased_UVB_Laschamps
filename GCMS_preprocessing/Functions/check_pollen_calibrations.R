
check_pollen_calibrations <- function(.x = gc_data_processed,
                                      grouping = c("calBatch", "date")){
  
  if(grouping== "date") grouping = "Date_factor"
  
  vanillic_acid_plot <- bind_rows(.x$polcal,
                                  .x$sample) %>% 
    ggplot(aes(x = Code, y = i196))+
    geom_point(aes(col = type))+
    scale_color_viridis_d(begin = 0.2, end = 0.6)+
    labs(color = "Sample Type")+
    theme_bw()+
    ylab("i196 (m/z)")+
    xlab("Sample Code")
  
  
  polCal_plot <- .x$polcal %>% 
    pivot_longer(cols = c("i161", "Ratio"), names_to = "Type", values_to = "p-CA" ) %>% 
    ggplot(aes(x = nGrains, y = `p-CA`))+
    geom_point(aes()) +
    geom_smooth(method = "lm", aes(col = .data[[grouping]]), se= FALSE) +
    #geom_smooth(method = "lm", aes(col = factor(Date)), se= FALSE) +
    facet_wrap(.~Type ,scales = "free_y")+
    scale_color_viridis_d(begin = 0, end = 0.8)+
    labs(color = grouping)+
    xlab("Number of Pollen Grains")+
    theme_bw()
  
  calStd_plot <- .x$cal %>% 
    pivot_longer(cols = c("i161", "Ratio"), names_to = "Type", values_to = "p-CA" ) %>% 
    ggplot(aes(x = pCA_g, y = `p-CA`))+
    geom_point(aes()) +
    geom_smooth(method = "lm", aes(col = .data[[grouping]]), se= FALSE) +
    facet_wrap(.~Type ,scales = "free_y")+
    scale_color_viridis_d(begin = 0, end = 0.8)+
    labs(color = grouping)+
    xlab("p-Coumaric Acid (g)")+
    theme_bw()
  
  
  plot_grid(vanillic_acid_plot, calStd_plot, polCal_plot, nrow = 3)
  
}