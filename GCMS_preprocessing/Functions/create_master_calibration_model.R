create_master_calib_model <- function(.x = gc_data_processed,
                                      method = "ratio"){
  
  .polcal <- .x$polcal %>% mutate(calibration_nest = "n1")
  .stdcal <- .x$cal %>% mutate(calibration_nest = "n1")
  
  if(method == "ratio") {
    .polcal$calData <- .polcal$Ratio
    .stdcal$calData <- .stdcal$Ratio
  } else {
    .polcal$calData <- .polcal$i161
    .stdcal$calData <- .stdcal$i161
  }
  
  .cal_model <- .stdcal %>% 
    nest(data = -calibration_nest) %>% 
    mutate(fit = map(data, ~ lm(pCA_g ~ calData, data = .)),
           results = map(fit, augment))
  
  ### Fit the calibration curve to the pollen samples
  calibrate_data <- function(.x, # the data used to calibrate
                             .y){ # the fitted model
    
    predict(.y, newdata= data.frame(calData = .x))
  }
  
  ratio_converted_pcG <- .polcal %>% 
    left_join(.cal_model, by = "calibration_nest") %>% 
    mutate(calibrated = map2_dbl(calData, fit, .f = calibrate_data )) %>% 
    dplyr::select(Date, Code, SampleID, i161, i196, Ratio, nGrains, calibrated)
  
  master_model <- lm(calibrated ~ nGrains, data = ratio_converted_pcG)
  
  # Plot calibration model
  master_model_plot <- ratio_converted_pcG %>% 
    ggplot(aes(x = nGrains, y = calibrated))+
    geom_point()+
    geom_smooth(method= "lm")+
    labs(x = "Number of Grains", y = "pCA (g)") +
    theme_bw()
  
  list(.master_model = master_model,
       .model_plot = master_model_plot)
  
}