chemical_calibration_function <- function(.x = gc_data_processed,
                                      method = "ratio"){
  
  .stdcal <- .x$cal %>% mutate(calibration_nest = "n1")
  .sample_to_cal <- .x$sample %>% mutate(calibration_nest = "n1")
  
  if(method == "ratio") {
    .stdcal$calData <- .stdcal$Ratio
    .sample_to_cal$calData <- .sample_to_cal$Ratio
  } else {
    .stdcal$calData <- .stdcal$i161
    .sample_to_cal$calData <- .sample_to_cal$i161
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
  
  calibrated_data <- .sample_to_cal %>% 
    left_join(.cal_model, by = "calibration_nest") %>% 
    mutate(calibrated = map2_dbl(calData, fit, .f = calibrate_data )) %>% 
    dplyr::select(Date, Code, SampleID, i161, i196, Ratio, nGrains, calibrated) %>% 
    mutate(calibrated_ngG = calibrated*(10^9)/100) # not nGrains here as have already corrected the i161 to rescale to 100
  
  

  
  # Plot calibration model
  master_model_plot <- .stdcal %>% 
    ggplot(aes(x = pCA_g, y = calData))+
    geom_point()+
    geom_smooth(method= "lm")+
    labs(x = "pCA (g)", y = "pCA: VA Ratio") +
    theme_bw()
  
  master_model_plot
  calibrated_data
}
