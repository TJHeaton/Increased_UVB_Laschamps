calibrate_sample_polcal <- function(.x = gc_data_processed, method = c("ratio, i161")){
  
  .polcal <- .x$polcal
  .sample <- .x$sample
  
  if(method == "ratio") {
    .polcal$calData <- .polcal$Ratio
    .sample$calData <- .sample$Ratio
  } else {
    .polcal$calData <- .polcal$i161
    .sample$calData <- .sample$i161
  }
  
  
  .cal_model <- .polcal %>% 
    nest(data = -calBatch) %>% 
    mutate(fit = map(data, ~ lm(nGrains ~ calData, data = .)),
           results = map(fit, augment))
  
  ### Fit the calibration curve to the pollen samples
  calibrate_data <- function(.x, # the data used to calibrate
                             .y){ # the fitted model
    
    predict(.y, newdata= data.frame(calData = .x))
  }
  
  .sample_cal <- .sample %>% 
    left_join(.cal_model, by = "calBatch") %>% 
    mutate(calibrated = map2_dbl(calData, fit, .f = calibrate_data )) %>% 
    dplyr::select(Date, Code, SampleID, i161, i196, Ratio, nGrains, calibrated)
  
  
  # Plot calibration model
  .cal_model_plot <- .cal_model %>% 
    unnest(results) %>% 
    ggplot(aes(x = nGrains, y = .fitted)) +
    geom_abline(intercept = 0, slope = 1, alpha = .5) +  # Line of perfect fit
    geom_point() +
    facet_grid(calBatch ~ .) +
    labs(x = "pCA", y = "Predicted Value") +
    theme_bw()
  
  list(.sample = .sample_cal,
       .polcal = .cal_model,
       .model_plot = .cal_model_plot)
}
