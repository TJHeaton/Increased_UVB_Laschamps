calibrate_grains_to_pCA <- function(.x = gc_data_calibrated,
                                    .y = master_model ){
  .x1 = .x$.sample
  .y1 = .y$.master_model
  
  .x1 <- .x1 %>% 
    mutate(calibrated = predict(.y1, newdata = data.frame(nGrains = .x1$calibrated))*(10^9))
  
  .x1
}
