prepare_gc_data <- function(
    .x = gc_data, # gc_data object
    .y = sample_info, # sample info object
    .z = standard_concentrations # standard concentrations sheet
    ){ 

  .x1 <- .x %>%   
    left_join(sample_info, by = "SampleID") %>% 
    mutate(Date_factor = as.factor(Date)) %>% 
    dplyr::select(Date, Date_factor, calBatch, Code, SampleID, std, i161,i196, Ratio, type, CoreCode, nGrains) %>%  
  filter(type %in% c("cal", "polcal", "sample")) 
  
  .cal <- .x1 %>% 
    filter(type == "cal") %>% 
    dplyr::select(-c("CoreCode", "nGrains")) %>% 
    left_join(.z, by = "SampleID")
  
  .polcal <- .x1 %>% 
    filter(type == "polcal")
  
  .sample <- .x1 %>% 
    filter(type == "sample")
  
  final_object <- list(cal = .cal, 
                       polcal = .polcal,
                       sample = .sample)
  final_object
}
