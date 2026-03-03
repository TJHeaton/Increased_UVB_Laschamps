library("tidyverse")
library("cowplot")
library("viridis")
library("broom")
library("here")

# Import functions from functions directory
sapply(list.files(here("Functions"), full.names = TRUE), source)
###################################
# Get data files

sample_info <- read_csv("data_greenhouse_2024/sample_info.csv")
sample_meta <- read_csv("data_greenhouse_2024/sample_meta.csv")
gc_log <- read_csv("data_greenhouse_2024/gc_log.csv")

# Fix anomalous i196 standard by dividing by two- assuming loaded twice
gc_log$i196[gc_log$i196 > 7*(10^6)] <- NA

# correct typo in datasheet
gc_log$SampleID[gc_log$Code == 3079] <- "J526"

standard_concentrations <- read_csv2(here("data_greenhouse_2024", "calibration_standards_sheet.csv"))


#####################################################
# Calibrating with the ratio 

gc_data_processed <- prepare_gc_data(.x = gc_log,.y = sample_info, .z = standard_concentrations )

# Rescale for nGrains being different sizes in the samples (we know this is a linear relationship based on pollen standards)
gc_data_processed$sample <- gc_data_processed$sample  %>% mutate(i161 = i161*150/nGrains, 
                                                                 Ratio = i161/i196)
gc_data_processed$sample$calBatch <-  gc_data_processed$sample$Date_factor
gc_data_processed$polcal$calBatch <-  gc_data_processed$polcal$Date_factor
gc_data_processed$cal$calBatch <-  gc_data_processed$cal$Date_factor

gc_data_processed$polcal$Ratio <- gc_data_processed$polcal$i161/gc_data_processed$polcal$i196
gc_data_processed$cal$Ratio <- gc_data_processed$cal$i161/gc_data_processed$cal$i196


### Calibrations 
check_pollen_calibrations(.x = gc_data_processed, grouping = "calBatch")



# Convert to quantities of pCA based on calibration curve using a master calibration model for pollen grains and pCA
master_model <- create_master_calib_model(.x = gc_data_processed, method= "ratio")

# Calibrate Sample results to number of pollen grains
gc_data_calibrated <- calibrate_sample_polcal(.x = gc_data_processed, method = "ratio") 
gc_data_calibrated <- calibrate_grains_to_pCA(.x = gc_data_calibrated, .y = master_model) 

# Align with metadata
calBatch <- dplyr::select(gc_data_processed$sample, c("SampleID", "calBatch"))

joined_data <- gc_data_calibrated  %>% 
  inner_join(sample_meta, by = "SampleID") %>% 
  left_join(calBatch, by = "SampleID") %>% 
  mutate(calibrated = calibrated/150) %>% 
  dplyr::select(Code, SampleID, calibrated, SampleCode, Frame, Treatment, TreeNumber, Sampling_Number)

write.csv(joined_data, file = here("output_data", "gc_cultivar_data.csv"))

