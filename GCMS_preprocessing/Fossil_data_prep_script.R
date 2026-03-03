library("tidyverse")
library("cowplot")
library("viridis")
library("broom")
library("here")

# Import functions from functions directory
sapply(list.files(here("Functions"), full.names = TRUE), source)
###################################
# Get data files

sample_info <- read_csv("data_fossil/sample_info.csv")
sample_meta <- read_csv("data_fossil/sample_meta.csv")
gc_log <- read_csv("data_fossil/gc_log.csv") 

standard_concentrations <- read_csv2(here("data_fossil", "calibration_standards_sheet.csv")) 

section_table <- tibble(Section = c("B-22", "C-02", "C-08", "C-13", "C-17"),
                        Time = c("65 ka BP", "01 ka BP", "08 ka BP", "41 ka BP", "53 ka BP"))


#####################################################
# Calibrating

gc_data_processed <- prepare_gc_data(.x = gc_log,.y = sample_info, .z = standard_concentrations )

# Rescale for nGrains being different sizes in the samples (we know this is a linear relationship based on pollen standards)
gc_data_processed$sample <- gc_data_processed$sample  %>% mutate(i161 = i161*100/nGrains, 
                                                                 Ratio = i161/i196)
gc_data_processed$sample$calBatch <-  gc_data_processed$sample$Date_factor
gc_data_processed$polcal$calBatch <-  gc_data_processed$polcal$Date_factor
gc_data_processed$cal$calBatch <-  gc_data_processed$cal$Date_factor

###  Check calibrations
#check_pollen_calibrations(.x = gc_data_processed, grouping = "calBatch")

# Convert to quantities of pCA based on calibration curve using a master calibration model for pollen grains and pCA
master_model <- create_master_calib_model(.x = gc_data_processed, method= "ratio")

# Calibrate Sample results to number of pollen grains
gc_data_calibrated <- calibrate_sample_polcal(.x = gc_data_processed, method = "ratio") 
gc_data_calibrated <- calibrate_grains_to_pCA(.x = gc_data_calibrated, .y = master_model) 


# Align with metadata
calBatch <- select(gc_data_processed$sample, c("SampleID", "calBatch"))

joined_data <- gc_data_calibrated  %>% 
  inner_join(sample_meta, by = "SampleID") %>% 
  left_join(calBatch, by = "SampleID") %>% 
  left_join(section_table, by = "Section") %>% 
  mutate(calibrated = calibrated/100) %>% 
  dplyr::select(Code, SampleID, calibrated, CoreCode, Section, Depth, Time)


write_csv(joined_data, file =here("output_data", "GCMS_pCA_day_dependent_calibration.csv"))



