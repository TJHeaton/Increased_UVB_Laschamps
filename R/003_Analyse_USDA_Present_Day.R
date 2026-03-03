# File containing calculation of present-day surface UV-B irradiances in Raleigh, N Carolina
# Data taken from USDA at https://uvb.nrel.colostate.edu/UVB/

# Note that some years (2020 and 2021) are heavily affected by COVID and have missing data
# Also 2024 has a lot of missing data in the period we are interested in
# These years are therefore not used
years_to_use <- c(2005:2019, 2022:2023)
location <- "NorthCarolina-Raleigh"
folder_name <- paste("data/USDA_UVB_Measurements/", location, "-uv-dailysums/", sep = "")

# Type of UV-B to use
type <- "Flint"

# Use 1st March - 31st May
start_day <- 60
end_day <- 151

UV_B_daily_doses <- c()
UV_B_yearly_average <- c()
for(year in years_to_use) {
  cat("Year = ", year, "\n")
  file_name <- paste(folder_name, location, "-", year, "-uv-dailysums.csv", sep = "")
  temp_uvb_data <- read.csv(file_name, skip = 9, header = TRUE)
  use_id <- which((temp_uvb_data$DOY >= start_day) & (temp_uvb_data$DOY <= end_day))
  temp_daily_doses <- temp_uvb_data[use_id, type]
  # Remove NAs/missing measurements (which are recorded as -9998 by USDA)
  below_zero_uvb <- which(temp_daily_doses < 0) # This identifies the missing measurements
  if(length(below_zero_uvb > 0)) {
    temp_daily_doses <- temp_daily_doses[-below_zero_uvb] # Removes the NAs/missing data
  }
  UV_B_yearly_average <- c(UV_B_yearly_average, mean(temp_daily_doses, na.rm = TRUE))
  UV_B_daily_doses <- c(UV_B_daily_doses, temp_daily_doses)
  cat("Minimum (Flint) Daily Dose = ", min(temp_daily_doses), "\n")
}

UV_B_yearly_average/1000

# Can calculate mean in two slightly different ways - both give the same result (14.2 kJ m-2 day-1)
mean(UV_B_daily_doses, na.rm = TRUE)/1000
mean(UV_B_yearly_average)/1000

cat("------------------- \n")
cat("Mean flux is ", round(mean(UV_B_yearly_average)/1000, 1), "kJ m-2 day-1 \n")




