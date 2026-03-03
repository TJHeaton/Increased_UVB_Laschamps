# Perform analysis of the individual replicates data

#### Create data frame containing all relevant informaiton on individual replicates
# Read in individual replicates
pca_individual_replicate_data <- read.csv(gcms_fname, header = TRUE)
pca_individual_replicate_data$Time <- as.factor(pca_individual_replicate_data$Time)

# Add information on time period
# Sort ordering so ordered by time exactly (i.e. interval and then depth)
list_depth <- strsplit(pca_individual_replicate_data$Depth, '-', perl = TRUE)
pca_individual_replicate_data$mid_sampled_depth <- sapply(list_depth, function(x) mean(as.double(x)))

list_ka_interval <- as.double(substr(pca_individual_replicate_data$Time, 1, 2))
dummy_order <- (1000 * list_ka_interval) + pca_individual_replicate_data$mid_sampled_depth
pca_individual_replicate_data <- pca_individual_replicate_data[order(dummy_order),]

# Add column with expected paleomag and interval and calage
look_up_in_full_ouput <- match(pca_individual_replicate_data$Depth, full_output$Depth)
pca_individual_replicate_data$expected_paleomag <- full_output$expected_paleomag[look_up_in_full_ouput]
pca_individual_replicate_data$interval <- full_output$interval[look_up_in_full_ouput]
pca_individual_replicate_data$calage <- full_output$calage[look_up_in_full_ouput]
pca_individual_replicate_data$plotcol <- full_output$plotcol[look_up_in_full_ouput]

# Add a column which specifies which interval (within each over time period) the replicate comes from
pca_individual_replicate_data  <- as.data.frame(pca_individual_replicate_data  %>%
                                                  group_by(Time) %>% mutate(order_within_period = match(Depth, unique(Depth))))


################## Now perform main pCA analysis

# Perform linear regression of p-CA
pca_linear_model_replicate_level <- lm(calibrated ~ expected_paleomag + calage, data = pca_individual_replicate_data)

# Print the output to the terminal
cat("\n------------------- \n")
cat("Summary of linear model of p-CA vs geodynamo and calendar age \n")
print(summary(pca_linear_model_replicate_level)) # Explains 40% of the variance
cat("\n95% confidence intervals for coefficients of explanatory variables \n")
print(confint(pca_linear_model_replicate_level))
cat("\n\n")
cat("------------------- \n \n")

### Create plots:
# SI Fig S4 - Linear model residuals
# Fig 1 - Main (two panel) figure showing pCA vs paleomag and pCA vs time period
source("R/Plot_pCA_Analysis_Day_Dependent.R")


####################################################################################
######### Perform analysis that compares the Laschamps vs other time periods

cat("\n\nAnalysis comparing p-CA in Laschamps vs other  time periods \n")

cat("Bartlett test of sample variances in p-CA amounts between all time period groupings \n")
# Perform Bartlett test of sample variance in each time period
print(bartlett.test(calibrated ~ Time, data = pca_individual_replicate_data))
# No real evidence (p = 0.30) of different variances between groups


cat("------------------- \n")
cat("Bartlett test of sample variances between individual groups \n")
cat("------------------- \n")
# Compare variance of each group against Laschamps 41 ka BP
for(period in c("01 ka BP", "08 ka BP", "53 ka BP", "65 ka BP")) {
  cat("Comparing variance of replicates from ", period, "vs Laschamps: \n")
  sub_replicates <- pca_individual_replicate_data[pca_individual_replicate_data$Time %in% c(period,"41 ka BP"),]
  print(bartlett.test(calibrated ~ Time, data = sub_replicates))
  cat("\n")
}
# Note some indication that variance of 8ka replicates are smaller than the Laschamps
# For Fieller's Test below, we treat the groups as same variance

########## Fieller's test to compare Laschamps vs other time periods
alpha <- (1 - conf_level)/2 # Uses conf_level from main analysis

# Find mean and se by time period
pCA_mean_by_period <- by(data = pca_individual_replicate_data[,c("calibrated")], pca_individual_replicate_data$Time, mean)
pCA_sd_by_period <- by(data = pca_individual_replicate_data[,c("calibrated")], pca_individual_replicate_data$Time, sd)
pCA_n_obs_by_period <- by(data = pca_individual_replicate_data[,c("calibrated")], pca_individual_replicate_data$Time, length)

# Be conservative in Fieller's test of pCA increase (by choosing large estimate of pooled variance)
# Exclude data from the 8ka period in calculating the pooled variabce
id_remove_conservative <- which(names(pCA_mean_by_period) %in% c("08 ka BP"))
conservative_pCA_pooled_sd <- sqrt(
  sum( (pCA_n_obs_by_period[-id_remove_conservative] - 1) * (pCA_sd_by_period[-id_remove_conservative]^2) ) /
    (sum(pCA_n_obs_by_period[-id_remove_conservative]) - length(pCA_n_obs_by_period[-id_remove_conservative]))
)
pCA_pooled_sd <- conservative_pCA_pooled_sd

cat("------------------- \n")
cat("Fieller test to analyse relative increases in p-CA in Laschamps vs each other time period \n")
cat("------------------- \n")

Fieller_test <- function(X, n, pooled_sd, compare_id) {
  X_1 <- X[compare_id[1]]
  X_2 <- X[compare_id[2]]
  n_1 <- n[compare_id[1]]
  n_2 <- n[compare_id[2]]
  m <- n_1 + n_2 - 2

  CI_ratio_numerator <- (X_1 * X_2) + c(1,-1) *
    qt(alpha, m) * pooled_sd *
    sqrt(
      (X_1^2/n_2) + (1/n_1)*(X_2^2 - qt(alpha, m)^2 * (pooled_sd^2 / n_2))
    )
  CI_ratio_denominator <- X_2^2 - qt(alpha, m)^2 * (pooled_sd^2 / n_2)
  CI_ratio <- CI_ratio_numerator/CI_ratio_denominator

  retlist <- list(mean_ratio = X_1/X_2, CI_ratio = CI_ratio)
  return(retlist)
}

# Now perform Fieller Test of Laschamps vs other periods
id_laschamps <- which(names(pCA_mean_by_period) == "41 ka BP")
id_1ka <- which(names(pCA_mean_by_period) == "01 ka BP")
id_8ka <- which(names(pCA_mean_by_period) == "08 ka BP")
id_53ka <- which(names(pCA_mean_by_period) == "53 ka BP")
id_65ka <- which(names(pCA_mean_by_period) == "65 ka BP")


# Series of test of Laschamps vs other periods
# Compare with 1ka
cat("Comparing pCA from 1ka vs Laschamps: \n")
print(Fieller_test(pCA_mean_by_period,
                   pCA_n_obs_by_period,
                   pCA_pooled_sd,
                   compare_id = c(id_laschamps, id_1ka)))
cat("------------------- \n")

# Compare with 8ka
cat("Comparing pCA from 8ka vs Laschamps: \n")
print(Fieller_test(pCA_mean_by_period,
                   pCA_n_obs_by_period,
                   pCA_pooled_sd,
                   compare_id = c(id_laschamps, id_8ka)))
cat("------------------- \n")

# Compare with 53ka
cat("Comparing pCA from 53ka vs Laschamps: \n")
print(Fieller_test(pCA_mean_by_period,
                   pCA_n_obs_by_period,
                   pCA_pooled_sd,
                   compare_id = c(id_laschamps, id_53ka)))
cat("------------------- \n")

# Compare with 65ka
cat("Comparing pCA from 65ka vs Laschamps: \n")
print(Fieller_test(pCA_mean_by_period,
                   pCA_n_obs_by_period,
                   pCA_pooled_sd,
                   compare_id = c(id_laschamps, id_65ka)))
cat("------------------- \n")

########################################
### Potential Extra Analysis: Testing p-CA in NGS vs the other periods
### We see weak evidence this is significantly different
########################################
extra_NGS_analysis <- FALSE
if(extra_NGS_analysis) {
  cat("------------------- \n")
  cat("Additional Extra: Fieller test to analyse relative increases in p-CA in NGS vs each other time period \n")
  cat("------------------- \n")


  cat("Comparing pCA from NGS 65ka vs 1ka: \n")
  print(Fieller_test(pCA_mean_by_period,
                     pCA_n_obs_by_period,
                     pCA_pooled_sd,
                     compare_id = c(id_65ka, id_1ka)))
  cat("------------------- \n")

  cat("Comparing pCA from NGS 65ka vs 8ka: \n")
  print(Fieller_test(pCA_mean_by_period,
                     pCA_n_obs_by_period,
                     pCA_pooled_sd,
                     compare_id = c(id_65ka, id_8ka)))
  cat("------------------- \n")

  cat("Comparing pCA from NGS 65ka vs 53ka: \n")
  print(Fieller_test(pCA_mean_by_period,
                     pCA_n_obs_by_period,
                     pCA_pooled_sd,
                     compare_id = c(id_65ka, id_53ka)))
  cat("------------------- \n")
}


