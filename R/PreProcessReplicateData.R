# Pre-process and combine the data and the paleomag
# This file reads in all the individual replicates/samples


# Read in paleomag data
paleomag <- read.csv(paleomag_fname, header = TRUE)
names(paleomag)[1] <- "kyr"
paleomag$ZAM2 <- 10 * paleomag$ADM.x10.22.Am.2.
paleomag$calage <- 1000 * paleomag$kyr

# Read in malformations data
malformations_data <- read.csv(malformations_fname, header = TRUE)
n_obs <- nrow(malformations_data)

# Just keep the combined data (and tidy up)
malformations_data <- subset(malformations_data, select = -c(sample_A_date_picked,
                                                             sample_B_date_picked))

names(malformations_data)[which(names(malformations_data) == "n_grains_combined")] <- "combined_n_grains"
names(malformations_data)[which(names(malformations_data) == "n_grains_malformed_combined")] <- "combined_n_grains_malformed"

# Individual replicates
pca_individual_replicate_data <- read.csv(gcms_fname, header = TRUE)

sample_level_chemistry <- pca_individual_replicate_data  %>%
  group_by(CoreCode, Time, Depth) %>%
  summarise(mean_gcms = mean(calibrated, na.rm= TRUE),
            sd_gcms = sd(calibrated, na.rm= TRUE),
            n_samples = n(),
            se = sd_gcms/sqrt(n_samples),
            upper = mean_gcms+ 2*se,
            lower= mean_gcms-2*se)

# Convert to dataframe
sample_level_chemistry <- as.data.frame(sample_level_chemistry)

# Sort ordering so ordered by time exactly (i.e., interval and then depth)
list_depth <- strsplit(sample_level_chemistry$Depth, '-', perl = TRUE)
sample_level_chemistry$mid_sampled_depth <- sapply(list_depth, function(x) mean(as.double(x)))


list_ka_interval <- as.double(substr(sample_level_chemistry$Time, 1, 2))
dummy_order <- (1000 * list_ka_interval) + sample_level_chemistry$mid_sampled_depth
sample_level_chemistry <- sample_level_chemistry[order(dummy_order),]

# Store number of normal grains (combined) for logistic regression
malformations_data$sample_A_n_grains_normal <- malformations_data$sample_A_n_grains - malformations_data$sample_A_n_grains_malformed
malformations_data$sample_B_n_grains_normal <- malformations_data$sample_B_n_grains - malformations_data$sample_B_n_grains_malformed
malformations_data$combined_n_grains_normal <- malformations_data$combined_n_grains - malformations_data$combined_n_grains_malformed

####################################################################
# Create a dataframe which will store the full output

full_output <- malformations_data
rm(malformations_data)

# Find expected paleomag for sampling design
full_output$expected_paleomag <- rep(NA, n_obs)
for(i in 1:n_obs) {
  calage <- full_output$calage[i]
  calsig <- full_output$calsig[i]

  # Now find the expected (integrated) value of paleomag bearing in mind calendar age uncertainty
  #  int [ paleomag(t) * phi(t, mu, sigma^2)] dt
  bound_mult <- 4

  cal_grid <- calage + seq(from = -bound_mult * calsig,
                           to = bound_mult * calsig,
                           by = 1)

  # Find probabilities on this grid
  cal_prob <- dnorm(cal_grid, mean = calage, sd = calsig)
  cal_prob <- cal_prob/sum(cal_prob)

  # Find paleomag on this grid
  cal_paleo_mag <- approx(x = paleomag$calage, y = paleomag$ZAM2, xout = cal_grid)$y

  # Integrate and store
  full_output$expected_paleomag[i] <- sum(cal_paleo_mag * cal_prob)
}

# Now combine with chemistry GC-MS

# Check sampling depth match across
chemistry_depth <- strsplit(sample_level_chemistry$Depth, '-', perl = TRUE)

if(!identical(full_output$Sample.top..cm.,
              sapply(chemistry_depth, function(x) mean(as.double(x[1]))))) {
  cat("Sample depths are not the same between chemistry and sampling design ")
}

full_output <- cbind(full_output, sample_level_chemistry)

########################################################
# Add estimate of malformation rate and CIs (for sample A, sample B, and combined)
full_output$sample_A_prob_malformed <- NA
full_output$sample_A_prob_malformed_lower_CI <- NA
full_output$sample_A_prob_malformed_upper_CI <- NA

full_output$sample_B_prob_malformed <- NA
full_output$sample_B_prob_malformed_lower_CI <- NA
full_output$sample_B_prob_malformed_upper_CI <- NA

full_output$combined_prob_malformed <- NA
full_output$combined_prob_malformed_lower_CI <- NA
full_output$combined_prob_malformed_upper_CI <- NA

# Create confidence intervals for each sample
for(i in 1:nrow(full_output)) {
  # For sample A
  if(!is.na(full_output$sample_A_n_grains[i])) {
    sample_A_binomial_test_results <- binom.test(full_output$sample_A_n_grains_malformed[i],
                                                 full_output$sample_A_n_grains[i],
                                                 conf.level = conf_level)
    full_output$sample_A_prob_malformed[i] <- sample_A_binomial_test_results$estimate
    full_output$sample_A_prob_malformed_lower_CI[i] <- sample_A_binomial_test_results$conf.int[1]
    full_output$sample_A_prob_malformed_upper_CI[i] <- sample_A_binomial_test_results$conf.int[2]
  } else {
    full_output$sample_A_prob_malformed[i] <- NA
    full_output$sample_A_prob_malformed_lower_CI[i] <- NA
    full_output$sample_A_prob_malformed_upper_CI[i] <- NA
  }
  # For sample B
  if(!is.na(full_output$sample_B_n_grains[i])) {
    sample_B_binomial_test_results <- binom.test(full_output$sample_B_n_grains_malformed[i],
                                                 full_output$sample_B_n_grains[i],
                                                 conf.level = conf_level)
    full_output$sample_B_prob_malformed[i] <- sample_B_binomial_test_results$estimate
    full_output$sample_B_prob_malformed_lower_CI[i] <- sample_B_binomial_test_results$conf.int[1]
    full_output$sample_B_prob_malformed_upper_CI[i] <- sample_B_binomial_test_results$conf.int[2]
  } else {
    full_output$sample_B_prob_malformed[i] <- NA
    full_output$sample_B_prob_malformed_lower_CI[i] <- NA
    full_output$sample_B_prob_malformed_upper_CI[i] <- NA
  }

  # For combined samples
  if(!is.na(full_output$combined_n_grains[i])) {
  combined_binomial_test_results <- binom.test(full_output$combined_n_grains_malformed[i],
                                               full_output$combined_n_grains[i],
                                               conf.level = conf_level)
  full_output$combined_prob_malformed[i] <- combined_binomial_test_results$estimate
  full_output$combined_prob_malformed_lower_CI[i] <- combined_binomial_test_results$conf.int[1]
  full_output$combined_prob_malformed_upper_CI[i] <- combined_binomial_test_results$conf.int[2]
  } else {
    full_output$combined_prob_malformed[i] <- NA
    full_output$combined_prob_malformed_lower_CI[i] <- NA
    full_output$combined_prob_malformed_upper_CI[i] <- NA
  }
}

# Tidy up the name order of the dataset
full_output <- full_output[, c("CoreCode",
                               "interval",
                               "calage",
                               "calsig",
                               "core",
                               "PSP.top..cm.",
                               "PSP.bottom..cm.",
                               "Sample.top..cm.",
                               "Sample.bottom..cm.",
                               "Taxa",
                               "sample_A_n_grains",
                               "sample_A_n_grains_malformed",
                               "sample_A_n_grains_normal",
                               "sample_A_prob_malformed",
                               "sample_A_prob_malformed_lower_CI",
                               "sample_A_prob_malformed_upper_CI",
                               "sample_B_n_grains",
                               "sample_B_n_grains_malformed",
                               "sample_B_n_grains_normal",
                               "sample_B_prob_malformed",
                               "sample_B_prob_malformed_lower_CI",
                               "sample_B_prob_malformed_upper_CI",
                               "combined_n_grains",
                               "combined_n_grains_malformed",
                               "combined_n_grains_normal",
                               "combined_prob_malformed",
                               "combined_prob_malformed_lower_CI",
                               "combined_prob_malformed_upper_CI",
                               "expected_paleomag",
                               "Time",
                               "Depth",
                               "mean_gcms",
                               "sd_gcms",
                               "n_samples",
                               "se",
                               "upper",
                               "lower",
                               "mid_sampled_depth")]







