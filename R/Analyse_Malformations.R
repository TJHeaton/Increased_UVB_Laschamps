# Perform logistic/binomial regression of grain malformations
full_output_sorted <- full_output[order(full_output$expected_paleomag),]
full_output_sorted$Laschamps <- (full_output_sorted$Time == "41 ka BP")


# Create separate row for each sample
pollen_malformations <- data.frame(n_grains_malformed = c(full_output_sorted$sample_A_n_grains_malformed,
                                                          full_output_sorted$sample_B_n_grains_malformed),
                                   n_grains_normal = c(full_output_sorted$sample_A_n_grains_normal,
                                                       full_output_sorted$sample_B_n_grains_normal),
                                   Time = rep(full_output_sorted$Time, 2),
                                   expected_paleomag = rep(full_output_sorted$expected_paleomag, 2),
                                   calage = rep(full_output_sorted$calage, 2),
                                   Laschamps = rep(full_output_sorted$Laschamps, 2))

# Remove rows with NA, i.e. not measured twice
pollen_malformations <- subset(pollen_malformations, !is.na(pollen_malformations$n_grains_normal))


# Change Time to be a factor variable in pollen malformations
pollen_malformations$Time <- as.factor(pollen_malformations$Time)
pollen_malformations$prop_malformed <- pollen_malformations$n_grains_malformed / (pollen_malformations$n_grains_malformed + pollen_malformations$n_grains_normal)


########
# Fit basic glm to the time period data (without overdispersion)
cat("Initial glm (binomial) without including overdispersion \n")
# Consider time period as factor variable
grain_malformation_logistic_model_group_all_periods <- glm(cbind(n_grains_malformed, n_grains_normal) ~ Time,
                                                           family = binomial,
                                                           data = pollen_malformations)
print(summary(grain_malformation_logistic_model_group_all_periods))
# Calculate overdispersion (it is highly significant - p-value of 3 x 10^-21)

cat("Calculate overdispersion - very highly significant (p-value of 3 x 10^-21) \n")
cat("This basic model is innapropriate, must account for this overdispersion \n \n")
malformations_overdispersion_group_all_periods <- overdispersion(grain_malformation_logistic_model_group_all_periods)
print(malformations_overdispersion_group_all_periods)
cat("\n")



#########################################################################
# Adjust accounting for overdispersion

# Fitting using quasi-likelihood approach of Williams (1982) where do not specify beta specifically
overdisperse_malformations_glm_all_periods <- quasibin(
  cbind(n_grains_malformed, n_grains_normal) ~ Time,
  data = pollen_malformations, phi = NULL, tol = 1e-8)
malformations_overdispersion_phi <- overdisperse_malformations_glm_all_periods@phi
malformations_df_all_periods <- overdisperse_malformations_glm_all_periods@fm$df.residual
malformations_pearson_resids_all_periods <- sum(residuals(overdisperse_malformations_glm_all_periods, type = "pearson")^2)

# Create confidence intervals for the different time periods
predict_grid <- data.frame(Time = levels(pollen_malformations$Time))

predict_vals <- predict(overdisperse_malformations_glm_all_periods, newdata = predict_grid, se = TRUE, type = "link")
predict_response <- predict(overdisperse_malformations_glm_all_periods, newdata = predict_grid, se = TRUE, type = "response")

# Create 95% CIs for values
round(coef(overdisperse_malformations_glm_all_periods) + q_val_tail * sqrt(diag(vcov(overdisperse_malformations_glm_all_periods))), digits = 2)
round(coef(overdisperse_malformations_glm_all_periods) - q_val_tail * sqrt(diag(vcov(overdisperse_malformations_glm_all_periods))), digits = 2)

summary_estimates <-  data.frame(predict_grid,
                                 link = predict_vals$fit,
                                 se_link = predict_vals$se.fit,
                                 lower_link = predict_vals$fit - q_val_tail * predict_vals$se.fit,
                                 upper_link = predict_vals$fit + q_val_tail * predict_vals$se.fit)

summary_estimates$odds <- exp(summary_estimates$link)
summary_estimates$lower_odds <- exp(summary_estimates$lower_link)
summary_estimates$upper_odds <- exp(summary_estimates$upper_link)

summary_estimates$prob <- summary_estimates$odds/ (1 + summary_estimates$odds)
summary_estimates$lower_prob <- summary_estimates$lower_odds/ (1 + summary_estimates$lower_odds)
summary_estimates$upper_prob <- summary_estimates$upper_odds/ (1 + summary_estimates$upper_odds)

# Create dataframe containing intervals for grain malformation by period
malformation_rate_output_by_period <- data.frame(period = summary_estimates$Time,
                                                 prob_malform = summary_estimates$prob,
                                                 prob_malform_upper_CI = summary_estimates$upper_prob,
                                                 prob_malform_lower_CI = summary_estimates$lower_prob)

###### Refit same model but using Laschamps as baseline (for providing p-values for the paper)
## Note this is exactly the same model as above (just with different baseline level)

cat("------------------- \n")
cat("Quasi-likelihood glm (binomial) approach of Williams (1982) including overdispersion \n")
cat("Relevelled data so that Laschamps forms the baseline for comparison \n")
cat("------------------- \n \n")
# Relevel model with Laschamps as baseline (to make comparison easier for manuscript)
pollen_malformations$TimeLaschampsBaseline <- relevel(pollen_malformations$Time,
                                                      ref = "41 ka BP")

# Fitting using quasi-likelihood approach of Williams (1982) where do not specify beta specifically
laschamps_baseline_overdisperse_malformations_glm_all_periods <- quasibin(
  cbind(n_grains_malformed, n_grains_normal) ~ TimeLaschampsBaseline,
  data = pollen_malformations, phi = NULL, tol = 1e-8)
print(laschamps_baseline_overdisperse_malformations_glm_all_periods)

# Create confidence intervals for the different time periods
laschamps_baseline_predict_grid <- data.frame(TimeLaschampsBaseline = levels(pollen_malformations$TimeLaschampsBaseline))

laschamps_baseline_predict_vals <- predict(laschamps_baseline_overdisperse_malformations_glm_all_periods,
                                           newdata = laschamps_baseline_predict_grid, se = TRUE, type = "link")
laschamps_baseline_predict_response <- predict(laschamps_baseline_overdisperse_malformations_glm_all_periods,
                                               newdata = laschamps_baseline_predict_grid, se = TRUE, type = "response")


# Create 95% CIs for values
round(coef(laschamps_baseline_overdisperse_malformations_glm_all_periods) + q_val_tail * sqrt(diag(vcov(laschamps_baseline_overdisperse_malformations_glm_all_periods))), digits = 2)
round(coef(laschamps_baseline_overdisperse_malformations_glm_all_periods) - q_val_tail * sqrt(diag(vcov(laschamps_baseline_overdisperse_malformations_glm_all_periods))), digits = 2)


###############################################################
# Now plot the malformation analysis in Figure 3
source("R/Plot_Malformation_Analysis.R")
