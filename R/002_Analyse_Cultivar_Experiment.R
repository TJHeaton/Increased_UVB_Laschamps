# File containing all R analysis of p-CA results from cultivar experiment
# You can run just this code and it will call everything else

library(tidyverse)
library(cowplot)

conf_level <- 0.95

cultivar <- read_csv("data/Cultivar_Experiment/gc_cultivar_data.csv")

# Exclude rows that have missing values (NA) in any of calibrated, TreeNumber, or Treatment
cultivar <- cultivar %>%
  drop_na(calibrated, TreeNumber, Treatment)

# Plot the p-CA values by tree
cultivar_plot_by_tree <- ggplot(cultivar, aes(x = as.factor(TreeNumber), fill = Treatment, y = calibrated)) +
  geom_boxplot(outlier.colour="red") +
  geom_jitter(shape=16, position=position_jitter(0.1))

cultivar_plot_by_tree
# Suggests that tree 52 (UV-B treatment) is an outlier
# This was confirmed in lab notes (very pale and pollen has underdeveloped sacci)
# Hence removed from further analysis

cat("------------------- \n \n")
cat("We identify tree 52 (UV-B treatment) as an outlier \n")
cat("This was checked against the lab notes and confirmed as containing underdeveloped sacci \n")
cat("and very pale in comparison to other cultivars. \n")
cat("Tree 52 (UV-B) samples were therefore removed from analysis \n \n")
cat("------------------- \n")


# Cultivars were sampled for pollen several times (on different occasions) denoted by sampling number
# Each of these samples were then measured three times on GC-MS (other than 2nd sampling of tree 57 which has only two GC-MS measurements)

# Find the mean values of each sampling (i.e. mean of three repeats measurements on each sampling)
# Same original pollen sample if share code, tree number, treatment and sampling number
cultivar_merge_repeats <- cultivar %>%
  group_by(SampleCode, TreeNumber, Treatment, Sampling_Number) %>%
  summarise(mean_calibrated = mean(calibrated))

# Use only the first sampling results
cultivar_first_sampling_only <- cultivar_merge_repeats %>%
  filter(Sampling_Number == "01")

# Remove tree 52 measurements(identified as not appropriate based on lab notes)
id_tree_52 <- which(cultivar_first_sampling_only$TreeNumber == "052")
cultivar_first_sampling_only_no_tree_52 <- cultivar_first_sampling_only[-id_tree_52,]

# Perform Bartlett test of sample variance in each treatment group (i.e., UV-B irradiated or control)
cat("Bartlett test to compare sample variance in each treatment group (i.e., UV-B irradiated or control) \n")
print(bartlett.test(mean_calibrated ~ Treatment, data = cultivar_first_sampling_only_no_tree_52))
# No evidence (p = 0.97) of different variances between groups
cat("No evidence (p = 0.97) of different variances between two treatment groups \n \n")
cat("------------------- \n \n")

cultivar_plot_by_treatment <- ggplot(cultivar_first_sampling_only_no_tree_52,
                                     aes(x = Treatment, y = mean_calibrated, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, position=position_dodge(0.9)) +
  geom_jitter(shape=16, position=position_jitter(0.1))
cultivar_plot_by_treatment <- cultivar_plot_by_treatment +
  ylab(expression(paste(italic("para"),"-Coumaric Acid (ng grain"^-1,")"))) +
  xlab("Treatment Group") +
  theme_bw(
    base_size = 9
  )  +
  scale_x_discrete(labels = c("Control Group", "UV-B Group")) +
  theme(legend.position = c(0.87, 0.9)) +
  scale_fill_discrete(name = "Treatment", labels = c(expression(paste("Control (0 kJ m"^-2, " day"[BE]^-1, ")")),
                                                     expression(paste("UV-B (16.8 kJ m"^-2, " day"[BE]^-1, ")"))))

cultivar_plot_by_treatment

ggsave("output/FigS6_Cultivar_Plot_by_Treatment.png",
       width = 8, height = 5)

########## Fieller's test to compare groups
alpha <- (1 - conf_level)/2 # Uses conf_level from main analysis

# Find mean, se AND n_obs by time period
cultivar_first_sampling_only_no_tree_52 %>%
  group_by(Treatment) %>%
  summarise(mean = mean(mean_calibrated),
            sd = sd(mean_calibrated),
            n = length(mean_calibrated))

cultivar_first_sampling_only_no_tree_52 <- as.data.frame(cultivar_first_sampling_only_no_tree_52)
pCA_mean_by_treatment <- by(data = cultivar_first_sampling_only_no_tree_52[,"mean_calibrated"], cultivar_first_sampling_only_no_tree_52$Treatment, mean)
pCA_sd_by_treatment <- by(data = cultivar_first_sampling_only_no_tree_52[,c("mean_calibrated")], cultivar_first_sampling_only_no_tree_52$Treatment, sd)
pCA_n_obs_by_treatment <- by(data = cultivar_first_sampling_only_no_tree_52[,c("mean_calibrated")], cultivar_first_sampling_only_no_tree_52$Treatment, length)


pCA_pooled_sd <- sqrt(
  sum( (pCA_n_obs_by_treatment - 1) * (pCA_sd_by_treatment^2) ) /
    (sum(pCA_n_obs_by_treatment) - length(pCA_n_obs_by_treatment))
)


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

# Now perform Fieller Test of Laschamps vs other treatments
id_control <- which(names(pCA_mean_by_treatment) == "CON")
id_treatment <- which(names(pCA_mean_by_treatment) == "UV")

# Print explaantion to screen
cat("Fieller test to analyse relative increases in p-CA in UV-B exposed group vs no UV-B controls \n")
cat("------------------- \n")

# Compare increase in treatment
cat("1st Measurement Only: Comparing pCA from control vs UVB treatment: \n")
print(Fieller_test(pCA_mean_by_treatment,
                   pCA_n_obs_by_treatment,
                   pCA_pooled_sd,
                   compare_id = c(id_treatment, id_control)))
cat("------------------- \n \n")



