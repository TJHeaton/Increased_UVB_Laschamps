# File containing all R analysis of UV-B pollen data from Lake Suigetsu
# You can run just this code and it will call everything else

#############################################################
# This code will perform the formal analysis and creates figures used in the paper

# Figure 2 - The GC-MS data (individual replicates):
#   Panel A will show against paleo-magnetic values
#   Panel B will show against time period group (with a box-plot)

# Figure 3 - The grain malformation data
# This will show against the time period group
# We will plot the individual replicates and also the CI showing the expected values

# We also create:
# Fig S3 - Plot of the design space

##############################################################


##### Select some overall parameters
# Choose confidence interval for plotting
conf_level <- 0.95
q_val_tail <- -qnorm((1-conf_level)/2)

#### Plotting adjustments to make the plots tidy and non-overlapping

# Select plot size
label_size <- 2.4




# Libraries to make plotting easier
library(ggplot2)
library(cowplot)
library(readr)
library(tidyverse)
library(colorBlindness) # For plotting colours
library(aod) # For including overdispersion in logistic regression of malformation counts

source("R/FigureLabel.R") # Functions to make nice labelling of panel plots


# Filenames of pCA, malformation and paleomag data
gcms_fname <- "data/Pollen/GC_MS_results/GCMS_pCA_day_dependent_calibration.csv"
malformations_fname <- "data/Pollen/Malformations/Suigetsu_Malformations.csv"
paleomag_fname <- "data/Paleomag/ADM_GGF100k.csv"


# Pre-process the data
source("R/PreProcessReplicateData.R") # Creates full_output which has grouped values

# Short function to calculate the overdispersion in the binomial modelling of the malformation count data
overdispersion <- function(model) {
  residuals_df <- df.residual(model)
  residuals_pearson <- residuals(model,type="pearson")
  pearson_chisq <- sum(residuals_pearson^2)
  pearson_ratio <- pearson_chisq/residuals_df
  p_value <- pchisq(pearson_chisq , df = residuals_df, lower.tail=FALSE)
  c(chi_sq = pearson_chisq, ratio = pearson_ratio, df = residuals_df, p = p_value)
}

############################################
# Create colour-blind friendly plotting palette
plot_palette <- palette.colors(palette = "Okabe-Ito")
plot_palette <- plot_palette[c(7, 2:5)]
full_output$plotcol <- plot_palette[as.integer(as.factor(full_output$Time))]
legend_txt <- levels(as.factor(full_output$Time))

# The greyscale to use for the points representing the individual replicates
replicates_greyscale_col <- grey(0.1, 0.6)

bespoke_adjust <- rep(1.7, 20)
nudge_y <- -0.0002

stop_overlap <- TRUE
if(stop_overlap) {
  bespoke_adjust[c(5, 8, 15)] <- -10
  nudge_y <- 0
}

#################################################################
####### Commence main analysis and plotting of figures
#################################################################


##################   Plot Fig S3 ##################################
####### Consider experimental space and plot it for the SI Fig S3

# Calculate correlation of paleomag and calage for sampled periods
cor(full_output$expected_paleomag, full_output$calage)

#Plot the design space
source("R/Plot_Design_Space.R")




###############################################################
# Analysis 1: Perform analysis of the GC-MS

source("R/Analyse_pCA_Data_Day_Dependent_Calibration.R")



###############################################################
# Analysis 2: Perform analysis of the grain malformations

cat("Analysis of grain malformation probabilities \n")
cat("------------------- \n")

source("R/Analyse_Malformations.R")


