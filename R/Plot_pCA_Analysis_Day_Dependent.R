############################################################
# Plotting of p-CA Analysis:
# Fig S3: Residuals of the linear model
# Fig 2: Main manuscript figure of p-CA vs paleomag strength and time period
############################################################

# Choose plotting region
fitted_xlim_plot <- c()
scale_location_legend_pos <- "bottomleft"
adjusted_location_legend_pos_x <- 0.47

# Choose where to put legend so doesn't cover data
tukey_adjusted_location_legend_pos_x <- -0.03
sca_loc_adjusted_location_legend_pos_x <- 0.5
lev_adjusted_location_legend_pos_y <- 2.6



TukeyPlot <- function() {
  plot(pca_linear_model_replicate_level,
       which = 1,
       pch = 19,
       col = pca_individual_replicate_data$plotcol,
       xlim = fitted_xlim_plot,
       labels.id = paste(as.character(pca_individual_replicate_data$calage), "yr BP"),
       sub.caption = expression(paste(italic("p"),"-CA"[i], " = ", beta[0], " + ", beta[1], "Paleomag"[i], " + ", beta[2], "CalAge"[i], " + ", epsilon[i])) )
  legend(x = adjusted_location_legend_pos_x, y = tukey_adjusted_location_legend_pos_x, pch = 19, legend = legend_txt, col = plot_palette, cex = 0.75)
  fig_label(LETTERS[1], cex = 2, region = "plot")
}

QQPlot <- function() {
  plot(pca_linear_model_replicate_level,
       which = 2,
       pch = 19,
       col = pca_individual_replicate_data$plotcol,
       labels.id = paste(as.character(pca_individual_replicate_data$calage), "yr BP"),
       sub.caption = expression(paste(italic("p"),"-CA"[i], " = ", beta[0], " + ", beta[1], "Paleomag"[i], " + ", beta[2], "CalAge"[i], " + ", epsilon[i])) )
  legend("bottomright", pch = 19, legend = legend_txt, col = plot_palette, cex = 0.75)
  fig_label(LETTERS[2], cex = 2, region = "plot")
}

ScaleLocationPlot <- function() {
  plot(pca_linear_model_replicate_level,
       which = 3,
       pch = 19,
       col = pca_individual_replicate_data$plotcol,
       xlim = fitted_xlim_plot,
       labels.id = paste(as.character(pca_individual_replicate_data$calage), "yr BP"),
       sub.caption = expression(paste(italic("p"),"-CA"[i], " = ", beta[0], " + ", beta[1], "Paleomag"[i], " + ", beta[2], "CalAge"[i], " + ", epsilon[i])) )
  legend(x = adjusted_location_legend_pos_x, y = sca_loc_adjusted_location_legend_pos_x, pch = 19, legend = legend_txt, col = plot_palette, cex = 0.75)
  fig_label(LETTERS[3], cex = 2, region = "plot")
}


LeveragePlot <- function() {
  plot(pca_linear_model_replicate_level,
       which = 5,
       pch = 19,
       col = pca_individual_replicate_data$plotcol,
       labels.id = paste(as.character(pca_individual_replicate_data$calage), "yr BP"),
       label.pos = c(4, 2),
       sub.caption = expression(paste(italic("p"),"-CA"[i], " = ", beta[0], " + ", beta[1], "Paleomag"[i], " + ", beta[2], "CalAge"[i], " + ", epsilon[i])) )
  legend(x = 0.002, y = lev_adjusted_location_legend_pos_y, pch = 19, legend = legend_txt, col = plot_palette, cex = 0.75)
  fig_label(LETTERS[4], cex = 2, region = "plot")
}



png(filename = paste("output/FigS4_UV_B_LinearModel_ResidualsFit_Replicate_Level.png", sep = ""),
    width = 1.2 * 8.59, height = 1.2 * 5.58,
    units = "in",
    res = 480)
par(mfrow = c(2, 2),
    mar = c(3, 3, 1.5, 1.5),
    mgp = c(2, 1, 0),
    oma = c(0, 0, 2, 0)) -> opar
TukeyPlot()
QQPlot()
ScaleLocationPlot()
LeveragePlot()
dev.off()




##################################################################
# Create Fig 2 plots of pCA for main manuscript:
# Panel A - pCA vs paleomag
# Panel B - pCA vs time period
##################################################################

#### Panel A - plot pCA vs paleomag (show all individual samples)
# Note that the smooth fitted is against the individual repicates (the correct thing to show)
gcms_lm_plot_legend <- full_output %>%
  ggplot(aes(x = expected_paleomag, y = mean_gcms)) +
  geom_errorbar(aes(ymax = upper, ymin = lower), colour = full_output$plotcol, width= 0.2)+
  geom_smooth(method = "lm",
              data = pca_individual_replicate_data,
              aes(x = expected_paleomag, y = calibrated),
              col= "red", alpha = 0.5) +
  geom_point(aes(fill = substr(Time, 1, 5)), size = 4, pch = 22) +
  theme_bw(base_size = 9) +
  xlab(expression(paste("Geomagnetic Axial Dipole Moment ( x 10"^21," Am"^2," / ZAm"^2, ")")))+
  ylab(expression(paste(italic("para"),"-Coumaric Acid (ng grain"^-1,")"))) +
  scale_fill_manual(values = as.character(plot_palette)) +
  labs(fill = "Time Period") +
  theme(
    legend.position = c(0.99, 0.99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(1, 1, 1, 1)
  )



# Add points showing individual replicates
gcms_lm_plot_legend <- gcms_lm_plot_legend +
  annotate("point",
           shape = 16,
           colour = replicates_greyscale_col,
           y = pca_individual_replicate_data$calibrated,
           x = pca_individual_replicate_data$expected_paleomag)


gcms_lm_plot_legend




########################################################################
#### Panel B - plot pCa vs Time Period (show mean pCA in each sample)
boxplot_width <- 0.8/2
big_text <- TRUE
if(big_text) {
  calage_textsize <- 2
  text_size_malformed <- 3
  text_size_subscript_malformed <- 3.5
} else {
  calage_textsize <- 1.6
  text_size_malformed <- 1.6
  text_size_subscript_malformed <- 1.6
}


n_obs_per_period <- 4
ticks <- seq(-boxplot_width, boxplot_width, length = n_obs_per_period + 1)
laschamps_ticks <- 3 + ticks
box_plot_locations <- -boxplot_width + ((1:n_obs_per_period) - 0.5) * (2 * boxplot_width) / n_obs_per_period

gcms_laschamps_ticks <- laschamps_ticks
gcms_laschamps_label_x_val <- 0.37


gcms_time_period_plot_all_replicates <- ggplot(pca_individual_replicate_data, aes(x = substr(Time, 1, 5), y = calibrated)) +
  geom_boxplot(aes(fill = substr(Time, 1, 5)), width = 2 * boxplot_width, outlier.shape = NA) +
  scale_fill_manual(values = as.character(plot_palette)) +
  annotate("point",
           shape = 16,
           colour = replicates_greyscale_col,
           y = pca_individual_replicate_data$calibrated,
           x = as.integer(pca_individual_replicate_data$Time) + box_plot_locations[pca_individual_replicate_data$order_within_period]) +
  geom_errorbarh(aes(xmin = gcms_laschamps_ticks[1], xmax = gcms_laschamps_ticks[2], y = gcms_laschamps_label_x_val),
                 width = 0.002) +
  geom_errorbarh(aes(xmin = gcms_laschamps_ticks[2], xmax = gcms_laschamps_ticks[3], y = gcms_laschamps_label_x_val),
                 width = 0.002) +
  geom_errorbarh(aes(xmin = gcms_laschamps_ticks[3], xmax = gcms_laschamps_ticks[4], y = gcms_laschamps_label_x_val),
                 width = 0.002) +
  geom_errorbarh(aes(xmin = gcms_laschamps_ticks[4], xmax = gcms_laschamps_ticks[5], y = gcms_laschamps_label_x_val),
                 width = 0.002) +
  labs(fill = "Time Period") +
  ylab(expression(paste(italic("para"),"-Coumaric Acid (ng grain"^-1,")"))) +
  xlab("Time Period") +
  theme_bw(
    base_size = 9
  ) +
  theme(legend.position="none") +
  annotate("text", label = as.character(seq(40700, 41300, length = n_obs_per_period + 1))[c(1,3,5)],
           x = laschamps_ticks[c(1,3,5)],
           y = rep(gcms_laschamps_label_x_val, n_obs_per_period + 1)[c(1,3,5)] - 0.004, size = calage_textsize, check_overlap = TRUE) +
  annotate("text", label = as.character(seq(40700, 41300, length = n_obs_per_period + 1))[c(2,4)],
           x = laschamps_ticks[c(2,4)],
           y = rep(gcms_laschamps_label_x_val, n_obs_per_period + 1)[c(2,4)] + 0.004, size = calage_textsize, check_overlap = TRUE)

# Also add intervals for the other periods
time_periods_start <- c(1000, 7700, 40700, 52700, 64750)
time_periods_end <- c(1600, 8300, 41300, 53300, 65350)
max_gcms_by_period <- by(data = pca_individual_replicate_data[,c("calibrated")], pca_individual_replicate_data$Time, max)

for(i in c(1,2,4,5)) {
  gcms_time_period_plot_all_replicates <- gcms_time_period_plot_all_replicates +
    annotate("text", label = as.character(seq(time_periods_start[i],time_periods_end[i], length = n_obs_per_period + 1))[c(1,3,5)],
             x = laschamps_ticks[c(1,3,5)] + (i-3),
             y = rep(max_gcms_by_period[i] + 0.01, n_obs_per_period + 1)[c(1,3,5)] - 0.004, size = calage_textsize, check_overlap = TRUE) +
    annotate("text", label = as.character(seq(time_periods_start[i],time_periods_end[i], length = n_obs_per_period + 1))[c(2,4)],
             x = laschamps_ticks[c(2,4)] + (i-3),
             y = rep(max_gcms_by_period[i] + 0.01, n_obs_per_period + 1)[c(2,4)] + 0.004, size = calage_textsize, check_overlap = TRUE) +
    annotate("errorbarh", xmin = gcms_laschamps_ticks[1] + (i-3), xmax = gcms_laschamps_ticks[2] + (i-3), y = max_gcms_by_period[i] + 0.01,
             width = 0.002) +
    annotate("errorbarh", xmin = gcms_laschamps_ticks[2] + (i-3), xmax = gcms_laschamps_ticks[3] + (i-3), y = max_gcms_by_period[i] + 0.01,
             width = 0.002) +
    annotate("errorbarh", xmin = gcms_laschamps_ticks[3] + (i-3), xmax = gcms_laschamps_ticks[4] + (i-3), y = max_gcms_by_period[i] + 0.01,
             width = 0.002) +
    annotate("errorbarh", xmin = gcms_laschamps_ticks[4] + (i-3), xmax = gcms_laschamps_ticks[5] + (i-3), y = max_gcms_by_period[i] + 0.01,
             width = 0.002)
}


#############################################################
# Create two panel plot of pCA data showing all replicates
plot_grid(gcms_lm_plot_legend,
          gcms_time_period_plot_all_replicates,
          ncol = 2, nrow= 1,
          labels = c("A", "B"))
ggsave(filename = paste("output/Fig2_UV_B_Analysis_pCA_All_Replicates.jpg", sep = ""),
       width = 8.59, height = 5.58, units = "in")

