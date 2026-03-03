############################################################
# Plotting of Grain Malformations Analysis:
# Fig 3: Main manuscript figure of grain malformations
############################################################


## Now create the plot
boxplot_width <- 0.85/2
calage_textsize <- 2.2
n_obs_per_period <- 4
ticks <- seq(-boxplot_width, boxplot_width, length = n_obs_per_period + 1)
box_plot_locations <- -boxplot_width + ((1:n_obs_per_period) - 0.5) * (2 * boxplot_width) / n_obs_per_period


# Plot every time period with smaller box plots
all_periods_malformation_laschamps_plot <- ggplot(malformation_rate_output_by_period, aes(x = substr(period, 1, 5), y = prob_malform)) +
  geom_crossbar(aes(ymax = prob_malform_upper_CI, ymin = prob_malform_lower_CI, fill = period),
                linewidth = 0.5, width = 2 * boxplot_width) +
  labs(fill = "Time Period") +
  scale_fill_manual(values = as.character(plot_palette)) +
  ylab("Probability of Grain Malformation") +
  xlab("Time Period") +
  theme_bw(
    base_size = 9
  ) +
  theme(legend.position="none")


# Add the individual observed proportions for each sample
all_periods_malformation_laschamps_plot <- all_periods_malformation_laschamps_plot +
  annotate("point", shape = 16, y = pollen_malformations$prop_malformed, x = as.integer(as.factor(pollen_malformations$Time)) + box_plot_locations)

max_malformation_by_period <- malformation_rate_output_by_period$prob_malform_upper_CI
cal_age_axis_adjust_malformation <- 0.0045

malformation_laschamps_axis_height <- malformation_rate_output_by_period$prob_malform_upper_CI[3] + cal_age_axis_adjust_malformation
i <- 3
# Now add calendar age bars for Laschamps
all_periods_malformation_laschamps_plot <- all_periods_malformation_laschamps_plot +
  annotate("text", label = as.character(seq(40700, 41300, length = n_obs_per_period + 1))[c(1,3,5)],
           x = gcms_laschamps_ticks[c(1,3,5)],
           y = rep(malformation_laschamps_axis_height, n_obs_per_period + 1)[c(1,3,5)] - 0.0008, size = calage_textsize, check_overlap = TRUE) +
  annotate("text", label = as.character(seq(40700, 41300, length = n_obs_per_period + 1))[c(2,4)],
           x = gcms_laschamps_ticks[c(2,4)],
           y = rep(malformation_laschamps_axis_height, n_obs_per_period + 1)[c(2,4)] + 0.0008, size = calage_textsize, check_overlap = TRUE) +
  annotate("errorbarh", xmin = gcms_laschamps_ticks[1] + (i-3), xmax = gcms_laschamps_ticks[2] + (i-3), y = malformation_laschamps_axis_height,
           width = 0.0004) +
  annotate("errorbarh", xmin = gcms_laschamps_ticks[2] + (i-3), xmax = gcms_laschamps_ticks[3] + (i-3), y = malformation_laschamps_axis_height,
           width = 0.0004) +
  annotate("errorbarh", xmin = gcms_laschamps_ticks[3] + (i-3), xmax = gcms_laschamps_ticks[4] + (i-3), y = malformation_laschamps_axis_height,
           width = 0.0004) +
  annotate("errorbarh", xmin = gcms_laschamps_ticks[4] + (i-3), xmax = gcms_laschamps_ticks[5] + (i-3), y = malformation_laschamps_axis_height,
           width = 0.0004)

poo_temp <- all_periods_malformation_laschamps_plot +
 geom_errorbarh(aes(xmin = gcms_laschamps_ticks[1] + (i-3), xmax = gcms_laschamps_ticks[2] + (i-3), y = malformation_laschamps_axis_height),
                width = 0.0004)

for(i in c(1,2,4,5)) {
  all_periods_malformation_laschamps_plot <- all_periods_malformation_laschamps_plot +
    annotate("text", label = as.character(seq(time_periods_start[i],time_periods_end[i], length = n_obs_per_period + 1))[c(1,3,5)],
             x = gcms_laschamps_ticks[c(1,3,5)] + (i-3),
             y = rep(max_malformation_by_period[i] + cal_age_axis_adjust_malformation, n_obs_per_period + 1)[c(1,3,5)] - 0.0008, size = calage_textsize, check_overlap = TRUE) +
    annotate("text", label = as.character(seq(time_periods_start[i],time_periods_end[i], length = n_obs_per_period + 1))[c(2,4)],
             x = gcms_laschamps_ticks[c(2,4)] + (i-3),
             y = rep(max_malformation_by_period[i] + cal_age_axis_adjust_malformation, n_obs_per_period + 1)[c(2,4)] + 0.0008, size = calage_textsize, check_overlap = TRUE) +
    annotate("errorbarh", xmin = gcms_laschamps_ticks[1] + (i-3), xmax = gcms_laschamps_ticks[2] + (i-3), y = max_malformation_by_period[i] + cal_age_axis_adjust_malformation,
             width = 0.0004) +
    annotate("errorbarh", xmin = gcms_laschamps_ticks[2] + (i-3), xmax = gcms_laschamps_ticks[3] + (i-3), y = max_malformation_by_period[i] + cal_age_axis_adjust_malformation,
             width = 0.0004) +
    annotate("errorbarh", xmin = gcms_laschamps_ticks[3] + (i-3), xmax = gcms_laschamps_ticks[4] + (i-3), y = max_malformation_by_period[i] + cal_age_axis_adjust_malformation,
             width = 0.0004) +
    annotate("errorbarh", xmin = gcms_laschamps_ticks[4] + (i-3), xmax = gcms_laschamps_ticks[5] + (i-3), y = max_malformation_by_period[i] + cal_age_axis_adjust_malformation,
             width = 0.0004)
}



# Add numbers of malformations on plot
height_malformed_text <- min(malformation_rate_output_by_period$prob_malform_lower_CI) - 0.001
height_malformed_adjust <- 0.0015
all_periods_malformation_laschamps_plot <- all_periods_malformation_laschamps_plot +
  annotate("text", label = "paste(n[malformed])", parse = TRUE,
           x = 0.9 + box_plot_locations[1] - (3 * boxplot_width) / n_obs_per_period ,
           y = height_malformed_text, size = text_size_subscript_malformed, hjust = "left") +
  annotate("text", label = "paste(n[grain])", parse = TRUE,
           x = 0.9 + box_plot_locations[1] - (3 * boxplot_width) / n_obs_per_period,
           y = height_malformed_text - height_malformed_adjust, size = text_size_subscript_malformed, hjust = "left") +
  expand_limits(x = 0.9 + box_plot_locations[1] - (3 * boxplot_width) / n_obs_per_period - 0.05)


n_malformation_by_period <- by(data = full_output[,c("combined_n_grains_malformed")], full_output$Time, sum)
n_grains_by_period <- by(data = full_output[,c("combined_n_grains")], full_output$Time, sum)


for(i in 1:5) {
  all_periods_malformation_laschamps_plot <- all_periods_malformation_laschamps_plot +
    annotate("text", label = bquote(.(n_malformation_by_period[i])),
             x = i,
             y = height_malformed_text, size = text_size_malformed) +
    annotate("text", label = bquote(.(n_grains_by_period[i])),
             x = i,
             y = height_malformed_text - height_malformed_adjust, size = text_size_malformed)
}

# Add legend in top right corner (just for clarity)
all_periods_malformation_laschamps_plot <- all_periods_malformation_laschamps_plot +
  labs(fill = "Time Period") +
  theme(
    legend.position = c(0.99, 0.99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(1, 1, 1, 1)
  )




# Create plot
all_periods_malformation_laschamps_plot
ggsave(filename = paste("output/Fig3_UV_B_Analysis_Malformations_by_Period.jpg", sep = ""),
       width = 8.59, height = 5.58, units = "in")

