# Create plots of the GGF100k paleomag record over the sampled time periods

paleomag <- read.csv("data/Paleomag/ADM_GGF100k.csv", header = TRUE)
paleomag$ZAM2 <- 10 * paleomag$ADM.x10.22.Am.2.
n_intervals <- 4
add_panel_label <- FALSE

source("R/FigureLabel.R")


# Plotting colours for sampled intervals
plot_interval_cols <- hcl.colors(4, alpha = 0.3)

#############################################
# Most recent time period ca 1000 cal yr BP

ylim_plot <- c(20, 110)

# 1st interval
intervals <- seq(from = 1, by = 0.15, length = n_intervals + 1)
xlim_plot <- rev(range(intervals) + c(-1,1))

png("output/Fig1_Panels/Detailed_Paleomag_Sampled_1ka_Period.png", width = 4, height = 4, units = "in", res = 480)

par(
  mgp = c(3, 0.6, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(3.1, 3, 3, 1.8) + 0.1,
  las = 1)


plot(paleomag$Age,
     paleomag$ZAM2,
     type = "l",
     col = "purple",
     xlim = xlim_plot, ylim = ylim_plot,
     ylab = "",
     xlab = "")

mtext(expression(paste("AMD (10"^21, " Am"^2, " / ZAm"^2, ")" )),side = 2, las = 0, line = 1.4, cex = 1)
mtext("Calendar Age (ka BP)", side = 1, line = 1.8, cex = 1)

for(i in 1:n_intervals) {
  polygon(x = c(rep(intervals[i], 2), rep(intervals[i+1], 2)),
          y = c(0, 120, 120, 0),
          col = plot_interval_cols[i])
}

inrange <- which(paleomag$Age < intervals[n_intervals +1] &
                   paleomag$Age > intervals[1])
mean_paleomag <- round(mean(paleomag$ZAM2[inrange]), 2)

mean_paleomag_string <- as.character(bquote(.(round(mean_paleomag, 1))))


# Add title and subtitle underneath
title(main = bquote("Core C-02 (1600 \u2013 1000 yr BP)"), line = 2)
title(main = bquote("Mean field strength" ~ .(round(mean_paleomag, 1)) ~ "ZAm"^2 ), line = 1)


lines(paleomag$Age,
      paleomag$ZAM2,
      col = "black", lwd = 2)

if(add_panel_label) {
  fig_label(LETTERS[2], cex = 2, region = "plot")
}

dev.off()



#########################################
# 2nd interval
intervals <- seq(from = 7.7, by = 0.15, length = n_intervals + 1)
xlim_plot <- rev(range(intervals) + c(-1,1))

png("output/Fig1_Panels/Detailed_Paleomag_Sampled_8ka_Period.png", width = 4, height = 4, units = "in", res = 480)

par(
  mgp = c(3, 0.6, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(3.1, 3, 3, 1.8) + 0.1,
  las = 1)

### Main plot
plot(paleomag$Age,
     paleomag$ZAM2,
     type = "l",
     col = "purple",
     xlim = xlim_plot, ylim = ylim_plot,
     ylab = "",
     xlab = "")

mtext(expression(paste("AMD (10"^21, " Am"^2, " / ZAm"^2, ")" )),side = 2, las = 0, line = 1.4, cex = 1)
mtext("Calendar Age (ka BP)", side = 1, line = 1.8, cex = 1)
#####


for(i in 1:n_intervals) {
  polygon(x = c(rep(intervals[i], 2), rep(intervals[i+1], 2)),
          y = c(0, 120, 120, 0),
          col = plot_interval_cols[i])
}


inrange <- which(paleomag$Age < intervals[n_intervals +1] &
                   paleomag$Age > intervals[1])
mean_paleomag <- round(mean(paleomag$ZAM2[inrange]), 2)

# Add title and subtitle underneath
title(main = bquote("Core C-08 (8300 \u2013 7700 yr BP)"), line = 2)
title(main = bquote("Mean field strength" ~ .(round(mean_paleomag, 1)) ~ "ZAm"^2 ), line = 1)



lines(paleomag$Age,
      paleomag$ZAM2,
      col = "black", lwd = 2)

if(add_panel_label) {
  fig_label(LETTERS[3], cex = 2, region = "plot")
}


dev.off()

#########################################
# 3rd interval
intervals <- seq(from = 40.7, by = 0.15, length = 5)
xlim_plot <- rev(range(intervals) + c(-1,1))

png("output/Fig1_Panels/Detailed_Paleomag_Sampled_41ka_Period.png", width = 4, height = 4, units = "in", res = 480)

par(
  mgp = c(3, 0.6, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(3.1, 3, 3, 1.8) + 0.1,
  las = 1)


### Main plot
plot(paleomag$Age,
     paleomag$ZAM2,
     type = "l",
     col = "purple",
     xlim = xlim_plot, ylim = ylim_plot,
     ylab = "",
     xlab = "")

mtext(expression(paste("AMD (10"^21, " Am"^2, " / ZAm"^2, ")" )),side = 2, las = 0, line = 1.4, cex = 1)
mtext("Calendar Age (ka BP)", side = 1, line = 1.8, cex = 1)
#####



for(i in 1:4) {
  polygon(x = c(rep(intervals[i], 2), rep(intervals[i+1], 2)),
          y = c(0, 120, 120, 0),
          col = plot_interval_cols[i])
}

inrange <- which(paleomag$Age < intervals[n_intervals +1] &
                   paleomag$Age > intervals[1])
mean_paleomag <- round(mean(paleomag$ZAM2[inrange]), 2)

# Add title and subtitle underneath
title(main = bquote("Core C-13 (41300 \u2013 40700 yr BP)"), line = 2)
title(main = bquote("Mean field strength" ~ .(round(mean_paleomag, 1)) ~ "ZAm"^2 ), line = 1)


lines(paleomag$Age,
      paleomag$ZAM2,
      col = "black", lwd = 2)

if(add_panel_label) {
  fig_label(LETTERS[4], cex = 2, region = "plot")
}

dev.off()

#########################################
# 4th interval
intervals <- seq(from = 52.7, by = 0.15, length = 5)
xlim_plot <- rev(range(intervals) + c(-1,1))

png("output/Fig1_Panels/Detailed_Paleomag_Sampled_53ka_Period.png", width = 4, height = 4, units = "in", res = 480)

par(
  mgp = c(3, 0.6, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(3.1, 3, 3, 1.8) + 0.1,
  las = 1)


### Main plot
plot(paleomag$Age,
     paleomag$ZAM2,
     type = "l",
     col = "purple",
     xlim = xlim_plot, ylim = ylim_plot,
     ylab = "",
     xlab = "")

mtext(expression(paste("AMD (10"^21, " Am"^2, " / ZAm"^2, ")" )),side = 2, las = 0, line = 1.4, cex = 1)
mtext("Calendar Age (ka BP)", side = 1, line = 1.8, cex = 1)
#####



for(i in 1:4) {
  polygon(x = c(rep(intervals[i], 2), rep(intervals[i+1], 2)),
          y = c(0, 120, 120, 0),
          col = plot_interval_cols[i])
}

inrange <- which(paleomag$Age < intervals[n_intervals +1] &
                   paleomag$Age > intervals[1])
mean_paleomag <- round(mean(paleomag$ZAM2[inrange]), 2)


# Add title and subtitle underneath
title(main = bquote("Core C-17 (53300 \u2013 52700 yr BP)"), line = 2)
title(main = bquote("Mean field strength" ~ .(round(mean_paleomag, 1)) ~ "ZAm"^2 ), line = 1)




lines(paleomag$Age,
      paleomag$ZAM2,
      col = "black", lwd = 2)

if(add_panel_label) {
  fig_label(LETTERS[5], cex = 2, region = "plot")
}


dev.off()

#########################################
# 5th interval
intervals <- seq(from = 64.75, by = 0.15, length = 5)
xlim_plot <- rev(range(intervals) + c(-1,1))

png("output/Fig1_Panels/Detailed_Paleomag_Sampled_65ka_Period.png", width = 4, height = 4, units = "in", res = 480)

par(
  mgp = c(3, 0.6, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(3.1, 3, 3, 1.8) + 0.1,
  las = 1)


### Main plot
plot(paleomag$Age,
     paleomag$ZAM2,
     type = "l",
     col = "purple",
     xlim = xlim_plot, ylim = ylim_plot,
     ylab = "",
     xlab = "")

mtext(expression(paste("AMD (10"^21, " Am"^2, " / ZAm"^2, ")" )),side = 2, las = 0, line = 1.4, cex = 1)
mtext("Calendar Age (ka BP)", side = 1, line = 1.8, cex = 1)
#####



for(i in 1:4) {
  polygon(x = c(rep(intervals[i], 2), rep(intervals[i+1], 2)),
          y = c(0, 120, 120, 0),
          col = plot_interval_cols[i])
}


inrange <- which(paleomag$Age < intervals[n_intervals +1] &
                   paleomag$Age > intervals[1])
mean_paleomag <- round(mean(paleomag$ZAM2[inrange]), 2)

# Add title and subtitle underneath
title(main = bquote("Core B-22 (65350 \u2013 64750 yr BP)"), line = 2)
title(main = bquote("Mean field strength" ~ .(round(mean_paleomag, 1)) ~ "ZAm"^2 ), line = 1)


lines(paleomag$Age,
      paleomag$ZAM2,
      col = "black", lwd = 2)

if(add_panel_label) {
  fig_label(LETTERS[6], cex = 2, region = "plot")
}


#####

dev.off()




