# SI Plot: Create plot of experimental space
design_space <- ggplot(full_output,
                       aes(x = expected_paleomag, y = calage)) +
  geom_path(aes(x = ZAM2, y = calage), data = paleomag[paleomag$calage < 66000,], colour = grey(0.15, 0.3)) +
  geom_point(aes(fill = substr(Time, 1, 5)), size = 4, pch = 22) +
  xlab(expression(paste("Geomagnetic Axial Dipole Moment ( x 10"^21," Am"^2," / ZAm"^2, ")")))+
  ylab("Calendar Age (yr BP)") +
  scale_fill_manual(values = as.character(plot_palette)) +
  labs(fill = "Time Period") +
  theme_bw(
    base_size = 9
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.05),
    legend.justification = c("left", "bottom"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6)
  )

design_space


# use function ggsave to export plot to jpeg...
ggsave("output/FigS3_UV_B_Analysis_Design_Space.jpg", width = 0.7 * 5.58, height = 0.7 * 5.58, units = "in")
