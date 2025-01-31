library(bayestestR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(svglite) # needed for SVG output
library(cowplot)
library(purrr)
library(timtamslamR)

## ============================================================
## Hard-code visualisation parameters
palette_green <- "#1B9E77"
palette_orange <- "#D95F02"
palette_purple <- "#7570B3"
scale_colour_vals <- c(palette_green, palette_orange, palette_purple)
scale_shape_vals <- c(15, 16, 17)
pointrange_size <- 0.7
legend_background_style <-
  element_rect(colour = "#363636", size = 0.25)
## ============================================================

output_png <- "out/manuscript/subsampling-experiment-combined-r0-timeseries.png"
ts_gg <- readRDS("out/subsampling-experiment/disaster-timeseries.rds")
r0_gg <- readRDS("out/subsampling-experiment/summary-plot-r0.rds")

com_gg <-
  cowplot::plot_grid(
             ts_gg + theme(legend.position.inside = c(0.2, 0.7),
                           legend.title = element_blank()),
             r0_gg + theme(legend.position = "none"),
             labels = c("A", "B"),
             nrow = 1, rel_widths = c(0.6, 0.4))

if (interactive()) {
  print(com_gg)
} else {
  plot_width <- 12 * 2
  plot_height <- 9
  ggsave(output_png, com_gg, dpi = 300,
         width = plot_width, height = plot_height, units = "cm")
  ggsave(filename = sub("png$", "svg", output_png), plot = com_gg,
         width = plot_width, height = plot_height, units = "cm")
  saveRDS(com_gg, sub("png$", "rds", output_png))
}
