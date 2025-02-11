library(bayestestR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(svglite) # needed for SVG output
library(cowplot)
library(purrr)
library(xml2)
library(timtamslamR)

## ============================================================ The
## current version of `cow_plot' needs the `get_legend' function
## updated. Here is a corrected version from the author as given here:
## https://github.com/wilkelab/cowplot/issues/202#issuecomment-1981765769

get_legend_35 <- function(plot, legend_number = 1) {
  # find all legend candidates
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  # find non-zero legends
  idx <- which(vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE))
  # return either the chosen or the first non-zero legend if it exists,
  # and otherwise the first element (which will be a zeroGrob)
  if (length(idx) >= legend_number) {
    return(legends[[idx[legend_number]]])
  } else if (length(idx) >= 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}
## ============================================================
## Hard-code visualisation parameters
palette_green <- "#1B9E77"
palette_orange <- "#D95F02"
palette_purple <- "#7570B3"
scale_colour_vals <- c(palette_green, palette_orange, palette_purple)
scale_shape_vals <- c(15, 16, 17)
pointrange_size <- 0.7
legend_background_style <-
  element_rect(colour = "#363636", linewidth = 0.25)
## ============================================================


output_png <- "out/manuscript/subsampling-experiment-combined-r0-timeseries.png"
ts_gg <- readRDS("out/subsampling-experiment/disaster-timeseries.rds")
r0_gg <- readRDS("out/subsampling-experiment/summary-plot-r0.rds")
hl_gg <- readRDS("out/subsampling-experiment/summary-plot-historysizes.rds")

no_lgnd <- theme(legend.position = "none")

## We use the `get_legend_35' function because `cow_plot::get_legend'
## is broken.
lgnd <- get_legend_35(ts_gg + theme(legend.background = element_rect(colour = "#ffffff")))
shifted_lgnd <- cowplot::ggdraw() +
  cowplot::draw_plot(lgnd, x = 0.1, y = -0.1,
                     width = 1, height = 1)

com_row_1 <-
  cowplot::plot_grid(
             ts_gg + no_lgnd,
             r0_gg + no_lgnd,
             nrow = 1,
             labels = c("A", "B"),
             label_size = 16,
             rel_widths = c(0.6, 0.4)
             )

com_row_2 <-
  cowplot::plot_grid(
             hl_gg + no_lgnd,
             shifted_lgnd,
             nrow = 1,
             labels = c("C"),
             rel_widths = c(0.7, 0.3)
           )

com_pre_gg <-
  cowplot::plot_grid(
             com_row_1,
             com_row_2,
             nrow = 2,
             rel_heights = c(0.5, 0.5)
           )
com_gg <- ggdraw() +
  draw_plot(com_pre_gg) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

if (interactive()) {
  print(com_gg)
} else {
  plot_width <- 12 * 2
  plot_height <- 9 * 2
  ggsave(output_png, com_gg, dpi = 300,
         width = plot_width, height = plot_height, units = "cm")
  ggsave(filename = sub("png$", "svg", output_png), plot = com_gg,
         width = plot_width, height = plot_height, units = "cm")
  saveRDS(com_gg, sub("png$", "rds", output_png))
}
