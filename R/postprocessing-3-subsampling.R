library(bayestestR)
library(dplyr)
## library(stringr)
library(reshape2)
library(ggplot2)
library(svglite) # needed for SVG output
## library(cowplot)
## library(jsonlite)
## library(xml2)
## library(lubridate)
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
  element_rect(colour = "#363636", linewidth = 0.25)
## ============================================================

read_r0_values <- function(beast_log) {
  beast_log |>
    read_beast2_log() |>
    select(starts_with("TTR0")) |>
    mutate(log_file = beast_log) |>
    melt(id.vars = "log_file", variable.name = "parameter", value.name = "value") |>
    group_by(log_file, parameter) |>
    summarise(median = median(value),
              lower = hdi(value, ci = 0.95)$CI_low,
              upper = hdi(value, ci = 0.95)$CI_high,
              .groups = "drop")
}

## ============================================================

output_png <- "out/subsampling-experiment/summary-plot-r0.png"

beast_logs <- c("out/timtam-2023-10-15-original.log",
                "out/timtam-2023-10-15-subsample-0p66.log",
                "out/timtam-2023-10-15-subsample-0p33.log")

## ============================================================

log_labels <- c("out/timtam-2023-10-15-original.log" = "Original",
                "out/timtam-2023-10-15-subsample-0p66.log" = "Subsample 66%",
                "out/timtam-2023-10-15-subsample-0p33.log" = "Subsample 33%")

parameter_labels <- c("TTR0.1" = "1",
                      "TTR0.2" = "2",
                      "TTR0.3" = "3",
                      "TTR0.4" = "4")


r0_post_df <- beast_logs |>
  map(read_r0_values) |>
  bind_rows() |>
  mutate(log_file = factor(log_file, levels = beast_logs))

## ============================================================

r0_gg <-
  ggplot(r0_post_df,
         aes(x = parameter,
                      y = median,
                      ymin = lower,
                      ymax = upper,
                      shape = log_file,
                      colour = log_file)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels = parameter_labels) +
  scale_colour_manual(name = "Time series",
                      labels = log_labels,
                      values = scale_colour_vals) +
  scale_shape_manual(name = "Time series",
                     labels = log_labels,
                     values = scale_shape_vals) +
  labs(x = "Date range",
       y = "Reproduction number") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.8),
        legend.background = legend_background_style)

if (interactive()) {
  print(r0_gg)
} else {
  plot_width <- 12
  plot_height <- 9
  ggsave(output_png,
         r0_gg,
         width = plot_width, height = plot_height, units = "cm")
  ggsave(filename = sub("png$", "svg", output_png),
         plot = r0_gg,
         width = plot_width, height = plot_height, units = "cm")
  ## Write the gg plot to a .rds file to make it easier to revive.
  saveRDS(r0_gg, sub("png$", "rds", output_png))
}
