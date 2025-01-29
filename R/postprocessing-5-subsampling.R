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

palette_green <- "#1B9E77"
palette_orange <- "#D95F02"
palette_purple <- "#7570B3"

read_prop_ts_values <- function(beast_log) {
  beast_log |>
    read_beast2_log() |>
    select(TTPropTS.2) |>
    mutate(log_file = beast_log) |>
    melt(id.vars = "log_file", variable.name = "parameter", value.name = "value") |>
    group_by(log_file, parameter) |>
    summarise(median = median(value),
              lower = hdi(value, ci = 0.95)$CI_low,
              upper = hdi(value, ci = 0.95)$CI_high) |>
    ungroup()
}

## ============================================================

output_png <- "out/subsampling-experiment/summary-plot-prop_ts.png"

beast_logs <- c("out/timtam-2023-10-15-original.log",
                "out/timtam-2023-10-15-subsample-0p66.log",
                "out/timtam-2023-10-15-subsample-0p33.log")

## ============================================================

log_labels <- c("out/timtam-2023-10-15-original.log" = "Original",
                "out/timtam-2023-10-15-subsample-0p66.log" = "Subsample 66%",
                "out/timtam-2023-10-15-subsample-0p33.log" = "Subsample 33%")

parameter_labels <- c("TTPropTS.2" = "2")


prop_ts_post_df <- beast_logs |>
  map(read_prop_ts_values) |>
  bind_rows() |>
  mutate(log_file = factor(log_file, levels = beast_logs))

## ============================================================

sensible_upper_y_lim <- round(prop_ts_post_df$upper * 1100) / 1000

prop_ts_gg <-
  ggplot() +
  geom_pointrange(data = prop_ts_post_df,
                  aes(x = parameter,
                      y = median,
                      ymin = lower,
                      ymax = upper,
                      colour = log_file),
                  position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels = parameter_labels) +
  scale_colour_manual(labels = log_labels,
                      values = c(palette_green, palette_orange, palette_purple)) +
  scale_y_continuous(limits = c(0, sensible_upper_y_lim)) +
  labs(x = "Date range",
       y = "Proportion of cases in time series",
       colour = "Time series") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2, 0.3),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

## print(prop_ts_gg)
plot_width <- 12
plot_height <- 9
ggsave(output_png, prop_ts_gg, width = plot_width, height = plot_height, units = "cm")
ggsave(filename = sub("png$", "svg", output_png), plot = prop_ts_gg, width = plot_width, height = plot_height, units = "cm")
