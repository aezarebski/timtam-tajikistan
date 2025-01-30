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

read_hs_values <- function(beast_log) {
  beast_log |>
    read_beast2_log() |>
    select(matches("HistorySize")) |>
    mutate(log_file = beast_log) |>
    melt(id.vars = "log_file", variable.name = "parameter", value.name = "value") |>
    group_by(log_file, parameter) |>
    summarise(median = median(value),
              lower = hdi(value, ci = 0.95)$CI_low,
              upper = hdi(value, ci = 0.95)$CI_high) |>
    ungroup()
}

## ============================================================

output_png <- "out/subsampling-experiment/summary-plot-historysizes.png"

beast_logs <- c("out/timtam-2023-10-15-original.log",
                "out/timtam-2023-10-15-subsample-0p66.log",
                "out/timtam-2023-10-15-subsample-0p33.log")

## ============================================================

log_labels <- c("out/timtam-2023-10-15-original.log" = "Original",
                "out/timtam-2023-10-15-subsample-0p66.log" = "Subsample 66%",
                "out/timtam-2023-10-15-subsample-0p33.log" = "Subsample 33%")

parameter_labels <- c("TTHistorySizes.1" = "1",
                      "TTHistorySizes.2" = "2",
                      "TTHistorySizes.3" = "3",
                      "TTHistorySizes.4" = "4",
                      "TTHistorySizes.5" = "5",
                      "TTHistorySizes.6" = "6",
                      "TTHistorySizes.7" = "7",
                      "TTHistorySizes.8" = "8",
                      "TTHistorySizes.9" = "9",
                      "TTHistorySizes.10" = "10",
                      "TTHistorySizes.11" = "11",
                      "TTHistorySizes.12" = "12",
                      "TTHistorySizes.13" = "13")


hs_post_df <- beast_logs |>
  map(read_hs_values) |>
  bind_rows() |>
  mutate(log_file = factor(log_file, levels = beast_logs))

## ============================================================

hs_gg <-
  ggplot() +
  geom_pointrange(data = hs_post_df,
                  aes(x = parameter,
                      y = median,
                      ymin = lower,
                      ymax = upper,
                      colour = log_file),
                  position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels = parameter_labels) +
  scale_y_log10(
    breaks = c(0, 10^(0:4)),
    labels = scales::label_log()
  ) +
  scale_colour_manual(labels = log_labels,
                      values = c(palette_green, palette_orange, palette_purple)) +
  labs(x = "Date",
       y = "History size parameter",
       colour = "Time series") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.background = element_rect(colour = "#eaeaea"))

## print(hs_gg)
plot_width <- 12
plot_height <- 9
ggsave(output_png, hs_gg, width = plot_width, height = plot_height, units = "cm")
ggsave(filename = sub("png$", "svg", output_png), plot = hs_gg, width = plot_width, height = plot_height, units = "cm")
