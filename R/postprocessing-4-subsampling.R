library(bayestestR)
library(dplyr)
## library(stringr)
library(reshape2)
library(ggplot2)
library(svglite) # needed for SVG output
## library(cowplot)
## library(jsonlite)
library(xml2)
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

read_hs_values <- function(beast_log, history_times) {
  tmp <- beast_log |>
    read_beast2_log() |>
    select(matches("HistorySize")) |>
    mutate(log_file = beast_log) |>
    melt(id.vars = "log_file", variable.name = "parameter", value.name = "value") |>
    group_by(log_file, parameter) |>
    summarise(median = median(value),
              lower = bayestestR::hdi(value, ci = 0.95, verbose = FALSE)[["CI_low"]],
              upper = bayestestR::hdi(value, ci = 0.95, verbose = FALSE)[["CI_high"]],
              ## There is a warning that gets thrown by the hdi
              ## function because the input are discrete. You can fuzz
              ## them with the code below to get the same results but
              ## without the wraning. Or you can just turn the warning
              ## off.
              ##
              ## lower = hdi(pmax(value + runif(length(value), min=-0.1, max=0.1), 0), ci = 0.95)[["CI_low"]],
              ## upper = hdi(pmax(value + runif(length(value), min=-0.1, max=0.1), 0), ci = 0.95)[["CI_high"]],
              .groups = "drop")
  tmp$history_times <- history_times
  return(tmp)
}

## ============================================================

beast_xml_orig <- "out/subsampling-experiment/xml/timtam-2023-10-15-original.xml"
output_png <- "out/subsampling-experiment/summary-plot-historysizes.png"

beast_logs <- c("out/timtam-2023-10-15-original.log",
                "out/timtam-2023-10-15-subsample-0p66.log",
                "out/timtam-2023-10-15-subsample-0p33.log")

## ============================================================
orig_xml <- "out/subsampling-experiment/xml/timtam-2023-10-15-original.xml"
stopifnot(file.exists(orig_xml))
orig_node <- xml2::read_xml(orig_xml)
r0_change_bkwd_times <-
  orig_node |>
  xml2::xml_find_first("//parameter[@name='r0ChangeTimes']") |>
  xml2::xml_text() |>
  strsplit(" ") |>
  purrr::pluck(1) |>
  as.numeric()


stopifnot(file.exists(beast_xml_orig))
beast_model <- xml2::as_list(xml2::read_xml(beast_xml_orig))
history_times_str <-
  beast_model$beast$run$distribution$distribution$distribution[5]$parameter[1][[1]]
history_times_num <-
  history_times_str |>
  strsplit(split = " ") |>
  pluck(1) |>
  as.numeric()

log_labels <- c(
  "out/timtam-2023-10-15-original.log" = "Original",
  "out/timtam-2023-10-15-subsample-0p66.log" = "Subsample 66%",
  "out/timtam-2023-10-15-subsample-0p33.log" = "Subsample 33%"
)


hs_post_df <- beast_logs |>
  map(read_hs_values, history_times = history_times_num) |>
  bind_rows() |>
  mutate(log_file = factor(log_file, levels = beast_logs))

## ============================================================

## This warning is fine, it is expected when there are zero counts in
## there. These will lock to the bottom of the plot, it is not a
## problem.
warning("DO NOT WORRY about a warning that scale_y_log10() produced infinite values.")

hs_gg <-
  ggplot(data = hs_post_df,
                  aes(x = history_times,
                      y = median,
                      ymin = lower,
                      ymax = upper,
                      shape = log_file,
                      colour = log_file)) +
  geom_vline(xintercept = r0_change_bkwd_times) +
  geom_pointrange(position = position_dodge(width = 10)) +
  scale_x_reverse() +
  scale_y_log10(
    breaks = 10^(0:4),
    labels = scales::label_log()
  ) +
  scale_colour_manual(name = "Time series",
                      labels = log_labels,
                      values = scale_colour_vals) +
  scale_shape_manual(name = "Time series",
                     labels = log_labels,
                     values = scale_shape_vals) +
  labs(x = "Backwards time (days)",
       y = "Hidden lineages") +
  theme_bw() +
  theme(legend.position = "none")


if (interactive()) {
  print(hs_gg)
} else {
  plot_width <- 12
  plot_height <- 9
  ggsave(output_png, hs_gg, width = plot_width, height = plot_height, units = "cm")
  ggsave(filename = sub("png$", "svg", output_png), plot = hs_gg, width = plot_width, height = plot_height, units = "cm")
  saveRDS(hs_gg, sub("png$", "rds", output_png))
}
