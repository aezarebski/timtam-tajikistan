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
library(furrr)
library(timtamslamR)


warning("DO NOT WORRY about a warning that there are random numbers without a seed.")
future::plan(future::multisession, workers = parallel::detectCores() - 2)

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

prev_summary_values <- function(prev_df, model_xml, history_times) {
  tmp <-
    prev_df |>
    dplyr::select(dplyr::matches("TTPrevalence")) |>
    mutate(xml_file = model_xml) |>
    melt(id.vars = "xml_file", variable.name = "prevalence", value.name = "value") |>
    group_by(xml_file, prevalence) |>
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
config <- as_list(read_xml("config.xml"))
input_present_rds <- config$files$results$intermediate$present[[1]]
rm(config)

output_png <- "out/subsampling-experiment/summary-plot-historysizes.png"

beast_xmls <- c(
  "out/subsampling-experiment/xml/timtam-2023-10-15-original.xml",
  "out/subsampling-experiment/xml/timtam-2023-10-15-subsample-0p33.xml",
  "out/subsampling-experiment/xml/timtam-2023-10-15-subsample-0p66.xml"
  )
beast_logs <- c("out/timtam-2023-10-15-original.log",
                "out/timtam-2023-10-15-subsample-0p66.log",
                "out/timtam-2023-10-15-subsample-0p33.log")
beast_trees <-
  gsub(beast_logs,
       pattern = ".log",
       replacement = "-timed-sequences.trees")

## ============================================================

present_val <- readRDS(input_present_rds)

## We will subset the samples to a managable amount by removing the
## first 10% as burn-in, and then sampling N uniformly selected points
## from the remaining samples.
##
## Here *N = 20* (see `subset_ixs')
##
num_to_burn <- 10001
subset_ixs <- floor(seq.int(from = 1, to = 90000, length.out = 500))

beast_models <-
  purrr::map(
           .x = beast_xmls,
           .f = \(x) timtamslamR::read_beast2_xml(x, present_val, "days")
         )
post_dfs <-
  furrr::future_map(
           .x = beast_logs,
           .f = \(x) timtamslamR::read_beast2_log(x, burn = num_to_burn)[subset_ixs,]
         )
tree_sets <-
  furrr::future_map(
           .x = beast_trees,
           .f = \(x) timtamslamR::read_beast2_trees(x, burn = num_to_burn)[subset_ixs]
         )
prev_dfs <-
  furrr::future_map(
           .x = seq_along(beast_xmls),
           .f = \(ix) timtamslamR::total_prevalence_log(
                                     post_dfs[[ix]],
                                     tree_sets[[ix]],
                                     beast_models[[ix]]
                                   )
         )


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

beast_xml_orig <- "out/subsampling-experiment/xml/timtam-2023-10-15-original.xml"
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

beast_xml_labels <- c(
  "out/subsampling-experiment/xml/timtam-2023-10-15-original.xml" = "Original",
  "out/subsampling-experiment/xml/timtam-2023-10-15-subsample-0p33.xml" = "Subsample 66%",
  "out/subsampling-experiment/xml/timtam-2023-10-15-subsample-0p66.xml" = "Subsample 33%"
  )

prev_post_df <-
  furrr::future_map2(
           .x = prev_dfs,
           .y = beast_xmls,
           .f = \(x, y) prev_summary_values(x, y, history_times = history_times_num)) |>
  bind_rows() |>
  mutate(xml_file = factor(xml_file, levels = beast_xmls))

## ============================================================

## This warning is fine, it is expected when there are zero counts in
## there. These will lock to the bottom of the plot, it is not a
## problem.
warning("DO NOT WORRY about a warning that scale_y_log10() produced infinite values.")

prev_gg <-
  ggplot(data = prev_post_df,
                  aes(x = history_times,
                      y = median,
                      ymin = lower,
                      ymax = upper,
                      shape = xml_file,
                      colour = xml_file)) +
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
       y = "Prevalence") +
  theme_bw() +
  theme(legend.position = "none")


if (interactive()) {
  print(prev_gg)
} else {
  plot_width <- 12
  plot_height <- 9
  ggsave(output_png, prev_gg, width = plot_width, height = plot_height, units = "cm")
  ggsave(filename = sub("png$", "svg", output_png), plot = prev_gg, width = plot_width, height = plot_height, units = "cm")
  saveRDS(prev_gg, sub("png$", "rds", output_png))
}
