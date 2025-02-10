library(dplyr)
library(ggpattern)
library(ggplot2)
library(reshape2)
library(lubridate)
library(timtamslamR)
library(xml2)
library(jsonlite)


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


config <- as_list(read_xml("config.xml"))

input_disaster_json <- config$files$results$subsampling$intermediate$disasterStringsJSON[[1]]
stopifnot(file.exists(input_disaster_json))
output_disaster_png <- config$files$results$subsampling$figures$disasterPlot[[1]]

disaster_data <- jsonlite::fromJSON(input_disaster_json)

plot_df <- data.frame(
  bwd_times = disaster_data$times$bwd_times,
  included_percent_100 = disaster_data$original$counts,
  included_percent_33 = disaster_data$subsamples$counts[[1]],
  included_percent_66 = disaster_data$subsamples$counts[[2]]
) |>
  melt(id.vars = "bwd_times",
       variable.name = "type",
       value.name = "counts") |>
  mutate(type = factor(type, levels = c("included_percent_100",
                                        "included_percent_66",
                                        "included_percent_33")))

log_labels <- c("included_percent_100" = "Original",
                "included_percent_66" = "Subsample 66%",
                "included_percent_33" = "Subsample 33%")

timeseries_gg <-
  ggplot(data = plot_df,
         mapping = aes(x = bwd_times,
                       y = counts,
                       shape = type,
                       colour = type)) +
  geom_point() +
  geom_line() +
  scale_x_reverse() +
  scale_colour_manual(name = "Time series",
                      labels = log_labels,
                      values = scale_colour_vals) +
  scale_shape_manual(name = "Time series",
                     labels = log_labels,
                     values = scale_shape_vals) +
  labs(x = "Backwards time",
       y = "Daily case count") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.background = legend_background_style)

if (interactive()) {
  print(timeseries_gg)
} else {
  ggsave(filename = output_disaster_png, plot = timeseries_gg,
         dpi = 300, width = 148 + 20, height = 105, units = "mm")
  ggsave(filename = sub("png$", "svg", output_disaster_png),
         plot = timeseries_gg,
         width = 148 + 20, height = 105, units = "mm")
  ## Write the gg plot to a .rds file to make it easier to revive.
  saveRDS(timeseries_gg, sub("png$", "rds", output_disaster_png))
}
