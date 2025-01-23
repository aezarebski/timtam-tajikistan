library(dplyr)
library(ggpattern)
library(ggplot2)
library(reshape2)
library(lubridate)
library(timtamslamR)
library(xml2)
library(jsonlite)


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
       value.name = "counts")

timeseries_gg <-
  ggplot(data = plot_df,
         mapping = aes(x = bwd_times,
                       y = counts,
                       colour = type)) +
  geom_point() +
  geom_line() +
  scale_x_reverse() +
  labs(x = "Backwards time",
       y = "Daily case count",
       colour = "Inclusion probability") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2, 0.8))

ggsave(filename = output_disaster_png, plot = timeseries_gg,
       dpi = 300, width = 148 + 20, height = 105, units = "mm")
