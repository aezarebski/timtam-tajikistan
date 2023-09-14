library(dplyr)
library(purrr)
library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(jsonlite)
library(magrittr)
library(xml2)
library(lubridate)
library(timtamslamR)

set.seed(1)

## Define a list of output files so they are all in one place.

output <- list(
  r_eff_png = "out/manuscript/parameter-r-eff.png",
  p_psi_png = "out/manuscript/parameter-p-psi.png",
  p_ts_png = "out/manuscript/parameter-p-ts.png"
)

## Define and extract relevant information for the input files after
## checking that they exist.

timtam_xml <- "xml/timtam-2023-09-13.xml"
stopifnot(file.exists(timtam_xml))
mcmc_config <- read_xml(timtam_xml)

timtam_log <-
  mcmc_config |>
  xml_find_first("//logger[@id='tracelog']") |>
  xml_attr("fileName")
## If there is a \code{$(filebase)} in the log file name, replace it
## with the the base name of the XML file specifying the analysis.
timtam_log <-
  timtam_log |>
  str_replace("\\$\\(filebase\\)",
              str_remove(timtam_xml, "\\.xml")) |>
  str_replace("xml/", "out/")
if (!file.exists(timtam_log)) {
  stop(sprintf("The log file %s is missing!", timtam_log))
}
num_to_burn <- 200
post_samples_df <- timtam_log |>
  read_beast2_log(burn = num_to_burn) |>
  select(TTR0.1, TTR0.2, TTPropPsi.2, TTPropTS.2) |>
  rename(r_eff_1 = TTR0.1, r_eff_2 = TTR0.2, p_psi = TTPropPsi.2, p_ts = TTPropTS.2) |>
  melt(id.vars = c())

is_r_eff_mask <- grepl(post_samples_df$variable, pattern = "r_eff_[0-9]+")
is_p_psi_mask <- grepl(post_samples_df$variable, pattern = "p_psi")
is_p_ts_mask <- grepl(post_samples_df$variable, pattern = "p_ts")

r_eff_gg <-
  ggplot() +
  geom_histogram(
    data = post_samples_df[is_r_eff_mask, ],
    mapping = aes(x = value, y = after_stat(density)),
    bins = 30
  ) +
  stat_function(
    fun = function(x) dnorm(x, mean = 2.0, sd = 2.0),
    geom = "line"
  ) +
  labs(title = "Reproduction number",
       y = "Posterior density") +
  facet_wrap(~variable, scales = "free") +
  theme_bw() +
  theme(axis.title.x = element_blank())

ggsave(filename = output$r_eff_png,
       plot = r_eff_gg,
       ## height = 14.8, width = 21.0, # A5
       height = 10.5, width = 21.0, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")

p_psi_gg <-
  ggplot() +
  geom_histogram(
    data = post_samples_df[is_p_psi_mask, ],
    mapping = aes(x = value, y = after_stat(density)),
    bins = 30
  ) +
  stat_function(
    fun = function(x) dbeta(x, shape1 = 2.0, shape2 = 301.0),
    geom = "line"
  ) +
  labs(title = "Proportion of infections sequenced",
       y = "Posterior density") +
  scale_x_continuous(
    limits = c(0, 0.01)
  ) +
  theme_bw() +
  theme(axis.title.x = element_blank())

ggsave(filename = output$p_psi_png,
       plot = p_psi_gg,
       ## height = 14.8, width = 21.0, # A5
       height = 10.5, width = 14.8, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")

p_ts_gg <- ggplot() +
  geom_histogram(
    data = post_samples_df[is_p_ts_mask, ],
    mapping = aes(x = value, y = after_stat(density)),
    bins = 30
  ) +
  stat_function(
    fun = function(x) dbeta(x, shape1 = 2.0, shape2 = 301.0),
    geom = "line"
  ) +
  labs(title = "Proportion of infections observed",
       subtitle = "Proportion of cases that are present in the time series and not sequenced",
       y = "Posterior density") +
  scale_x_continuous(
    limits = c(0, 0.03)
  ) +
  theme_bw() +
  theme(axis.title.x = element_blank())

ggsave(filename = output$p_ts_png,
       plot = p_ts_gg,
       ## height = 14.8, width = 21.0, # A5
       height = 10.5, width = 14.8, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")