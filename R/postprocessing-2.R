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
library(xtable)
library(timtamslamR)

set.seed(1)

## Define all the files that are either used or created by this script
## see the configuration XML for these values and a short description
## of what the files contain.

config_xml <- "config.xml"
stopifnot(file.exists(config_xml))
config <- as_list(read_xml(config_xml))

output <- list(
  r_eff_png = config$files$results$figures$manuscript$posterior$r[[1]],
  p_psi_png = config$files$results$figures$manuscript$posterior$propPsi[[1]],
  p_ts_png = config$files$results$figures$manuscript$posterior$propTs[[1]],
  sigma_png = config$files$results$figures$manuscript$posterior$sigma[[1]],
  summary_tex = config$files$results$tables$manuscript$summary[[1]]
)
timtam_xml <- config$files$results$intermediate$beastXML[[1]]
stopifnot(file.exists(timtam_xml))

## Define a label for the y-axis in the prior posterior plots to keep
## things DRY.

Y_AXIS_TITLE <- "Prior/posterior density"

## Define and extract relevant information for the input files after
## checking that they exist.

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
  str_replace("xml/", "")
if (!file.exists(timtam_log)) {
  stop(sprintf("The log file %s is missing!", timtam_log))
}

## We select and rename a subset of the variables to produce a data
## frame that is easier to work with. Note that if you want to plot
## extra variables this is where the changes should start.

num_to_burn <- 200
post_samples_df <-
  timtam_log |>
  read_beast2_log(burn = num_to_burn) |>
  select(starts_with("TTR0"),
         TTPropPsi.2,
         TTPropTS.2,
         TTNetRemovalRate) |>
  rename(p_psi = TTPropPsi.2,
         p_ts = TTPropTS.2,
         sigma = TTNetRemovalRate) |>
  rename_with(.fn = \(n) str_replace(string = n,
                                     pattern = "TTR0.",
                                     replacement = "r_eff_"),
              .cols = starts_with("TTR0")) |>
  melt(id.vars = c())

num_r_eff_vars <-
  post_samples_df$variable |>
  unique() |>
  str_detect(pattern = "r_eff") |>
  sum()

is_r_eff_mask <- grepl(post_samples_df$variable, pattern = "r_eff_[0-9]+")
is_p_psi_mask <- grepl(post_samples_df$variable, pattern = "p_psi")
is_p_ts_mask <- grepl(post_samples_df$variable, pattern = "p_ts")
is_sigma_mask <- grepl(post_samples_df$variable, pattern = "sigma")

## Because we want to be able to see a table of the parameter
## estimates we will compute the summary statistics and use the xtable
## package to write this to a latex table.
##
## The strange sequence of letters used in the call to \code{xtable}
## is so that it knows the type of variables involved and how they
## should be displayed.

stop("USE HDPI FOR INTERVALS")
sink(file = output$summary_tex)
post_samples_df |>
  group_by(variable) |>
  summarise(est = median(value),
            cri_lower = quantile(value, probs = 0.025),
            cri_upper = quantile(value, probs = 0.975)) |>
  xtable(math.style.exponents = TRUE,
         digits = 2,
         display = c("s", "s", "E", "E", "E")) |>
  print(include.rownames = FALSE)
sink()

## Because we want to have nice labels for the facets of this plot, we
## need to construct a labeller object to give to the
## \code{facet_wrap} function. We put them all in a single column so
## the figures sit nicely in a two column layout.

r_eff_cols <- paste0("r_eff_", 1:num_r_eff_vars)
r_eff_labels <- c("Before 20 April",
                  "20 April - 4 May",
                  "4 May - 18 May",
                  "After 18 May")
labelling_list <- setNames(r_eff_labels, r_eff_cols)
facet_labeller <- labeller(variable = labelling_list)
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
  labs(y = Y_AXIS_TITLE,
       x = "Effective reproduction number") +
  facet_wrap(~variable, scales = "free",
             ncol = 1,
             labeller = facet_labeller) +
  theme_bw() +
  theme()

ggsave(filename = output$r_eff_png,
       plot = r_eff_gg,
       ## height = 14.8, width = 21.0, # A5
       ## height = 10.5, width = 21.0, # A6
       height = 4 * 7.4 - 5, width = 10.5, # A7
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
  labs(x = "Proportion of infections sequenced",
       y =  Y_AXIS_TITLE) +
  scale_x_continuous(
    limits = c(0, 0.01)
  ) +
  theme_bw() +
  theme()

ggsave(filename = output$p_psi_png,
       plot = p_psi_gg,
       ## height = 14.8, width = 21.0, # A5
       height = 10.5, width = 14.8, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")

p_ts_gg <-
  ggplot() +
  geom_histogram(
    data = post_samples_df[is_p_ts_mask, ],
    mapping = aes(x = value, y = after_stat(density)),
    bins = 30
  ) +
  stat_function(
    fun = function(x) dbeta(x, shape1 = 2.0, shape2 = 301.0),
    geom = "line"
  ) +
  labs(x = "Proportion of infections observed (in time series)",
       y = Y_AXIS_TITLE) +
  scale_x_continuous(
    limits = c(0, 0.03)
  ) +
  theme_bw() +
  theme()

ggsave(filename = output$p_ts_png,
       plot = p_ts_gg,
       ## height = 14.8, width = 21.0, # A5
       height = 10.5, width = 14.8, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")

sigma_gg <-
  ggplot() +
  geom_histogram(
    data = post_samples_df[is_sigma_mask, ],
    mapping = aes(x = value, y = after_stat(density)),
    breaks = seq(from = 0.1,
                 to = max(post_samples_df[is_sigma_mask, ]$value) * 1.01,
                 by = 0.002)
  ) +
  stat_function(
    fun = function(x) dunif(x, min = 0.1, max = 1.0),
    geom = "line"
  ) +
  labs(x = "Net becoming uninfectious rate",
       y = Y_AXIS_TITLE) +
  theme_bw() +
  theme()

ggsave(filename = output$sigma_png,
       plot = sigma_gg,
       ## height = 14.8, width = 21.0, # A5
       height = 10.5, width = 14.8, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")
