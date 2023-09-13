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

## Add a helper function that is helpful for smoothing out the
## lubridate interface.

#' @title Convert numeric days to a lubridate period
#'
#' @param x A numeric vector of days.
#'
#' @return A lubridate period object.
#'
days_to_period <- function(x) {
  seconds(round(x * 24 * 60 * 60))
}

## Define a list of output files so they are all in one place.

output <- list(
  params_png = "out/demo-reff.png",
  prev_png = "out/demo-prevalence.png",
  combined_png = "out/combined-plot.png"
)

## Define and extract relevant information for the input files after
## checking that they exist.

timtam_xml <- "xml/timtam-2023-09-07.xml"
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
stopifnot(file.exists(timtam_log))
timtam_tree_id <-
  mcmc_config |>
  xml_find_first("//data[@id]/child::sequence/..") |>
  xml_attr("id")
timtam_tree_log <-
  mcmc_config |>
  xml_find_first("//logger[@mode='tree']") |>
  xml_attr("fileName") |>
  str_replace("\\$\\(filebase\\)",
              str_remove(timtam_xml, "\\.xml")) |>
  str_replace("\\$\\(tree\\)",
              timtam_tree_id) |>
  str_replace("xml/", "out/")
stopifnot(file.exists(timtam_tree_log))

origin_time <-
  mcmc_config |>
  xml_find_first("//parameter[@name='originTime']") |>
  xml_text() |>
  as.numeric()
r_change_times <-
  mcmc_config |>
  xml_find_first("//parameter[@name='r0ChangeTimes']") |>
  xml_text() |>
  str_split(" ") |>
  unlist() |>
  as.numeric()
num_to_burn <- 200
timtam_post <- timtam_log |>
  read_beast2_log(burn = num_to_burn)

## START STEP FUNCTION SUMMARY
##
## Summarising the step function is a bit fiddly, so feel free to skip
## this part of the code on a first reading. All you need to know is
## that we are building a dataframe that contains a summary of the
## estimated values of the basic reproduction number through time
## based on the parameters of the model.

#' Compute a posterior summary of the samples of a piece-wise constant
#' function.
#'
#' @param ts the times at which there is a step.
#' @param vss the posterior samples of the value of the function at
#'   the step.
#'
step_function_cri <- function(ts, vss) {
  .ts <- c(head(ts, 1), rep(tail(head(ts, -1), -1), each = 2), tail(ts, 1))
  n <- length(.ts)
  result <- data.frame(
    t = .ts,
    lower_bound = rep(NA, n),
    median = rep(NA, n),
    upper_bound = rep(NA, n)
  )
  for (ix in seq.int(length(vss))) {
    summ <- as.numeric(quantile(x = vss[[ix]], probs = c(0.025, 0.5, 0.975)))
    .ix <- 2 * ix - 1
    result[.ix, 2:4] <- summ
    result[.ix + 1, 2:4] <- summ
  }
  return(result)
}

buffer_after_last_seq <- 10 # days

step_fun_times <- origin_time - c(origin_time, r_change_times, - buffer_after_last_seq)
eval(parse(text = paste0(c("step_fun_samples <- list(", paste0(sprintf("timtam_post$TTR0.%d", 1:(length(r_change_times) + 1)), collapse = ","), ")"), collapse = "")))
estimate_cri <- step_function_cri(step_fun_times, step_fun_samples)

## END STEP FUNCTION SUMMARY

gg_r_eff <-
  ggplot() +
  geom_ribbon(
    data = estimate_cri,
    mapping = aes(x = t, ymin = lower_bound, ymax = upper_bound),
    alpha = 0.5
  ) +
  geom_line(
    data = estimate_cri,
    mapping = aes(x = t, y = median)
  ) +
  geom_vline(
    xintercept = origin_time,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 1.0,
    linetype = "dotted"
  ) +
  labs(
    x = "Days since outbreak began",
    y = "Reproduction number"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13)
  )

ggsave(
  filename = output$params_png,
  plot = gg_r_eff,
  height = 10.5, width = 14.8,
  units = "cm"
)

my_present <- readRDS("out/present.Rds")
my_units <- "days"


my_beast2_log <- read_beast2_log(timtam_log, burn = num_to_burn)
my_beast2_trees <- read_beast2_trees(timtam_tree_log, burn = num_to_burn)
my_beast2_model <- read_beast2_xml(timtam_xml, my_present, my_units)

## Prevalence calculation can take a while so we will subsample the
## posterior samples first to keep things managable.

log_len <- nrow(my_beast2_log)
subsample_ixs <-
  round(seq(from = 1, to = nrow(my_beast2_log), length.out = 200))
prev_helper <- total_prevalence_log(
  my_beast2_log[subsample_ixs,],
  my_beast2_trees[subsample_ixs],
  my_beast2_model
)
prev_times  <- my_beast2_model$hist_times

## We can make a diagnostic plot of the prevalence through time
## estimates to double check that the MCMC has converged.

trace_plot <-
  ggplot() +
  geom_line(
    data = melt(prev_helper,
                id.vars = "Sample",
                value.name = "Prevalence",
                variable.name = "Time"),
    mapping = aes(x = Sample, y = Prevalence, colour = Time)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

## The following is just some boring munging to get the prevalence
## estimates into a sensible format for plotting.

bkwd_hist_times <- as.numeric(my_beast2_model$present$date_time - my_beast2_model$hist_times) / (60 * 60 * 24)

time_df <- data.frame(
  point = sprintf("size%d", seq_along(prev_times)),
  time = origin_time - bkwd_hist_times
)

summary_vec <- function(prev_samples) {
  as.numeric(quantile(prev_samples, probs = c(0.025, 0.5, 0.975)))
}

## Use some meta-programming to get the prevalence estimates into a
## suitable data frame because I am lazy.

prev_df <- data.frame(type = c("lower", "mid", "upper"))
for (ix in seq_along(prev_times)) {
  eval(parse(text = paste0(c("prev_df$size", ix, " <- summary_vec(prev_helper$TTPrevalence.", ix, ")"), collapse = "")))
}

prev_df <-
  prev_df |>
  melt(id.vars = "type",
       variable.name = "point",
       value.name = "size") |>
  mutate(point = as.character(point)) |>
  full_join(time_df, by = "point") |>
  select(time, type, size) |>
  dcast(time ~ type, value.var = c("size"))

## Make a nice plot that demonstrates all of the prevalence results.

prev_fig <-
  ggplot() +
  geom_linerange(
    data = prev_df,
    mapping = aes(x = time, ymin = lower, ymax = upper),
    linewidth = 3,
    alpha = 0.5
  ) +
  geom_point(
    data = prev_df,
    mapping = aes(x = time, y = mid),
    size = 2
  ) +
  geom_vline(
    xintercept = origin_time,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    limits = c(0, origin_time + buffer_after_last_seq)
  ) +
  labs(
    x = "Days after introduction",
    y = "Prevalence of infection"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13)
  )

ggsave(
  filename = output$prev_png,
  plot = prev_fig,
  height = 10.5, width = 14.8,
  units = "cm"
)

example_plot <- plot_grid((gg_r_eff +
                           theme(axis.text.x = element_blank(),
                                 axis.title.x = element_blank())),
                          prev_fig,
                          ncol = 1)

ggsave(filename = output$combined_png,
       plot = example_plot,
       width = 21.0, height = 14.8,
       ## height = 14.8, width = 21.0, # A5
       ## height = 10.5, width = 14.8, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")
