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


config <- as_list(read_xml("config.xml"))

## Define all the files that are either used or created by this script
## see the configuration XML for these values and a short description
## of what the files contain.

input_present_rds <- config$files$results$intermediate$present[[1]]
timtam_xml <- config$files$results$intermediate$beastXML[[1]]
output <- list(
  params_png = config$files$results$figures$posteriorR[[1]],
  prev_png = config$files$results$figures$posteriorPrev[[1]],
  combined_png = config$files$results$figures$manuscript$combinedParameters[[1]],
  combined_2_png = config$files$results$figures$manuscript$combinedEverything[[1]]
)
data_plot_rds <- config$files$results$figures$dataPlotRds[[1]]

stopifnot(file.exists(data_plot_rds))
stopifnot(file.exists(input_present_rds))
stopifnot(file.exists(timtam_xml))

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
  str_replace("xml/", "")
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

my_present <- readRDS(input_present_rds)
my_units <- "days"

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
estimate_cri <-
  step_function_cri(step_fun_times, step_fun_samples) |>
  mutate(date = my_present$date + (t - origin_time))

## END STEP FUNCTION SUMMARY

gg_r_eff <-
  ggplot() +
  geom_ribbon(
    data = estimate_cri,
    mapping = aes(x = date, ymin = lower_bound, ymax = upper_bound),
    alpha = 0.5
  ) +
  geom_line(
    data = estimate_cri,
    mapping = aes(x = date, y = median)
  ) +
  geom_hline(
    yintercept = 1.0,
    linetype = "dotted"
  ) +
  scale_x_date(
    breaks = c(ymd("2010-01-01"),
               ymd("2010-03-01"),
               ymd("2010-05-01"),
               ymd("2010-07-01")),
    expand = c(0, 20),
    limits = range(estimate_cri$date),
    date_labels = "%b %d",
    name = "Days since outbreak began"
  ) +
  labs(
    y = "Reproduction number"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.title.x = element_blank()
  )

ggsave(
  filename = output$params_png,
  plot = gg_r_eff,
  height = 10.5, width = 14.8,
  units = "cm"
)

## In addition to the figure showing the effective reproduction number
## through time, it would be nice to have an analogous figure that
## shows these values in comparison to the previous estimates from Li
## /et al/.

library(jsonlite)

li_r0_estimates <- read_json("data/li-estimates.json", simplifyVector = TRUE)

tjk_weighted_value <- function(v_0to5, v_over5) {
  tjk_0to5_prop <- li_r0_estimates$tajikistanPop2013$age0to5 / li_r0_estimates$tajikistanPop2013$total
  v_0to5 * tjk_0to5_prop + v_over5 * (1 - tjk_0to5_prop)
}

li_r0_est_df <-
  data.frame(
  age_group = c("adult", "child", "weighted_average"),
  point_est = c(li_r0_estimates$r0$adult$pointEst,
                li_r0_estimates$r0$child$pointEst,
                tjk_weighted_value(li_r0_estimates$r0$child$pointEst,
                                   li_r0_estimates$r0$adult$pointEst)),
  cri_upper = c(li_r0_estimates$r0$adult$credInt[1],
                li_r0_estimates$r0$child$credInt[1],
                tjk_weighted_value(li_r0_estimates$r0$child$credInt[1],
                                   li_r0_estimates$r0$adult$credInt[1])),
  cri_lower = c(li_r0_estimates$r0$adult$credInt[2],
                li_r0_estimates$r0$child$credInt[2],
                tjk_weighted_value(li_r0_estimates$r0$child$credInt[2],
                                   li_r0_estimates$r0$adult$credInt[2]))
  )

gg_r_eff_comparison <-
  ggplot() +
  geom_ribbon(
    data = estimate_cri,
    mapping = aes(x = date, ymin = lower_bound, ymax = upper_bound),
    alpha = 0.5
  ) +
  geom_line(
    data = estimate_cri,
    mapping = aes(x = date, y = median)
  ) +
  geom_hline(
    yintercept = 1.0,
    linetype = "dotted"
  ) +
  geom_errorbar(
    data = li_r0_est_df,
    mapping = aes(x = ymd("2009-09-15"), ymin = cri_lower, ymax = cri_upper, colour = age_group),
    width = 0.1
  ) +
  geom_point(
    data = li_r0_est_df,
    mapping = aes(x = ymd("2009-09-15"), y = point_est, colour = age_group)
  ) +
  scale_x_date(
    breaks = c(ymd("2010-01-01"),
               ymd("2010-03-01"),
               ymd("2010-05-01"),
               ymd("2010-07-01")),
    expand = c(0, 20),
    limits = range(estimate_cri$date) + c(-40, 0),
    date_labels = "%b %d",
    name = "Days since outbreak began"
  ) +
  scale_colour_manual(
    values = c("adult" = "#7fc97f",
               "child" = "#beaed4",
               "weighted_average" = "#fdc086"),
    labels = c("Adult (over 5)", "Child (5 and under)", "Weighted average")
  ) +
  labs(
    y = "Reproduction number",
    colour = "Population"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.8)
  )

ggsave(
  filename = output$params_comparison_png,
  plot = gg_r_eff_comparison,
  height = 10.5, width = 14.8,
  units = "cm"
)

## Prevalence is a parameter we are particularly interested in, so we
## should definitely make a nice plot to display out estimates of
## this.

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

#' The following is just some boring munging to get the prevalence
#' estimates into a sensible format for plotting. We use the
#' \code{difftime} function to create the time difference explicitly
#' so that we can specify the units. This way it prevents R from
#' trying to be clever and potentially measuring the difference in the
#' wrong units.

bkwd_difftime <-
  difftime(my_beast2_model$present$date_time,
           my_beast2_model$hist_times,
           units = "secs")
bkwd_hist_times <- as.numeric(bkwd_difftime) / (60 * 60 * 24)

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
  dcast(time ~ type, value.var = c("size")) |>
  mutate(date = my_present$date + (time - origin_time))

## Make a nice plot that demonstrates all of the prevalence results.

prev_fig <-
  ggplot() +
  geom_linerange(
    data = prev_df,
    mapping = aes(x = date, ymin = lower, ymax = upper),
    linewidth = 3,
    alpha = 0.5
  ) +
  geom_point(
    data = prev_df,
    mapping = aes(x = date, y = mid),
    size = 2
  ) +
  scale_x_date(
    breaks = c(ymd("2010-01-01"),
               ymd("2010-03-01"),
               ymd("2010-05-01"),
               ymd("2010-07-01")),
    date_labels = "%b %d",
    expand = c(0, 20),
    limits = range(estimate_cri$date)
  ) +
  labs(
    x = NULL,
    y = "Prevalence of infection"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 13),
    axis.title.x = element_blank()
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
                          align = "v", axis = "l",
                          ncol = 1)

ggsave(filename = output$combined_png,
       plot = example_plot,
       width = 21.0, height = 14.8,
       units = "cm")

data_gg <- readRDS(data_plot_rds) +
  ## my_present$date - origin_time
  geom_point(
    data = data.frame(x = my_present$date - origin_time,
                      y = 2),
    mapping = aes(x = x, y = y),
    shape = 6
  ) +
  geom_text(
    data = data.frame(x = my_present$date - origin_time,
                      y = 2,
                      label = "Origin (15 Oct 2009)"),
    mapping = aes(x = x, y = y, label = label),
    hjust = 0,
    vjust = -1,
    angle = 30,
    size = 3
  ) +
  scale_x_date(
    breaks = c(ymd("2010-01-01"),
               ymd("2010-03-01"),
               ymd("2010-05-01"),
               ymd("2010-07-01")),
    date_labels = "%b %d",
    expand = c(0, 20),
    limits = range(estimate_cri$date)
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 100, 20),
    expand = c(0, 2)
  ) +
  theme(legend.position = c(0.3, 0.6),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

theme_tweak <-
  theme(
    axis.title.y = element_text(size = 11)
  )

example_plot_2 <-
  plot_grid(data_gg,
            prev_fig + theme(axis.text.x = element_blank()) + theme_tweak,
            gg_r_eff + theme_tweak,
            align = "v", axis = "l",
            rel_heights = c(1, 0.7, 0.7),
            ncol = 1,
            labels = c("A", "B", "C"),
            hjust = -0.1, vjust = 1.1)

## hjust more negative moves it to the right
## vjust more positive moves it down

ggsave(filename = output$combined_2_png,
       plot = example_plot_2,
       height = 21.0, width = 21.0,
       units = "cm")
