library(dplyr)
library(ggpattern)
library(ggplot2)
library(lubridate)

## Read in the data from WPD and munged into a nice data.frame. By
## default WPD labels the columns of a bar-chart "BAR0", "BAR1", ...
## and the weeks are labelled 1, 2, ... so these are parsed and used
## to create a column of the first day of each week. Based on the
## timing of the vaccination rounds, it looks like week 18 covers the
## dates 4--8 May, so we have shifted the epidemiological weeks to
## match this.
##
## To account for the non-integer values exported by WPD we just round
## to the nearest integer. The total case count matches the reported
## number of cases so it seems like this has been accurate.

epi_week_start_date <- ymd("2010-01-01") - days(3)

who_df <-
  read.csv("./data/time-series.csv", header = FALSE) |>
  rename(week_zrd = V1, cases = V2) |>
  mutate(week_zrd = as.integer(gsub(x = week_zrd,
                                    pattern = "Bar",
                                    replacement = "")) + 1,
         week_date = epi_week_start_date + week_zrd * days(7),
         cases = round(cases))

## It is useful to have the dates of the vaccination rounds in a data
## frame so we can plot them later.

interval_df <-
  data.frame(
    start_date = ymd(c("2010-05-04")),
    end_date = ymd(c("2010-05-08")),
    label = c("mOPV1 Round 1")
  )

## It is useful to have a plot which displays this data along with the
## various times at which interventions where enacted.

outbreak_gg <-
  ggplot() +
  geom_col_pattern(
    data = who_df,
    mapping = aes(x = week_date,
                  y = cases),
    fill = "white", colour = "black",
    pattern_spacing = 0.015, pattern_angle = 45
  ) +
  geom_vline(xintercept = ymd("2010-05-04")) +
  scale_pattern_manual(
    values = c("stripe")
  ) +
  labs(x = NULL, y = "Confirmed cases") +
  theme_bw()
