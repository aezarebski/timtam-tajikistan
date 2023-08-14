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
    start_date = ymd(c("2010-05-04",
                       "2010-05-18",
                       "2010-06-01",
                       "2010-06-15",
                       "2010-10-04",
                       "2010-11-08",
                       "2010-09-13")),
    end_date = ymd(c("2010-05-08",
                     "2010-05-22",
                     "2010-06-05",
                     "2010-06-19",
                     "2010-10-08",
                     "2010-11-12",
                     "2010-09-17")),
    label = c("Round 1",
              "Round 2",
              "Round 3",
              "Round 4",
              "Round 5",
              "Round 6",
              "Mop-up"),
    y = c(80, 80, 80, 80, 80, 80, 80)
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
  geom_point(
    data = interval_df,
    mapping = aes(x = start_date, y = y),
    shape = 6
  ) +
  geom_text(
    data = interval_df,
    mapping = aes(x = start_date, y = y, label = label),
    hjust = 0,
    vjust = -1,
    angle = 30,
    size = 3
  ) +
  scale_pattern_manual(
    values = c("stripe")
  ) +
  scale_x_date(
    date_breaks = "2 month",
    date_labels = "%b",
    limits = c(ymd("2010-01-01"), ymd("2010-08-01")),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 100, 20),
    expand = c(0, 2)
  ) +
  labs(x = NULL, y = "Confirmed cases") +
  theme_bw()


ggsave(filename = "out/manuscript/data-plot.png",
       plot = outbreak_gg,
       height = 0.7 * 14.8, width = 21.0, # A5
       ## height = 10.5, width = 14.8, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")
