library(ape)
library(dplyr)
library(ggpattern)
library(ggplot2)
library(reshape2)
library(stringr)
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

## 4 May is a Tuesday, \code{weekdays(ymd("2010-05-04")) == "Tuesday"}
## so we can use this to add columns for the Sunday and Saturday of
## each week to the data frame above. The weekend columns are useful
## for plotting the sequence dates later.

who_df$week_start_date <- who_df$week_date - days(2)
who_df$week_end_date <- who_df$week_date + days(4)

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

## And if we want to include the sequenced samples as a bar with
## different cross hatching, we will need to extract the weekly
## sequence counts from the data as well.

li_nexus <- "data/li-alignment.nexus"
stopifnot(file.exists(li_nexus))
seqs <- ape::read.nexus.data(file = li_nexus)
seq_dates <-
  seqs |>
  names() |>
  str_extract("[0-9]{4}-[0-9]{2}-[0-9]{2}") |>
  purrr::discard(.p = is.na) |>
  ymd()

seq_df <-
  who_df |>
  select(week_date, week_start_date, week_end_date) |>
  mutate(seq_count = 0)
for (ix in seq_len(nrow(seq_df))) {
  week_start <- seq_df[ix, "week_start_date"]
  week_end <- seq_df[ix, "week_end_date"]
  seq_df[ix, "seq_count"] <-
    sum(seq_dates >= week_start & seq_dates <= week_end)
}

## It appears that the sequences don't match up nicely with the number
## of cases. Quite possibly there is a difference between the date
## associated with the case and the date associated with the sequence.
## To avoid double counting cases, I have used the cases minus the
## number of sequences associated with that week, and for the two
## weeks in which this is negative (a value of -1) I have set the
## value to 0. This should only make a minor change to the data and
## should have a far smaller effect than over counting the cases.

plt_df <-
  who_df |>
  inner_join(seq_df, by = c("week_date")) |>
  mutate(cases_minus_seqs = pmax(0, cases - seq_count)) |>
  select(week_date, seq_count, cases_minus_seqs) |>
  melt(id.vars = c("week_date"),
       variable.name = "data_type",
       value.name = "count") |>
  mutate(data_type = factor(data_type,
                            levels = c("cases_minus_seqs", "seq_count"),
                            labels = c("Cases", "Sequences")))

## It is useful to have a plot which displays this data along with the
## various times at which interventions where enacted.

outbreak_gg <-
  ggplot() +
  geom_col_pattern(
    data = plt_df,
    mapping = aes(x = week_date,
                  y = count,
                  pattern = data_type),
    fill = "white",
    pattern_spacing = 0.015, pattern_angle = 45
  ) +
  scale_pattern_manual(
    values = c("stripe", "crosshatch")
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
  scale_x_date(
    date_breaks = "2 month",
    date_labels = "%b %d",
    limits = c(ymd("2010-01-01"), ymd("2010-08-01")),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 100, 20),
    expand = c(0, 2)
  ) +
  labs(x = NULL, y = NULL, pattern = NULL, pattern_angle = NULL) +
  theme_bw() +
  theme(
    legend.position = c(0.2, 0.3)
  )

ggsave(filename = "out/manuscript/data-plot.png",
       plot = outbreak_gg,
       height = 0.7 * 14.8, width = 21.0, # A5
       ## height = 10.5, width = 14.8, # A6
       ## height = 7.4, width = 10.5, # A7
       units = "cm")
