library(ape)
library(dplyr)
library(ggpattern)
library(ggplot2)
library(reshape2)
library(stringr)
library(lubridate)
library(timtamslamR)

output_fasta <- "out/timed-sequences.fasta"
if (file.exists(output_fasta)) {
  stop("Output file already exists. Please delete it and try again.")
}

li_nexus <- "data/li-alignment.nexus"
stopifnot(file.exists(li_nexus))
seqs <- read.nexus.data(file = li_nexus) |> as.DNAbin()

## Filter for only the sequences from Tajiikistan.

seqs <- seqs[grep("^TJK", names(seqs))]

## Uniformly distribute the sequences across the day on which they are
## are dated so that we can model them as a point process rather than
## a time series. Write the results to a new file which we can use for
## BEAUti. There are plotting functions commented out below which can
## be used to check this has produced sensible output.

## plot_dates(seqs)
timed_seqs <- rename_dates_to_times_a(seqs)
## plot_times(timed_seqs)
write_fasta(timed_seqs, output_fasta)

p <- get_present(seqs, timed_seqs)

## Align the time series to the "present" (i.e. the most recent
## sequence). This is a little bit fiddly because the last week
## included doesn't have any cases but is probably important to
## include as a zero.

z1 <-
  read.csv("out/who_df.csv") |>
  mutate(count = cases,
         week_start = as.Date(week_start_date),
         week_end = as.Date(week_end_date)) |>
  select(count, week_start, week_end) |>
  spread_across_days()
z2 <- rename_time_series(p, z1)
write.csv(z2, "out/who_df_timestamped.csv",
          row.names=FALSE)

## Finally, it is most helpful to have the disaster data formatted as
## a string that we can easily copy and paste into the XML file.

sink("out/disaster-strings.txt")
print("Here are the disaster sizes:\n")
paste(z2$count, sep = "", collapse = " ")
print("Here are the backward-times of the disasters:\n")
paste(z2$bwd_times, sep = "", collapse = " ")
sink()
