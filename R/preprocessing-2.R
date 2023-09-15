library(ape)
library(dplyr)
library(ggpattern)
library(ggplot2)
library(reshape2)
library(stringr)
library(lubridate)
library(timtamslamR)
library(xml2)

config <- as_list(read_xml("config.xml"))

## Define all the files that are either used or created by this script
## see the configuration XML for these values and a short description
## of what the files contain.

input_cases_csv <- config$files$results$intermediate$timeSeries[[1]]
li_nexus <- config$files$data$sequences[[1]]
output_fasta <- config$files$results$intermediate$sequences[[1]]
output_present_rds <- config$files$results$intermediate$present[[1]]
output_disaster_txt <- config$files$results$intermediate$disasterStrings[[1]]
session_rdata_out <-
  sprintf("out/preprocessing-2-workspace-%s.RData", Sys.Date())

stopifnot(file.exists(input_cases_csv))
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
saveRDS(p, output_present_rds)

## Align the time series to the "present" (i.e. the most recent
## sequence). This is a little bit fiddly because the last week
## included doesn't have any cases but is probably important to
## include as a zero.

z1 <- input_cases_csv |>
  read.csv() |>
  filter(data_type == "Cases") |>
  mutate(week_start = as.Date(week_start_date),
         week_end = as.Date(week_end_date)) |>
  select(count, week_start, week_end) |>
  spread_across_days()

z2 <- rename_time_series(p, z1)

## Finally, it is most helpful to have the disaster data formatted as
## a string that we can easily copy and paste into the XML file. We
## remove the date prior to the first of February because had there
## been surveillence at that point it probably would have captured one
## of the cases as a sequence.

z3 <- z2[z2$date >= ymd("2010-02-01"),]
sink(output_disaster_txt)
print("Here are the disaster sizes:\n")
paste(z3$count, sep = "", collapse = " ")
print("Here are the backward-times of the disasters:\n")
paste(z3$bwd_times, sep = "", collapse = " ")
sink()

save.image(file = session_rdata_out)
