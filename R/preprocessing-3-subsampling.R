library(ape)
library(dplyr)
library(ggpattern)
library(ggplot2)
library(reshape2)
library(stringr)
library(lubridate)
library(timtamslamR)
library(xml2)
library(jsonlite)


config <- as_list(read_xml("config.xml"))

## Define all the files that are either used or created by this script
## see the configuration XML for these values and a short description
## of what the files contain.

input_cases_csv <- config$files$results$intermediate$timeSeries[[1]]
stopifnot(file.exists(input_cases_csv))
li_nexus <- config$files$data$sequences[[1]]
stopifnot(file.exists(li_nexus))

output_present_rds <- config$files$results$subsampling$intermediate$present[[1]]
stopifnot(file.exists(dirname(output_present_rds)))
output_disaster_txt <- config$files$results$subsampling$intermediate$disasterStrings[[1]]
stopifnot(file.exists(dirname(output_disaster_txt)))
output_disaster_json <- config$files$results$subsampling$intermediate$disasterStringsJSON[[1]]

## ===================================================================
## Since the subsampling is stochastic, it is useful to set the seed
## so it is reproducible. We also need to pull out the various
## subsampling probabilities form the configuration. NOTE that these
## probabilities reflext the probability that an observed case will be
## included in the sub-sample of the data considered. So full data is
## probability 1.0, removing half is 0.5 and removing all of it is
## 0.0.

prng_seed <-
  config$files$results$subsampling |>
  attr("prngSeed") |>
  as.integer()
set.seed(prng_seed)

ss_inclusion_probs <-
  config$files$results$subsampling |>
  attr("inclusionProbs") |>
  strsplit(" ") |>
  unlist() |>
  as.numeric()
## ===================================================================

seqs <- read.nexus.data(file = li_nexus) |> as.DNAbin()

## Filter for only the sequences from Tajiikistan.

seqs <- seqs[grep("^TJK", names(seqs))]

## Uniformly distribute the sequences across the day on which they are
## are dated so that we can model them as a point process rather than
## a time series. Write the results to a new file which we can use for
## BEAUti. There are plotting functions commented out below which can
## be used to check this has produced sensible output.

#### ## plot_dates(seqs)
timed_seqs <- rename_dates_to_times_a(seqs)
#### ## plot_times(timed_seqs)
#### write_fasta(timed_seqs, output_fasta)

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
##
## We also include a final line to be printed with times for the
## change in R0 which are a bit easier to interpret relative to the
## date of the start of vaccination.
##
## The initial values for the history sizes are also printed. They may
## look a little odd at first, but this is just a Gaussian
## distribution with a peak at 60 days. Nothing too magical, it just
## initialises them at some reasonable values.

space_sep_str <- \(v) paste(v, sep = "", collapse = " ")

z3 <- z2[z2$date >= ymd("2010-02-01"),]

sink(output_disaster_txt)
print("Here are the disaster sizes:\n")
disaster_str <- space_sep_str(z3$count)
print(disaster_str)
print("Here are the backward-times of the disasters:\n")
space_sep_str(z3$bwd_times)

bwd_hist_times <- 21 * ( 12:0 ) - 1 / 24
init_hist_sizes <- ceiling(20000 * exp(-0.001 * (bwd_hist_times - 60)^2) + 1)
## Here is a little plot that shows what these values look like in real terms.
## plot(bwd_hist_times, init_hist_sizes,
##      xlab = "Backwards time (days)",
##      ylab = "History size",
##      main = "Rough initial values for history size parameter",
##      xlim = rev(range(bwd_hist_times)))
print("Here are the backward-times of the history size estimates:\n")
space_sep_str(bwd_hist_times)
print("Here are some initial values to use for the history sizes:\n")
space_sep_str(init_hist_sizes)

print("Here are the backward-times of the parameter change times:\n")
space_sep_str(c(61.0, 153.0) + 1 / 3)
print("In particular, here are the values for the R0 change times:\n")
space_sep_str(61.0 + 1/3 + c(14, 0, -14))
sink()

## ===================================================================
## We need the subsampled time series for the MCMC XML. The following
## generates that subsample and stores it all in a JSON file that is
## convenient for programmatic use down the line.

result <-
  list(original = list(inclusion_prob = 1.0,
                       str = disaster_str,
                       counts = z3$count),
       times = list(str = space_sep_str(z3$bwd_times),
                    bwd_times = z3$bwd_times),
       subsamples = list())
for (ss_p_ix in seq_along(ss_inclusion_probs)) {
  ss_p <- ss_inclusion_probs[[ss_p_ix]]
  new_counts <- rbinom(n = length(z3$count),
                       size = z3$count,
                       prob = ss_p)
  result$subsamples[[ss_p_ix]] <- list(
    inclusion_prob = ss_p,
    str = space_sep_str(new_counts),
    counts = new_counts
  )
}

jsonlite::write_json(result,
                     output_disaster_json,
                     pretty = TRUE,
                     auto_unbox = TRUE)
## ===================================================================
