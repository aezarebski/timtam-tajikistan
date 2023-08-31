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
y <- rename_dates_to_times_a(seqs)
## plot_times(y)
write_fasta(y, output_fasta)
