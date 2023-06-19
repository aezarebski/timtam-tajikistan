## [[file:scratch.org::*Packages][Packages:1]]
library(rentrez)
library(ape)
library(xml2)
library(purrr)
library(stringr)
library(lubridate)

fetch_as_xml <- function(...) {
  temp_file <- tempfile()
  on.exit(file.remove(temp_file))

  entrez_obj <- entrez_fetch(
    ...,
    rettype = "xml",
    parsed = TRUE
  )
  XML::saveXML(entrez_obj, file = temp_file)
  return(xml2::read_xml(temp_file))
}
## Packages:1 ends here

## [[file:scratch.org::*Define files and accession numbers][Define files and accession numbers:1]]
fasta_file <- "tajikistan-poliomyelitis.fasta"
entrez_xml <- "tajikistan-poliomyelitis.xml"

tajikistan_accession_numbers <-
  seq.int(from = 880365, to = 880521)
accession_numbers <- c(
  seq.int(from = 800662, to = 800683),
  seq.int(from = 812248, to = 812257),
  tajikistan_accession_numbers
)
## Define files and accession numbers:1 ends here

## [[file:scratch.org::*Download the data and extract relevant information][Download the data and extract relevant information:1]]
seqs_xml <- fetch_as_xml(
  db = "nucleotide",
  id = sprintf("KC%d", accession_numbers)
)
write_xml(x = seqs_xml, file = entrez_xml)
## Download the data and extract relevant information:1 ends here

## [[file:scratch.org::*Download the data and extract relevant information][Download the data and extract relevant information:2]]
qualifiers <- xml_find_all(
  seqs_xml, "//GBQualifier[GBQualifier_name[text()='collection_date']]"
)

collection_dates <-
  xml_find_first(qualifiers, "./GBQualifier_value") |>
  xml_text()

accession_texts <-
  xml_find_all(seqs_xml, "//GBSeq/GBSeq_primary-accession") |>
  xml_text()

sequences <-
  xml_find_all(seqs_xml, "//GBSeq_sequence") |>
  xml_text()
## Download the data and extract relevant information:2 ends here

## [[file:scratch.org::*Format as FASTA][Format as FASTA:1]]
seqs_dnabin <-
  sequences |>
  map(str_split, "") |>
  map(unlist) |>
  as.DNAbin()

date_mask <- str_detect(collection_dates, "[0-9]+-[A-Za-z]+-[0-9]+")

seqs_dnabin <- seqs_dnabin[date_mask]
new_names <- format(
  dmy(collection_dates[date_mask]),
  "%d-%m-%Y"
)
names(seqs_dnabin) <- str_c(accession_texts[date_mask], new_names, sep = "_")

write.dna(seqs_dnabin, file = fasta_file, format = "fasta")
## Format as FASTA:1 ends here

## [[file:scratch.org::*Checking for temporal signal in the sequences][Checking for temporal signal in the sequences:1]]
alignment_dna <- read.dna("test-alignment.fasta", format = "fasta")
dist_matrix <- dist.dna(alignment_dna, model = "K80", pairwise.deletion = TRUE)
nj_phylo <- nj(dist_matrix)
write.nexus(nj_phylo, file = "test-nj-tree.nexus")
## Checking for temporal signal in the sequences:1 ends here

## [[file:scratch.org::*Check the results from TempEst][Check the results from TempEst:1]]
tempest_df <- read.table("tempest-data.txt", header = TRUE, sep = "\t")

tempest_gg <-
  ggplot(data = tempest_df,
	 mapping = aes(x = date, y = distance)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_continuous(limits = c(2009.7,2011)) +
  scale_y_continuous(limits = c(0, 0.015)) +
  coord_cartesian(ylim = c(0, 0.015)) +
  theme_bw()
print(tempest_gg)
## Check the results from TempEst:1 ends here
