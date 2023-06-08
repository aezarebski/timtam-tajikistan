## [[file:scratch.org::*Retrieving and organising the data][Retrieving and organising the data:1]]
library(rentrez)
library(ape)
library(xml2)
## Retrieving and organising the data:1 ends here

## [[file:scratch.org::*Retrieving and organising the data][Retrieving and organising the data:2]]
fasta_file <- "tajikistan-poliomyelitis.fasta"
metadata_file <- "tajikistan-poliomyelitis-metadata.xml"

tajikistan_accession_numbers <-
  seq.int(from = 880365, to = 880521)
accession_numbers <- c(
  seq.int(from = 800662, to = 800683),
  seq.int(from = 812248, to = 812257),
  tajikistan_accession_numbers
)
## Retrieving and organising the data:2 ends here

## [[file:scratch.org::*Retrieving and organising the data][Retrieving and organising the data:3]]
seqs <- entrez_fetch(
  db = "nucleotide",
  id = sprintf("KC%d", accession_numbers),
  rettype = "fasta"
)

writeLines(seqs, fasta_file)
## Retrieving and organising the data:3 ends here

## [[file:scratch.org::*Retrieving and organising the data][Retrieving and organising the data:4]]
foo <- entrez_fetch(
  db = "nucleotide",
  id = sprintf("KC%d", tajikistan_accession_numbers),
  rettype = "xml",
  parsed = TRUE
)
XML::saveXML(foo, file = metadata_file)
foo <- read_xml(metadata_file)
## Retrieving and organising the data:4 ends here

## [[file:scratch.org::*Retrieving and organising the data][Retrieving and organising the data:5]]
qualifiers <- xml_find_all(
  foo, "//GBQualifier[GBQualifier_name[text()='collection_date']]"
)
collection_dates <-
  xml_text(
    xml_find_first(qualifiers, "./GBQualifier_value")
  )

primary_accessions <-
  xml_find_all(foo, "//GBSeq/GBSeq_primary-accession")
accession_texts <- xml_text(primary_accessions)
## Retrieving and organising the data:5 ends here
