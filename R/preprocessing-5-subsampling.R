library(dplyr)
library(stringr)
library(lubridate)
library(timtamslamR)
library(xml2)
library(jsonlite)


config <- as_list(read_xml("config.xml"))

input_disaster_json <- config$files$results$subsampling$intermediate$disasterStringsJSON[[1]]
stopifnot(file.exists(input_disaster_json))
disaster_data <- jsonlite::fromJSON(input_disaster_json)
mcmc_xml <- config$files$results$intermediate$beastXML[[1]]
out_dir <- config$files$results$subsampling$intermediate$newXMLDir[[1]]
stopifnot(file.exists(out_dir))

mcmc_xml_filepath <- function(inclusion_prob, out_dir, xml_basename) {
  ip_str <-
    inclusion_prob |>
    as.character() |>
    stringr::str_replace(pattern = "\\.", replacement = "p")
  tmp_basename <-
    stringr::str_replace(xml_basename,
                         pattern = "\\.xml",
                         replacement = stringr::str_interp(string = "-subsample-${ip_str}.xml"))
  file.path(out_dir, tmp_basename)
}

## Create a backup for the MCMC XML in the new output directory to use
## as a point of comparison. This way all the XML will be formatted in
## the same way making it easier to file diff them.

original_mcmc <- xml2::read_xml(mcmc_xml)
new_orig_mcmc_path <-
  stringr::str_replace(
             basename(mcmc_xml),
             pattern = "\\.xml",
             replacement = "-original.xml") |>
  (\(x) file.path(out_dir, x))()

xml2::write_xml(
  x = original_mcmc,
  file = new_orig_mcmc_path
)

## As a sanity check, we can assert that the disaster time series
## found in the XML should be the same as the original values.

original_str_a <- disaster_data$original$str
original_str_b <- xml2::xml_text(xml2::xml_find_first(original_mcmc, xpath = "//disasterSizes"))
stopifnot(original_str_a == original_str_b)

## Iterate over the inclusion probabilities and disaster time series
## strings to create new XML files with the appropriate values.

ix <- 1
for (ix in seq_along(disaster_data$subsamples$inclusion_prob)) {
  new_mcmc <- xml2::as_xml_document(as.character(original_mcmc))
  ip_val <- disaster_data$subsamples$inclusion_prob[ix]
  new_d_str <- disaster_data$subsamples$str[ix]
  new_d_node <- xml2::xml_find_first(new_mcmc, xpath = "//disasterSizes")
  xml2::xml_text(new_d_node) <- new_d_str
  new_xml_path <- mcmc_xml_filepath(ip_val, out_dir, basename(mcmc_xml))
  xml2::write_xml(new_mcmc, new_xml_path)
}
