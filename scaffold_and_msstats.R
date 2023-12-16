library(data.table)
library(MSstats)
library(tidyverse)
library(MSstatsConvert)

peptides <- data.table::fread(
  "data/Vika_20231116_Sh690-Sh695_SYNCRIP-FLAG-IP-CP_PeptideReport.xls",
  check.names = TRUE,
  skip = 73,
  fill = TRUE,
  sep = "\t"
)

proteins <- data.table::fread(
  "data/Vika_20231116_Sh690-Sh695_SYNCRIP-FLAG-IP-CP_ProteinReport.xls",
  check.names = TRUE,
  skip = 73,
  fill = TRUE,
  sep="\t"
)

peptides <- peptides %>%
  mutate(Total.TIC = as.numeric(Total.TIC))

# test_df <- peptides %>% 
#   dplyr::filter(
#     is.na(Total.TIC)
#   )

test_df <- peptides %>% 
  dplyr::mutate(
    Total.TIC.NA = case_when(is.na(Total.TIC) ~ TRUE, TRUE ~ FALSE),
    Total.TIC = case_when(is.na(Total.TIC) ~ min(Total.TIC, na.rm = TRUE), TRUE ~ Total.TIC)
  )

test_df %>% 
  ggplot(aes(x = log2(Total.TIC), fill = Total.TIC.NA)) +
  geom_histogram()

peptides %>% 
  dplyr::mutate(
    Total.TIC.NA = case_when(is.na(Total.TIC) ~ TRUE, TRUE ~ FALSE),
    Total.TIC = case_when(is.na(Total.TIC) ~ 1)
  ) %>% 
  ggplot(aes(x = log2(Total.TIC), fill = Total.TIC.NA)) +
  geom_histogram()

annotation <- data.table::fread(
  "data/annotation.csv",
  check.names = TRUE
)

split_string <- function(input_string) {
  # Split the input string by comma and remove extra spaces
  split_values <- strsplit(input_string, ",\\s*")[[1]]
  return(split_values)
}

req_col <- split_string(
  "Protein.accession.numbers, Peptide.sequence, MS.MS.sample.name, Total.TIC, Assigned, Other.Proteins"
)

pep_ms <- peptides %>% 
  dplyr::select(
    all_of(req_col)
  )

pep_ms_test <- pep_ms %>% 
  dplyr::group_by(Peptide.sequence, MS.MS.sample.name) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n > 1)

names(pep_ms) <- c("ProteinName", "PeptideSequence", "Run", "Intensity")

pep_ms <- pep_ms %>% 
  dplyr::mutate(
    PrecursorCharge = NA,
    FragmentIon = NA,
    ProductCharge = NA
  )

prot_ms <- proteins %>% 
  dplyr::select(
    Protein.accession.numbers, MS.MS.sample.name
  )

names(prot_ms) <- c("ProteinName", "Run")

prot_ms_cleaned <- MSstatsConvert::MSstatsClean(prot_ms, protein_id_col = "ProteinName")

input <- MSstats::MaxQtoMSstatsFormat(
  evidence = pep_ms,
  annotation = annotation,
  proteinGroups = prot_ms
)

peptides %>% 
  ggplot(aes(x = log2(Total.TIC))) +
  geom_histogram()

input <- merge(
  pep_ms,
  annotation,
  by = "Run"
) %>% 
  dplyr::select(!c("Raw.file"))

missed <- input %>% 
  dplyr::filter(
    grepl("_|-|>|<", PeptideSequence)
  )

quant_msstats <- dataProcess(
  raw = input,
  normalization = 'equalizeMedians',
  summaryMethod = 'TMP',
  censoredInt = "NA", 
  # cutoffCensored = "minFeature",
  MBimpute = TRUE,
  maxQuantileforCensored = 0.999,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = "logs/MSstats.log"
)
