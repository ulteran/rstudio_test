library(data.table)
library(tidyverse)

mq <- fread('data/Reports_for_Vika_20230815_MaxQuant.xlsx')

df <- fread('data/ResultMulti_df.csv')

df <- left_join(
  df,
  select(mq, Protein.IDs, Gene.simbol),
  by = c("Protein" = "Protein.IDs")
) %>% distinct()

df <- left_join(
  cells_comparison_results,
  select(mq, Protein.IDs, Gene.simbol),
  by = c("Protein" = "Protein.IDs")
) %>% distinct()

library(eulerr)

set.seed(1)

s2 <- c(A = 1, B = 2)

plot(venn(s2))
plot(euler(s2), quantities = TRUE)

# Assuming seq_column represents the column name containing sequences in both data frames

# For MaxQuant data frame
unique_sequences_maxquant <- unique(evidence$Sequence)

# For Scaffold data frame
unique_sequences_scaffold <- unique(peptides$Peptide.sequence)
num_unique_sequences_scaffold <- length(unique_sequences_scaffold)

# Creating a data frame to use with eulerr package
library(eulerr)

euler_data <- list(
  MaxQuant = unique_sequences_maxquant,
  Scaffold = unique_sequences_scaffold
)

plot(euler(euler_data), quantities = TRUE)

pep_ms <- peptides %>% 
  dplyr::select(
    all_of(req_col)
  )

pep_ms_test <- pep_ms %>% 
  dplyr::group_by(Peptide.sequence, MS.MS.sample.name) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(
    n > 1,
    Peptide.sequence == "IWHHTFYNELR"
  )

pep_ms_test_new <- pep_ms_test %>%
  group_by(Peptide.sequence, MS.MS.sample.name) %>%
  mutate(
    Protein.accession.numbers = case_when(
      Assigned == TRUE ~ paste(Protein.accession.numbers, Other.Proteins, sep = ","),
      TRUE ~ Protein.accession.numbers
    ) 
  ) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate()

