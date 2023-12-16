library(data.table)
library(MSstats)
library(tidyverse)
library(MSstatsConvert)
library(stringr)
library(plotly)

evidence <- data.table::fread(
  "data/evidence.txt",
  check.names = TRUE
)

mq_peptides <- as.list(evidence$Sequence) %>% 
  dplyr::distinct()

proteins <- data.table::fread(
  "data/proteinGroups.txt",
  check.names = TRUE
)

test_df <- evidence %>% 
  dplyr::filter(Potential.contaminant == "+") %>% 
  dplyr::select(Proteins, Raw.file, Intensity, Leading.razor.protein) %>%
  merge(
    dplyr::select(proteins, Protein.IDs, Fasta.headers, Majority.protein.IDs),
    by.x = "Proteins", by.y = "Protein.IDs"
  ) %>% 
  dplyr::mutate(
    Fasta.headers = case_when(Fasta.headers = "" ~ Majority.protein.IDs)
  )

test_df %>% 
  ggplot(aes(x = Intensity, y = Raw.file, fill = Leading.razor.protein)) +
  geom_bar(stat = "identity")

test_df <- proteins %>% 
  dplyr::filter(Potential.contaminant == "+")

evidence <- evidence %>% 
  dplyr::filter(Potential.contaminant == "+")

test_df_long <- test_df %>% 
  dplyr::select(Protein.IDs, starts_with("LFQ")) %>% 
  left_join(
    dplyr::select(evidence, Proteins, Leading.razor.protein),
    by = join_by(Protein.IDs == Proteins)
  ) %>% 
  distinct() %>% 
  tidyr::pivot_longer(starts_with("LFQ"), names_to = "sample", values_to = "LFQ_intensity") %>% 
  dplyr::mutate(
    Proteins.ID = str_extract(Protein.IDs, "^[^;]+"),
    Leading.razor.protein = case_when(is.na(Leading.razor.protein) ~ Proteins.ID, TRUE ~ Leading.razor.protein)
  )

g <- test_df_long %>% 
  ggplot(aes(x = log2(LFQ_intensity), y = sample, fill = reorder(Leading.razor.protein, LFQ_intensity))) +
  geom_bar(stat = "identity")

ggplotly(g)

proteins %>% 
  dplyr::select(Protein.IDs, starts_with("LFQ"), Potential.contaminant) %>% 
  tidyr::pivot_longer(starts_with("LFQ"), names_to = "sample", values_to = "LFQ_intensity") %>% 
  ggplot(aes(x = LFQ_intensity, y = sample, fill = reorder(Potential.contaminant, LFQ_intensity),
             color = reorder(Potential.contaminant, LFQ_intensity))) +
  geom_bar(stat = "identity")

