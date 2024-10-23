


library(tidyverse)

root_dir <- "/g/data/te53/sj2852/blgrp/analyses/flye_assemblies/NEW_FLYE/hapdup/output_files"

sub_dirs <- list.dirs(root_dir, full.names = T, recursive = F)

all_dataframes <- list()

for (folder in sub_dirs) {
  
  hp1_file <- file.path(folder, "phased_blocks_hp1.bed")
  hp2_file <- file.path(folder, "phased_blocks_hp2.bed")
  
  if (file.exists(hp1_file)) {
    hp1_df <- read_tsv(hp1_file, col_names = FALSE, show_col_types = FALSE) %>%
      mutate(hp = "hp1", folder_name = basename(folder))  # Add 'hp' and 'folder_name' columns
    all_dataframes <- append(all_dataframes, list(hp1_df))
  }
  
  if (file.exists(hp2_file)) {
    hp2_df <- read_tsv(hp2_file, col_names = FALSE, show_col_types = FALSE) %>%
      mutate(hp = "hp2", folder_name = basename(folder))  # Add 'hp' and 'folder_name' columns
    all_dataframes <- append(all_dataframes, list(hp2_df))
  }
}

full_df <- bind_rows(all_dataframes)

#Add in Acnestry

full_df <- full_df %>%
  mutate(ancestry = case_when(
    str_starts(folder_name, "P") ~ "NCIG",
    TRUE ~ "KGP" 
  ))

KGP_ancestry <- read_tsv("/g/data/te53/sj2852/blgrp/metadata/ancestry_files/merged_ancestry.txt")

KGP_ancestry <- KGP_ancestry %>% dplyr::rename(folder_name = V1)

full_df <- full_df %>%
  left_join(KGP_ancestry %>%
              select(folder_name, Description), by = "folder_name") %>%
  mutate(ancestry = ifelse(!is.na(Description), Description, ancestry)) %>%
  select(-Description)


#Count unique contigs

count_unique_contigs <- full_df %>%
  group_by(folder_name, hp, ancestry) %>% 
  summarize(unique_contig_count = n_distinct(X1), .groups = "drop")

count_unique_contigs_nohap <- count_unique_contigs %>%
  group_by(folder_name, ancestry) %>%
  summarize(total_contigs_per_sample = sum(unique_contig_count), .groups = "drop")


count_unique_contigs <- count_unique_contigs %>%
  mutate(folder_name = factor(folder_name, 
                              levels = unique(folder_name[order(ancestry)])))

count_unique_contigs$folder_name <- factor(count_unique_contigs$folder_name, 
                                           levels = unique(count_unique_contigs$folder_name))

colour_palette <- c("AFR" = "gold2", "EUR" = "deepskyblue", "AMR" = "red3", "EAS" = "green2", "SAS" = "purple3", "NCIG" = "darkorange3")


ggplot(count_unique_contigs, aes(x = folder_name, y = unique_contig_count, fill = ancestry)) +
  geom_bar(stat = "identity", position = "dodge") + 
  labs(
    x = "Sample ID", 
    y = "Number of Contigs", 
    title = "Contig Count by Sample and Ancestry"
  ) +
  scale_fill_manual(values = colour_palette) +
  geom_hline(yintercept = 45, color = "black", linetype = "solid", size = 1) +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(count_unique_contigs_nohap, aes(x = folder_name, y = total_contigs_per_sample, fill = ancestry)) +
  geom_bar(stat = "identity", position = "dodge") + 
  labs(
    x = "Sample ID", 
    y = "Number of Contigs", 
    title = "Contig Count by Sample and Ancestry"
  ) +
  scale_fill_manual(values = colour_palette) +
  geom_hline(yintercept = 90, color = "black", linetype = "solid", size = 1) +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


#Types of Contigs

unique_type_counts <- full_df %>%
  mutate(X4_clean = str_extract(X4, "^MissingConcordancy|^Discordancy|^ContigEnd")) %>%  # Extract only core terms
  dplyr::filter(!is.na(X4_clean)) %>%  # Remove rows where X4_clean is NA (i.e., not matching the core terms)
  group_by(folder_name, ancestry, X4_clean) %>%  # Group by folder_name, ancestry, and cleaned X4
  summarize(count = n(), .groups = "drop") %>%  # Count occurrences
  pivot_wider(
    names_from = X4_clean,        # Make each unique cleaned X4 value a column
    values_from = count,          # Fill with the count of each unique X4 value
    values_fill = list(count = 0)  # Fill missing values with 0
  )
