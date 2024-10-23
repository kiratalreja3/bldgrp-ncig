#Variant Call Stats from PGGB

library(patchwork)


#Number of SNPs per Region
snp_dir <- "/g/data/te53/sj2852/blgrp/analyses/pangenome/analysis-v4/statistics/SNPs"

snp_files <- list.files(path = snp_dir, pattern = "*.txt", full.names = T)

dataframes_list <- list()

for (file in snp_files) {
  
  df <- read.table(file, header = TRUE)  
  
  df$SourceFile <- basename(file) 
  
  df$gene <- sub("\\..*", "", basename(file))
  
  dataframes_list[[file]] <- df
}

merged_df <- bind_rows(dataframes_list)

df_chm13 <- merged_df %>% dplyr::filter(grepl("chm13\\.vcf$", Sample))
df_grch38 <- merged_df %>% dplyr::filter(grepl("grch38\\.vcf$", Sample))

snp_summary <- merged_df %>%
  group_by(gene, sample_type = ifelse(grepl("chm13\\.vcf$", .data$Sample), "chm13", "grch38")) %>%
  summarise(number_of_SNPs = sum(number_of_SNPs, na.rm = TRUE), .groups = 'drop')

blgrp_regions_chm13 <- read_delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_allregionsref.bed", col_names = T)

blgrp_regions_chm13 <- blgrp_regions_chm13 %>%
  mutate(region_length = reg_end - reg_str)

blgrp_regions_grch38 <- read_delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/grch38/blgrp_coords_grch38_allregionsref.bed", col_names = T)

blgrp_regions_grch38 <- blgrp_regions_grch38 %>%
  mutate(region_length = reg_end - reg_str)

#Normalise # of SNPs

df_chm13 <- blgrp_regions_chm13 %>%
  mutate(source = "chm13") %>% 
  select(gene, region_length, source)

df_grch38 <- blgrp_regions_grch38 %>%
  mutate(source = "grch38") %>% 
  select(gene, region_length, source)

region_lengths <- bind_rows(df_chm13, df_grch38)

snp_summary <- region_lengths %>%
  inner_join(snp_summary, by = c("gene" = "gene", "source" = "sample_type")) %>%
  select(gene, region_length, everything())

snp_summary$normalised_snps <- snp_summary$number_of_SNPs / snp_summary$region_length


ggplot(snp_summary, aes(x = gene, y = normalised_snps, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Gene", y = "Number of SNPs (Normalised)", title = "Number of SNPs per Gene Compared to CHM13 and GRCH38") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Small and Large Variants (Reference Free)

var_dir <- "/g/data/te53/sj2852/blgrp/analyses/pangenome/vgdonstruct/vcf_parsing"

files <- list.files(var_dir, pattern = "_merged_vars_chr_cnt.txt$", recursive = TRUE, full.names = TRUE)

var_df_og <- files %>%
  map_dfr(~ read_delim(.x, delim = "\t"))


var_df <- var_df_og %>%
  mutate(smallvar_1to49 = indel_1to2 + indel_3to12 + indel_13to20 + indel_21to49)

var_df <- var_df %>%
  mutate(largevar_50up = indel_50to1K + indel_1K1to10K + indel_10K1to100K + indel_100K1to1M + indel_over1M1)

var_df <- var_df %>%
  select(-indel_1to2, -indel_3to12, -indel_13to20, -indel_21to49, -indel_50to1K, -indel_1K1to10K, -indel_10K1to100K, -indel_100K1to1M, -indel_over1M1)

var_df <- var_df %>%
  mutate(gene = str_extract(chr_mrg, "^[^:]+"))



blgrp_regions_chm13 <- read_delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_allregionsref.bed", col_names = T)

blgrp_regions_chm13 <- blgrp_regions_chm13 %>%
  mutate(region_length = `3_anchor_str` - `5_anchor_end`)

var_df <- var_df %>%
  left_join(blgrp_regions_chm13 %>% select(gene, region_length), by = "gene")

var_df$smallvar_norm <- var_df$smallvar_1to49 / var_df$region_length

var_df$largevar_norm <- var_df$largevar_50up / var_df$region_length

small_var<- ggplot(var_df, aes(x = gene, y = smallvar_norm)) + 
  geom_bar(stat = "identity", fill = "black") + 
  labs(x = "Gene", y = "Number of Small Variants <49bp (Normalised)", title = "Number of Small Variants Per Blood Group Region") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

large_var <- ggplot(var_df, aes(x = gene, y = largevar_norm)) + 
  geom_bar(stat = "identity", fill = "black") + 
  labs(x = "Gene", y = "Number of Large Variants >50bp (Normalised)", title = "Number of Large Variants Per Blood Group Region") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

small_var + large_var

#New Version - Stacked Barchart

var_df_largevar <- var_df_og %>%
  select(1, 7:11)

var_df_largevar <- var_df_largevar %>%
  mutate(gene = str_extract(chr_mrg, "^[^:]+"))

var_df_largevar <- var_df_largevar %>%
  left_join(blgrp_regions_chm13 %>% select(gene, region_length), by = "gene")

var_df_largevar <- var_df_largevar %>%
  mutate(across(2:6, ~ . / var_df_largevar[[8]]))

largevar_proportions <- var_df_largevar %>%
  mutate(row_total = rowSums(select(., 2:6), na.rm = TRUE)) %>%
  mutate(across(2:6, ~ . / row_total)) %>%
  select(-row_total)

largevar_long <- largevar_proportions %>%
  pivot_longer(cols = 2:6, names_to = "variable", values_to = "value")

custom_colours <- c("orange", "limegreen", "gold", "dodgerblue1", "firebrick")

prop <- ggplot(largevar_long, aes(x = gene, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene", y = "Proportion of Total Indels", fill = "Indel Size") +
  scale_fill_manual(values = custom_colours) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

large_var + prop


var_df_smallvar <- var_df_og %>%
  select(1, 3:6)

var_df_smallvar <- var_df_smallvar %>%
  mutate(gene = str_extract(chr_mrg, "^[^:]+"))

var_df_smallvar <- var_df_smallvar %>%
  left_join(blgrp_regions_chm13 %>% select(gene, region_length), by = "gene")

var_df_smallvar <- var_df_smallvar %>%
  mutate(across(2:5, ~ . / var_df_smallvar[[7]]))

smallvar_proportions <- var_df_smallvar %>%
  mutate(row_total = rowSums(select(., 2:5), na.rm = TRUE)) %>%
  mutate(across(2:5, ~ . / row_total)) %>%
  select(-row_total)

smallvar_long <- smallvar_proportions %>%
  pivot_longer(cols = 2:5, names_to = "variable", values_to = "value")

custom_colours_2 <- c("gold", "dodgerblue1", "firebrick", "limegreen")

prop_small <- ggplot(smallvar_long, aes(x = gene, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene", y = "Proportion of Total Indels", fill = "Indel Size") +
  scale_fill_manual(values = custom_colours_2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

small_var + prop_small
