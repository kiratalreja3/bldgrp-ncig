#plotting locityper results

#this code plots the qv_locityper figure

.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))

install.packages("ggplot")
install.packages("tidyverse")
install.packages("jsonlite")

library(ggplot2)
library(tidyverse)
library(jsonlite)
library(tidyr)
library(dplyr)
library(purrr)
library(patchwork)

#load in files
json_files_NCIG <- list.files(path = "/g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/json_NCIG", pattern = "*.json", full.names = T)

json_files_KGP <- list.files(path="/g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/json_KGP", pattern = "*.json", full.names = T)


read_json_file <- function(file) {
  sample_name <- sub(".*/(.*)_res\\.json$", "\\1", basename(file))
  json_data <- fromJSON(file, flatten = TRUE)  
  json_data$options <- NULL
  json_data$sample_name <- sample_name
  return(json_data)
}

json_data_list_NCIG <- lapply(json_files_NCIG, read_json_file)
json_data_list_KGP <- lapply(json_files_KGP, read_json_file)

NCIG_data <- do.call(rbind, json_data_list_NCIG)
KGP_data <- do.call(rbind, json_data_list_KGP)

#make sure to check if the NCIG or KGP_data file is a dataframe or matrix! convert to dataframe if it is!

NCIG_data <- as.data.frame(NCIG_data)
KGP_data <- as.data.frame(KGP_data)

NCIG_data$quality <- as.numeric(NCIG_data$quality)
KGP_data$quality <- as.numeric(KGP_data$quality)

#edit some unwanted columns

NCIG_data <- subset(NCIG_data, select = -sample_name)
KGP_data <- subset(KGP_data, select = -sample_name)

NCIG_data <- NCIG_data %>%
  dplyr::rename(sample_name = warnings)

KGP_data <- KGP_data %>%
  dplyr::rename(sample_name = warnings)

if (is.list(NCIG_data$locus)) {
  NCIG_data$locus <- sapply(NCIG_data$locus, as.character)
}

if (is.list(KGP_data$locus)) {
  KGP_data$locus <- sapply(KGP_data$locus, as.character)
}

NCIG_data$locus <- as.factor(NCIG_data$locus)
KGP_data$locus <- as.factor(KGP_data$locus)

NCIG_data$quality <- as.factor(NCIG_data$quality)
KGP_data$quality <- as.factor(KGP_data$quality)

str(NCIG_data)
str(KGP_data)


#plotting QV scores

NCIG_data$quality_bin <- cut(NCIG_data$quality, 
                                 breaks = c(-Inf, 17, 23, 33, Inf), 
                                 labels = c("1", "2", "3", "4"), 
                                 right = F)

KGP_data$quality_bin <- cut(KGP_data$quality, 
                            breaks = c(-Inf, 17, 23, 33, Inf), 
                            labels = c("1", "2", "3", "4"), 
                            right = F)
                            

NCIG_data$quality_bin <- as.factor(NCIG_data$quality_bin)
KGP_data$quality_bin <- as.factor(KGP_data$quality_bin)

NCIG_data$locus <- factor(NCIG_data$locus, levels = c("REFERENCE", sort(Biostrings::setdiff(unique(NCIG_data$locus), "REFERENCE"))))
KGP_data$locus <- factor(KGP_data$locus, levels = c("REFERENCE", sort(Biostrings::setdiff(unique(KGP_data$locus), "REFERENCE"))))


ggplot(NCIG_data, aes(x = locus, fill = quality_bin)) +
  geom_bar(position = "stack")+
  scale_fill_manual(values = c("1" = "red3", "2" = "yellow2", "3" = "lightgreen", "4" = "turquoise3"), 
                    labels = c("0-17", "17-23", "23-33", "33+"))+
  labs(x = "Blood Group Region", y = "Sample Count (NCIG)", fill = "Accuracy (QV)") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 30))+
  theme(axis.text = element_text(angle = 45, hjust = 1))

ggplot(KGP_data, aes(x = locus, fill = quality_bin)) +
  geom_bar(position = "stack")+
  scale_fill_manual(values = c("1" = "red3", "2" = "yellow2", "3" = "lightgreen", "4" = "turquoise3"), 
                    labels = c("0-17", "17-23", "23-33", "33+"))+
  labs(x = "Blood Group Region", y = "Sample Count (KGP)", fill = "Accuracy (QV)") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, by = 30))+
  theme(axis.text = element_text(angle = 45, hjust = 1))


#plotting Number of Unexplained Reads


NCIG_data$unexpl_reads <- as.numeric(NCIG_data$unexpl_reads)
KGP_data$unexpl_reads <- as.numeric(KGP_data$unexpl_reads)

NCIG_data$locus <- factor(NCIG_data$locus, levels = c("REFERENCE", sort(Biostrings::setdiff(unique(NCIG_data$locus), "REFERENCE"))))
KGP_data$locus <- factor(KGP_data$locus, levels = c("REFERENCE", sort(Biostrings::setdiff(unique(KGP_data$locus), "REFERENCE"))))

KGP_data$unexplained_bin <- cut(KGP_data$unexpl_reads, 
                                 breaks = c(-Inf, 10, 100, 1000, Inf), 
                                 labels = c("1", "2", "3"), 
                                 right = F)

NCIG_data$unexplained_bin <- cut(NCIG_data$unexpl_reads, 
                                     breaks = c(-Inf, 100, 1000, Inf), 
                                     labels = c("1", "2", "3"), 
                                     right = F)

ggplot(KGP_data, aes(x = locus, fill = unexplained_bin)) +
  geom_bar(position = "stack")+
  scale_fill_manual(values = c("1" = "lightgreen", "2" = "yellow2", "3" = "red3"), 
                    labels = c("<100", "101-1000", "1000+"))+
  labs(x = "Blood Group Region", y = "Sample Count (KGP)", fill = "Number of Unexplained Reads") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 178, by = 30))+
  theme(axis.text = element_text(angle = 45, hjust = 1))

ggplot(NCIG_data, aes(x = locus, fill = unexplained_bin)) +
  geom_bar(position = "stack")+
  scale_fill_manual(values = c("1" = "lightgreen", "2" = "yellow2", "3" = "red3"), 
                    labels = c("<100", "101-1000", "1000+"))+
  labs(x = "Blood Group Region", y = "Sample Count (NCIG)", fill = "Number of Unexplained Reads") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 30))+
  theme(axis.text = element_text(angle = 45, hjust = 1))

#plotting Weight Distribution

NCIG_data$weight_dist <- unlist(NCIG_data$weight_dist)
KGP_data$weight_dist <- unlist(KGP_data$weight_dist)

NCIG_data$weight_dist <- as.numeric(as.character(NCIG_data$weight_dist))
KGP_data$weight_dist <- as.numeric(as.character(KGP_data$weight_dist))


NCIG_data$weight_dist_bin <- cut(NCIG_data$weight_dist, 
                                     breaks = c(-Inf, 1, 50, Inf), 
                                     labels = c("1", "2", "3"), 
                                     right = F)

KGP_data$weight_dist_bin <- cut(KGP_data$weight_dist, 
                                 breaks = c(-Inf, 1, 50, Inf), 
                                 labels = c("1", "2", "3"), 
                                 right = F)

KGP_data <- KGP_data[-4110, ]
NCIG_data <- NCIG_data[-5247, ]

ggplot(NCIG_data, aes(x = locus, fill = weight_dist_bin)) +
  geom_bar(position = "stack")+
  scale_fill_manual(values = c("1" = "lightgreen", "2" = "yellow2", "3" = "red3"), 
                    labels = c("<1", "1-50", "51+"))+
  labs(x = "Blood Group Region", y = "Sample Count (NCIG)", fill = "Weight Distribution of Genotype Prediction") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 30))+
  theme(axis.text = element_text(angle = 45, hjust = 1))

ggplot(KGP_data, aes(x = locus, fill = weight_dist_bin)) +
  geom_bar(position = "stack")+
  scale_fill_manual(values = c("1" = "lightgreen", "2" = "yellow2", "3" = "red3"), 
                    labels = c("<1", "1-50", "51+"))+
  labs(x = "Blood Group Region", y = "Sample Count (KGP)", fill = "Weight Distribution of Genotype Prediction") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 178, by = 30))+
  theme(axis.text = element_text(angle = 45, hjust = 1))


setwd("/g/data/te53/sj2852/blgrp/data/figures")


#Counting the Number of Quality Scores Per Sample/Region


quality_NCIG <- NCIG_data %>%
  group_by(locus) %>%
  summarise(
    less_than_20 = sum(quality <= 20, na.rm = T), 
    greater_than_20 = sum(quality > 20, na.rm = T)
  )

quality_KGP <- KGP_data %>%
  group_by(locus) %>%
  summarise(
    less_than_20 = sum(quality <= 20, na.rm = T), 
    greater_than_20 = sum(quality > 20, na.rm = T)
  )

quality_KGP$less_than_20 <- quality_KGP$less_than_20 *-1
quality_KGP$greater_than_20 <- quality_KGP$greater_than_20 *-1

quality_KGP$Group <- "KGP"
quality_NCIG$Group <- "NCIG"

combined_df <- rbind(quality_KGP, quality_NCIG)

df_long <- combined_df %>%
  pivot_longer(cols = c("less_than_20", "greater_than_20"),
               names_to = "Category",
               values_to = "Count")

df_long1 <- df_long %>%
  mutate(locus = factor(locus, levels = df_long %>%
                          dplyr::filter(Group == "NCIG", Category == "greater_than_20") %>%
                          arrange(desc(Count)) %>%
                          pull(locus)))

ggplot(df_long1, aes(x = locus, y = Count, fill = interaction(Category, Group))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = c("less_than_20.KGP" = "darkred",
                               "greater_than_20.KGP" = "limegreen",
                               "less_than_20.NCIG" = "darkred",
                               "greater_than_20.NCIG" = "limegreen"
                               )) +
  labs(y = "KGP Samples                   NCIG Samples", x = "Blood Group Genes", 
       title = "QV Score Distribution Across Blood Group Genes for Two Datasets", fill = "Samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.key = element_rect(fill = "transparent", color = "transparent"), 
        legend.key.size = unit(0.6, "cm"),
        legend.box.background = element_blank()) + 
  guides(fill = guide_legend(title = "QV Score", override.aes = list(fill = c("limegreen", "darkred"), 
  
                                                                                                                                        color = c("limegreen", "darkred")))) 
#Plotting by Ancestry

KGP_ancestry <- read_delim("/g/data/te53/sj2852/blgrp/metadata/ancestry_files/merged_ancestry.txt", col_names = F)

KGP_data <- KGP_data %>%
  mutate(sample_name_trimmed = str_extract(sample_name, "^[^_]+"))

KGP_data <- KGP_data %>%
  left_join(KGP_ancestry, by = c("sample_name_trimmed" = "X1"))

#AFR

KGP_quality_AFR <- KGP_data %>%
  dplyr::filter(X2 == "AFR") %>%  
  group_by(locus) %>%  
  summarise(
    less_than_20 = sum(quality <= 20, na.rm = TRUE), 
    greater_than_20 = sum(quality > 20, na.rm = TRUE)  
  )

AFR_long <- KGP_quality_AFR %>%
  pivot_longer(cols = c(greater_than_20, less_than_20), names_to = "category", values_to = "count")

ggplot(AFR_long, aes(x = locus, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(x = "Locus", y = "Count", fill = "Category") +
  theme_minimal() +
  scale_fill_manual(values = c("greater_than_20" = "limegreen", "less_than_20" = "darkred")) +
  labs(y = "Samples (AFR)", x = "Blood Group Genes", 
       title = "QV Score Distribution Across Blood Group Genes for KGP-AFR Samples", fill = "QV Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#AMR

KGP_quality_AMR <- KGP_data %>%
  dplyr::filter(X2 == "AMR") %>%  
  group_by(locus) %>%  
  summarise(
    less_than_20 = sum(quality <= 20, na.rm = TRUE),  
    greater_than_20 = sum(quality > 20, na.rm = TRUE)  
  )

AMR_long <- KGP_quality_AMR %>%
  pivot_longer(cols = c(greater_than_20, less_than_20), names_to = "category", values_to = "count")

ggplot(AMR_long, aes(x = locus, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(x = "Locus", y = "Count", fill = "Category") +
  theme_minimal() +
  scale_fill_manual(values = c("greater_than_20" = "limegreen", "less_than_20" = "darkred")) +
  labs(y = "Samples (AMR)", x = "Blood Group Genes", 
       title = "QV Score Distribution Across Blood Group Genes for KGP-AMR Samples", fill = "QV Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#EUR

KGP_quality_EUR <- KGP_data %>%
  dplyr::filter(X2 == "EUR") %>%  
  group_by(locus) %>%  
  summarise(
    less_than_20 = sum(quality <= 20, na.rm = TRUE),  
    greater_than_20 = sum(quality > 20, na.rm = TRUE)  
  )

EUR_long <- KGP_quality_EUR %>%
  pivot_longer(cols = c(greater_than_20, less_than_20), names_to = "category", values_to = "count")

ggplot(EUR_long, aes(x = locus, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(x = "Locus", y = "Count", fill = "Category") +
  theme_minimal() +
  scale_fill_manual(values = c("greater_than_20" = "limegreen", "less_than_20" = "darkred")) +
  labs(y = "Samples (EUR)", x = "Blood Group Genes", 
       title = "QV Score Distribution Across Blood Group Genes for KGP-EUR Samples", fill = "QV Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#EAS

KGP_quality_EAS <- KGP_data %>%
  dplyr::filter(X2 == "EAS") %>%  
  group_by(locus) %>%  
  summarise(
    less_than_20 = sum(quality <= 20, na.rm = TRUE),  
    greater_than_20 = sum(quality > 20, na.rm = TRUE)  
  )

EAS_long <- KGP_quality_EAS %>%
  pivot_longer(cols = c(greater_than_20, less_than_20), names_to = "category", values_to = "count")

ggplot(EAS_long, aes(x = locus, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(x = "Locus", y = "Count", fill = "Category") +
  theme_minimal() +
  scale_fill_manual(values = c("greater_than_20" = "limegreen", "less_than_20" = "darkred")) +
  labs(y = "Samples (EAS)", x = "Blood Group Genes", 
       title = "QV Score Distribution Across Blood Group Genes for KGP-EAS Samples", fill = "QV Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#SAS

KGP_quality_SAS <- KGP_data %>%
  dplyr::filter(X2 == "SAS") %>%  
  group_by(locus) %>%  
  summarise(
    less_than_20 = sum(quality <= 20, na.rm = TRUE),  
    greater_than_20 = sum(quality > 20, na.rm = TRUE)  
  )

SAS_long <- KGP_quality_SAS %>%
  pivot_longer(cols = c(greater_than_20, less_than_20), names_to = "category", values_to = "count")

ggplot(SAS_long, aes(x = locus, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(x = "Locus", y = "Count", fill = "Category") +
  theme_minimal() +
  scale_fill_manual(values = c("greater_than_20" = "limegreen", "less_than_20" = "darkred")) +
  labs(y = "Samples (SAS)", x = "Blood Group Genes", 
       title = "QV Score Distribution Across Blood Group Genes for KGP-SAS Samples", fill = "QV Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Investigate <20

AFR_selected <- KGP_quality_AFR %>% select(locus, less_than_20)
AMR_selected <- KGP_quality_AMR %>% select(locus, less_than_20)
EAS_selected <- KGP_quality_EAS %>% select(locus, less_than_20)
EUR_selected <- KGP_quality_EUR %>% select(locus, less_than_20)
SAS_selected <- KGP_quality_SAS %>% select(locus, less_than_20)
NCIG_selected <- quality_NCIG %>% select(locus, less_than_20)


AFR_selected <- AFR_selected %>% dplyr::rename(AFR = less_than_20)
AMR_selected <- AMR_selected %>% dplyr::rename(AMR = less_than_20)
EAS_selected <- EAS_selected %>% dplyr::rename(EAS = less_than_20)
EUR_selected <- EUR_selected %>% dplyr::rename(EUR = less_than_20)
SAS_selected <- SAS_selected %>% dplyr::rename(SAS = less_than_20)
NCIG_selected <- NCIG_selected %>% dplyr::rename(NCIG = less_than_20)

merged_ancestry_df <- purrr::reduce(list(AFR_selected, AMR_selected, EAS_selected, EUR_selected, SAS_selected, NCIG_selected), full_join, by = "locus")

merged_ancestry_df <- merged_ancestry_df %>%
  mutate(
    AFR = AFR / 64,
    AMR = AMR / 20,
    EAS = EAS / 28,
    EUR = EUR / 30,
    SAS = SAS / 36,
    NCIG = NCIG / 142
  )


merged_ancestry_1 <- merged_ancestry_df %>% dplyr::slice(1:22)    
merged_ancestry_2 <- merged_ancestry_df %>% dplyr::slice(23:n())

ma1_long <- merged_ancestry_1 %>%
  pivot_longer(cols = c(AFR, AMR, EAS, EUR, SAS, NCIG),   
               names_to = "category",              
               values_to = "value") 

ma2_long <- merged_ancestry_2 %>%
  pivot_longer(cols = c(AFR, AMR, EAS, EUR, SAS, NCIG),   
               names_to = "category",              
               values_to = "value")

bar_colors <- c("AFR" = "orange",
                "AMR" = "red3",
                "EAS" = "green2",
                "EUR" = "deepskyblue",
                "SAS" = "purple3",
                "NCIG" = "darkorange3")

ma1 <- ggplot(ma1_long, aes(x = locus, y = value, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +  
  scale_fill_manual(values = bar_colors) +           
  labs(x = NULL, y = NULL, fill = "Ancestral Superpopulation") +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ma2 <- ggplot(ma2_long, aes(x = locus, y = value, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +  
  scale_fill_manual(values = bar_colors) +           
  labs(x = NULL, y = NULL, fill = "Ancestral Superpopulation") +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot <- ma1 / ma2 +
  labs(x = "Blood Group Gene", y = "                                                                                      Number of Samples QV <20")

combined_plot


#Plot Total Number of Haplotypes

combined_df <- rbind(NCIG_data, KGP_data)

combined_df <- combined_df %>%
  separate(genotype, into = c("hap1", "hap2"), sep = ",", fill = "right")

hap_counts <- combined_df %>%
  mutate(sample_group = ifelse(grepl("^P", sample_name), "NCIG", "KGP")) %>%
  group_by(locus, sample_group) %>%
  summarise(unique_haplotypes = n_distinct(na.omit(c(hap1, hap2)))) %>%
  ungroup()

hap_counts_neg <- hap_counts %>%
  mutate(unique_haplotypes = if_else(sample_group == "NCIG", -unique_haplotypes, unique_haplotypes))

custom_colours <- c("KGP" = "skyblue", "NCIG" = "orange2")

ggplot(hap_counts_neg, aes(x = locus, y = unique_haplotypes, fill = sample_group)) +
  geom_bar(stat = "identity") +
  labs(x = "Locus", y = "Count of Unique Haplotypes", title = "Distinct Haplotypes by Locus") +
  theme_minimal() +
  scale_fill_manual(values = custom_colours) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

NCIG_group <- combined_df %>%
  dplyr::filter(grepl("^P", sample_name)) %>%
  select(locus, hap1, hap2)

KGP_group <- combined_df %>%
  dplyr::filter(!grepl("^P", sample_name)) %>%
  select(locus, hap1, hap2)

count_unique_haplotypes <- function(locus_name) {
  haplotypes_NCIG <- unique(c(NCIG_group %>% dplyr::filter(locus == locus_name) %>% select(hap1, hap2) %>% unlist()))
  haplotypes_KGP <- unique(c(KGP_group %>% dplyr::filter(locus == locus_name) %>% select(hap1, hap2) %>% unlist()))
  
  unique_in_NCIG <- base::setdiff(haplotypes_NCIG, haplotypes_KGP)
  unique_in_KGP <- base::setdiff(haplotypes_KGP, haplotypes_NCIG)
  
  return(data.frame(locus = locus_name, unique_in_NCIG = length(unique_in_NCIG), unique_in_KGP = length(unique_in_KGP)))
}

unique_loci <- unique(combined_df$locus)

result_df <- bind_rows(lapply(unique_loci, count_unique_haplotypes))
