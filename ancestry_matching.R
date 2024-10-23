#this code matched 178 KGP-ONT samples to their ancestry and created the files 
#/g/data/te53/sj2852/blgrp/metadata/ancestry_files/merged_ancestry.txt
#/g/data/te53/sj2852/blgrp/metadata/ancestry_files/merged_ancestry.csv

library(tidyverse)
library(dplyr)

meta_dataframe <- read.csv("/g/data/te53/sj2852/blgrp/metadata/ancestry_files/coreill_catalog_export_sample_ancestry.csv", header = T)

basenames <- read.table("/g/data/te53/sj2852/blgrp/metadata/ancestry_files/sample_list.txt")

basenames <- basenames %>%
  dplyr::rename(ID = V1)


sample_ancestry <- basenames %>%
  inner_join(meta_dataframe, by = "ID")

sample_ancestry <- sample_ancestry %>%
  distinct(ID, .keep_all = T)

sample_ancestry <- sample_ancestry[, -c(3, 4, 5, 6, 7, 9)]
  
sample_ancestry <- sample_ancestry %>%
  mutate(
    Description = case_when(
      Description %in% c("YORUBA IN IBADAN, NIGERIA", 
                         "LUHYA IN WEBUYE, KENYA",
                         "AFRICAN CARIBBEAN IN BARBADOS",
                         "GAMBIAN IN WESTERN DIVISION - MANDINKA",
                         "ESAN IN NIGERIA",
                         "MENDE IN SIERRA LEONE") ~ "AFR",
      Description %in% c("BRITISH FROM ENGLAND AND SCOTLAND",
                         "FINNISH IN FINLAND",
                         "IBERIAN POPULATIONS IN SPAIN") ~ "EUR",
      Description %in% c("HAN CHINESE SOUTH, CHINA",
                         "CHINESE DAI IN XISHUANGBANNA, CHINA", 
                         "KINH IN HO CHI MINH CITY, VIETNAM") ~ "EAS",
      Description %in% c("PUERTO RICAN IN PUERTO RICO", 
                         "PERUVIAN IN LIMA, PERU", 
                         "COLOMBIAN IN MEDELLIN, COLOMBIA") ~ "AMR",
      Description %in% c("PUNJABI IN LAHORE, PAKISTAN",
                         "BENGALI IN BANGLADESH", 
                         "SRI LANKAN TAMIL IN THE UK", 
                         "INDIAN TELUGU IN THE UK") ~ "SAS",
      TRUE ~ Description
    )
  ) %>%
  distinct(ID, .keep_all = TRUE)


write.csv(sample_ancestry, "/g/data/te53/sj2852/blgrp/metadata/ancestry_files/merged_ancestry.csv", row.names = FALSE)

write.table(sample_ancestry, "/g/data/te53/sj2852/blgrp/metadata/ancestry_files/merged_ancestry.txt", sep = "\t", row.names = F, quote = F)

#Add in Assembly Alignment Details

chm13 <- read.table("/g/data/te53/sj2852/blgrp/metadata/ancestry_files/chm13_samples.txt", header = F)

chm13$alignment <- "chm13"

grch38 <- read.table("/g/data/te53/sj2852/blgrp/metadata/ancestry_files/grch38_samples.txt", header = F)

grch38$alignment <- "grch38"

alignment_combined <- rbind(chm13, grch38)

ancestry_alignment <- left_join(sample_ancestry, alignment_combined, by = c("ID" = "V1"))

#Plot Ancestry

colour_mapping <- c("AFR" = "darkorange2", "EUR" = "deepskyblue", "AMR" = "red3", "EAS" = "green2", "SAS" = "purple3")

colour_mapping_new <-   c("AFR.chm13" = "darkorange2",    
"AFR.grch38" = "orange",         
"EUR.chm13" = "deepskyblue",     
"EUR.grch38" = "skyblue",        
"AMR.chm13" = "red3",            
"AMR.grch38" = "red",            
"EAS.chm13" = "green2",          
"EAS.grch38" = "lightgreen",    
"SAS.chm13" = "purple3",         
"SAS.grch38" = "mediumpurple", 
"NCIG.chm13" = "darkorange3"
)

sample_counts <- sample_ancestry %>%
  group_by(Description) %>%
  summarise(Count = n())

alignment_counts <- ancestry_alignment %>%
  group_by(Description, alignment) %>%
  summarise(Count = n(), .groups = "drop")

new_rows <- data.frame(
  Description = c("NCIG", "NCIG"),
  alignment = c("chm13", "grch38"),
  Count = c(142, 0)
)

alignment_counts <- rbind(alignment_counts, new_rows)

alignment_counts$Description <- factor(alignment_counts$Description, 
                                       levels = c("NCIG", "AFR", "AMR", "EUR", "EAS", "SAS" ))

alignment_counts_filtered <- alignment_counts %>% dplyr::filter(Count != 0)

ggplot(alignment_counts_filtered, aes(x = Description, y = Count, fill = interaction(Description, alignment))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colour_mapping_new) +
  theme_bw() +
  coord_flip() +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14), 
        legend.position = "none") +
  geom_text(aes(label = alignment), 
            position = position_stack(vjust = 0.5), 
            size = 4, 
            color = "black") +
  labs(title = "Ancestry of 320 Samples (KGP-ONT and NCIG)",
       x = "Superpopulation",
       y = "Count")


