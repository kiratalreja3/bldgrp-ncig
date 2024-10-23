.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))

install.packages("ggplot")
install.packages("dyplr")
install.packages("readr")
install.packages("tidyverse")

library(tidyverse)
library(ggplot2)
library(readr)
library(dyplr)


setwd("/g/data/te53/sj2852/blgrp/data/KGPdata")

chm13 <- read_tsv("/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/chm13_combined_assembly_stats_longform_transposed.tsv", col_names = T)

grch38 <- read_tsv("/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/grch38_combined_assembly_stats_longform_transposed.tsv", col_names = T)

NCIG <- read_tsv("/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/NCIG_combined_assembly_stats_longform_transposed.tsv", col_names = T)



column_SID <- c("Sample_ID",
                "NA", 
                "raw total sequences",
                "filtered sequences",
                "sequences",
                "is sorted",
                "1st fragments",
                "last fragments",
                "reads mapped",
                "reads mapped and paired",
                "reads unmapped",
                "reads properly paired",
                "reads paired",
                "reads duplicated",
                "reads MQ0",
                "reads QC failed",
                "non-primary alignments",
                "supplementary alignments",
                "total length",
                "total first fragment length",
                "total last fragment length",
                "bases mapped",
                "bases mapped (cigar)",
                "bases trimmed",
                "bases duplicated",
                "mismatches",
                "error_rate",
                "average_length",
                "average first fragment length",
                "average last fragment length",
                "maximum length",
                "maximum first fragment length",
                "maximum last fragment length",
                "average_quality",
                "insert size average",
                "insert size standard deviation",
                "inward oriented pairs",
                "outward oriented pairs",
                "pairs with other orientation",
                "pairs on different chromosomes")

colnames(chm13) <- column_SID
colnames(grch38) <- column_SID
colnames(NCIG) <- column_SID

#Merge Dataframes Together

mergedKGP <- bind_rows(chm13, grch38, NCIG)

#Cleaning

mergedKGP <- mergedKGP %>%
  mutate(Sample_ID = gsub("\\.blgrp$", "", Sample_ID))

mergedKGP <- mergedKGP %>%
  mutate(Group = ifelse(startsWith(Sample_ID, "P"), "NCIG", "KGP"))

#Violin Plot of Average Read Length
ggplot(mergedKGP, aes(x = Group, y = average_length))+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, outlier.shape = NA, colour = "black", aes(fill = Group))+
  scale_fill_manual(values = c("NCIG" = "orange3", "KGP" = "lightblue"))+
  labs(x = "", 
      y= "Average Read Length")+
  theme_minimal()

#BoxPlot
mean_value_length <- mean(mergedKGP$average_length, na.rm = T)

ggplot(mergedKGP, aes(x = Sample_ID, y = average_length))+
  geom_col(fill = "navyblue", colour = "black")+
  geom_hline(yintercept = mean_value, colour = "red", linetype = "dashed", linewidth = 1)+
  labs(x = "Sample ID", 
       y = "Average Read Length") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  


#Violin Plot of Average Quality
ggplot(mergedKGP, aes(x = Group, y = average_quality))+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, outlier.shape = NA, colour = "black", aes(fill = Group))+
  scale_fill_manual(values = c("NCIG" = "orange3", "KGP" = "lightblue"))+
  labs(x = "", 
       y= "Average Read Quality")+
  theme_minimal()

#Boxplot
mean_value_qual <- mean(mergedKGP$average_quality, na.rm = T)

ggplot(mergedKGP, aes(x = Sample_ID, y = average_quality))+
  geom_col(fill = "navyblue", colour = "black")+
  geom_hline(yintercept = mean_value_qual, colour = "red", linetype = "dashed", linewidth = 1)+
  labs(x = "Sample ID", 
       y = "Average Read Quality") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))



#Violin Plot Error Rate

ggplot(mergedKGP, aes(x = Group, y = error_rate))+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, outlier.shape = NA, colour = "black", aes(fill = Group))+
  scale_fill_manual(values = c("NCIG" = "orange3", "KGP" = "lightblue"))+
  labs(x = "", 
       y= "Average Error Rate Per Sample")+
  theme_minimal()

#Boxplot
mean_value_depth <- mean(merged_long$Average_Coverage, na.rm = T)

ggplot(merged_long, aes(x = sample_id, y = Average_Coverage))+
  geom_col(fill = "navyblue", colour = "black")+
  geom_hline(yintercept = mean_value_depth, colour = "red", linetype = "dashed", linewidth = 1)+
  labs(x = "Sample ID", 
       y = "Average Read Depth") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))


