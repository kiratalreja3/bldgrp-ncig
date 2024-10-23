#Mapping Quality

libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))

install.packages("ggplot")
install.packages("dyplr")
install.packages("readr")
install.packages("tidyverse")
install.packages("pafr")

library(tidyverse)
library(ggplot2)

#Mapping Quality Per Sample

sample_directory <- ("/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/mapping_quality/sample_level")

file_list <- list.files(sample_directory, full.names = T)

#Empty list to store dataframes
data_list <- list()

#read each file into dataframe
for (file in file_list) {
  df <- read.table(file, header = F, col.names = c("Count", "Mapping_Quality"))
  
  df$Sample <- gsub("_mapping_quality.txt", "", basename(file))
  
  data_list[[file]] <- df

}

all_data <- do.call(rbind, data_list)


average_mapping_quality_sample <- all_data %>%
  group_by(Sample) %>%
  summarize(AverageMQ = sum(Mapping_Quality * Count) / sum(Count))

average_mapping_quality_sample <- average_mapping_quality_sample %>%
  mutate(Group = ifelse(startsWith(Sample, "P"), "NCIG", "KGP"))


#Plot 
ggplot(average_mapping_quality_sample, aes(x = Group, y = AverageMQ))+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, outlier.shape = NA, colour = "black", aes(fill = Group))+
  scale_fill_manual(values = c("NCIG" = "orange3", "KGP" = "lightblue")) + 
  labs(x = "Samples", 
       y= "Average Mapping Quality")+
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))


#combine all dataframes into single dataframe
all_data <- bind_rows(data_list)

#Mapping Quality Per Region - done in vs code





