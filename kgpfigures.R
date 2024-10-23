.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))
install.packages("ggplot")
install.packages("dyplr")

library(tidyverse)
library(ggplot2)


setwd("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords")

chm13 <- read_tsv("/g/data/te53/sj2852/blgrp/metadata/chm13_combined_assembly_stats_transposed_sorted.tsv", col_names = F)
grch38 <- read_tsv("/g/data/te53/sj2852/blgrp/metadata/grch38_combined_assembly_stats_longform_transposed_sorted.tsv", col_names = F)

chm13 <- within(chm13, rm(X2))
grch38 <- within(grch38, rm(X2))


#Combine two stats files together for one big metadata file
combined <- rbind(chm13, grch38)

combined_size_sorted <- combined %>%
  arrange(desc(X9))

combined_size_sorted$X1 <- factor(combined_size_sorted$X1, levels = combined_size_sorted$X1)


#Plot Number of Mapped Reads Per Sample

ggplot(combined_size_sorted, aes(x = X1, y = X9))+
  geom_col(width = 0.5)+
  labs(title = "# Mapped Reads Per Sample", 
       x = "Number of Mapped Reads", 
       y = "KGP Samples")+
  geom_hline(yintercept = 13945.57, color = "red", linetype = "dotted", size = 1) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

mean(combined_size_sorted$X9)
