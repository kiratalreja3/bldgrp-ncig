# Plotting Closeness and Betweenness

library(tidyverse)
library(ggplot2)


dataframe <- read.csv("/g/data/te53/sj2852/blgrp/analyses/pangenome/analysis-v4/communities/centrality_measures_summary.csv", header = T)

betweenness <- ggplot(dataframe, aes(x = gene, y = mean_betweenness)) +
  geom_bar(stat = "identity") +
  scale_y_log10() + 
  theme_minimal() +
  xlab("Gene") +
  ylab("Mean Betweenness") +
  ggtitle("Mean Betweenness per Pangenome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

closeness <- ggplot(dataframe, aes(x = gene, y = mean_closeness)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Gene") +
  ylab("Mean Closeness") +
  ggtitle("Mean Closeness per Pangenome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

closeness + betweenness
