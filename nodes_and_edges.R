#Counting Edges and Nodes

install.packages("yaml")

library(tidyverse)
library(yaml)
library(ggplot2)

yaml_dir <- "/g/data/te53/sj2852/blgrp/analyses/pangenome/analysis-v4/statistics/nodes_edges/stats_files"

yaml_files <- list.files(path = yaml_dir, pattern = "*.yaml", full.names = T)

data_list <- list()

for (file in yaml_files) {
 
  yaml_data <- read_yaml(file)
  
  
  extracted_data <- list(
    length = yaml_data$length,
    nodes = yaml_data$nodes,
    edges = yaml_data$edges,
    paths = yaml_data$paths,
    steps = yaml_data$steps
  )
  
  
  extracted_data$gene <- sub("\\..*", "", basename(file))
  
  
  data_list[[file]] <- as.data.frame(t(unlist(extracted_data)))
}

final_df <- bind_rows(data_list)

final_df$nodes <- as.numeric(final_df$nodes)
final_df$edges <- as.numeric(final_df$edges)
final_df$length <- as.numeric(final_df$length)
final_df$paths <- as.numeric(final_df$paths)

#Gene Lengths

blgrp_regions <- read_delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_allregionsref.bed", col_names = T)

blgrp_regions <- blgrp_regions %>%
  mutate(region_length = reg_end - reg_str)

final_df <- final_df %>%
  left_join(blgrp_regions %>% select(gene, region_length), by = "gene")

final_df$normalized_nodes <- final_df$nodes / final_df$region_length
final_df$normalized_nodes1 <- final_df$nodes / final_df$length

final_df$normalised_edges <- final_df$edges / final_df$region_length
final_df$normalised_edges1 <- final_df$edges / final_df$length

final_df$node_to_edge <- final_df$edges / final_df$nodes

final_df$path_to_node <- final_df$paths / final_df$nodes

final_df$node_density <- final_df$nodes / final_df$length

#Plot Number of Nodes Per Sample

ggplot(final_df, aes(x = gene, y = normalized_nodes1)) + 
  geom_bar(stat = "identity", fill = "darkgrey") + 
  labs(x = "Gene", y = "Number of Nodes (Normalised)", title = "Number of Nodes per Pangenome") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


#Plot Number of Edges Per Sample

ggplot(final_df, aes(x = gene, y = normalised_edges1)) + 
  geom_bar(stat = "identity", fill = "darkgrey") + 
  labs(x = "Gene", y = "Number of Edges (Normalised)", title = "Number of Edges per Pangenome") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


#Edge to Node Ratio

ggplot(final_df, aes(x = gene, y = node_to_edge)) + 
  geom_bar(stat = "identity", fill = "darkgrey") + 
  labs(x = "Gene", y = "Edge to Node Ratio", title = "Edge to Node Ratio per Pangenome") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


#Path to Node Ratio

ggplot(final_df, aes(x = gene, y = path_to_node)) + 
  geom_bar(stat = "identity", fill = "darkgrey") + 
  labs(x = "Gene", y = "Path to Node Ratio", title = "Path to Node Ratio per Pangenome") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
