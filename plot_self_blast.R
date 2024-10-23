#Plotting Self Blast Loop

library(tidyverse)
library(ggplot2)


#Function to apply to all

plot_self_blast_all <- function(directory_path, filter_identity = 95, filter_threshold = 100, title = "Self-BLAST Comparison") {
  
  files <- list.files(directory_path, full.names = TRUE)
  
  for (file_path in files) {
    
    input_table <- read.table(file_path, header = FALSE, sep = "\t")
    
    filtered_table <- input_table %>%
      dplyr::filter(V3 != filter_threshold, V3 > filter_identity)
    
    file_name <- basename(file_path)
    title <- gsub("_alignment\\.tbl$", "", file_name)
    
    # Plot the filtered table
    plot <- filtered_table %>%
      ggplot() +
      geom_segment(aes(x = V7, y = V9, xend = V8, yend = V10)) +
      labs(
        title = title,                 # Use the title parameter
        x = "Query Sequence Position",  # X-axis label
        y = "Subject Sequence Position" # Y-axis label
      )
    
    output_file <- file.path(directory_path, paste0(title, "_plot.png"))
    ggsave(output_file, plot, width =5, height = 4) 
  }
  
}

write.table(NCIG_data, "/g/data/te53/sj2852/blgrp/analyses/locityper_kgp/NEW/NCIG_data_tbl.txt", sep = "\t", row.names = F)

