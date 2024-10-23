#this code is to identify the locations of duplication to determine the best position for unique anchors


library(tidyverse)
library(ggplot2)

#Determine Gene Location Relative to 100kb coordinates

hundredkb_coords <- read.table("/g/data/te53/sj2852/blgrp/analyses/anchors/full_script_test/full_region.bed", sep = "\t", header = F)

hundredkb_coords <- hundredkb_coords %>%
  dplyr::rename(
    chr = V1, 
    flank_start = V2, 
    flank_end = V3, 
    gene_region = V4
  )

gene_coords <- read_delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/chm13/bloodcoords_liftover_pre_mergeflank.bed", delim = "\t", col_names = T)

gene_coords <- gene_coords %>% select(-...5) %>% dplyr::slice(-50)

gene_coords <- gene_coords %>%
  mutate(chr = paste0("chr", chr)) %>%  
  select(chr, everything())

merged_df <- hundredkb_coords %>%
  inner_join(gene_coords, by = "chr", relationship=many-to-many) %>% 
  dplyr::filter(start > flank_start & end < flank_end) %>%  
  select(chr, flank_start, flank_end, gene_region, start, end, gene) 

merged_df[7, "start"] <- 206761475
merged_df[7, "end"] <- 206888507

merged_df[7, "flank_start"] <- 206661475
merged_df[7, "flank_end"] <- 206908262


merged_df <- merged_df %>%
  mutate(region_size = flank_end - flank_start)

merged_df <- merged_df %>%
  mutate(
    scaling_factor = region_size / (flank_end - flank_start), 
    
    scaled_start = 1 + (scaling_factor * (start - flank_start)), 
    scaled_end = 1 + (scaling_factor * (end - flank_start))
  )

merged_df <- merged_df %>%
  mutate(coordinate = paste(chr, flank_start, flank_end, sep = "_")) 

#Function to apply to all

plot_self_blast_all <- function(directory_path, filter_identity = 95, filter_threshold = 100, title = "Self-BLAST Comparison") {
  
  files <- list.files(directory_path, full.names = TRUE)
  
  for (file_path in files) {
  
  input_table <- read.table(file_path, header = FALSE, sep = "\t")
  
  filtered_table <- input_table %>%
    filter(V3 != filter_threshold, V3 > filter_identity)
  
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

input_table <- read.table("/g/data/te53/sj2852/blgrp/analyses/anchors/full_script_test/blast_tables/chr1_25003561-25371735_alignment.tbl", header = FALSE, sep = "\t")

filtered_table <- input_table %>%
  dplyr::filter(V3 != 100, V3 > 95)

filtered_table %>%
  ggplot() +
  geom_segment(aes(x = V7, y = V9, xend = V8, yend = V10)) +
  
  {if (nrow(filtered_table) > 0) {
    geom_segment(aes(x = 100000, y = 100000, 
                     xend = merged_df$scaled_start, 
                     yend = merged_df$scaled_end, 
                     colour = "Custom Segment"), linewidth = 1)
  } else {
    NULL
  }
  }+ 
  scale_colour_manual(
    name = "Segments", 
    values = c("Custom Segment" = "red")
  ) +
  
  labs(
    title = "RHD_RHCE_Self_Blast",                 # Use the title parameter
    x = "Query Sequence Position",  # X-axis label
    y = "Subject Sequence Position" # Y-axis label
  )}

  

#Individual

#RHD/RHCE
#blast alignment of 1MB region to itself in VS code

RHD_RHCE_table <- read.table("/g/data/te53/sj2852/blgrp/analyses/anchors/RHD_RHCE_anchor_identifying/rhd.vs.rhd.100kb.tbl", header = F, sep = "\t")

RHD_RHCE_table <- RHD_RHCE_table %>% 
  dplyr::filter(V3!=100 & V3>75 & V4>100) 


RHD_RHCE_table %>% 
  ggplot() + 
  geom_segment(aes(x=V7,y=V9,xend=V8,yend=V10)) +
  geom_segment(aes(x = 100000, y = 100000, xend = 164978, yend = 164978, 
                   color = "RHD Gene"), linewidth = 1) +
  geom_segment(aes(x = 192786, y = 192786, xend = 267774, yend = 267774, 
                   color = "RHCE Gene"), linewidth = 1) +
  geom_segment(aes(x = 98260, y = 98449, xend = 46031, yend = 45831, 
                   colour = "Custom"), linewidth = 1)+
  scale_color_manual(
    name = "Segments", 
    values = c("RHD Gene" = "red", "RHCE Gene" = "red3", "Custom" = "blue")
  ) +
  
  labs(
    title = "Self-BLAST Comparison: RHD_RHCE Region", 
    x = "Query Sequence Position", 
    y = "Subject Sequence Position"
  )
  





RHD_RHCE_filtered <- RHD_RHCE_table %>% 
  dplyr::filter(V3 >=90)
  # for identifying two big alignments 5' of the RHD gene start


#CR1

CR1_table <- read.table("/g/data/te53/sj2852/blgrp/analyses/anchors/CR1_anchor_identifying/cr1.vs.cr1.tbl", header = F, sep = "\t")

CR1_table <- CR1_table %>%
  filter(V3!=100 & V4>100 & V3>75)

CR1_table %>%
  ggplot() +
  geom_segment(aes(x=V7, y=V9, xend=V8, yend=V10)) +
  geom_segment(aes(x = 100000, y = 100000 , xend = 227032, yend = 227032, 
                   color = "CR1 Gene"), linewidth = 1) +
  scale_color_manual(
    name = "Segments", 
    values = c("CR1 Gene" = "red")) +
  labs(
    title = "Self-BLAST Comparison: CR1 Region", 
    x = "Query Sequence Position", 
    y = "Subject Sequence Position")


+ xlim(c(95000,110000))


CR1_filtered <- CR1_table %>% 
  filter(V3!=100 & V7<206761475 & V3>98 & V4 >200)


#GYPA_GYPB

GYPA_GYPB_table_new <- read.table("/g/data/te53/sj2852/blgrp/analyses/anchors/GYPA_GYPB_anchor_identifying/300kb_5_100kb_3/gypa.vs.gypb.300_100kb.tbl", header = F, sep = "\t")

GYPA_GYPB_table_new <- GYPA_GYPB_table_new %>%
  filter(V3!=100 & V4>100 & V3>75 & V8>452514)

GYPA_GYPB_table_new %>%
  ggplot() +
  geom_segment(aes(x=V7, y=V9, xend=V8, yend=V10)) +
  geom_segment(aes(x = 414097, y = 414097, xend = 452514, yend = 452514, 
                   color = "GYPA Gene"), linewidth = 1) +
  geom_segment(aes(x = 300000, y = 300000, xend = 333669, yend = 333669, 
                   color = "GYPB Gene"), linewidth = 1) +
  scale_color_manual(
    name = "Segments", 
    values = c("GYPB Gene" = "red3", "GYPA Gene" = "red"))+
  labs(
    title = "Self-BLAST Comparison: GYPA_GYPB Region", 
    x = "Query Sequence Position", 
    y = "Subject Sequence Position"
  )


GYPA_GYPB_filtered_5 <- GYPA_GYPB_table %>% 
  filter(V3!=100 & V7<147310608 & V3>97 & V4 >200)


#PIEZO1

PIEZO1_table <- read.table("/g/data/te53/sj2852/blgrp/tmp/PIEZO1_test/self_blast.tbl", header = F, sep = "\t")

PIEZO1_filtered <- PIEZO1_table %>% 
  dplyr::filter(V3!=100 & V3>90)


title <- "PIEZO"

PIEZO1_filtered %>%
  ggplot() +
  geom_segment(aes(x = V7, y = V9, xend = V8, yend = V10)) +
  annotate("segment", x = 300000, y = 300000 , xend = 380742, yend = 380742, 
                   color = "red", linewidth = 1) +
  labs(
    title = title,                 # Use the title parameter
    x = "Query Sequence Position",  # X-axis label
    y = "Subject Sequence Position" # Y-axis label
  )
