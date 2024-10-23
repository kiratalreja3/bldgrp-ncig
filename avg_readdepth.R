#Plot of Average Read Depth Per Blood Group Gene Per Sample

#this code creates the figure avg_readdepth_per_region (/g/data/te53/sj2852/blgrp/data/figures/avg_readdepth_per_region.png)

library(tidyverse)
library(ggplot2)


#load in read depth files
input_folder_chm13 <- "/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/bedcov_chm13"
input_folder_grch38 <- "/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/bedcov_grch38"
input_folder_NCIG <- "/g/data/te53/sj2852/blgrp/metadata/assembly_stats/NEW/bedcov_NCIG"

bed_files_chm13 <- list.files(path = input_folder_chm13, pattern = "\\.bed$", full.names = T)
bed_files_grch38 <- list.files(path = input_folder_grch38, pattern = "\\.bed$", full.names = T)
bed_files_NCIG <- list.files(path = input_folder_NCIG, pattern = "\\.bed$", full.names = T)


read_bed_file <- function(file_path) {
  bed_data <- read.table(file_path, header = F)
  colnames(bed_data) <- c("chr", "chr_str", "chr_end", "read_depth")
  file_name <- basename(file_path)
  sample_name <- str_split(file_name, "\\.", simplify = TRUE)[1]
  
  bed_data <- bed_data %>%
    mutate(file_name = sample_name)
  
  return(bed_data)
           
}

chm13 <- lapply(bed_files_chm13, read_bed_file) %>%
  bind_rows()

grch38 <- lapply(bed_files_grch38, read_bed_file) %>%
  bind_rows()

NCIG <- lapply(bed_files_NCIG, read_bed_file) %>%
  bind_rows()

#load in blood group gene files 

bloodcoords_chm13 <- read.table("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchorflank_named.bed", header = F)

bloodcoords_chm13 <- bloodcoords_chm13 %>%
  set_names(c("chr", "reg_str", "reg_end", "gene"))


bloodcoords_grch38 <- read.table("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/grch38/blgrp_coords_grch38_anchorflank.bed", header = F)

bloodcoords_grch38$gene <- bloodcoords_chm13[, 4]

bloodcoords_grch38 <- bloodcoords_grch38 %>%
  set_names(c("chr", "reg_str", "reg_end", "gene"))


chm13_joined <- chm13 %>%
  left_join(bloodcoords_chm13, by = c("chr" = "chr", "chr_str" = "reg_str", "chr_end" = "reg_end"))

grch38_joined <-grch38 %>%
  left_join(bloodcoords_grch38, by = c("chr" = "chr", "chr_str" = "reg_str", "chr_end" = "reg_end"))

NCIG_joined <- NCIG %>%
  left_join(bloodcoords_chm13, by = c("chr" = "chr", "chr_str" = "reg_str", "chr_end" = "reg_end"))


#merge the two files together

KGP_data <- bind_rows(chm13_joined, grch38_joined)


#Average Read Depth Per Blood Group Region and Per Sample

avr_readdepth_region_KGP <- KGP_data %>%
  group_by(gene) %>%
  summarise(avg_read_depth = mean(read_depth, na.rm = T))

avr_readdepth_region_NCIG <- NCIG_joined %>%
  group_by(gene) %>%
  summarise(avg_read_depth = mean(read_depth, na.rm = T))


avr_readdepth_KGP <- KGP_data %>%
  group_by(file_name) %>%
  summarise(avg_read_depth = mean(read_depth, na.rm = T))

avr_readdepth_NCIG <- NCIG_joined %>%
  group_by(file_name) %>%
  summarise(avg_read_depth = mean(read_depth, na.rm = T))

avr_readdepth_region_KGP$source <- "KGP"
avr_readdepth_region_NCIG$source <- "NCIG"

KGP_data$source <- "KGP"
NCIG_joined$source <- "NCIG"

avr_readdepth_KGP$source <- "KGP"
avr_readdepth_NCIG$source <- "NCIG"

  
combined_data <- bind_rows(avr_readdepth_region_KGP, avr_readdepth_region_NCIG)

combined_data2 <- bind_rows(avr_readdepth_KGP, avr_readdepth_NCIG)


combined_data_noaverage <- bind_rows(KGP_data, NCIG_joined)

combined_data_outlier <- combined_data_noaverage %>%
  group_by(gene, source) %>%
  mutate(outlier_flag = ifelse(read_depth %in% 
                                 boxplot.stats(read_depth)$out, 
                               "outlier", "normal"))


#Plot Average Read Depth Per Blood Group Region


#Boxplot with Outliers
ggplot(combined_data_outlier %>% dplyr::filter(gene != "REFERENCE"), aes(x = gene, y = read_depth, fill = source)) +
  geom_boxplot(outlier.colour = NA, outlier.shape = 16, outlier.size = 2) +
  geom_point(data = subset(combined_data_outlier, outlier_flag == "outlier" & gene != "REFERENCE"), 
             aes(x = gene, y = read_depth, color = source), 
             size = 1) +
  labs(title = "Read Coverage per Blood Group Region", 
       x = "Blood Group Region", 
       y = "Coverage") +
  theme_minimal()+
  theme(axis.text = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("skyblue2", "orange3"))+
  scale_color_manual(values = c("skyblue2", "orange3"))

#Calculate Average Per Region

combined_data <- combined_data %>%
  group_by(gene) %>%
  summarise(avg_read_depth = mean(avg_read_depth, na.rm = T))



