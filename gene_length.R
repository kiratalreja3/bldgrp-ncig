#Plotting Gene Lengths



library(karyoploteR)
library(GenomicRanges)

#this code generates the figure called gene_length_barchart (/g/data/te53/sj2852/blgrp/data/figures/gene_length_barchart.png)


bloodcoords <- read.table("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchorflank_withname.bed", header = F)

name_data <- read.delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/OLD_COORDS/chm13/bloodcoords_flank_withregion_chm13.bed", header = F, sep = "\t")



column_4 <- name_data$V4



bloodcoords <- bloodcoords %>%
  mutate(chr = V1, flank_str = V2, flank_end = V3, gene = column_4) %>%
  mutate(
    f_anchor_str = flank_str + 1000, 
    f_anchor_end = flank_str + 1200, 
    t_anchor_str = flank_end - 1200, 
    t_anchor_end = flank_end - 1000
  )

bloodcoords <- subset(bloodcoords, gene != "REFERENCE")

region_data <- read.delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_region_withname.bed", header = F, sep = "\t")

region_str <- region_data$V2

region_end <- region_data$V3

bloodcoords <- bloodcoords %>%
  mutate(reg_str = region_str, reg_end = region_end) %>%
  mutate(overall_length = flank_end - flank_str,
         region_length = reg_end - reg_str)



transform_coord <- function(x, flank_str, flank_end, overall_length) {
  scale <- (x - flank_str) / (flank_end - flank_str) * overall_length
  scale[flank_end == flank_str] <- 0  # Handle cases where start == end
  return(scale)
}

output_dataframe <- bloodcoords %>%
  mutate(
    f_anchor_str_trans = transform_coord(f_anchor_str, flank_str, flank_end, overall_length),
    f_anchor_end_trans = transform_coord(f_anchor_end, flank_str, flank_end, overall_length),
    t_anchor_str_trans = transform_coord(t_anchor_str, flank_str, flank_end, overall_length),
    t_anchor_end_trans = transform_coord(t_anchor_end, flank_str, flank_end, overall_length),
    reg_str_trans = transform_coord(reg_str, flank_str, flank_end, overall_length),
    reg_end_trans = transform_coord(reg_end, flank_str, flank_end, overall_length)
  ) %>%
  select(gene, f_anchor_str_trans, f_anchor_end_trans, t_anchor_str_trans, t_anchor_end_trans, reg_str_trans, reg_end_trans)



#Testing KaryoploteR
f_anchor_ranges <- GRanges(
  seqnames = output_dataframe$gene, 
  ranges = IRanges(
    start = output_dataframe$f_anchor_str_trans,
    end = output_dataframe$f_anchor_end_trans
  )
)


t_anchor_ranges <- GRanges(
  seqnames = output_dataframe$gene,  
  ranges = IRanges(
    start = output_dataframe$t_anchor_str_trans,
    end = output_dataframe$t_anchor_end_trans
  )
)


reg_ranges <- GRanges(
  seqnames = output_dataframe$gene,  
  ranges = IRanges(
    start = output_dataframe$reg_str_trans,
    end = output_dataframe$reg_end_trans
  )
)



custom_genome <- bloodcoords %>%
  select(gene, overall_length) %>%
  dplyr::rename(chr = gene, length = overall_length)


gr_custom_genome <- GRanges(
  seqnames = custom_genome$chr,
  ranges = IRanges(start = 1, end = custom_genome$length)
)



kp <- plotKaryotype(genome = custom_genome)


kpPlotRegions(kp, f_anchor_ranges, col = "red", border = NA, lwd = 4, r0 = -0.1, r1 = -0.35)

kpPlotRegions(kp, t_anchor_ranges, col = "red", border = NA, lwd = 4, r0 = -0.1, r1 = -0.35)

kpPlotRegions(kp, reg_ranges, col = "black", border = NA, lwd = 4, r0 = -0.1, r1 = -0.35)



#Plotting Gene Length


df_combined <- output_dataframe %>%
  left_join(bloodcoords %>% select(gene, overall_length, region_length), by = "gene")


df_combined <- df_combined %>%
  mutate(
    gap_5_2_region = reg_str_trans - f_anchor_end_trans,
    gap_3_2_region = t_anchor_str_trans - reg_end_trans,
    end_to_t = overall_length - t_anchor_end_trans, 
    `5_anchor` = 600,  # New column with constant value
    `3_anchor` = 600   # New column with constant value
  )

df_long <- df_combined %>%
  pivot_longer(cols = c(f_anchor_str_trans, `5_anchor`, gap_5_2_region, 
                        region_length, gap_3_2_region, `3_anchor`, end_to_t), 
               names_to = "category", 
               values_to = "value") %>%
  mutate(category = factor(category, levels = c("f_anchor_str_trans", 
                                                "5_anchor", 
                                                "gap_5_2_region", 
                                                "region_length", 
                                                "gap_3_2_region", 
                                                "3_anchor", 
                                                "end_to_t")))  # Set factor levels for stacking order

ggplot(df_long, aes(x = gene, y = value, fill = category)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene", y = "Overall Length", title = "Stacked Bar Chart of Overall Length by Gene") +
  scale_fill_manual(values = c(
    "f_anchor_str_trans" = "grey", 
    "5_anchor" = "red", 
    "gap_5_2_region" = "grey", 
    "region_length" = "black", 
    "gap_3_2_region" = "grey", 
    "3_anchor" = "red", 
    "end_to_t" = "grey"
  )) +  # Customize colors if needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")













ggplot(bloodcoords, aes(x = gene, y = overall_length)) +
  geom_bar(stat = "identity")+
  labs(title = "Blood Group Region Length",
       x = "Gene Region",
       y = "Length (bp)")+
  theme_minimal()+
  theme(axis.text = element_text(angle = 45, hjust = 1))
  




ggplot(bloodcoords_long, aes(x = region, y = length, fill = segment))+
  geom_bar(stat = "identity")+
  labs(title = "Blood Group Region Length",
       x = "Region",
       y = "Length (bp)") +
  scale_fill_manual(values = c("gene_length" = "skyblue", 
                               "str_flank" = "lightcoral", 
                               "end_flank" = "lightcoral"), 
                    labels = c("5' Flanking Region", "Gene Length", "3' Flanking Region")) +
  theme_minimal()+
  theme(axis.text = element_text(angle = 45, hjust = 1))
