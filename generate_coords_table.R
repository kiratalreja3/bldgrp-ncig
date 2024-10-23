
#Generate Table with all Coordinates for CHM13 and GRCH38

library(tidyverse)

#CHM13

anchorflank <- read.delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchorflank_named.bed", header = F)

gene_coords <- read.delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/OLD_COORDS/chm13/bloodcoords_noflank_withregionnoheader_chm13.bed", header = F)


full_df <- anchorflank %>%
  mutate(
    gene = V4, 
    chr = V1, 
    reg_str = V2, 
    `5_anchor_str` = V2 + 1000,
    `5_anchor_end` = V2 + 1200,
    `3_anchor_str` = V3 - 1200,
    `3_anchor_end` = V3 - 1000,
    reg_end = V3
  ) %>%
  select(gene, chr, reg_str, `5_anchor_str`, `5_anchor_end`, `3_anchor_str`, `3_anchor_end`, reg_end)
  )

full_df <- full_df %>%
  mutate(
    gene_str = gene_coords$V2,
    gene_end = gene_coords$V3
  )

full_df <- full_df %>%
  select(
    gene, chr, reg_str, `5_anchor_str`, `5_anchor_end`, 
    `gene_str`, `gene_end`, `3_anchor_str`, `3_anchor_end`, reg_end
  )


write_delim(full_df, "/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_allregionsref.bed", col_names = T)


anchor_only_df <- full_df %>%
  select(gene, chr, `5_anchor_str`, `5_anchor_end`, `3_anchor_str`, `3_anchor_end`)


write_delim(anchor_only_df, "/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/chm13/blgrp_coords_chm13_anchoronlyref.bed", col_names = T)

#GRCH38

anchorflank <- read.delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/grch38/blgrp_coords_grch38_anchorflank_named.bed", header = F)

gene_coords <- read.delim("/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/OLD_COORDS/grch38/bloodcoords_noflank_withregionnoheader_grch38.bed", header = F)


full_df <- anchorflank %>%
  mutate(
    gene = V4, 
    chr = V1, 
    reg_str = V2, 
    `5_anchor_str` = V2 + 1000,
    `5_anchor_end` = V2 + 1200,
    `3_anchor_str` = V3 - 1200,
    `3_anchor_end` = V3 - 1000,
    reg_end = V3
  ) %>%
  select(gene, chr, reg_str, `5_anchor_str`, `5_anchor_end`, `3_anchor_str`, `3_anchor_end`, reg_end)
)

full_df <- full_df %>%
  mutate(
    gene_str = gene_coords$V2,
    gene_end = gene_coords$V3
  )

full_df <- full_df %>%
  select(
    gene, chr, reg_str, `5_anchor_str`, `5_anchor_end`, 
    `gene_str`, `gene_end`, `3_anchor_str`, `3_anchor_end`, reg_end
  )


write_delim(full_df, "/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/grch38/blgrp_coords_cgrch38_allregionsref.bed", col_names = T)


anchor_only_df <- full_df %>%
  select(gene, chr, `5_anchor_str`, `5_anchor_end`, `3_anchor_str`, `3_anchor_end`)


write_delim(anchor_only_df, "/g/data/te53/sj2852/blgrp/metadata/blgrp_coords/NEW_COORDS/grch38/blgrp_coords_grch38_anchoronlyref.bed", col_names = T)

