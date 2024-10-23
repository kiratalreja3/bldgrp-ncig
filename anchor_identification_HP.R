
#Anchor Identification - Hardip's Code

library(tidyverse)
library(GenomicRanges)

alnfiles <- list.files("/g/data/te53/sj2852/blgrp/tmp/full_coordinate_test/take_2/aln", ".tbl", full.names = T)

clean <- list()

for (i in 1:length(alnfiles)) {
  print(alnfiles[i])
  ##pad length
  pad <- 3e5
  
  ## read the alignment file
  aln <- read_delim(alnfiles[i], delim = "\t", col_names = F)
  colnames(aln) <- c("query", "target", "pid", "alnlen", "mm", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bitscore")
  ## get the genomic coordinates for the range
  rchr <- aln %>% pull(query) %>% unique() %>% str_split_i(":",1)
  rstart <- aln %>% pull(query) %>% unique() %>% str_split_i(":",2) %>% str_split_i("-",1) %>% as.numeric()
  rend <- aln %>% pull(query) %>% unique() %>% str_split_i(":",2) %>% str_split_i("-",2) %>% as.numeric()
  ## region length
  rlen <- rend - rstart + 1
  
  ## gene length, it may be off by few base pairs, adjust the following as required
  glen <- rlen - pad * 2
  
  ## gene coordinates within the region
  genecoords <- data.frame(seqid = "thisseq", start = pad, end = pad + glen)
  
  ## get regions that align to non-self within the region
  queryalncoords <- aln %>% 
    dplyr::filter(alnlen != rlen & alnlen > 150 & pid > 90) %>% 
    mutate(seqid = "thisseq") %>% 
    select(seqid, start = qstart, end = qend) 
  
  ## target coordinates need adjustment as start can be higher than end
  targetalncoords <- aln %>% 
    dplyr::filter(alnlen != rlen & alnlen > 150 & pid > 90) %>% 
    mutate(seqid = "thisseq") %>% select(seqid, tstart, tend) %>%
    mutate(start = case_when(tstart > tend ~ tend, TRUE ~ tstart), end = case_when(tend < tstart ~ tstart, TRUE ~ tend)) %>%
    select(seqid, start, end)
  
  ## collect all alignment coordinates as a single data frame
  alncoords <- bind_rows(queryalncoords,targetalncoords)
  alncoords <- distinct(alncoords)
  ## add in the gene coordinates as we don't want to include this in the flank
  alncoords <- bind_rows(genecoords, alncoords)
  
  ## make GRanges
  alnranges <- makeGRangesFromDataFrame(alncoords)
  ## collapse overlapping alignments and alignments within 200bp of each other as single range
  alnranges <- GenomicRanges::reduce(alnranges, min.gapwidth = 200)
  ## identify the range that overlaps with gene coordinates to ensure if there is overlapping alignments near ends of the gene
  generange <- makeGRangesFromDataFrame(genecoords)
  olap <- findOverlaps(generange, alnranges)
  regioncomplex <- alnranges[olap@to]
  
  ## get the nearest five/three prime alignments
  five <- data.frame(alnranges[follow(regioncomplex, alnranges)])
  three <- data.frame(alnranges[precede(regioncomplex, alnranges)])
  
  x <- data.frame(id = aln %>% pull(query) %>% unique(), 
                  rchr = rchr, rstart = rstart, rend = rend, 
                  fstart = rstart + five$end + 1, fend = rstart + data.frame(regioncomplex)$start - 1,
                  tstart = rstart + data.frame(regioncomplex)$end + 1, tend = rstart + three$start - 1
  )
  clean[[i]] <- x    
}

clean <- bind_rows(clean)

write_delim(clean, "/g/data/te53/sj2852/blgrp/tmp/rerunning_new_coords/NEW_COORDS/new_blgrp_coords", delim = "\t", col_names = T)


