#Author: Sarah Jackson

.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))
install.packages("ggplot")
install.packages("dyplr")

library(tidyverse)
library(ggplot2)

setwd("/g/data/te53/t2t2024/analyses/t2t_figure/QV_vs_mappingrate")

mappingqv <- read.table("/g/data/te53/t2t2024/analyses/t2t_figure/QV_vs_mappingrate/mappingqv.txt", header = T)

mappingrate.primary <- read.table("/g/data/te53/t2t2024/analyses/t2t_figure/QV_vs_mappingrate/mappingrate.primary.txt", header = T, sep = "\t")

#Process MappingQV

sliced_mapping_qv <- mappingqv %>%
  group_by(donor) %>%
  slice(c(1, 3)) %>%
  ungroup() %>%
  mutate(hap = str_replace(hap, "^h", "")) %>%
  select(-percent, -tech)
  
sliced_mapping_qv$hap <- as.integer(sliced_mapping_qv$hap)

#Process MappingRate.Primary

sliced_mappingrate.primary <- mappingrate.primary %>%
  group_by(Group, tech) %>%
  slice(c(1,2)) %>%
  ungroup() %>%
  rename(donor = Group) %>%
  select(-Types, -total_count) %>%
  filter(!donor %in% c("C8180020", "G5581483", "L1576200", "N003302")) %>% #these samples do not have Illumina data
  mutate(hap = as.integer(hap))


#Merge dataframes

merged <- left_join(sliced_mappingrate.primary, sliced_mapping_qv, by = c("donor", "hap"))

#Plot Data, x= average qv, y=average mapping

qv_vs_mapping <- ggplot()+
  geom_point(data = merged, aes(x = qv, y = mean_percentage, colour = tech), size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(52, 80), breaks = seq(50, 80, by = 2)) +
  ylim(94, 100) +
  theme(text = element_text(size=12), 
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size = 12), 
        axis.line = element_line(linewidth = 0.5))+
  labs(x = "Quality Value", y = "Mapping Rate (%)", colour = "Tech")


#Save
ggsave("/g/data/te53/t2t2024/analyses/t2t_figure/QV_vs_mappingrate/qv_vs_mapping.pdf", plot = qv_vs_mapping, width = 3, height = 4)
