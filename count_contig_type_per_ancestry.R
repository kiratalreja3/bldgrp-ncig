#Contig Type per Region

complete <- read_delim("/g/data/te53/sj2852/blgrp/tmp/communities.complete.txt", col_names =  F)

single <- read_delim("/g/data/te53/sj2852/blgrp/tmp/communities.single.txt", col_names = F)

complete <- complete %>%
  mutate(X3 = sub("#.*", "", X3))

ancestry <- read_delim("/g/data/te53/sj2852/blgrp/metadata/ancestry_files/merged_ancestry.txt", col_names = T)

complete <- complete %>%
  left_join(ancestry %>% select(V1, Description), by = c("X3" = "V1"))

count_complete <- complete %>%
  group_by(X1, Description) %>%
  summarise(count = n(), .groups = 'drop') 

count_complete <- count_complete %>%
  mutate(Description = ifelse(is.na(Description), "NCIG", Description))

count_NCIG <- count_complete %>%
  filter(Description == "NCIG") %>%  # Filter for NCIG values
  group_by(X1) %>%
  summarise(count = sum(count), .groups = 'drop') %>%
  mutate(Description = "NCIG") %>%
  mutate(normalized_count = count / 284) 

count_KGP <- count_complete %>%
  filter(Description != "NCIG") %>%
  group_by(X1) %>%
  summarise(count = sum(count), .groups = 'drop') %>%
  mutate(Description = "KGP") %>%
  mutate(normalized_count = count / 356)

# Combine the counts
final_count <- bind_rows(count_NCIG, count_KGP)

new_row <- data.frame(X1 = "PIEZO1", count = 0, Description = "NCIG", normalized_count = 0)

final_count <- rbind(final_count, new_row)



# Create the bar plot
complete_graph <- ggplot(final_count, aes(x = X1, y = count, fill = Description)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) + 
  labs(x = "Gene", y = "Count", fill = "Source", title = "Number of Contigs Containing Both 5' and 3' Anchors") +
  scale_fill_manual(values = c("NCIG" = "orange2", "KGP" = "skyblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Single Anchors

single <- single %>%
  mutate(X3 = sub("#.*", "", X3))

ancestry <- read_delim("/g/data/te53/sj2852/blgrp/metadata/ancestry_files/merged_ancestry.txt", col_names = T)

single <- single %>%
  left_join(ancestry %>% select(V1, Description), by = c("X3" = "V1"))

single <- single %>%
  mutate(Description = ifelse(is.na(Description), "NCIG", "KGP"))

single$X2 <- as.character(single$X2)

count_single <- single %>%
  group_by(X1, Description, X2) %>%
  summarise(count = n(), .groups = 'drop') 

summarized_data <- count_single%>%
  group_by(X1, X2, Description) %>%
  summarise(total_count = sum(count), .groups = 'drop')

single_contigs <- ggplot(count_single, aes(x = X1, y = count, fill = Description)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ X2, nrow = 1, labeller = as_labeller(c('3' = "3' Anchors", '5' = "5' Anchors"))) + 
  labs(x = "Gene", y = "Count", title = "Number of Contigs Containing only a Single Anchor per Sample Group", 
       fill = "Source") +
  theme_minimal() +
  theme(strip.text = element_text(size = 8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("KGP" = "skyblue", "NCIG" = "orange2")) 


complete_graph / single_contigs
