library(tidyverse)
library(data.table)
library(plotly)

## Plot Species level relative abundance from MetaPhlAn Result ##### 
merged_abundance_table <- read_table2("MetaPhAn/merged_abundnace_table_species.txt")
abundance_tbl_PRJNA255893 <- read_table2("MetaPhAn/final_abundance_table_PRJNA255893.txt")
abundance_tble_PRJNA385349 <- read_table2("MetaPhAn/PRJNA385349_species_table.txt")

merged_abundance_table_clean <- merged_abundance_table %>%
  pivot_longer(cols = 2:7, 
               names_to = "sample", 
               values_to = "relative_abundance")

abundance_tbl_PRJNA385349_clean <- abundance_tble_PRJNA385349 %>%
  select(!NCBI_tax_id) %>%
  pivot_longer(cols = 2:14, 
               names_to = "sample", 
               values_to = "relative_abundance")

abundance_tbl_PRJNA255893_clean <- abundance_tbl_PRJNA255893 %>%
  pivot_longer(cols = 2, 
               names_to = "sample", 
               values_to = "relative_abundance")



abundance_table <- rbind(merged_abundance_table_clean, abundance_tbl_PRJNA255893_clean)

#### Stacked Bar Chart
plot <- ggplot(abundance_tbl_PRJNA255893_clean) + 
  aes(x = sample, fill = clade_name, y = relative_abundance) + 
  geom_bar(position = "fill", stat = "identity") + 
  labs(x = "Sample", y = "Relative Abundance", fill = "Species") + 
  theme_classic() +
  coord_flip()

pdf("spp_abundance_tbl_PRJNA255893.pdf")  
print(plot)
dev.off()

#### Bubble Plot
plot2 <- ggplot(abundance_tbl_PRJNA385349_clean, aes(x = sample, y = clade_name)) + 
  geom_tile(aes(fill = relative_abundance), colour = "white") +
  scale_fill_gradient(low = "grey", high = "green") + 
  theme_classic()

pdf("spp_abundance_hm.pdf")
print(plot2)
dev.off()


