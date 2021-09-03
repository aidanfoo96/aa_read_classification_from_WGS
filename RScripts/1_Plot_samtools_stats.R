library(tidyverse)
library(ggthemes)
#import CSVs  
setwd("../..")
align_stats <- read_csv("samtools_stats_alignment_statistics.csv")

align_stats_clean <- align_stats %>%
  pivot_longer(cols = c(total_seq, reads_mapped, reads_unmapped), 
               names_to = "assembly_type", 
               values_to = "total_read_number") %>%
  select(PRJN, SRR, assembly_type, total_read_number)

# Plot
tol_col_pal <- c("#BBCCEE", "#CCEEFF", "#CCDDAA", "#EEEEBB", "#FFCCCC", "#DDDDDD")

align_stats_clean_plot <- ggplot(align_stats_clean) + 
  aes(x = SRR, fill = assembly_type, y = total_read_number) + 
  geom_bar(position = "dodge", stat = "identity") + 
  labs(x = "Sample", y = "Read Count", fill  = "Alignment Type") + 
  theme_minimal() + 
  coord_flip() +
  scale_fill_manual(values = tol_col_pal)

pdf("samtools_statistics2.pdf")
print(align_stats_clean_plot)
dev.off()

fig <- ggplotly(align_stats_clean_plot)
