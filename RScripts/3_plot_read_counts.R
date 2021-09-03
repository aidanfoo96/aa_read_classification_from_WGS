# WRANGLE SPECIES ABUNDANCE DATA
library(tidyverse)

## Plot PRJNA385 METAPHLAN #####
setwd("../")
ID = c("SRR5519646", "SRR5562853", 
       "SRR5562854", "SRR5562855", "SRR5562856",
       "SRR5562858", "SRR5562859", "SRR5562860",
       "SRR5562861",
       "SRR5562862", "SRR5562863", "SRR5562864", 
       "SRR5562865", "SRR5562866", "SRR5562867", 
       "SRR5562868", "SRR5562869", "SRR5562870", 
      "SRR5562871", "SRR5562872", 
       "SRR5562873", "SRR5562875", "SRR5562876", 
       "SRR5562877", "SRR5562878",  "SRR5562878")
columns <- c("Spp", "Relative_Abundance", "Coverage", "Est_Num_Reads")
 
filesnames <- list.files(pattern = ".txt")
filelist <- lapply(filesnames, read.table)
filelist2 <- lapply(filelist, setNames, columns)
filelist3 <- mapply('[<-', filelist2, "SampleID", value = ID, SIMPLIFY = F)

PRJNA385349_metaphl_abun <- as_tibble(do.call("rbind", filelist3))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_per_sample_reads_metphln <- function(filtered_read_table, filter_thresh){
  filtered_read_table2 <- filtered_read_table %>%
    group_by(Spp, SampleID) %>%
    summarise(total_read_count = sum(Est_Num_Reads)) %>%
    filter(total_read_count > filter_thresh) %>%
    filter(Spp != "UNKNOWN")
  
  total_reads <- filtered_read_table2 %>%
    group_by(SampleID) %>%
    summarise(total_reads = sum(total_read_count))
  
  normalised_join <- filtered_read_table2 %>%
    left_join(x = total_reads, by = "SampleID") %>%
    mutate(normalised_read_count = total_read_count/total_reads *100)
  
  normalised_join_plot <- normalised_join %>%
    ggplot() +
    aes(x = SampleID, y = Spp) + 
    geom_tile(aes(fill = normalised_read_count), colour = "black") +
    scale_fill_gradientn(colours = c("yellow", "orange", "blue"), 
                         values = scales::rescale(c(0, 5, 10, 20, 40, 60, 80))) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

    return(normalised_join_plot) 
  
}

PRJNA385349_metaphl_abun_per_samp <- plot_per_sample_reads_metphln(PRJNA385349_metaphl_abun, 0)

pdf("PRJNA385349_metaphl_abun_per_samp.pdf")
print(PRJNA385349_metaphl_abun_per_samp)
dev.off()

PRJNA385349_metaphl_ov_plot <- PRJNA385349_metaphl_abun %>%
  filter(Spp != "UNKNOWN") %>%
  filter(Est_Num_Reads > 500) %>%
  group_by(Spp) %>%
  summarise(total_read_count = sum(Est_Num_Reads)) %>%
  ggplot() +
  aes(x = reorder(Spp, total_read_count), y = total_read_count) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_minimal()

pdf("PRJNA385349_metaphl_ov_plot.pdf")
print(PRJNA385349_metaphl_ov_plot)
dev.off()





## Plot Summaries
grouped_ab_plot <- ggplot(grouped_ab_tbl) +
  aes(x = reorder(Spp, total_read_count), y = total_read_count) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  theme_classic()

pdf("grouped_ab_plot.png")
print(grouped_ab_plot)
dev.off()

grouped_ab_ng <- ggplot(grouped_ab_tbl_ng) +
  aes(x = reorder(Spp, total_read_count), y = total_read_count) +
  geom_bar(stat = "identity") +
  coord_flip() 


## Plot PRJNA419379 METAPHLAN #####
ID = c("SRR6313629", "SRR6313630", "SRR6313631", "SRR6313632")
setwd(dir = "PRJNA419379_metaphlan_abun_stats/")
columns_419 <- c("Spp", "Relative_Abundance", "Coverage", "Est_Num_Reads")
filesnames <- list.files(pattern = ".txt")
filelist <- lapply(filesnames, read.table)
filelist2 <- lapply(filelist, setNames, columns)
filelist3 <- mapply('[<-', filelist2, "SampleID", value = ID, SIMPLIFY = F)

PRJNA419379_metaphl_abun <- as_tibble(do.call("rbind", filelist3))


plot_per_sample_reads_metphln <- function(filtered_read_table, filter_thresh){
  filtered_read_table2 <- filtered_read_table %>%
    group_by(Spp, SampleID) %>%
    summarise(total_read_count = sum(Est_Num_Reads)) %>%
    filter(total_read_count > filter_thresh) %>%
    filter(Spp != "UNKNOWN")
  
  total_reads <- filtered_read_table2 %>%
    group_by(SampleID) %>%
    summarise(total_reads = sum(total_read_count))
  
  normalised_join <- filtered_read_table2 %>%
    left_join(x = total_reads, by = "SampleID") %>%
    mutate(normalised_read_count = total_read_count/total_reads *100)
  
  normalised_join_plot <- normalised_join %>%
    ggplot() +
    aes(x = SampleID, y = Spp) + 
    geom_tile(aes(fill = normalised_read_count), colour = "black") +
    scale_fill_gradientn(colours = c("yellow", "orange", "blue"), 
                         values = scales::rescale(c(0, 5, 10, 20, 40, 60, 80))) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  
  return(normalised_join_plot) 
  
}

PRJNA419379_metaphl_abun_per_samp <- plot_per_sample_reads_metphln(PRJNA419379_metaphl_abun, 0)

pdf("PRJNA419379_metaphl_abun_per_samp.pdf")
print(PRJNA419379_metaphl_abun_per_samp)
dev.off()

PRJNA419379_metaphl_ov_plot <- PRJNA419379_metaphl_abun %>%
  filter(Spp != "UNKNOWN") %>%
  filter(Spp != "unclassified") %>%
  group_by(Spp) %>%
  summarise(total_read_count = sum(Est_Num_Reads)) %>%
  ggplot() +
  aes(x = reorder(Spp, total_read_count), y = total_read_count) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_minimal()

pdf("PRJNA419379_metaphl_ov_plot.pdf")
print(PRJNA419379_metaphl_ov_plot)
dev.off()



## Plot PRJNA255893 KRAKEN #####
setwd(dir = "../PRJNA255_metaphln/")
ID = c("SRR1523064", "SRR1523066", "SRR1523067", "SRR1523068")
columns <- c("Spp", "Relative_Abundance", "Coverage", "Est_Num_Reads")
filesnames <- list.files(pattern = ".txt")
filelist <- lapply(filesnames, read.table)
filelist2 <- lapply(filelist, setNames, columns)
filelist3 <- mapply('[<-', filelist2, "SampleID", value = ID, SIMPLIFY = F)

PRJNA255893_metaphl_abun <- as_tibble(do.call("rbind", filelist3))

PRJNA255893_metaphl_abun_per_samp <- plot_per_sample_reads_metphln(PRJNA255893_metaphl_abun, 0)

pdf("PRJNA255893_metaphl_abun_per_samp.pdf")
print(PRJNA255893_metaphl_abun_per_samp)
dev.off()

PRJNA255893_metaphl_abun_ov_plot <- PRJNA255893_metaphl_abun %>%
  filter(Spp != "UNKNOWN") %>%
  filter(Spp != "unclassified") %>%
  group_by(Spp) %>%
  summarise(total_read_count = sum(Est_Num_Reads)) %>%
  ggplot() +
  aes(x = reorder(Spp, total_read_count), y = total_read_count) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_minimal()

pdf("PRJNA255893_metaphl_abun_ov_plot.pdf")
print(PRJNA255893_metaphl_abun_ov_plot)
dev.off()

PRJNA385349_metaphl_abun
PRJNA419379_metaphl_abun
PRJNA255893_metaphl_abun

metaphln_abun_tbl <- bind_rows(PRJNA385349_metaphl_abun,
                               PRJNA419379_metaphl_abun, 
                               PRJNA255893_metaphl_abun)

metaphln_abun_tbl_unknown <- metaphln_abun_tbl %>%
  filter(Spp == "UNKNOWN") %>%
  group_by(SampleID) %>%
  summarise(unknown_reads = sqrt(sum(Est_Num_Reads)))

metaphln_abun_tbl_known <- metaphln_abun_tbl %>%
  dplyr::filter(Spp != "UNKNOWN") %>%
  group_by(SampleID) %>%
  dplyr::summarise(known_reads = sqrt(sum(Est_Num_Reads)))

metaphln_abun_tbl_joined <- metaphln_abun_tbl_known %>%
  left_join(metaphln_abun_tbl_unknown, by = "SampleID") %>%
  pivot_longer(cols = c("known_reads", "unknown_reads"), names_to = "read_type", values_to = "no_reads") %>%
  ggplot() +
  aes(x = SampleID, y = no_reads, fill = read_type) + 
  geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() + 
  theme_minimal()

pdf("metaphln_abun_tbl_joined.pdf")
print(metaphln_abun_tbl_joined)
dev.off()

PRJNA255893_summaries <- PRJNA255893_metaphl_abun %>%
  group_by(Spp) %>%
  summarise(sum(Est_Num_Reads))

PRJNA419379_summaries <- PRJNA419379_metaphl_abun %>%
  group_by(Spp) %>%
  summarise(sum(Est_Num_Reads))

PRJNA385349_summaries <- PRJNA385349_metaphl_abun %>%
  group_by(Spp) %>%
  summarise(sum_reads = sum(Est_Num_Reads)) %>%
  filter(sum_reads > 2000) %>%
  arrange(sum_reads)
