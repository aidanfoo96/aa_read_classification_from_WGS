## Load packages
library(tidyverse)
columns = c("Spp", "read_count")
## Plot multiple reads from KRAKEN2 (PRJNA385349) #####
##
ID = c("SRR5519646", "SRR5562853", 
       "SRR5562854", "SRR5562855", "SRR5562861", 
       "SRR5562862", "SRR5562863", "SRR5562864", 
       "SRR5562865", "SRR5562866", "SRR5562867", 
       "SRR5562868", "SRR5562869",  
       "SRR5562871", "SRR5562872", "SRR5562873", 
       "SRR5562875", "SRR5562876", "SRR5562877", 
       "SRR5562878", "SRR5562879")

setwd("PRJNA385349_KRAKN_abundance/")
filesnames <- list.files(pattern = ".txt")
filelist <- lapply(filesnames, read_tsv)
filelist2 <- lapply(filelist, setNames, columns)
filelist3 <- mapply('[<-', filelist2, "SampleID", value = ID, SIMPLIFY = F)

PRJNA385349_kraken_abun <- bind_rows(filelist3)

clean_output <- function(kraken_table){
  kraken_table_clean <- kraken_table %>%
    separate(Spp, 
             c("d", "domain", "k", "kingdom", "p", "phylum", 
               "c", "class", "o", "order", "f", "familiy", "g", "genus",
               "species")) %>%
    select(domain, kingdom, phylum, class, order, familiy, genus, species, read_count, SampleID) %>%
    unite("genus_spp", genus:species)
  return(kraken_table_clean)
}

PRJNA385349_kraken_abun_clean <- clean_output(PRJNA385349_kraken_abun)

#### Plot bacterial read abundance against different samples
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

PRJNA385349_kraken_abun_clean %>%
  group_by(SampleID) %>%
  summarise(reads = sum(read_count))

PRJNA385349_kraken_prop_plot <- PRJNA385349_kraken_abun_clean %>%
  filter(domain =="Bacteria") %>%
  filter(read_count > 2000) %>%
  filter(SampleID != "SRR551970") %>%
  filter(genus_spp != "NA_NA") %>%
  filter(genus_spp != "symbiont_NA") %>%
  filter(genus_spp != "viridis_NA") %>%
  filter(genus_spp != "NIES_3974") %>%
  ggplot() + 
  aes(x = reorder(SampleID, read_count), y = read_count, fill = genus_spp) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_hue(l = 70, c = 60) +
  theme_classic() + 
  coord_flip()

pdf("PRJNA385349_kraken_prop_plot.pdf")
print(PRJNA385349_kraken_prop_plot)
dev.off()

plot_per_sample_reads <- function(filtered_read_table, filter_thresh){
  filtered_read_table2 <- filtered_read_table %>%
    filter(domain == "Bacteria") %>%
    filter(genus_spp != "NIES_3974") %>%
    filter(SampleID != "SRR5562870") %>%
    filter(genus_spp != "symbiont_NA") %>%
    filter(genus_spp != "viridis_NA") %>%
    filter(genus_spp != "NA_NA") %>%
    group_by(genus_spp, SampleID) %>%
    summarise(total_read_count = sum(read_count)) %>%
    filter(total_read_count > filter_thresh)
  
  total_reads <- filtered_read_table2 %>%
    group_by(SampleID) %>%
    summarise(total_reads = sum(total_read_count))
  
  normalised_join <- filtered_read_table2 %>%
    left_join(x = total_reads, by = "SampleID") %>%
    mutate(normalised_read_count = total_read_count/total_reads *100)
  
  normalised_join_plot <- normalised_join %>%
    ggplot() +
    aes(x = SampleID, y = genus_spp) + 
    geom_tile(aes(fill = normalised_read_count), colour = "black") +
    scale_fill_gradientn(colours = c("yellow", "orange", "blue"), 
                         values = scales::rescale(c(0, 5, 10, 20, 40, 60, 80))) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  return(normalised_join_plot) 
  
}

PRJNA385349_kraken_abun_per_samp <- plot_per_sample_reads(PRJNA385349_kraken_abun_clean, 2000)

pdf("PRJNA385349_kraken_abun_per_samp.pdf")
print(PRJNA385349_kraken_abun_per_samp)
dev.off()

#### Plot overall read abundance accross samples
PRJNA385349_kraken_ov_plot <- PRJNA385349_kraken_abun_clean %>%
  filter(domain == "Bacteria") %>%
  filter(genus_spp != "NA_NA") %>%
  filter(genus_spp != "NIES_3974") %>%
  filter(genus_spp != "symbiont_NA") %>%
  filter(genus_spp != "viridis_NA") %>%
  filter(genus_spp != "Sedis_g") %>%
  group_by(genus_spp) %>%
  summarise(total_read_count = sum(read_count)) %>%
  filter(total_read_count > 20000) %>%
  ggplot() +
  aes(x = reorder(genus_spp, total_read_count), y = total_read_count) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_minimal()
  
pdf("PRJNA385349_kraken_ov_plot.pdf")
print(PRJNA385349_kraken_ov_plot)
dev.off()

## Plot multiple reads from KRAKEN2 (PRJNA419) #####
ID = c("SRR6313629", "SRR6313630", "SRR6313631", "SRR6313632")

setwd("../PRJNA419379_metaphlan_abun_stats_krak/ver2/")
filesnames <- list.files(pattern = ".txt")
filelist <- lapply(filesnames, read_tsv)
filelist2 <- lapply(filelist, setNames, columns)
filelist3 <- mapply('[<-', filelist2, "SampleID", value = ID, SIMPLIFY = F)

PRJNA419379_kraken_abun <- bind_rows(filelist3)

clean_output <- function(kraken_table){
  kraken_table_clean <- kraken_table %>%
    separate(Spp, 
             c("d", "domain", "k", "kingdom", "p", "phylum", 
               "c", "class", "o", "order", "f", "familiy", "g", "genus",
               "species")) %>%
    select(domain, kingdom, phylum, class, order, familiy, genus, species, read_count, SampleID) %>%
    unite("genus_spp", genus:species)
  return(kraken_table_clean)
}

PRJNA419379_kraken_abun_clean <- clean_output(PRJNA419379_kraken_abun)

### Plot this!
#### Plot bacterial read abundance against different samples
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

PRJNA419379_kraken_prop_plot <- PRJNA419379_kraken_abun_clean %>%
  filter(read_count > 2000) %>%
  filter(domain == "Bacteria") %>%
  filter(genus_spp != "NA_NA") %>%
  filter(genus_spp != "NIES_3974") %>%
  filter(genus_spp != "viridis_NA") %>%
  ggplot() + 
  aes(x = reorder(SampleID, read_count), y = read_count, fill = genus_spp) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_hue(l = 70, c = 60) +
  theme_classic() + 
  coord_flip()

pdf("PRJNA419379_kraken_prop_plot.pdf")
print(PRJNA419379_kraken_prop_plot)
dev.off()


PRJNA419379_kraken_abun_per_samp <- plot_per_sample_reads(PRJNA419379_kraken_abun_clean, 2000)

pdf("PRJNA419379_kraken_abun_per_samp.pdf")
print(PRJNA419379_kraken_abun_per_samp)
dev.off()

#### Plot overall read abundance accross samples
PRJNA419379_kraken_overall <- PRJNA419379_kraken_abun_clean %>%
  filter(domain == "Bacteria") %>%
  filter(genus_spp != "NIES_3974") %>%
  filter(genus_spp != "symbiont_NA") %>%
  filter(genus_spp != "viridis_NA") %>%
  filter(genus_spp != "NA_NA") %>%
  group_by(genus_spp) %>%
  summarise(average_read_count = sum(read_count)) %>%
  filter(average_read_count > 5000) %>%
  ggplot() +
  aes(x = reorder(genus_spp, average_read_count), y = average_read_count) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_minimal()

pdf("PRJNA419379_kraken_overall.pdf")
print(PRJNA419379_kraken_overall)
dev.off()

## Plot multiple reads from KRAKEN2 (PRJNA255) #####
ID = c("SRR1523064", "SRR1523066", "SRR1523067", "SRR1523068")

setwd("../../PRJNA255893_abun_stats/version2/")
filesnames <- list.files(pattern = ".txt")
filelist <- lapply(filesnames, read_tsv)
filelist2 <- lapply(filelist, setNames, columns)
filelist3 <- mapply('[<-', filelist2, "SampleID", value = ID, SIMPLIFY = F)

PRJNA255893_kraken_abun <- bind_rows(filelist3)

clean_output <- function(kraken_table){
  kraken_table_clean <- kraken_table %>%
    separate(Spp, 
             c("d", "domain", "k", "kingdom", "p", "phylum", 
               "c", "class", "o", "order", "f", "familiy", "g", "genus",
               "species")) %>%
    select(domain, kingdom, phylum, class, order, familiy, genus, species, read_count, SampleID) %>%
    unite("genus_spp", genus:species)
  return(kraken_table_clean)
}

PRJNA255893_kraken_abun_clean <- clean_output(PRJNA255893_kraken_abun)

### Plot this!

#### Plot bacterial read abundance against different samples
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

PRJNA255893_prop_plot <- PRJNA255893_kraken_abun_clean %>%
  filter(read_count > 2000) %>%
  filter(domain == "Bacteria") %>%
  filter(genus_spp != "NA_NA") %>%
  filter(genus_spp != "NIES_3974") %>%
  ggplot() + 
  aes(x = reorder(SampleID, read_count), y = read_count, fill = genus_spp) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_hue(l = 70, c = 60) +
  theme_classic() + 
  coord_flip()

pdf("PRJNA255893_prop_plot.pdf")
print(PRJNA255893_prop_plot)
dev.off()

PRJNA255893_kraken_abun_per_samp <- plot_per_sample_reads(PRJNA255893_kraken_abun_clean, 2000)

pdf("PRJNA255893_kraken_abun_per_samp.pdf")
print(PRJNA255893_kraken_abun_per_samp)
dev.off()

#### Plot overall read abundance accross samples
exclude = c("NIES_3974",
            "symbiont_NA",
            "viridis_NA")
PRJNA255893_ab_plot <- PRJNA255893_kraken_abun_clean %>%
  filter(domain == "Bacteria") %>%
  filter(genus_spp != "NIES_3974") %>%
  filter(genus_spp != "symbiont_NA") %>%
  filter(genus_spp != "viridis_NA") %>%
  filter(genus_spp != "NA_NA") %>%
  filter(genus_spp != "PSE_93") %>%
  filter(genus_spp != "NIES_3708") %>%
  group_by(genus_spp) %>%
  summarise(average_read_count = sum(read_count)) %>%
  filter(average_read_count > 4000) %>%
  ggplot() +
  aes(x = reorder(genus_spp, average_read_count), y = average_read_count) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_minimal()

pdf("PRJNA255893_ab_plot.pdf")
print(PRJNA255893_ab_plot)
dev.off()


## Plot overall abundances found

PRJNA255893_tot_reads <- PRJNA255893_kraken_abun_clean %>%
  group_by(SampleID) %>%
  summarise(total_reads = sum(read_count))

PRJNA255893_kraken_abun_clean_row_add <- PRJNA255893_kraken_abun_clean %>%
  left_join(x = PRJNA255893_tot_reads, by = "SampleID") %>%
  mutate(ProjectID = "PRJNA255893") %>%
  mutate(total_read_norm = read_count / total_reads)

PRJNA419379_tot_reads <- PRJNA419379_kraken_abun_clean %>%
  group_by(SampleID) %>%
  summarise(total_reads = sum(read_count))

PRJNA419379_kraken_abun_clean_row_add <- PRJNA419379_kraken_abun_clean %>%
  left_join(x = PRJNA419379_tot_reads, by = "SampleID") %>%
  mutate(ProjectID = "PRJNA419379") %>%
  mutate(total_read_norm = read_count / total_reads)


PRJNA385349_tot_reads <- PRJNA385349_kraken_abun_clean %>%
  group_by(SampleID) %>%
  summarise(total_reads = sum(read_count))

PRJNA385349_kraken_abun_clean_row_add <- PRJNA385349_kraken_abun_clean %>%
  left_join(x = PRJNA385349_tot_reads, by = "SampleID") %>%
  mutate(ProjectID = "PRJNA385349") %>%
  mutate(total_read_norm = read_count / total_reads)

sample_abundance_table <- rbind(PRJNA255893_kraken_abun_clean_row_add,
                                PRJNA419379_kraken_abun_clean_row_add,
                                PRJNA385349_kraken_abun_clean_row_add)

total_read_count <- sample_abundance_table %>%
  group_by(ProjectID) %>%
  summarise(total_reads = sum(read_count))

read_aln_statistics_project <- read_aln_staistics_clean %>%
  filter(read_type == "reads_mapped") %>%
  group_by(PRJN) %>%
  summarise(total_read = sum(read_stat)) %>%
  rename("ProjectID" = PRJN)


sample_abundance_class_plot <- sample_abundance_table %>%
  filter(domain == "Bacteria") %>%
  filter(read_count > 2000) %>%
  group_by(phylum, SampleID) %>%
  summarise(total_read_norm = sum(read_count)) %>%
  ggplot() + 
  aes(x = SampleID, y = total_read_norm, fill = phylum) +
  geom_bar(stat = "identity", position = "fill") + 
  coord_flip()

pdf("sample_abundance_class_plot.pdf")
print(sample_abundance_class_plot)
dev.off()

sample_abundance_table2 <- sample_abundance_table %>%
  group_by(domain, SampleID) %>%
  summarise(total_read_norm = sum(total_read_norm)) %>%
  ggplot() + 
  aes(x = SampleID, y = total_read_norm, fill = domain) +
  geom_bar(stat = "identity", position = "stack") + 
  coord_flip()

pdf("sample_abundance_table2.pdf")
print(sample_abundance_table2)
dev.off()

sample_abundance_table %>%
  filter(domain == "Viruses") %>%
  filter(read_count > 500) %>%
  group_by(SampleID, order) %>%
  summarise(total_read_norm = sum(total_read_norm)) %>%
  ggplot() + 
  aes(x = SampleID, y = total_read_norm, fill = order) + 
  geom_bar(stat = "identity", position = "stack") + 
  coord_flip()


## Plot Sample Heat Maps #####

PRJNA255893_kraken_abun_clean %>%
  filter(domain == "Bacteria") %>%
  filter(genus_spp != "NIES_3974") %>%
  filter(genus_spp != "symbiont_NA") %>%
  filter(genus_spp != "viridis_NA") %>%
  filter(genus_spp != "NA_NA") %>%
  group_by(genus_spp, SampleID) %>%
  summarise(average_read_count = sum(read_count)) %>%
  filter(average_read_count > 1500) %>%
  ggplot() +
  aes(x = SampleID, y = genus_spp) + 
  geom_tile(aes(fill = average_read_count), colour = "black") +
  scale_fill_gradient(low = "grey", high = "green") + 
  theme_classic()

PRJNA385349_kraken_abun_clean2 <- PRJNA385349_kraken_abun_clean %>%
  filter(domain == "Bacteria") %>%
  filter(genus_spp != "NIES_3974") %>%
  filter(SampleID != "SRR551970") %>%
  filter(genus_spp != "symbiont_NA") %>%
  filter(genus_spp != "viridis_NA") %>%
  filter(genus_spp != "NA_NA") %>%
  group_by(genus_spp, SampleID) %>%
  summarise(total_read_count = sum(read_count)) %>%
  filter(total_read_count > 2000)

PRJNA385349_total_reads <- PRJNA385349_kraken_abun_clean2 %>%
  group_by(SampleID) %>%
  summarise(total_reads = sum(total_read_count))

PRJNA385349_normalised_join <- PRJNA385349_kraken_abun_clean2 %>%
  left_join(x = PRJNA385349_total_reads, by = "SampleID") %>%
  mutate(normalised_read_count = total_read_count/total_reads * 100)

PRJNA385349_kraken_abun_per_samp <- PRJNA385349_normalised_join %>%
  mutate(normalised_read_count2 = ifelse(
    normalised_read_count > 90, 90, ifelse(normalised_read_count < 15, 15, normalised_read_count))) %>%
  ggplot() +
  aes(x = SampleID, y = genus_spp) + 
  geom_tile(aes(fill = normalised_read_count), colour = "black") +
  scale_fill_gradientn(colours = c("yellow", "orange", "blue"), 
                       values = scales::rescale(c(0, 5, 10, 20, 40, 60, 80))) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
## Plot KRAKEN mapping rate

read_aln_statistics_clean <- align_stats %>%
  select(SRR, total_seq, reads_mapped, reads_unmapped) %>%
  rename("SampleID" = SRR)

kraken_mapped <- bind_rows(PRJNA385349_kraken_abun_clean, 
                           PRJNA419379_kraken_abun_clean, 
                           PRJNA255893_kraken_abun_clean)

kraken_mapped_total_read_count <- kraken_mapped %>%
  group_by(SampleID) %>%
  summarise(total_read_count = sum(read_count))

kraken_mapped_taxon_count <- kraken_mapped %>%
  group_by(SampleID, domain) %>%
  summarise(total_read_count = sum(read_count))

overall_kraken_statistic <- read_aln_statistics_clean %>%
  left_join(kraken_mapped_total_read_count, by = "SampleID") %>%
  pivot_longer(cols = c("reads_unmapped", "total_read_count"), names_to = "read_type", values_to = "no.reads") %>%
  select(SampleID, read_type, no.reads) %>%
  filter(no.reads != "NA") %>%
  mutate(norm_no.reads = sqrt(no.reads))

overall_kraken_statistic_plot <- overall_kraken_statistic %>%
  ggplot() + 
  aes(x = SampleID, y = norm_no.reads, fill = read_type) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_calc() + 
  coord_flip()

pdf("overall_kraken_statistic_plot.pdf")
print(overall_kraken_statistic_plot)
dev.off()
