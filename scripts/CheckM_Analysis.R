##Looking at the CheckM Output for MiniOmics Co-Assemblies. 
##January 4th 2022
##Laura Schaerer


##This first part is getting the data into a usable format

#packages
library(tidyverse)

#load files
CM_L1 <- read.table(file = '/home/lgschaer/old/MiniOmics/Metagenomics/Laura1/Laura1/CheckM/storage/bin_stats_ext.tsv', sep = '\t', header = FALSE)
CM_E2 <- read.table(file = '/home/lgschaer/old/MiniOmics/Metagenomics/Emma2/Emma2/CheckM/storage/bin_stats_ext.tsv', sep = '\t', header = FALSE)

#try to get un usable format
head(CM_L1)

Summary <- CM_L1 %>%
  full_join(CM_E2) %>%
  mutate(BinID = V1) %>%
  separate(V1, into = c("Enrichment", "Excess"), sep = "metabat") %>%
  dplyr::select(BinID, Enrichment, V2) %>%
  separate(V2, into = c("Marker_Genes", "Genomes", "Markers", "Marker_Sets", "Ambiguous_Bases", "Scaffolds", "Contigs", "Predicted_Genes"), sep = ", #") %>%
  separate(Marker_Sets, into = c("Marker_Sets", "Completeness"), sep = ", Com") %>%
  separate(Completeness, into = c("Completeness", "Contamination", "GC", "GC_std", "Genome_Size"), sep = "\\,") %>%
  mutate(
    Genomes = as.double(gsub("[^0-9.-]", "", Genomes)),
    Markers = as.double(gsub("[^0-9.-]", "", Markers)),
    Completeness = as.double(gsub("[^0-9.-]", "", Completeness)),
    Contamination = as.double(gsub("[^0-9.-]", "", Contamination)),
    GC = as.double(gsub("[^0-9.-]", "", GC)),
    GC_std = as.double(gsub("[^0-9.-]", "", GC_std)),
    Genome_Size = as.double(gsub("[^0-9.-]", "", Genome_Size)),
    Ambiguous_Bases = as.double(gsub("[^0-9.-]", "", Ambiguous_Bases)),
    Scaffolds = as.double(gsub("[^0-9.-]", "", Scaffolds)),
    Quality = ifelse(Completeness < 50 & Contamination < 10, "Low", "Contaminated"),
    Quality = ifelse(Completeness >= 50 & Contamination < 10, "Medium", Quality),
    Quality = ifelse(Completeness > 90 & Contamination < 5, "High", Quality)
  )
head(Summary)
View(Summary)

#write_csv(Summary, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_Metagenome_Quality.csv")


#Make sure dimensions add up
dim(Summary)
dim(CM_E2)
dim(CM_L1)

expected <- length(CM_E2$V1)+length(CM_L1$V1)
expected

##This portion of the script is to make summary plots to look at the quality of our metagenomes

#Summary of completeness and contamination
ggplot(Summary)+
  geom_histogram(aes(x = Completeness), fill = "blue", color = "black", alpha = 0.4)+
  geom_histogram(aes(x = Contamination), fill = "red", color = "black", alpha = 0.4)+
  facet_grid(cols = vars(Enrichment))+
  xlab("Blue = Completeness, Red = Contamination")+
  theme_bw()

#Summary of bin quality
Summary$Quality <- factor(Summary$Quality, levels = c("Contaminated", "Low", "Medium", "High"))

ggplot(Summary)+
  geom_bar(aes(x = Quality), fill = "blue", color = "black")+
  facet_grid(cols = vars(Enrichment))+
  theme_bw()

Summary_Stats <- Summary %>%
  group_by(Enrichment, Quality) %>%
  summarise(
    Count = n()
  )
View(Summary_Stats)
