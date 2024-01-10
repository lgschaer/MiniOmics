##
##

#load packages
library(tidyverse)

#load files
cat_L1 <- read.table(file = '/home/lgschaer/old/MiniOmics/Metagenomics/Laura1/Laura1/out.BAT.bin2classification_w_tax.txt', sep = '*', header = FALSE)
cat_E2 <- read.table(file = '/home/lgschaer/old/MiniOmics/Metagenomics/Emma2/Emma2/out.BAT.bin2classification_w_tax.txt', sep = ',', header = FALSE)

#clean up files, make them usable
taxa_clean_pre <- cat_L1 %>%
  full_join(cat_E2) %>%
  separate(V1, into = c("BinID", "col2", "No_ORFs", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16", "col17"), sep = "\\t", fill = "right") %>%
  pivot_longer(cols=c(col8, col9, col10, col11, col12, col13, col14, col15, col16, col17), names_to = "throw_away", values_to = "tax_info") %>%
  filter(!is.na(tax_info) & BinID != "") %>%
  separate(BinID, into = c("BinID", "not_needed"), sep = ".f") %>%
  dplyr::select(-c(throw_away, col2, col4, col5, col6, col7, not_needed)) %>%
  separate(tax_info, into = c("Classification", "tax_info"), sep = " \\(") %>%
  separate(tax_info, into = c("Tax_Level", "Confidence"), sep = "\\): ")
head(taxa_clean_pre)

taxa_clean_summary <- taxa_clean_pre %>%
  group_by(BinID, Tax_Level) %>%
  summarise(
    count = n()
  ) %>% 
  filter(count>1)
View(taxa_clean_summary)

taxa_clean <- taxa_clean_pre %>%
  mutate(should_I_keep = ifelse(BinID %in% taxa_clean_summary$BinID & Tax_Level %in% taxa_clean_summary$Tax_Level, 
                                "sus", "ok"),
         Tax_Level2 = Tax_Level,
         Confidence2 = Confidence) %>%
  unite(Classification_Confidence, Classification, Confidence2, Tax_Level2, sep = "_") %>%
  filter(should_I_keep=="ok") %>%
  dplyr::select(-c(Confidence, should_I_keep)) %>%
  pivot_wider(names_from = Tax_Level, values_from = Classification_Confidence) %>%
  mutate(
    Best_Classification = ifelse(!is.na(superkingdom), superkingdom, "fill later"),
    Best_Classification = ifelse(!is.na(phylum), phylum, Best_Classification),
    Best_Classification = ifelse(!is.na(subphylum), subphylum, Best_Classification),
    Best_Classification = ifelse(!is.na(clade), clade, Best_Classification),
    Best_Classification = ifelse(!is.na(class), class, Best_Classification),
    Best_Classification = ifelse(!is.na(order), order, Best_Classification),
    Best_Classification = ifelse(!is.na(family), family, Best_Classification),
    Best_Classification = ifelse(!is.na(genus), genus, Best_Classification),
    Best_Classification = ifelse(!is.na(species), species, Best_Classification)
  )
head(taxa_clean)
View(taxa_clean)

taxa_summary <- taxa_clean %>%
  full_join(Summary) %>%
  dplyr::select(BinID, phylum, class, order, family, genus, species, Quality) %>%
  filter(Quality=="High")%>%
  separate(phylum, into = c("phylum", "phylum_confidence", "discard"), sep = "_") %>%
  separate(class, into = c("class", "class_confidence", "discard"), sep = "_") %>%
  separate(order, into = c("order", "order_confidence", "discard"), sep = "_") %>%
  separate(family, into = c("family", "family_confidence", "discard"), sep = "_") %>%
  separate(genus, into = c("genus", "genus_confidence", "discard"), sep = "_") %>%
  separate(species, into = c("species", "species_confidence", "discard"), sep = "_") %>%
  dplyr::select(BinID, phylum, phylum_confidence, class, class_confidence, order, order_confidence, family, family_confidence, genus, genus_confidence, species, species_confidence) %>%
  group_by(genus) %>%
  summarise(
    count = n()
  ) %>%
  arrange(-count)
View(taxa_summary)

#write_csv(taxa_clean, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_Taxa_Assignments.csv")


#make sure dimensions make sense
dim(cat_L1)
dim(cat_E2)
dim(taxa_clean)

length(cat_L1$V1)+length(cat_E2$V1)


Taxa_Quality_Summary <- taxa_clean %>%
  full_join(Summary) %>%
  dplyr::select(BinID, Enrichment, Best_Classification, Completeness, Contamination, Quality) %>%
  mutate(
    Best_Classification2 = Best_Classification,
  ) %>%
  separate(Best_Classification2, into = c("Classification", "Confidence", "Tax_Level"), sep = "_") %>%
  mutate(Confidence = as.double(Confidence))
head(Taxa_Quality_Summary)

#write_csv(Taxa_Quality_Summary, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_Taxa_Quality_Summary.csv")

ggplot(Taxa_Quality_Summary, aes(y = Classification, x = Quality))+
  facet_grid(cols = vars(Enrichment))+
  geom_tile(aes(fill=Completeness))
