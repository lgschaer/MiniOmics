#All Files Transcriptomics Analysis

#packages needed
library(tidyverse)
library(csv)

#read tsv files from prokka analysis of contigs.fa
L1_genes <- read.csv(file = "/home/lgschaer/old/MiniOmics/Metagenomics/Laura1/Laura1/prokka/Laura1.contigs/Laura1.contigs.tsv", sep = '\t') 
head(L1_genes)

E2_genes <- read.csv(file = "/home/lgschaer/old/MiniOmics/Metagenomics/Emma2/Emma2/prokka/Emma2.contigs/Emma2.contigs.tsv", sep = '\t') 
head(E2_genes)

gene_list <- L1_genes %>%
  full_join(E2_genes) %>%
  #filter(product != "hypothetical protein") %>%
  mutate(temp = locus_tag) %>%
  separate(temp, into = c("BinID", "geneID"), sep = "_") %>%
  mutate(BinID2 = BinID) %>%
  separate(BinID2, into = c("Enrichment", "trash"), sep = ".c") %>%
  dplyr::select(-c(trash, BinID))
head(gene_list)
dim(gene_list)

#save the filtered genes as a csv
#write_csv(gene_list, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_contig_gene_list.csv")



###### Get a complete list of known aromatic-degrading genes from KEGG

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("KEGGREST")

library("KEGGREST")
#Manually get pathway IDs from: https://www.genome.jp/pathway/ec00362

test <- keggGet(c("ec00624",   #Polycyclic aromatic hydrocarbon degradation
                  "ec00362",   #Benzoate degradation
                  "ec00627",   #Aminobenzoate degradation
                  "ec00020"))  #TCA cycle

e1 <- as.data.frame(test[[1]]$ENZYME) %>% 
  mutate(EC_number = test[[1]]$ENZYME, 
         General_Pathway = test[[1]]$NAME)

e2 <- as.data.frame(test[[2]]$ENZYME) %>% 
  mutate(EC_number = test[[2]]$ENZYME, 
         General_Pathway = test[[2]]$NAME)

e3 <- as.data.frame(test[[3]]$ENZYME) %>% 
  mutate(EC_number = test[[3]]$ENZYME, 
         General_Pathway = test[[3]]$NAME)

e4 <- as.data.frame(test[[4]]$ENZYME) %>% 
  mutate(EC_number = test[[4]]$ENZYME, 
         General_Pathway = test[[4]]$NAME)

head(e1)

interesting_enzymes <- e1 %>%
  full_join(e2, by = c("EC_number", "General_Pathway")) %>%
  full_join(e3, by = c("EC_number", "General_Pathway")) %>%
  full_join(e4, by = c("EC_number", "General_Pathway")) %>%
  dplyr::select(c(EC_number, General_Pathway)) %>%
  distinct()
head(interesting_enzymes)


####Load transcriptomics files

## Load GTF Files

L1_TPA_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/L1_TPA.gtf", header = FALSE)
L1_DCPET_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/L1_DCPET.gtf", header = FALSE)
L1_HDPE_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/L1_HDPE.gtf", header = FALSE)
L1_TA_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/L1_TA.gtf", header = FALSE)
L1_EG_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/L1_EG.gtf", header = FALSE)
E2_TPA_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/E2_TPA.gtf", header = FALSE)
E2_DCPET_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/E2_DCPET.gtf", header = FALSE)
E2_HDPE_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/E2_HDPE.gtf", header = FALSE)
E2_TA_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/E2_TA.gtf", header = FALSE)
E2_EG_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/E2_EG.gtf", header = FALSE)

head(L1_TPA_gtf)

clean_gtfs <- L1_TPA_gtf %>%
  rbind(L1_DCPET_gtf) %>%
  rbind(L1_HDPE_gtf) %>%
  rbind(L1_TA_gtf) %>%
  rbind(L1_EG_gtf) %>%
  rbind(E2_TPA_gtf) %>%
  rbind(E2_DCPET_gtf) %>%
  rbind(E2_HDPE_gtf) %>%
  rbind(E2_TA_gtf) %>%
  rbind(E2_EG_gtf) %>%
  mutate(contigName = V1) %>%
  separate(V2, into = c("extra", "prodigal"), sep = ":") %>%
  separate(V9, into = c("transcript_id", "gene_id"), sep = "; ") %>%
  separate(gene_id, into = c("excess", "locus_tag"), sep = " ") %>%
  dplyr::select(contigName, locus_tag, prodigal) %>%
  filter(!is.na(locus_tag)) %>%
  distinct()
head(clean_gtfs)

#check dimensions - did this to make sure the number of rows made sense running only the rbind lines of code
#expected <- length(L1_TPA_gtf$V1)+length(L1_DCPET_gtf$V1)+length(L1_HDPE_gtf$V1)+length(L1_TA_gtf$V1)+length(L1_EG_gtf$V1)+length(E2_TPA_gtf$V1)+length(E2_DCPET_gtf$V1)+length(E2_HDPE_gtf$V1)+length(E2_TA_gtf$V1)+length(E2_EG_gtf$V1)
#expected
#length(clean_gtfs$contigName)


E2_BinID_Key <- read.csv("/home/lgschaer/old/MiniOmics/Metagenomics/Emma2/Emma2/Metabat2_Bin_Cat_Output/out.BAT.ORF2LCA.txt", header = TRUE, sep = '\t')
L1_BinID_Key <- read.csv("/home/lgschaer/old/MiniOmics/Metagenomics/Laura1/Laura1/Metabat2_Bin_Cat_Output/out.BAT.ORF2LCA.txt", header = TRUE, sep = '\t')
head(L1_BinID_Key)

#checking dimensions
#length(E2_BinID_Key$bin)+length(L1_BinID_Key$bin)
#length(BinID_Key$bin)

BinID_Key <- L1_BinID_Key %>%
  full_join(E2_BinID_Key) %>%
  separate(X..ORF, into = c("BinID", "contigName"), sep = ".fa_") %>%
  separate(contigName, into = c("contigNameA", "contigNameB", "contigNumber"), sep = "_") %>%
  unite(contigName, contigNameA, contigNameB, sep = "_") %>%
  separate(bin, into = c("Enrichment", "Bin_Number"), sep = "metabat.") %>%
  separate(Bin_Number, into = c("Bin_Number", "Trash", sep = ".")) %>%
  dplyr::select(BinID, contigName, Enrichment, Bin_Number) %>%
  distinct()
head(BinID_Key)
dim(BinID_Key)

#load taxa and quality data (cleaned up in CheckM and CAT analysis scripts)
Taxa_Quality_Summary <- read.csv("/home/lgschaer/old/MiniOmics/clean_data/L1_E2_Taxa_Quality_Summary.csv")
head(Taxa_Quality_Summary)

sum(duplicated(Taxa_Quality_Summary$BinID))

#load count data for each sample
L1_TPA_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/counts/L1_TPA_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, TPA = V2, .keep="none")
L1_DCPET_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/counts/L1_DCPET_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, DCPET = V2, .keep="none")
L1_HDPE_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/counts/L1_HDPE_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, HDPE = V2, .keep="none")
L1_TA_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/counts/L1_TA_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, TA = V2, .keep="none")
L1_EG_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/counts/L1_EG_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, EG = V2, .keep="none")
E2_TPA_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/counts/E2_TPA_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, TPA = V2, .keep="none")
E2_DCPET_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/counts/E2_DCPET_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, DCPET = V2, .keep="none")
E2_HDPE_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/counts/E2_HDPE_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, HDPE = V2, .keep="none")
E2_TA_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/counts/E2_TA_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, TA = V2, .keep="none")
E2_EG_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Emma2/counts/E2_EG_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, EG = V2, .keep="none")


full_count_data_w_quality_taxa <- L1_TPA_counts %>%
  full_join(L1_DCPET_counts) %>%
  full_join(L1_HDPE_counts) %>%
  full_join(L1_TA_counts) %>%
  full_join(L1_EG_counts) %>%
  full_join(E2_TPA_counts) %>%  
  full_join(E2_DCPET_counts) %>%
  full_join(E2_HDPE_counts) %>%
  full_join(E2_TA_counts) %>%
  full_join(E2_EG_counts) %>%
  group_by(locus_tag) %>%
  filter(locus_tag != "__too_low_aQual" 
         & locus_tag != "__not_aligned" 
         & locus_tag != "__no_feature" 
         & locus_tag != "__ambiguous") %>%
  left_join(gene_list, by = "locus_tag") %>%
  left_join(interesting_enzymes) %>%
  left_join(clean_gtfs) %>%
  left_join(BinID_Key) %>% ##
  left_join(Taxa_Quality_Summary) %>%
  mutate(
    General_Pathway = ifelse(is.na(General_Pathway), "Unknown_Pathway", General_Pathway)
  ) %>%
  filter(!is.na(BinID)) %>%
  distinct()
head(full_count_data_w_quality_taxa)
dim(full_count_data_w_quality_taxa)
View(full_count_data_w_quality_taxa)
#write_csv(full_count_data_w_quality_taxa, file = "/home/lgschaer/old/MiniOmics/clean_data/full_count_data_w_quality_taxa.csv")

sum(duplicated(full_count_data_w_quality_taxa$locus_tag))

#####Plotting genes

known_genes_only <- full_count_data_w_quality_taxa %>%
  filter(product == "Mono(2-hydroxyethyl) terephthalate hydrolase"|
           gene == "tphA1I"|gene == "tphBI_1"|gene == "tphA3I"|
           gene == "tphA1II"|gene == "tphBI"|gene == "tphA2I"|
           gene == "tphA3I"|gene == "tphA3II"|gene == "tphA2I_1"|
           gene == "pcaG"|gene == "pcaG_1"|gene == "pcaG_2"|
           gene == "pcaH"|gene == "pcaH_1"|gene == "pcaH_2"|gene=="tcbC"|
           gene == "alkB1"|gene == "alkB2"|gene=="alkB2") %>%
  pivot_longer(cols = c("TPA", "DCPET", "TA", "EG", "HDPE"), names_to = "Substrate", values_to = "Counts") %>%
  filter(!is.na(Counts) & Counts != 0) %>%
  mutate(
    Num_Transcripts = ifelse(Counts <= 100, "1 to 100",
                             ifelse(Counts > 100 & Counts <= 500, "100 to 500",
                                    ifelse(Counts > 500 & Counts <= 1000, "500 to 1,000", 
                                           ifelse(Counts > 1000, "> 1,000", "other"))))
  )
View(known_genes_only)

unique(known_genes_only$product)

order <- c("Mono(2-hydroxyethyl) terephthalate hydrolase",
           "Terephthalate 1,2-dioxygenase, terminal oxygenase component subunit alpha 1", 
           "Terephthalate 1,2-dioxygenase, terminal oxygenase component subunit beta 2",
           "1,2-dihydroxy-3,5-cyclohexadiene-1,4-dicarboxylate dehydrogenase",
           "Protocatechuate 3,4-dioxygenase alpha chain",
           "Protocatechuate 3,4-dioxygenase beta chain")

labels <- c("Mono(2-hydroxyethyl)\nterephthalate hydrolase",
           "Terephthalate\n1,2-dioxygenase\nterminal oxygenase\ncomponent subunit alpha 1", 
           "Terephthalate\n1,2-dioxygenase\nterminal oxygenase\ncomponent subunit beta 2",
           "1,2-dihydroxy\n-3,5-cyclohexadiene\n-1,4-dicarboxylate\ndehydrogenase",
           "Protocatechuate\n3,4-dioxygenase\nalpha chain",
           "Protocatechuate\n3,4-dioxygenase\nbeta chain")

#known_genes_only$product <- reorder(known_genes_only$product, order)


ggplot(known_genes_only, aes(x = product, y = Best_Classification))+
  geom_tile(aes(fill = Num_Transcripts), color = "black") +
  facet_grid(rows = vars(Substrate), cols = vars(Enrichment), scales = "free_y") +
  scale_fill_manual(values = colors) +
  scale_x_discrete(limits = order, labels = labels) +
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
        #axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 9),
        axis.title.x = element_blank())


filt_count_data <- full_count_data_w_quality_taxa %>%
  filter(Quality == "High") %>%
  #filter(General_Pathway != "Unknown_Pathway") %>%
  pivot_longer(cols = c("TPA", "DCPET", "HDPE", "TA", "EG"), names_to = "Substrate", values_to = "Counts") %>%
  filter(Counts != 0) %>%
  mutate(
    Num_Transcripts = ifelse(Counts >= 1 & Counts <= 100, "1 to 100",
                             ifelse(Counts > 100 & Counts <= 500, "100 to 500",
                                    ifelse(Counts > 500 & Counts <= 1000, "500 to 1,000", 
                                           ifelse(Counts > 1000, "> 1,000", "other"))))
  ) %>%
  dplyr::select(-c("ftype", "length_bp", "COG", "geneID", "prodigal"))
head(filt_count_data)
#View(filt_count_data)

filt_count_data$Num_Transcripts <- factor(filt_count_data$Num_Transcripts, levels = c("1 to 100", "100 to 500", "500 to 1,000", "> 1,000"))

colors <- c("dodgerblue", "firebrick", "darkgreen", "orange", "purple4")

ggplot(filt_count_data, aes(x = gene, y = Best_Classification))+
  facet_grid(rows = vars(Enrichment), scales = "free_y")+
  #geom_tile(aes(fill = Num_Transcripts)) +
  geom_jitter(aes(fill = Substrate, shape = Substrate), color = "black", size = 3) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(legend.position = "bottom",
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust =  0.5, size = 10),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        legend.text = element_text(size = 18),
        axis.title.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=20)))

filt_count_data2 <- filter(filt_count_data, General_Pathway != "Unknown_Pathway")

ggplot(filt_count_data2, aes(x = product, y = Best_Classification))+
  facet_grid(rows = vars(Substrate, Enrichment), scales = "free")+
  geom_tile(aes(fill = Num_Transcripts), color = "black")+
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 35, vjust = 1, hjust =  1, size = 6),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        axis.title.x = element_blank())

#####stopped here

ggplot(filt_L1_TPA2, aes(x = product, y = Best_Classification))+
  geom_tile(aes(fill = Num_Transcripts)) +
  #facet_grid(rows = vars(General_Pathway)) +
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust =  0.5, size = 10),
        #axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        axis.title.x = element_blank())

##### TPA pathway genes only

#####Plotting genes involved in TPA breakdown with taxonomy classification

TPA_genes <- filt_L1_TPA %>%
  filter(product == "Mono(2-hydroxyethyl) terephthalate hydrolase"|
           gene == "tphA1I"|gene == "tphBI_1"|gene == "tphA3I"|
           gene == "tphA1II"|gene == "tphBI"|gene == "tphA2I"|
           gene == "tphA3I"|gene == "tphA3II"|gene == "tphA2I_1"|
           gene == "pcaG"|gene == "pcaG_1"|gene == "pcaG_2"|
           gene == "pcaH"|gene == "pcaH_1"|gene == "pcaH_2"|gene=="tcbC"|
           gene == "alkB1"|gene == "alkB2"|gene=="alkB2") #%>%
#  mutate(
#   gene = ifelse(product == "Mono(2-hydroxyethyl) terephthalate hydrolase", "no_name", gene),
#  P1 = ifelse(gene == "alkB1"|gene == "alkB2"|gene=="alkB2", "Alkane", gene),
# P2 = ifelse(P1 == "pcaG"|P1 == "pcaG_1"|P1 == "pcaG_2"|
#              P1 == "pcaH"|P1 == "pcaH_1"|P1 == "pcaH_2", "Protocatechuate", P1),
#    P3 = ifelse(P2 == "no_name", "Cleave Ethylene Glycol", P2),
#   P4 = ifelse(P3 == "tcbC", "Catechol", P3),
#  Pathway = ifelse(P4 == "tphA1I"|P4 == "tphBI_1"|P4 == "tphA3I"|
#                    P4 == "tphA1II"|P4 == "tphBI"|P4 == "tphA2I"|
#                   P4 == "tphA3I"|P4 == "tphA3II"|P4 == "tphA2I_1", "Terephthalate", P4) 
#  ) #%>%
#unite(Taxa_Name, Enrichment, Taxa_Name, sep = "_", remove = FALSE) %>%
#select(c("Bin_ID", "Enrichment", "Bin_Number", "Pathway", "product", "gene", "Taxa_Name", "Completeness", "Contamination"))
head(TPA_genes)
#View(TPA_genes)

colors <- c("firebrick", "orange", "green", "dodgerblue", "purple")

ggplot(TPA_genes, aes(x = gene, y = Classification))+
  geom_tile(aes(fill = gene), color = "black") +
  #  facet_grid(rows = vars(Enrichment), shrink = TRUE)+
  scale_fill_manual(values = colors) +
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(), #angle = 0, size = 10
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        axis.title.x = element_blank())


#High quality only

hq_only <- Taxa_Quality_Summary %>% filter(Quality == "High")

HQ_Aromatic_Genes <- hq_only %>%
  inner_join(genes_filt) %>%
  inner_join(interesting_enzymes, by = "EC_number") %>%
  dplyr::select(c("BinID", "Enrichment", "product", "gene", "Classification", "Confidence", "Tax_Level", "Quality", "EC_number", "General_Pathway")) %>%
  filter(!is.na(Classification))%>%
  filter(!is.na(General_Pathway))
head(HQ_Aromatic_Genes)

write_csv(HQ_Aromatic_Genes, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_high_quality_bins_only_aromatic_genes.csv")

any(is.na(HQ_Aromatic_Genes$General_Pathway))

colors <- c("firebrick", "orange", "green", "dodgerblue", "purple")

heat<-ggplot(HQ_Aromatic_Genes, aes(x = product, y = Classification))+
  geom_tile(aes(fill = General_Pathway), color = "black") +
  facet_grid(rows = vars(Enrichment), scales = "free_y")+
  scale_fill_manual(values = colors) +
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text.x = element_blank(), #angle = 0, size = 10
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        axis.title.x = element_blank())
heat



MG_Example <- read.csv("/home/lgschaer/old/MiniOmics/clean_data/L1_E2_genes_unfiltered.csv")
head(MG_Example)


MG_Example2 <- MG_Example %>%
  separate(locus_tag, into =c("Enrichment", "tag"), sep = "metabat.") %>%
  dplyr::select(c("product", "gene", "EC_number", "COG", "length_bp", "Enrichment"))
head(MG_Example2)

write_csv(MG_Example2, "/home/lgschaer/old/MiniOmics/example_files/metagenomics_example.csv")
