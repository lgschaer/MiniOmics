#Contigs Prokka Analysis

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
  dplyr::select(-trash)
head(gene_list)
dim(gene_list)

#save the filtered genes as a csv
#write_csv(gene_list, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_contig_gene_list.csv")

MetaCyc <- read.table("/home/lgschaer/old/MiniOmics/MetaCyc_PET_to_TCA_full_pathways_LGS_03022022.txt", header = TRUE, sep = "\t", dec = "_")
head(MetaCyc)

enzyme_list <- MetaCyc

###### Trying to get more complete list of genes

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")

library("KEGGREST")
# Pathway IDs from: https://www.genome.jp/pathway/ec00362

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
  dplyr::select(c(EC_number, General_Pathway))
head(interesting_enzymes)


####Load transcriptomics files
L1_TPA_gtf <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/L1_TPA.gtf", header = FALSE) %>%
  mutate(contigName = V1) %>%
  separate(V2, into = c("extra", "prodigal"), sep = ":") %>%
  separate(V9, into = c("transcript_id", "gene_id"), sep = "; ") %>%
  separate(gene_id, into = c("excess", "locus_tag"), sep = " ") %>%
  dplyr::select(contigName, locus_tag, prodigal) %>%
  filter(!is.na(locus_tag)) %>%
  distinct()
head(L1_TPA_gtf)
dim(L1_TPA_gtf)

BinID_Key <- read.csv("/home/lgschaer/old/MiniOmics/Metagenomics/Laura1/Laura1/Metabat2_Bin_Cat_Output/out.BAT.ORF2LCA.txt", header = TRUE, sep = '\t') %>%
  separate(X..ORF, into = c("BinID", "contigName"), sep = ".fa_") %>%
  separate(contigName, into = c("contigNameA", "contigNameB", "contigNumber"), sep = "_") %>%
  unite(contigName, contigNameA, contigNameB, sep = "_") %>%
  separate(bin, into = c("Enrichment", "Bin_Number"), sep = "metabat.") %>%
  separate(Bin_Number, into = c("Bin_Number", "Trash", sep = ".")) %>%
  dplyr::select(BinID, contigName, Enrichment, Bin_Number) %>%
  distinct()
head(BinID_Key)
dim(BinID_Key)

Taxa_Quality_Summary <- read.csv("/home/lgschaer/old/MiniOmics/clean_data/L1_E2_Taxa_Quality_Summary.csv")
head(Taxa_Quality_Summary)

L1_TPA_genes_plus_counts <- read.delim("/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/counts/L1_TPA_count.txt", header = FALSE) %>%
  mutate(locus_tag = V1, 
         L1_TPA = V2) %>%
  dplyr::select(-c("V1", "V2")) %>%
  filter(L1_TPA > 0) %>%
  left_join(L1_genes, by = "locus_tag") %>%
  left_join(interesting_enzymes, by = "EC_number") %>%
  left_join(L1_TPA_gtf, by = "locus_tag") %>%
  left_join(BinID_Key, by = "contigName") %>%
  left_join(Taxa_Quality_Summary, by = c("BinID", "Enrichment")) %>%
  mutate(
    General_Pathway = ifelse(is.na(General_Pathway), "Unknown_Pathway", General_Pathway)
  ) %>%
  filter(!is.na(BinID))
head(L1_TPA_genes_plus_counts)
dim(L1_TPA_genes_plus_counts)
#View(L1_TPA_genes_plus_counts)

#####Plotting genes

filt_L1_TPA <- L1_TPA_genes_plus_counts %>%
  filter(Quality == "High") %>%
  mutate(
    Num_Transcripts = ifelse(L1_TPA >= 1 & L1_TPA <= 100, "1 to 100",
                            ifelse(L1_TPA > 100 & L1_TPA <= 500, "100 to 500",
                                   ifelse(L1_TPA > 500 & L1_TPA <= 1000, "500 to 1,000", 
                                          ifelse(L1_TPA > 1000, "> 1,000", "other"))))
  ) 
  #filter(General_Pathway != "Unknown_Pathway")
head(filt_L1_TPA)

range(filt_L1_TPA$L1_TPA)

filt_L1_TPA$Num_Transcripts <- factor(filt_L1_TPA$Num_Transcripts, levels = c("1 to 100", "100 to 500", "500 to 1,000", "> 1,000"))

colors <- c("dodgerblue", "firebrick", "darkgreen", "orange")

ggplot(filt_L1_TPA, aes(x = gene, y = Best_Classification))+
  geom_tile(aes(fill = Num_Transcripts)) +
  scale_fill_manual(values = colors)+
  theme_classic()+
  theme(legend.position = "bottom",
       # axis.text.x = element_text(angle = 90, hjust = 1, vjust =  0.5, size = 10),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        axis.title.x = element_blank())

ggplot(filt_L1_TPA, aes(x = L1_TPA, y = product))+
  facet_grid(cols = vars(General_Pathway), scales = "free")+
  geom_boxplot() +
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust =  1, size = 6),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        axis.title.x = element_blank())

unique(filt_L1_TPA$General_Pathway)

filt_L1_TPA2 <- filt_L1_TPA %>%
  filter(General_Pathway != "Unknown_Pathway")
dim(filt_L1_TPA)
dim(filt_L1_TPA2)

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

ggplot(TPA_genes, aes(x = product, y = Classification))+
  geom_tile(aes(fill = product), color = "black") +
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



