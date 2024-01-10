#Prokka Analysis

#packages needed
library(tidyverse)
library(csv)
library(readr)

#move all tsv files to be analyzed from the prokka directory to the same directory, set that directory as wd

setwd("/home/lgschaer/old/MiniOmics/Metagenomics/prokka_output/")
getwd()

file_names <- dir("./") #where you have your files
file_names

genes <- do.call(rbind,lapply(file_names,read_tsv))
head(genes)
dim(genes)

#write_csv(genes, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_genes_unfiltered.csv")

genes_filt <- genes %>%
  filter(product != "hypothetical protein") %>%
  separate(locus_tag, into = c("BinID", "geneID"), sep = "_") %>%
  mutate(BinID2 = BinID) %>%
  separate(BinID2, into = c("Enrichment", "Bin_Number"), sep = "metabat.")
head(genes_filt)
View(genes_filt)

#save the filtered genes as a csv
#write_csv(genes_filt, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_genes_filtered.csv")

ggplot(genes_filt, aes(x = gene, y = BinID))+
  geom_tile(aes(fill=Enrichment)) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
    #axis.text.x = element_text(angle = 90)
  )


#I got this from here: https://metacyc.org/group?id=biocyc17-58703-3855250609# 
MetaCyc <- read.table("/home/lgschaer/old/MiniOmics/MetaCyc_PET_to_TCA_full_pathways_w_genes_LGS_03022022.txt", header = TRUE, sep = "\t", dec = "_")
head(MetaCyc)

MetaCyc2 <- MetaCyc %>%
  select(Names, Enzymes.of.pathway, Genes.of.pathway) %>%
  separate(Names, into = c("MetaCyc_Pathway_Name", "Names2", "Names3", "Names4", "Names5", "Names6", "Names7"), sep = " // | / ") %>%
  separate(Enzymes.of.pathway, into = c("Enz1", "Enz2", "Enz3", "Enz4", "Enz5", "Enz6", "Enz7", "Enz8", "Enz9", "Enz10",
                                        "Enz11", "Enz12", "Enz13", "Enz14", "Enz15", "Enz16", "Enz17", "Enz18", "Enz19"), sep = " // | / ") %>%
  separate(Genes.of.pathway, into = c("Gen1", "Gen2", "Gen3", "Gen4", "Gen5", "Gen6", "Gen7", "Gen8", "Gen9", "Gen10", "Gen11", 
                                      "Gen12", "Gen13", "Gen14", "Gen15", "Gen16", "Gen17", "Gen18", "Gen19", "Gen20", "Gen21", "Gen22"), sep = " // | / ") %>%
  pivot_longer(cols = c("Enz1", "Enz2", "Enz3", "Enz4", "Enz5", "Enz6", "Enz7", "Enz8", "Enz9", "Enz10", "Enz11", "Enz12", "Enz13", "Enz14", "Enz15", "Enz16", "Enz17", "Enz18", "Enz19"), 
               names_to = "Enzyme_No", values_to = "product") %>%
  pivot_longer(cols = c("Gen1", "Gen2", "Gen3", "Gen4", "Gen5", "Gen6", "Gen7", "Gen8", "Gen9", "Gen10", "Gen11", 
                        "Gen12", "Gen13", "Gen14", "Gen15", "Gen16", "Gen17", "Gen18", "Gen19", "Gen20", "Gen21", "Gen22"), 
               names_to = "Gene_No", values_to = "MetaCyc_gene") %>%
  mutate(
    product = sub("</i>", "", product),
    product = sub("<i>", "", product)
 ) %>%
 # separate(product, into = c("descriptor", "product"), sep = " of ", fill = "right") %>%
  dplyr::select(-c("Names2", "Names3", "Names4", "Names5", "Names6", "Names7", "Enzyme_No", "Gene_No", "MetaCyc_gene")) %>%#, "excess")) %>%
  filter(!is.na(product)) %>% #|!is.na(MetaCyc_gene)) 
  group_by(product) %>%
  unique()
#View(MetaCyc2)
head(MetaCyc2)
#View(MetaCyc2)

#write_csv(MetaCyc2, "/home/lgschaer/old/MiniOmics/MetaCyc_files/MetaCyc_PET_to_TCA_enzyme_list_03042022_LGS.csv")
MetaCyc3 <- as.csv("/home/lgschaer/old/MiniOmics/MetaCyc_files/edited_MetaCyc_PET_to_TCA_enzyme_list_03042022_LGS.csv")
head(MetaCyc3)


head(genes_filt)

ssm <- keggGet(c("ec00500"))   #Starch and sucrose metabolism

ssm_e1 <- as.data.frame(ssm[[1]]$ENZYME) %>% 
  mutate(EC_number = ssm[[1]]$ENZYME, 
         General_Pathway = ssm[[1]]$NAME)
head(ssm_e1)

cellulose_genes <- genes_filt %>%
  inner_join(Taxa_Quality_Summary) %>%
  left_join(ssm_e1) %>%
  filter(grepl("Cellulose.", product)|grepl(".cellulose.", product)|grepl(".cellulose", product)|!is.na(General_Pathway))
head(cellulose_genes)
dim(cellulose_genes)

write.csv(cellulose_genes, "/home/lgschaer/old/MiniOmics/clean_data/cellulose_genes_only_03222022_LGS.csv")

TPA_complete_pathway <- genes_filt %>%
  mutate(
    prokka_product = product
  ) %>%
  left_join(MetaCyc3, by = c("prokka_product")) %>%
  filter(!is.na(MetaCyc_Pathway_Name))
View(TPA_complete_pathway)
dim(TPA_complete_pathway)

ggplot(TPA_complete_pathway, aes(x = product, y = BinID))+
  geom_tile(aes(fill = MetaCyc_Pathway_Name))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, size = 11),
        axis.text.y = element_text())

#Making a pathway figure with pathview

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")

library(pathview)
data(gse16873.d)

head(gse16873.d[,1])
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873")
pv.out

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



#####Plotting genes involved in TPA breakdown with taxonomy classification

TPA_genes <- genes_filt %>%
  full_join(Taxa_Quality_Summary) %>%
  filter(product == "Mono(2-hydroxyethyl) terephthalate hydrolase"|
           gene == "tphA1I"|gene == "tphBI_1"|gene == "tphA3I"|
           gene == "tphA1II"|gene == "tphBI"|gene == "tphA2I"|
           gene == "tphA3I"|gene == "tphA3II"|gene == "tphA2I_1"|
           gene == "pcaG"|gene == "pcaG_1"|gene == "pcaG_2"|
           gene == "pcaH"|gene == "pcaH_1"|gene == "pcaH_2"|gene=="tcbC"|
           gene == "alkB1"|gene == "alkB2"|gene=="alkB2") %>%
  mutate(
    gene = ifelse(product == "Mono(2-hydroxyethyl) terephthalate hydrolase", "no_name", gene),
    P1 = ifelse(gene == "alkB1"|gene == "alkB2"|gene=="alkB2", "Alkane", gene),
    P2 = ifelse(P1 == "pcaG"|P1 == "pcaG_1"|P1 == "pcaG_2"|
                  P1 == "pcaH"|P1 == "pcaH_1"|P1 == "pcaH_2", "Protocatechuate", P1),
    P3 = ifelse(P2 == "no_name", "Cleave Ethylene Glycol", P2),
    P4 = ifelse(P3 == "tcbC", "Catechol", P3),
    Pathway = ifelse(P4 == "tphA1I"|P4 == "tphBI_1"|P4 == "tphA3I"|
                       P4 == "tphA1II"|P4 == "tphBI"|P4 == "tphA2I"|
                       P4 == "tphA3I"|P4 == "tphA3II"|P4 == "tphA2I_1", "Terephthalate", P4) 
  ) #%>%
  #unite(Taxa_Name, Enrichment, Taxa_Name, sep = "_", remove = FALSE) %>%
  #select(c("Bin_ID", "Enrichment", "Bin_Number", "Pathway", "product", "gene", "Taxa_Name", "Completeness", "Contamination"))
head(TPA_genes)
#View(TPA_genes)


names_of_genes <- TPA_genes %>%
  dplyr::select(c("gene", "product")) %>%
  group_by(gene, product) %>%
  summarise(
    count = n()
  )
View(names_of_genes)

#TPA_genes$gene <- factor(TPA_genes$gene, c("no_name","tphA1I","tphBI_1","tphA3I","tphA1II","tphBI","tphA2I",
#                "tphA3II","tphA2I_1","tcbC","pcaG","pcaG_1","pcaG_2","pcaH",
#                "pcaH_1","pcaH_2","alkB1","alkB2"))

#making the dendogram - clustered by only TPA genes
TPA_genes_matrix <- as_tibble(TPA_genes) %>%
  dplyr::select(c("Classification", "gene")) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = gene, values_from = n, values_fn = sum, values_fill = 0) %>%
  column_to_rownames(var = "Classification")
head(TPA_genes_matrix)

library(ggdendro)

hc <- hclust(dist(TPA_genes_matrix), "ave")
dend<-ggdendrogram(hc, rotate = TRUE, size = 2)+
  theme(axis.text.x = element_text(angle = 0, size = 11),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank())

colors <- c("firebrick", "orange", "green", "dodgerblue", "purple")

ddata <- dendro_data(hc)
tax_labels <- ddata$labels

#enrichment_colors <- ifelse(TPA_genes$Enrichment=="Emma2","purple", 
#                           ifelse(TPA_genes$Enrichment=="Laura1","red",
#                                 ifelse(TPA_genes$Enrichment=="Laura2","orange",
#                                       ifelse(TPA_genes$Enrichment=="Lindsay2","green",
#                                             ifelse(TPA_genes$Enrichment=="Steve2","grey", "black")))))
#enrichment_colors

heat<-ggplot(TPA_genes, aes(x = gene, y = Classification))+
  geom_tile(aes(fill = Pathway), color = "black") +
  facet_grid(rows = vars(Enrichment), space = "free", scales = "free")+
  scale_fill_manual(values = colors) +
  scale_y_discrete(                                                  #change x-axis labels
    limits = tax_labels$label,
    breaks = tax_labels$label,
    labels = tax_labels$label
  )+ 
#  scale_x_discrete(                                                  #change x-axis labels
 #   limits = gene_order,
  #  breaks = gene_order,
   # labels = gene_order
  #)+ 
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0, size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        axis.title.x = element_blank())
heat


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

head(genes)
head(Taxa_Quality_Summary)

Aromatic_Genes <- genes_filt %>%
  inner_join(Taxa_Quality_Summary) %>%
  inner_join(interesting_enzymes, by = "EC_number") %>%
  dplyr::select(c("BinID", "Enrichment", "product", "gene", "Classification", "Confidence", "Tax_Level", "Quality", "EC_number", "General_Pathway")) %>%
  filter(!is.na(General_Pathway)) %>%
  dplyr::select(EC_number, everything())
head(Aromatic_Genes)

write_csv(Aromatic_Genes, file = "/home/lgschaer/old/MiniOmics/clean_data/aromatic_genes_EC_first_03022022.csv")



colors <- c("firebrick", "orange", "green", "dodgerblue", "purple")

heat<-ggplot(Aromatic_Genes, aes(x = product, y = Classification))+
  geom_tile(aes(fill = General_Pathway), color = "black") +
  facet_grid(rows = vars(Enrichment), shrink = TRUE)+
  scale_fill_manual(values = colors) +
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text.x = element_blank(), #angle = 0, size = 10
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, size = 11),
        axis.title.x = element_blank())
heat


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



