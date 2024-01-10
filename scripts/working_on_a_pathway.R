
#packages needed
library(tidyverse)
library(csv)
library(readr)

#move all tsv files to be analyzed from the prokka directory to the same directory, set that directory as wd

#setwd("/home/lgschaer/old/MiniOmics/Metagenomics/prokka_output/")
#getwd()

#file_names <- dir("./") #where you have your files
#file_names

#genes <- do.call(rbind,lapply(file_names,read_tsv))
#head(genes)
#dim(genes)

#write_csv(genes, file = "/home/lgschaer/old/MiniOmics/clean_data/L1_E2_genes_unfiltered.csv")
genes<- as.csv("/home/lgschaer/old/MiniOmics/clean_data/L1_E2_genes_unfiltered.csv")
head(genes)

genes_filt <- genes %>%
  filter(product != "hypothetical protein") %>%
  separate(locus_tag, into = c("BinID", "geneID"), sep = "_") %>%
  mutate(BinID2 = BinID) %>%
  separate(BinID2, into = c("Enrichment", "Bin_Number"), sep = "metabat.")
head(genes_filt)
dim(genes_filt)

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
MetaCyc <- read.table("/home/lgschaer/old/MiniOmics/MetaCyc_files/MetaCyc_PET_to_TCA_full_pathways_w_genes_LGS_03022022.txt", header = TRUE, sep = "\t", dec = "_")
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


dim(genes_filt)

TPA_complete_pathway <- genes_filt %>%
  mutate(
    prokka_product = product
  ) %>%
  left_join(MetaCyc3, by = c("prokka_product")) %>%
  filter(!is.na(MetaCyc_Pathway_Name))
head(TPA_complete_pathway)
dim(TPA_complete_pathway)

ggplot(TPA_complete_pathway, aes(x = product, y = BinID))+
  geom_tile(aes(fill = MetaCyc_Pathway_Name))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, size = 11),
        axis.text.y = element_text())

#load data
full_count_data_w_quality_taxa <- read.csv(file = "/home/lgschaer/old/MiniOmics/clean_data/full_count_data_w_quality_taxa.csv", sep = ',') %>%
  mutate(
    BinID = as.character(BinID)
  ) %>%
  dplyr::select(c(locus_tag, TPA, DCPET, HDPE, TA, EG, product, Enrichment, General_Pathway, BinID, Best_Classification, Completeness, Contamination, Quality))
head(full_count_data_w_quality_taxa)
dim(full_count_data_w_quality_taxa)


TPA_complete_pathway_w_counts2 <- full_count_data_w_quality_taxa %>% 
  filter(TPA>0|DCPET>0|TA>0|HDPE>0|EG>0) %>%
  group_by(BinID, product, Enrichment, Best_Classification) %>%
  summarise(
    TPA = sum(TPA),
    DCPET = sum(DCPET), 
    TA = sum(TA),
    EG = sum(EG),
    HDPE = sum(HDPE)
  ) %>%
  filter(!is.na(TPA)|!is.na(TA)|!is.na(DCPET|!is.na(EG)|!is.na(HDPE))) %>%
  mutate(
    TPA = ifelse(is.na(TPA), 0, TPA),
    DCPET = ifelse(is.na(DCPET), 0, DCPET),
    HDPE = ifelse(is.na(HDPE), 0, HDPE),
    TA = ifelse(is.na(TA), 0, TA),
    EG = ifelse(is.na(EG), 0, EG)
  ) %>%
  separate(Best_Classification, into = c("Taxa", "Similarity", "Taxa_Level"), sep = "_")%>%
  pivot_longer(cols = c(TPA, DCPET, HDPE, EG, TA), names_to = "Substrate", values_to = "Counts") %>%
  filter(Counts != 0) %>%
  group_by(BinID, Substrate) %>%
  mutate(expressed_genes_per_bin = sum(Counts))
head(TPA_complete_pathway_w_counts2)
#View(TPA_complete_pathway_w_counts2)

dim(TPA_complete_pathway_w_counts)


unique(TPA_complete_pathway_w_counts$MetaCyc_Pathway_Name)

pathway_names <- c(
  "3-oxoadipate degradation"   ="3-Oxoadipate\nDegradation" ,                           
  "TCA cycle I (prokaryotic)"  ="TCA Cycle\n(prokaryotic)",                             
  "terephthalate degradation"  ="Terephthalate\nDegradation",                              
  "protocatechuate degradation II (ortho-cleavage pathway)" ="PCA Degradation\n(ortho-cleavage pathway)", 
  "protocatechuate degradation I (meta-cleavage pathway)"   ="PCA Degradation\n(meta-cleavage pathway)" , 
  "2-hydroxypenta-2,4-dienoate degradation"  ="2-Hydroxypenta\n2,4-Dienoate\nDegradation" ,               
  "glycolate and glyoxylate degradation I"   ="Glyoxylate\nDegradation I" ,                
  "glycolate and glyoxylate degradation II" ="Glyoxylate\nDegradation II" 
)

unique(TPA_complete_pathway_w_counts$product)

product_names <- c(
  "3-ketoacyl-CoA thiolase" =  "3-ketoacyl-CoA\nthiolase",                                                   
  "Isocitrate dehydrogenase [NADP]" = "Isocitrate\ndehydrogenase",                                          
  "Succinate dehydrogenase flavoprotein subunit"  = "Succinate\ndehydrogenase",                             
  "Terephthalate 1,2-dioxygenase, terminal oxygenase component subunit alpha 1" = "TPA\n1,2-dioxygenase\nalpha subunit" ,
  "Terephthalate 1,2-dioxygenase, terminal oxygenase component subunit beta 1"  = "TPA\n1,2-dioxygenase\nbeta subunit" ,
  "Quinone oxidoreductase 2"  = "Quinone\noxidoreductase",                                                      
  "3-oxoadipate enol-lactonase 2" = "3-oxoadipate\nenol-lactonase",                                             
  "4-carboxy-4-hydroxy-2-oxoadipate aldolase" = "4-carboxy-4-hydroxy\n2-oxoadipate\naldolase",                                 
  "Protocatechuate 4,5-dioxygenase alpha chain" =  "PCA\n4,5-dioxygenase\nalpha chain",                              
  "Protocatechuate 3,4-dioxygenase alpha chain"  =  "PCA\n3,4-dioxygenase\nalpha chain",                             
  "Citrate synthase"  =  "Citrate\nsynthase",                                                        
  "Acetaldehyde dehydrogenase"  =  "Acetaldehyde\ndehydrogenase",                                              
  "Aconitate/2-methylaconitate hydratase"  = "Aconitate\n2-methylaconitate\nhydratase",                                    
  "Glycerate 2-kinase" = "Glycerate\n2-kinase",                                                        
  "Ureidoglycolate dehydrogenase (NAD(+))" = "Ureidoglycolate\ndehydrogenase",                                    
  "4-carboxy-2-hydroxymuconate-6-semialdehyde dehydrogenase" = "4-carboxy-2\nhydroxymuconate\n6-semialdehyde\ndehydrogenase"
)

length(unique(TPA_complete_pathway_w_counts_filt$Taxa))
length(taxa_colors)

cl <- colors(distinct = TRUE)
set.seed(11031993) # to set random generator seed
pre_taxa_colors <- sample(cl, 45)
pre_taxa_colors

taxa_colors <- c(
  "firebrick4",           "orchid1",              "green",                "cyan",
  "cyan",                 "darkgreen",            "palegoldenrod",        "tan4",
  "darkblue",             "darkcyan",               
  "dodgerblue",           "firebrick",            "magenta",
  "limegreen",              "chocolate4",           "skyblue2",             "plum1",           
  "plum1",                "orange3",              "palevioletred2",      
  "mediumpurple3",        "cadetblue4",           "grey47",                "darkorange4",              
  "purple3",              "darkolivegreen",                  
  "tomato1",              "blue4",                "pink4",                          
  "chartreuse4",          "yellow",              "yellow",               "orange",               
  "gold3",                "hotpink4",                 
  "hotpink4",             "goldenrod4",           "goldenrod4",           "steelblue",           
  "deepskyblue",          "orchid4"
)

head(TPA_complete_pathway_w_counts)

TPA_complete_pathway_w_counts_filt <- TPA_complete_pathway_w_counts %>%
  filter(Counts > 0) %>%
  mutate(Taxa = ifelse(is.na(Taxa), "Unclassified", Taxa)) %>%
  group_by(Enrichment, product, Taxa, General_Pathway) %>%
  summarise(Counts = sum(Counts))
head(TPA_complete_pathway_w_counts_filt)

ggplot(TPA_complete_pathway_w_counts_filt, aes(x = product, y = Counts, fill = Taxa))+
  geom_col(position = "fill", color = "black")+
  ylab("Percent of TPA Pathway Counts")+
  #scale_y_log10()+
  scale_x_discrete(labels = product_names)+
  #scale_fill_manual(values = taxa_colors)+
  facet_grid(cols = vars(General_Pathway), rows = vars(Enrichment), scales = "free", space = "free", 
             labeller = labeller(General_Pathway = pathway_names))+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5), 
    legend.position = "bottom"
  )
  
TPA_complete_pathway_w_counts3 <- TPA_complete_pathway %>%
  left_join(TPA_complete_pathway_w_counts2, by = c("BinID", "product", "Enrichment")) %>%
  filter(!is.na(Substrate)) %>%
  group_by(BinID, Substrate) %>%
  mutate(
    percent_expressed_gene_per_bin = Counts/expressed_genes_per_bin
  )
head(TPA_complete_pathway_w_counts3)
#View(TPA_complete_pathway_w_counts3)


#write_csv(TPA_complete_pathway_w_counts, "/home/lgschaer/old/MiniOmics/MetaCyc_files/TPA_complete_pathway_w_counts_03072022_LGS.csv")


colors <- c("blue", "firebrick", "dodgerblue", "green", "lightgoldenrod", "purple4",
            "orange", "magenta", "green4", "lightblue", "tan4", "yellow", "lightpink",
            "purple1", "red", "grey")

#counts log transformed
ggplot(TPA_complete_pathway_w_counts3, aes(y = Taxa, x = Counts))+
  facet_grid(cols = vars(Substrate), rows = vars(Enrichment), scales = "free")+
  geom_col(aes(fill = metacyc_product), color = "black")+
  scale_fill_manual(values = colors)+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, size = 11),
        axis.text.y = element_text())

#counts as a percentage of the total
ggplot(TPA_complete_pathway_w_counts3, aes(y = Taxa, x = percent_expressed_gene_per_bin))+
  facet_grid(cols = vars(Substrate), rows = vars(Enrichment))+
  geom_col(aes(fill = metacyc_product), color = "black")+
  scale_fill_manual(values = colors)+
  xlab("Count of Gene/Total Expressed Genes")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, size = 11),
        axis.text.y = element_text())

#counts as a percentage of the total, log transformed
ggplot(TPA_complete_pathway_w_counts3, aes(y = Taxa, x = percent_expressed_gene_per_bin))+
  facet_grid(cols = vars(Substrate), rows = vars(Enrichment))+
  geom_col(aes(fill = metacyc_product), color = "black")+
  scale_fill_manual(values = colors)+
  scale_x_log10()+
  xlab("Log10(Count of Gene/Total Expressed Genes)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, size = 11),
        axis.text.y = element_text())


#making alternate versions of above plot

##Laura1 Only
taxa_to_include <- c("Ancylobacter",                      "Phyllobacteriaceae",               
                     "Rhodospirillales",                  "Alcaligenaceae",                    "Pseudomonas",                      
                     "Parapedobacter",                   
                     "Pelagibacterium",                   "Betaproteobacteria",                "Achromobacter",                    
                     "Hyphomicrobiaceae",                 "Luteimonas",                        "Comamonadaceae",                   
                     "Alphaproteobacteria",               "Rhizobiales",                      
                     "Hyphomonas",                        "Proteobacteria",                    "Shinella",                         
                     "Rhodococcus",                       
                     "Brevundimonas",                     "Sphingobacteriaceae",               "Hydrogenophaga",                   
                     "Actinobacteria",                   
                     "Chelatococcus",                     "Acidobacteria",                     "Sphingobacterium"                
                     )

TPA_complete_pathway_w_counts_Laura1 <- TPA_complete_pathway_w_counts3 %>%
  filter(Enrichment=="Laura1") %>%
  filter(!is.na(Taxa)) %>%
  filter(Taxa %in% taxa_to_include) %>%
  filter(Substrate != "HDPE") %>%
  mutate(Substrate = as.factor(Substrate, c("DCPET", "TPA", "EG", "TA")))

unique(TPA_complete_pathway_w_counts_Laura1$Taxa)

product_labels <- c(
  "LigC" = "LigC",                                                                
  "subunit of 3-oxoadipate enol-lactone hydrolase" = "3-oxoadipate enol-lactone hydrolase",                       
  "LigK" = "LigK",                                                                   
  "citrate synthase" = "citrate synthase",                                                   
  "succinate dehydrogenase" ="succinate dehydrogenase",                                              
  "acetaldehyde dehydrogenase" = "acetaldehyde dehydrogenase",                                           
  "isocitrate dehydrogenase (NADP-dependent)" = "isocitrate dehydrogenase",                           
  "beta-ketoadipyl CoA thiolase" = "beta-ketoadipyl CoA thiolase",                                        
  "malate:quinone oxidoreductase"  = "quinone oxidoreductase",                                       
  "glycolate dehydrogenase"  = "glycolate dehydrogenase",                                             
  "glycerate 2-kinase 2"   =   "glycerate 2-kinase 2",                                             
  "bifunctional aconitate hydratase B and 2-methylisocitrate dehydratase" = "2-methylisocitrate dehydratase",
  "protocatechuate 3,4-dioxygenase"    = "PCA 3,4-dioxygenase",                                   
  "terephthalate 1,2-dioxygenase oxygenase component" = "TPA 1,2-dioxygenase",                    
  "protocatechuate 4,5-dioxygenase"  = "PCA 4,5-dioxygenase")

ggplot(TPA_complete_pathway_w_counts_Laura1, aes(y = Taxa, x = Counts))+
  facet_grid(cols = vars(Substrate), rows = vars(Enrichment), scales = "free")+
  geom_col(aes(fill = metacyc_product), color = "black")+
  scale_fill_manual(values = colors, labels = product_labels)+
  scale_x_log10()+
  xlab("Log of HTSeq Counts")+
  scale_y_discrete(limits = rev)+
  theme_linedraw()+
  theme(axis.text.x = element_text(size = 45, angle = 0, color = "black"),
        axis.title.x = element_text(size = 50, face = "bold"),
        axis.title.y =  element_blank(),
        axis.text.y = element_text(size = 55, color = "black"),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 75),
        legend.title = element_blank(),
        legend.key.height = unit(15, "mm"),
        legend.key.width = unit(15, "mm"),
        legend.text = element_text(size = 50))

#Only Rhodococcus, Hydrogenophaga, Pseudomonas, Chelatococcus, and Alcaligenaceae
TPA_complete_pathway_w_counts_specific_taxa <- TPA_complete_pathway_w_counts %>%
  filter(Enrichment=="Laura1") %>%
  filter(Taxa=="Rhodococcus"|Taxa=="Hydrogenophaga"|Taxa=="Pseudomonas"|Taxa=="Chelatococcus"|Taxa=="Alcaligenaceae")

ggplot(TPA_complete_pathway_w_counts_specific_taxa, aes(y = Taxa, x = Counts))+
  facet_grid(cols = vars(Substrate), rows = vars(Enrichment), scales = "free")+
  geom_col(aes(fill = metacyc_product), color = "black")+
  scale_fill_manual(values = colors)+
  scale_x_log10()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 0, size = 12, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = rel(1.5)),
        #legend.title = element_blank(),
        legend.text = element_text(size = 11))

#make a bar chart with complete gene counts
TPA_complete_pathway_w_counts <- full_count_data_w_quality_taxa %>%
  filter(TPA>0|DCPET>0|TA>0|HDPE>0|EG>0) %>%
  group_by(BinID, product, Enrichment, Best_Classification, General_Pathway) %>%
  summarise(
    TPA = sum(TPA),
    DCPET = sum(DCPET), 
    TA = sum(TA),
    EG = sum(EG),
    HDPE = sum(HDPE)
  ) %>%
  filter(!is.na(TPA)|!is.na(TA)|!is.na(DCPET|!is.na(EG)|!is.na(HDPE))) %>%
  mutate(
    TPA = ifelse(is.na(TPA), 0, TPA),
    DCPET = ifelse(is.na(DCPET), 0, DCPET),
    HDPE = ifelse(is.na(HDPE), 0, HDPE),
    TA = ifelse(is.na(TA), 0, TA),
    EG = ifelse(is.na(EG), 0, EG)
  ) %>%
  separate(Best_Classification, into = c("Taxa", "Similarity", "Taxa_Level"), sep = "_")%>%
  pivot_longer(cols = c(TPA, DCPET, HDPE, EG, TA), names_to = "Substrate", values_to = "Counts")
head(TPA_complete_pathway_w_counts)

#write_csv(TPA_complete_pathway_w_counts, "/home/lgschaer/old/MiniOmics/MetaCyc_files/TPA_complete_pathway_w_counts_03072022_LGS.csv")


#colors <- c("blue", "firebrick", "dodgerblue", "green", "lightgoldenrod", "purple4",
 #           "orange", "magenta", "green4", "lightblue", "tan4", "yellow", "lightpink",
  #          "purple1", "red", "grey")

ggplot(TPA_complete_pathway_w_counts, aes(y = Taxa, x = Counts))+
  facet_grid(cols = vars(Substrate), rows = vars(Enrichment), scales = "free")+
  geom_col(fill = "grey", color = "black")+
  scale_x_log10()+
  #scale_fill_manual(values = colors)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, size = 11),
        axis.text.y = element_text())

