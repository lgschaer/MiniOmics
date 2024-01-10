#load packages
library(csv)
library(tidyverse)
library(phyloseq)
library(DESeq2)


######### Section 1: Prepare the large data frame with transcriptomics data for phyloseq.
#                    Divide the data into three tables:
#                    1. A metadata table with information about each bin: enrichment conditions, taxa classification, etc
#                    2. A count table with bins names + substrate as rows and genes as columns
#                    3. A gene table (replaces the taxa table from a traditional ps object) that describes each gene, 
#                       expected product, EC number, etc.
#load data
full_count_data_w_quality_taxa <- read.csv(file = "/home/lgschaer/old/MiniOmics/clean_data/full_count_data_w_quality_taxa.csv", sep = ',') 
head(full_count_data_w_quality_taxa)
#View(full_count_data_w_quality_taxa)

sum(duplicated(full_count_data_w_quality_taxa$geneID))

#metadata
sdata <- full_count_data_w_quality_taxa %>%
  filter(product != "hypothetical protein" & gene != "") %>%
  dplyr::select(BinID, TPA, DCPET, HDPE, TA, EG, Enrichment, Best_Classification, Confidence, Completeness, Contamination, Quality) %>%
  pivot_longer(cols = c("TPA", "DCPET", "HDPE", "TA", "EG"), names_to = "Substrate", values_to = "Counts") %>%
  filter(Counts > 0) %>%
  unite(BinID_Substrate, BinID, Substrate, sep = "_", remove = FALSE) %>%
  mutate(BinID_Substrate2 = BinID_Substrate) %>%
  select(-c("Counts")) %>%
  distinct()
dim(sdata) 
head(sdata)
any(duplicated(sdata$BinID_Substrate))

sdata3 <- sdata %>%
  column_to_rownames(var = "BinID_Substrate2")
head(sdata3)

#load count table
sequence_table <- full_count_data_w_quality_taxa %>%
  filter(product != "hypothetical protein" & gene != "") %>%
  dplyr::select(BinID, TPA, DCPET, HDPE, TA, EG, product) %>%
  pivot_longer(cols = c("TPA", "DCPET", "HDPE", "TA", "EG"), names_to = "Substrate", values_to = "Counts") %>%
  #mutate(
   # gene_prefix = "Gene",
    #geneID = as.character(geneID)
  #)%>%
  #unite(geneID, gene_prefix, geneID) %>%
  filter(Counts > 0) %>%
  unite(BinID_Substrate, BinID, Substrate, sep = "_") %>%
  group_by(BinID_Substrate, product) %>%
  summarise(
    Counts = sum(Counts)
  ) %>%
  pivot_wider(names_from = "product", values_from = "Counts", values_fill = 0) %>%
  column_to_rownames(var = "BinID_Substrate")
#head(sequence_table)
#View(sequence_table)
sequence_table[1:5,1:5]

any(is.na(sequence_table$gene))
any(duplicated(sequence_table$newcol))

# remove column names
colnames(sequence_table) <- NULL
sequence_table[1:5,1:5]

# convert to data frame
sequence_table <- as.data.frame(sequence_table)
class(sequence_table)

#dim(seq_table)
dim(sdata)
dim(sequence_table)

#rename
subset_all <- sequence_table

#remove column names from the sequence table
colnames(subset_all) <- NULL
subset_all[1:5,1:5]

#load "taxa" table (gene information)
taxa <- full_count_data_w_quality_taxa %>%
  dplyr::select(product, EC_number, contigName) %>%
  #mutate(
   # gene_prefix = "Gene",
    #geneID = as.character(geneID)
  #)%>%
  #unite(geneID, gene_prefix, geneID) #%>%
distinct()# %>%
#mutate(geneID2 = geneID) %>%
# column_to_rownames(var = "geneID")
head(taxa)
dim(taxa)


####### Section 2: making a phyloseq object: bring the tables together into one phyloseq object

#make phyloseq object
samdata = sample_data(sdata3)
seqtab = otu_table(subset_all, taxa_are_rows = FALSE)
taxtab = tax_table(taxa)

sample_names(samdata)
sample_names(seqtab)
sample_names(taxtab)

taxa_names(seqtab)
taxa_names(taxtab)

#combine all components into a phyloseq object
ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
ps

head(sample_data(ps))

#save the phyloseq object as an RDS object
#write_rds(ps, "/home/lgschaer/old/MiniOmics/clean_data/ps_transcripts_grouped_by_product_02012022.rds")

########### Section 3: Now we can do some analysis!

head(sdata3)

sample_colors <- c("purple4", "green", "orange", "dodgerblue", "firebrick")

#violin plot
ps %>%                                                              #phyloseq object
  plot_richness(
    x = "Substrate",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Substrate), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+   #set fill colors
  # scale_x_discrete(                                                  #change x-axis labels
  # breaks = sample_types)+                   
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

### Section 3: Normalize the data

#samplesover1000_all <- subset_samples(ps, sample_sums(ps) > 100)
#any(taxa_sums(samplesover1000_all) == 0)
#sum(taxa_sums(samplesover1000_all) == 0)
#prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
#any(taxa_sums(prune_samplesover1000_all) == 0)

#for reproducible data
#set.seed(81)
#rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))
#pps<- rarefy_samplesover1000_all
#pps

#PCoA plot
#Making a PCoA plot
#ordination
all_pcoa <- ordinate(
  physeq = ps, 
  method = "PCoA", 
  distance = "bray"
)

head(sdata)

colors <- c("purple", "orange", "blue", "red", "white")

#plot
plot_ordination(
  physeq = ps,                                                         #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Substrate, shape = Enrichment), size = 5) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command


##gene abundance plot
head(sample_data(ps))
head(tax_table(ps))

# Convert phyloseq object into a data frame with relative abundance counts
abundance_by_contig <- ps %>%
  tax_glom(taxrank = c("ta1", "ta2")) %>%                     # agglomerate at phylum level
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(ta1) 
head(abundance_by_contig)
dim(abundance_by_contig)

abundance_by_contig_filt <- abundance_by_contig %>%
  filter(Abundance > 0) %>%
  dplyr::mutate(gene = ifelse(Abundance < 0.5, "< 50%", ta1))
dim(abundance_by_contig_filt)

hist(abundance_by_contig_filt$Abundance)

# Save phylum colors
gene_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "blue",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "blue",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "blue",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

# Create gene plot!
head(abundance_by_contig_filt)

ggplot(abundance_by_contig_filt)+
  geom_col(mapping = aes(x = Substrate, y = Abundance, fill = gene), position = "fill", show.legend = TRUE, color = "black")+
  #facet_grid(rows = vars(Environment), cols = vars(Transfer))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = gene_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))


### DESeq2: Differential Abundance Analysis

#install DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

#packages used
library(tidyverse)
library(DESeq2)
library(phyloseq)


#phyloseq object
ps

#subset phyloseq object to make sample comparison categories, Media_Carbon = DCPET_BH or TPA_BH
head(sample_data(ps))
#unique(sample_info$Substrate)


#my_data[["RNA"]]@counts <- as.matrix(my_data[["RNA"]]@counts)+1
ps@otu_table <- as.matrix(ps@otu_table)+1

#phyloseq to Deseq but for the whole phyloseq
ps_deseq = phyloseq_to_deseq2(ps, ~Substrate)
ps_deseq2 <- ps_deseq
ps_deseq3 <- ps_deseq
ps_deseq4 <- ps_deseq
ps_deseq$type<-relevel(ps_deseq$Substrate, "TPA")
ps_deseq2$type<-relevel(ps_deseq2$Substrate, "TA")
ps_deseq3$type<-relevel(ps_deseq3$Substrate, "DCPET")
ps_deseq4$type<-relevel(ps_deseq4$Substrate, "EG")


#run DESeq
dds<-DESeq(ps_deseq)
dds2<-DESeq(ps_deseq2)
dds3<-DESeq(ps_deseq3)
dds4<-DESeq(ps_deseq4)

resultsNames(dds)  # This shows you all of the comparisons that have been performed
resultsNames(dds2) 
resultsNames(dds3)
resultsNames(dds4)

col_fec <- results(dds, contrast=c("type","colon","feces"))
small_fec <- results(dds, contrast=c("type","small","feces"))
small_col <- results(dds2, contrast = c("type", "small", "colon"))


#This command takes the significance table and appends the taxonomy assignments from the phyloseq object
sigdsqtax_ta_tpa = cbind(as(sigdsq_ta_tpa, "data.frame"), as(tax_table(justbacteria)[rownames(sigdsq_ta_tpa), ], "matrix"))
sigdsqtax_ta_c = cbind(as(sigdsq_ta_c, "data.frame"), as(tax_table(justbacteria)[rownames(sigdsq_ta_c), ], "matrix"))
sigdsqtax_tpa_c = cbind(as(sigdsq_tpa_c, "data.frame"), as(tax_table(justbacteria)[rownames(sigdsq_tpa_c), ], "matrix"))








#This command creates a new file of just the ASVs that are significant
#remove hash marks in this block of code to get a df with just the significant ASVs.
sigdsq_L1_E2_results_b = L1_E2_results_b[which(L1_E2_results_b$padj < alpha), ]
head(sigdsq_L1_E2_results_b)

#This command takes the significance table and appends the taxonomy assignments from the phyloseq object
sigdsqtax_L1_E2 = cbind(as(sigdsq_L1_E2_results_b, "data.frame"), as(tax_table(pps)[rownames(sigdsq_L1_E2_results_b), ], "matrix"))

#This command outputs the signficance table as a .csv
write.csv(sigdsqtax_L1_E2,"/home/lgschaer/old/MiniOmics/clean_data/dsq_out/L1_E2_enriched_and_nonenriched_w_tax.csv")

#make a volcano plot.

#add column with labels for significantly enriched ASVs.
sigdsqtax_L1_E2$threshold = as.factor(abs(sigdsqtax_L1_E2$log2FoldChange) > 2 & sigdsqtax_L1_E2$padj < alpha)

#rename columns enriched or not enriched
library(tidyverse)

L1_E2 <- sigdsqtax_L1_E2 %>%
  mutate(
    Comparison = "Laura1 vs. Emma2",
    threshold = as.factor(ifelse(log2FoldChange > 2 & padj < 0.01 | log2FoldChange < -2 & padj < 0.01, "Enriched", "Not_Enriched")),
    ta2 = as.character(ta2),
    Product = ifelse(threshold == "Enriched", ta2, "Not_Enriched"))#,
#Where_Enriched = ifelse(log2FoldChange > 2 & padj < 0.01, "Terephthalate", 
# ifelse(log2FoldChange < -2 & padj < 0.01, "Terephthalamide", "Not_Enriched"))) %>%
#filter(threshold == "Enriched" | threshold == "Not_Enriched")
head(L1_E2)

######

#enriched <- ta_tpa %>% 
# full_join(ta_c) %>%
#full_join(tpa_c) %>%
#  mutate(
#   Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified", Enriched_Genus)
#)
#head(enriched)

#dim(ta_tpa)
#dim(ta_c)
#dim(tpa_c)

#dim(enriched)

#write.csv(enriched,"/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/dsq_out/enriched_and_notenriched_all3comparisons.csv")

#This is a color palette that I like
cbbPalette <- c("red3","lightblue","darkgray", "black","#009E73","#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colors11 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   "tan1", "purple1",
  "grey77",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "orange",  "darkblue",      "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue", "green4", 
  "grey50", "#009E73","#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "white"
) 

#volcano plot:
PlotA <- ggplot(data=L1_E2, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Product), color = "black", size=4, shape = 21) +
  facet_grid(cols = vars(Comparison))+
  # xlim(c(-30, 30)) + 
  # ylim(c(0, 40)) +
  scale_fill_manual(values=colors11) +
  #scale_shape_manual(values=c(21,22,23)) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  #ggtitle("Differentially Abundant Taxa")+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom",
        title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))
PlotA

##### Additional DESeq comparisons for substrate

#phyloseq object
pps

#subset phyloseq object to make sample comparison categories, Media_Carbon = DCPET_BH or TPA_BH
head(sample_data(pps))
#unique(sample_info$Substrate)

TPA_DCPET <- subset_samples(pps, Substrate == "TPA" | Substrate == "DCPET")

##This command takes a phyloseq object (phyloseq_all) and perpares it for DESeq analysis by the category "sampletype"
dsq_TPA_DCPET = phyloseq_to_deseq2(TPA_DCPET, ~Substrate)

##This command sets the first sample in the TYPE category to be one that OIL as its class
###(type = sampletype, oil = openwater??)
dsq_TPA_DCPET$Substrate <- relevel(dsq_TPA_DCPET$Substrate, "TPA")

##This command converts the column data in the DESeq object as a data frame.
as.data.frame(colData(dsq_TPA_DCPET))


## this block of commands will calculate the geometric means for the count data in your samples
#in order to prepare the data for analysis and normalization in DESeq
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dsq_TPA_DCPET), 1, gm_mean)
dsq_TPA_DCPET = estimateSizeFactors(dsq_TPA_DCPET, geoMeans = geoMeans)

##This command runs the DESeq analysis (SLOW)
TPA_DCPET_results = DESeq(dsq_TPA_DCPET, fitType="local")

##This command outputs the results from DESeq
TPA_DCPET_results_b = results(TPA_DCPET_results, cooksCutoff = FALSE)



#This command sets the significance cut off to be 0.01
alpha = 0.01

#This command creates a new file of just the ASVs that are significant
#remove hash marks in this block of code to get a df with just the significant ASVs.
sigdsq_TPA_DCPET_results_b = TPA_DCPET_results_b[which(TPA_DCPET_results_b$padj < alpha), ]
head(sigdsq_TPA_DCPET_results_b)

#This command takes the significance table and appends the taxonomy assignments from the phyloseq object
sigdsqtax_TPA_DCPET = cbind(as(sigdsq_TPA_DCPET_results_b, "data.frame"), as(tax_table(pps)[rownames(sigdsq_TPA_DCPET_results_b), ], "matrix"))

#This command outputs the signficance table as a .csv
write.csv(sigdsqtax_TPA_DCPET,"/home/lgschaer/old/MiniOmics/clean_data/dsq_out/L1_E2_enriched_and_nonenriched_w_tax.csv")

#make a volcano plot.

#add column with labels for significantly enriched ASVs.
sigdsqtax_TPA_DCPET$threshold = as.factor(abs(sigdsqtax_TPA_DCPET$log2FoldChange) > 2 & sigdsqtax_TPA_DCPET$padj < alpha)

#rename columns enriched or not enriched
library(tidyverse)

TPA_DCPET <- sigdsqtax_TPA_DCPET %>%
  mutate(
    Comparison = "DCPET vs. HDPE",
    threshold = as.factor(ifelse(log2FoldChange > 2 & padj < 0.01 | log2FoldChange < -2 & padj < 0.01, "Enriched", "Not_Enriched")),
    ta2 = as.character(ta2),
    Product = ifelse(threshold == "Enriched", ta2, "Not_Enriched"))#,
#Where_Enriched = ifelse(log2FoldChange > 2 & padj < 0.01, "Terephthalate", 
# ifelse(log2FoldChange < -2 & padj < 0.01, "Terephthalamide", "Not_Enriched"))) %>%
#filter(threshold == "Enriched" | threshold == "Not_Enriched")
head(TPA_DCPET)

######

#enriched <- ta_tpa %>% 
# full_join(ta_c) %>%
#full_join(tpa_c) %>%
#  mutate(
#   Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified", Enriched_Genus)
#)
#head(enriched)

#dim(ta_tpa)
#dim(ta_c)
#dim(tpa_c)

#dim(enriched)

#write.csv(enriched,"/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/dsq_out/enriched_and_notenriched_all3comparisons.csv")

#This is a color palette that I like
cbbPalette <- c("red3","lightblue","darkgray", "black","#009E73","#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colors11 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   "tan1", "purple1",
  "grey77",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "orange",  "darkblue",      "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue", "green4", 
  "grey50", "#009E73","#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "white"
) 

#volcano plot:
PlotA <- ggplot(data=TPA_DCPET, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Product), color = "black", size=4, shape = 21) +
  facet_grid(cols = vars(Comparison))+
  # xlim(c(-30, 30)) + 
  # ylim(c(0, 40)) +
  scale_fill_manual(values=colors11) +
  #scale_shape_manual(values=c(21,22,23)) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  #ggtitle("Differentially Abundant Taxa")+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom",
        title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))
PlotA

#ggarrange(PlotA,PlotB, align = "hv", common.legend = FALSE)

head(ta_tpa)

ta_tpa2 <- ta_tpa %>%
  group_by(Where_Enriched)%>%
  mutate(
    N_total = n(),
    Class = as.character(Class),
    Class = ifelse(is.na(Class), "Unclassified", Class)
  ) %>%
  group_by(Class) %>%
  mutate(
    N = n()/N_total,
    Class = ifelse(N<0.01, "< 1%", Class)
  )
head(ta_tpa2)

colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "coral1",        "blue",
  "grey47",  "cyan",         "yellow",    "darkgreen",   "palegoldenrod",     "tan4",
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1",    "purple4",
  "white",   "dodgerblue",   "firebrick", "yellowgreen", "magenta"
)  

ggplot(ta_tpa2, aes(x = Where_Enriched))+
  geom_bar(position = "fill", aes(fill = Class), color = "black", show.legend = TRUE)+
  scale_fill_manual(values = colors10)+
  xlab(NULL)+
  ylab("Percent of ASVs")+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 15, face = "bold", color = "black"))


