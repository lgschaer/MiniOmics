#Making a phylogenetic tree from metagenome data

#Load and install packages
#install.packages("adegenet")
library(adegenet)
library(tidyverse)
library(csv)

#load fasta file
dna <- fasta2DNAbin(file="/home/lgschaer/old/MiniOmics/Metagenomics/Laura1/Laura1/Contig_Bins/Laura1.contigs.fa")
head(dna)
