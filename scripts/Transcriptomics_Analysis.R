#load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")
BiocManager::install("readr")

library(tximport)
library(GenomicFeatures)
library(readr)
library(tidyverse) 
library(csv)
library(DESeq2)

samples <- read.csv("/home/lgschaer/old/MiniOmics/clean_data/samples.csv", header = TRUE)
head(samples)

dir <- system.file("extdata", package = "tximportData")
files <- file.path(dir, "/home/lgschaer/old/MiniOmics/Transcriptomics/Laura1/quants", samples$sample, "quant.sf")
names(files) <- paste0(samples$sample)
txi.inf.rep <- tximport(files, type = "salmon", txOut = TRUE)
names(txi.inf.rep)
head(txi.inf.rep$counts)

sampleTable <- data.frame(condition = factor(rep(c("DCPET", "EG", "TPA", "TA", "HDPE"), each = 1)))
rownames(sampleTable) <- colnames(txi.inf.rep$counts)
head(sampleTable)

dds <- DESeqDataSetFromTximport(txi = txi.inf.rep, sampleTable, ~condition)
dds

#this command gives an error, https://support.bioconductor.org/p/134505/ suggests that it is because
#there is only one sample in each category so it can't do a statistical comparison.
dds_DESeq <- DESeq(dds, test="Wald", fitType="parametric")

