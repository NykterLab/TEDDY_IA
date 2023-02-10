# R script for geneSelection (last updated 10.02.2023)
# Author: Elaheh Moradi, University of Easatern Finland, Kuopio ,Finland (elaheh.moradi@uef.fi)
#
# Requirements:
# R version 4.1.1 or later
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
#
# Usage: 
# Set working directory to directory containing R script
# source('geneselection.R')
# sel_genes <- geneSelection(data_, condition, pair, hla)
#
# Parameters: 
# data_: RNA-seq data (only protein coding genes and after filtering for low count genes)
# condition: case, control
# pair: sample pair
# hla: HLA information

# Function for selection of differentially expressed genes
geneSelection = function(data_, condition, pair, hla){
  
  library(DESeq2)
  coldata<- data.frame(condition = condition, pair= pair, hla= hla)
  rownames(coldata)<- colnames(data_)
  coldata$condition <- factor(coldata$condition)
  coldata$hla <- factor(coldata$hla)
  coldata$pair <- factor(coldata$pair)
  dds <- DESeqDataSetFromMatrix(countData = data_, colData = coldata, design = ~ pair+hla+condition)
  dds <- DESeq(dds)
  res<- results(dds)
  sel_genes= rownames(res)[which(res$padj< 0.05 & abs(res$log2FoldChange) > 0.5)]
  return(sel_genes)
}