# R script for deconvolution (last updated 10.02.2023)
# Author: Karoliina Salenius karoliina.salenius@tuni.fi
#
# install.packages("BiocManager")
# BiocManager::install("preprocessCore")
#
# Download regressionAnalysis_fun.R from https://github.com/NykterLab/GBM_immune
# Set working directory to directory containing R script
# source('regressionAnalysis_fun.R')
#
#
# data = DESeq2 normalized RNAseq counts from data to be deconvoluted
# reference_data =  DESeq2 normalized RNAseq counts from reference  cell types + median control sample taken from data
#
#
# Quantile normalize the sample counts and reference counts together
library(preprocessCore)
all_counts <- cbind(data, reference_data)
quantile_counts <- normalize.quantiles(as.matrix(all_counts), copy=FALSE)
quantile_counts_log2 <- log2(quantile_counts+1)

# Separate the sample counts from the reference
data_normalized = quantile_counts_log2[,1:ncol(data)]
reference_normalized=quantile_counts_log2[,(ncol(data)+1):ncol(quantile_counts_log2)]

# Run deconvolution with alpha 0.25
proportions <- regressionAnalysis(mixture = data_normalized, reference_cells = reference_normalized, alpha = 0.25)
