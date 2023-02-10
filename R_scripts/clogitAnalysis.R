# R script for Conditional logistic regression for integrative analysis of virome exposures, gene expression, and cell type proportions (last updated 10.02.2023)
# Author: Elaheh Moradi, University of Easatern Finland, Kuopio ,Finland (elaheh.moradi@uef.fi)
#
# Requirements:
# R version 4.1.1 or later
# install.packages("survival")

# Usage: 
# Set working directory to directory containing R script
# source('clogitAnalysis.R')
# fit <- clogitAnalysis(data_, condition, pair, hla)
#
# Parameters: 
# gene: The selected gene
# condition: case, control
# pair: sample pair
# hla: HLA information
#snps: SNPs information
#virus: Virus variable 


# Function for Conditional logistic regression for integrative analysis of virome exposures, gene expression, and cell type proportions

clogitAnalysis= function(gene, pair, condition, hla, snps, virus){
  library(survival)
  table_ <- data.frame(pair,condition, gene,hla,snps,virus)
  fit_ <- clogit(condition ~ gene+strata(pair)+ hla+snps+virus, method="exact", data=table_)
  fit = summary(fit_)
  return (fit)
}