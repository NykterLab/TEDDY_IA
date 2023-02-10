# R script for Temporal filtering of differentially expressed genes (last updated 10.02.2023)
# Author: Elaheh Moradi, University of Easatern Finland, Kuopio ,Finland (elaheh.moradi@uef.fi)
#
# Requirements:
# R version 4.1.1 or later

#
# Usage: 
# Set working directory to directory containing R script
# source('TemporalFiltering_geneSelection.R')
# sel_genes <- TemporalFiltering_geneSelection(Xdata, sel_genes, sample_ids)
#
# Parameters: 
#sel_genes: selected genes for analaysis 
#Xdata= DESeq2 normalized gene count data
#sample_ids: a list of case and control sample_ids for different time points (0,3,6,9,12)

TemporalFiltering_geneSelection=function(Xdata, sel_genes, sample_ids){
  
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  
  Xdata <-t(apply(Xdata[sel_genes,],1,normalize))
  
  #Calculating the median for all genes (selected genes) across case and controls, for different time points
  median_case<- {}
  median_control<- {}
  for (i in 1: length(sample_ids)){
    ids= sample_ids[[i]]
    case_sample_id <-ids[,"case_sample_id"]
    ctl_sample_id<- ids[,"control_sample_id"]
    
    Xcase <-Xdata[,as.character(case_sample_id)]
    Xcontrol <-Xdata[,as.character(ctl_sample_id)]
    
    median_case<- cbind(median_case, apply(Xcase, FUN=median, 1))
    median_control<- cbind(median_control,apply(Xcontrol, FUN=median, 1))
  }
  c_names<- c("Baseline","3m", "6m","9m","12m")
  colnames(median_case)<- c_names
  colnames(median_control)<-c_names
  rownames(median_case)<- rownames(Xdata)
  rownames(median_control)<- rownames(Xdata)
  
  ###Calculating the slopes between the median of a gene for each time point, for case and control groups
  max_slop_case <- vector()
  max_slop_control <- vector()
  median_slop_case<- vector()
  median_slop_control<- vector()

  mat_slop_case<- {}
  mat_slop_control <-{}
  for(i in 1:dim(Xdata)[1]){
    
    gene<- median_case[i,]
    SS_case<- diff(gene)
    mat_slop_case<- rbind(mat_slop_case, SS_case)
    max_slop_case<- c(max_slop_case, max(abs(SS_case)))
    median_slop_case<- c(median_slop_case, median(abs(SS_case)))
    
    
    gene<- median_control[i,]
    SS_control<- diff(gene)
    mat_slop_control<- rbind(mat_slop_control, SS_control)
    max_slop_control<- c(max_slop_control, max(abs(SS_control)))
    median_slop_control<- c(median_slop_control, median(abs(SS_control)))
    
  }
  
  #Selecting genes with higer slope in cases and ower slope in controls
  ind <- which(normalize(max_slop_case)> 0.15 & median_slop_control<(0.5*median_slop_case) )
  sel_genes <- rownames(median_case)[ind]
  return(sel_genes)
}
