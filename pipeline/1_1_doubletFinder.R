#' Following the tutorial for merged dataset
#' https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets
#'
#' Load packages
setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/")
set.seed(42)
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(RCurl)
  library(cowplot)
  library(viridis)
  library(qs)
  library(magrittr)
  library(Matrix)
  library(purrr)
  library(reshape2)
  library(S4Vectors)
  library(tibble)
  library(ComplexHeatmap)
  library(pheatmap)
  library(scales)
  library(ggplot2)
  library(parallel)
  library(DoubletFinder)
})
#'
#' Run doublet finder on the filtered and merged data
#'
load("data/SeuratOut/filt0.9_seurat.RData")
pop.split <- SplitObject(filtered_seurat, split.by = "sample") 
#'
#' Chen et al, 2021 used DoubletFinder tool with following criteria: 
#' number of artificial doublets (pN) of 0.25. 
#' For optimal neighborhood size (pK), the function “paramSweep_v3” was done 
#' PCs =1:20 and the maximum of pK value was selected as optimal pK 
#' nExp, pN, and pK were set to “223, 0.25, 0.005” for wood
#'
#'
#' Loop through samples to find doublets
for (i in 1:length(pop.split)) {
  # print the sample we are on
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  pop.sample <- NormalizeData(pop.split[[i]])
  pop.sample <- FindVariableFeatures(pop.sample)
  pop.sample <- ScaleData(pop.sample)
  pop.sample <- RunPCA(pop.sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- pop.sample[["pca"]]@stdev
  sum.stdv <- sum(pop.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  pop.sample <- RunUMAP(pop.sample, dims = 1:min.pc)
  pop.sample <- FindNeighbors(object = pop.sample, dims = 1:min.pc)              
  pop.sample <- FindClusters(object = pop.sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(pop.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  nExp <- round(optimal.pk * nrow(pop.sample@meta.data))
  
  # run DoubletFinder
  pop.sample <- doubletFinder_v3(seu = pop.sample, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp)
  metadata <- pop.sample@meta.data
  colnames(metadata)[13] <- "doublet_finder"
  pop.sample@meta.data <- metadata 
  
  # subset
  pop.singlets <- subset(pop.sample, doublet_finder == "Singlet")
  pop.split[[i]] <- pop.singlets
  remove(pop.singlets)
}
#'
pop.singlets <- merge(x = pop.split[[1]], y = pop.split[[2]],
                        project = "singletSNR")

#'
#'
metadata_singlet <- pop.singlets@meta.data
dim(metadata_singlet)
UMIGene <- pop.singlets[[c("nUMI", "nGene")]]
VlnPlot(pop.singlets, features = c("nGene", "nUMI"), ncol = 2)
#'
#'
#'
seurat_phase <- NormalizeData(pop.singlets, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)
#'
#' go to seuratSNRIII.R and continue from step 2.4