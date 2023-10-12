#' Tutorials from the following sources were followed:
# For single dataset:
# https://rpubs.com/kenneditodd/doublet_finder_example
# For merged dataset:
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

# Load packages
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
  library(here)
  library(parallel)
  library(DoubletFinder)
})

# setwd("~/mnt/picea/home/schoudhary/shruti/SNRIII/data/SeuratOut/")

#' PART 1: run doublet finder on the filtered and merged data

load("/mnt/picea/home/schoudhary/shruti/SNRIII/data/SeuratOut/filt0.9_seurat.RData")
table(filtered_seurat$sample)

pop.split <- SplitObject(filtered_seurat, split.by = "sample") 
pop.split <- SplitObject(filt0.8_seurat, split.by = "sample") 

# Chen et al, 2021, used DoubletFinder tool with following criteria: 
# number of artificial doublets (pN) of 0.25. 
# For optimal neighborhood size (pK), the function “paramSweep_v3” was executed 
# using “PCs =1:20” and the maximum of pK value was selected as optimal pK 
# For number of expected real doublets (nExp), the doublet formation rate was assumed
# as 7.5% and nExp value was then adjusted according to homotypic doublet proportion.
# nExp, pN, and pK were set to “223, 0.25, 0.005” for wood

# loop through samples to find doublets
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
  # nExp <- round(ncol(pop.sample@meta.data) * 0.08) 
  ## Assuming 7.5% doublet formation rate - tailor for your dataset

  # run DoubletFinder
  pop.sample <- doubletFinder_v3(seu = pop.sample, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp)
  metadata <- pop.sample@meta.data
  colnames(metadata)[13] <- "doublet_finder"
  pop.sample@meta.data <- metadata 
  
  # subset and save
  pop.singlets <- subset(pop.sample, doublet_finder == "Singlet")
  pop.split[[i]] <- pop.singlets
  remove(pop.singlets)
}

pop.singlets <- merge(x = pop.split[[1]], y = pop.split[[2]],
                        project = "singletSNR")

pop.singlets
table(pop.singlets$sample)

# ctrl   kno 
# 11048 12672 

# kclmtcp knomtcp 
#  

metadata_singlet <- pop.singlets@meta.data

#' Re-assess QC if you needed
#' 
#' Visualize number of cell counts per sample and per cluster etc.: not all plots are useful
#' Plot cell =barcode = orig.ident; nCount_RNA = nUMI; nFeature_RNA = nGene
metadata_singlet %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells")
#'
metadata_singlet %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000) +
  ggtitle("UMI per Cell")
#'
metadata_singlet %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500)+
  ggtitle("Genes per Cell")
#'
#' via boxplot
metadata_singlet %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells vs Number of Genes")
#'
#' Overall complexity of the gene expression from genes detected per UMI
metadata_singlet %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.85) +
  ggtitle("Genes Per UMI")
#'

dim(metadata_singlet)

UMIGene <- pop.singlets[[c("nUMI", "nGene")]]

VlnPlot(pop.singlets, features = c("nGene", "nUMI"), ncol = 2)

# mean number of counts for each cell
counts_per_cell <- Matrix::colSums(pop.singlets, slot = 'counts')
counts_per_gene <- Matrix::rowSums(pop.singlets, slot = 'counts')

# mean number of counts for each cell you can run:
mean_counts_per_cell <- Matrix::colMeans(pop.singlets, slot = 'counts')
# 
head(AggregateExpression(pop.singlets))

seurat_phase <- NormalizeData(pop.singlets, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)
#' 
#' go to seuratSNRIII.R and continue from step 2.5
#'  
#' PART 2: If you want to run on individual sample (not merged)
pop.split <- SplitObject(filtered_seurat, split.by = "sample")

# pop.sample <- NormalizeData(pop.split$knomtcp)
pop.sample <- NormalizeData(pop.split$kclmtcp)

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

# 20 -knomtcp

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
pop.sample <- doubletFinder_v3(seu = pop.sample, PCs = 1:min.pc, 
                               pK = optimal.pk, nExp = nExp)
DF.name = colnames(pop.sample@meta.data)[grepl("DF.classification", colnames(pop.sample@meta.data))]

kclmtcp= pop.sample[, pop.sample@meta.data[, DF.name] == "Singlet"]
# 17855

knomtcp= pop.sample[, pop.sample@meta.data[, DF.name] == "Singlet"]
# 11131

# pop.singlets <- merge(x = kclmtcp, y = knomtcp, project = "singletSNRIII")
# table(pop.singlets$sample)
# kclmtcp knomtcp 
# 12122    9974 

# save singlets merged data

# seurat_phase <- NormalizeData(singlets)