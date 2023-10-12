#' ---
#' title: "Single Cell data analysis for quick nitrate repsonse in poplar"
#' author: "Shruti"
#' #' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' 
#' Important info:
#' 
#' The analysis is for single cell RNA-Seq data for short term nitrate response 
#' (SNRIII) in xylem of hybrid aspen (T89). Two pools of protoplast suspension - 
#' kcl (control) and kno (KNO3 treated) were generated from 5 individual plants per treatment.
#' 
#' The data was generated using 10x standard kit targeting 10000 cells 
#' where cDNA was amplified for 12 cycles at step2.2 while
#' sample index PCR for library prep for 11 (kno) and 11 (kcl) cycles at step 3.5 
#' and sequenced on NovaSeq6000 with dual index kit TT set A well E5 and E6; 
#' 28nt(Read1)-10nt(Index1)-10nt(Index2) -90nt(Read2) setup using 'NovaSeqXp' workflow in 'SP' mode flowcell. 
#'
#' The Raw chromium data is in: /mnt/picea/projects/aspseq/htuominen/SingleCellSeqSNRIII/raw/
#' It was linked by me into respective sample folders: 
#' 1. /mnt/picea/home/schoudhary/shruti/SNRIII/data/kcl for P27752_1001_S1
#' 2. /mnt/picea/home/schoudhary/shruti/SNRIII/data/kno for P27752_1002_S2
#'
#' Cellranger is run by Shruti with defaults- see /mnt/picea/home/schoudhary/shruti/SNRIII/pipeline/runCellRangerCount.sh
#' and the output, for some reason, always comes in pipeline folder which was then 
#' moved to: /mnt/picea/home/schoudhary/shruti/SNRIII/data/CellrangerCount
#' 
#' Tutorials from the following sources were followed:
#' https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
#' https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html
#' https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md
#' https://github.com/Woodformation1136/SingleCell
#' http://barcwiki.wi.mit.edu/wiki/SOP/scRNA-seq/Slingshot
#' https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html#slingshot-map-pseudotime
#' https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/advanced-exercises.html
#' https://github.com/cellgeni/notebooks/blob/master/notebooks/new-10kPBMC-Scanpy.ipynb
#' https://bioconductor.org/books/3.17/OSCA.basic/quality-control.html#quality-control-motivation
#' 
#' 1. Load packages

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
  library(data.table)
})

#' 1. Load cellranger output (.h5 matrix) for control (kcl) and treated (kno3)
# setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/data/SeuratOut/")
ctrl1.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/SNRIII/data/CellRangerCount/kcl/outs/filtered_feature_bc_matrix.h5")
# 20539 x 37184

kno1.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/SNRIII/data/CellRangerCount/kno/P27752_1002/outs/filtered_feature_bc_matrix.h5")
# 14217 x 37184 

#' If you want to process the data from the previous runs as well (SNRII) and general LT
#' ctrl.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/single_cell_analysis_poplar/data/SNR1/cellranger_count/P26151_1005/outs/filtered_feature_bc_matrix.h5")
#' kno.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/single_cell_analysis_poplar/data/SNR1/cellranger_count/P26151_1006/outs/filtered_feature_bc_matrix.h5")
#' n.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/single_cell_analysis_poplar/data/LT1expUBQGFP/cellranger_count/P24151/outs/filtered_feature_bc_matrix.h5")

#' Initialize Seurat object with the non-normalized or raw data and filter:
cdata1 <- CreateSeuratObject(counts = ctrl1.data, project = "kcl1", 
                             min.cells = 3, min.features = 200)
# 20539 x 27687

tdata1 <- CreateSeuratObject(counts = kno1.data, project = "kno1", 
                             min.cells = 3, min.features = 200)
# 14217 x 26963

#' 2. Standard Pre-processing and plots
#'
#' 2.1. Merge the two data sets together into a single Seurat object. 
#' This way it will be easier to run QC steps for all samples/ groups together
#' and to compare their data quality.
#'  
merged_seurat <- merge(x = cdata1, y = tdata1, add.cell.id = c("ctrl", "kno"))
# 34756 x 28217

# Retrieve specific values from the metadata
# https://satijalab.org/seurat/articles/essential_commands.html

# Remove the elements not needed after merging: 
rm(cdata1,tdata1,ctrl1.data,kno1.data)
# garbage collect to free up memory
gc()

# 2.2. Generating quality metrics
# No. of genes detected per UMI: this metric with give us an idea of 
# data complexity (more genes detected per UMI, more complex our data)
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / 
  log10(merged_seurat$nCount_RNA)

#' Compute percent mito and  chloroplast ratio
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT")
merged_seurat[["percent.cp"]] <- PercentageFeatureSet(merged_seurat, pattern = "^CP")
#' 
#' To add more info to metadata slot in Seurat object for QC metrics, 
#' like cell IDs, condition etc., use the $ operator. But we extract dataframe 
#' into a separate variable instead to avoid affecting the original merged_seurat object.
#' Create a metadata frame by extracting metadata slot from the seurat object
metadata <- merged_seurat@meta.data

#'Add cell IDs to metadata
metadata$cells <- rownames(metadata)

#'Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
#'
#' Get sample names for each cell based on cell prefix:
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl1"
metadata$sample[which(str_detect(metadata$cells, "^kno1_"))] <- "kno1"

metadata$sample[which(str_detect(metadata$cells, "^kclmtcp_"))] <- "kclmtcp"
metadata$sample[which(str_detect(metadata$cells, "^knomtcp"))] <- "knomtcp"

#'
#' Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

#'
#'2.3. Assessing quality metrics
#' Cell counts:
#' Visualize number of cell counts per sample and per cluster etc.: not all plots are useful
#' Plot cell =barcode = orig.ident; nCount_RNA = nUMI; nFeature_RNA = nGene
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells Per Sample")
#'
#' Plot UMI (transcripts) per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000) +
  ggtitle("UMI Per Cell")
#'
#' Plot Genes per cell
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500)+
  ggtitle("Genes Per Cell")
#'
#' Plot genes per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
#'
#'Plot correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

#' Overall complexity of the gene expression can be obtained from genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.85) +
  ggtitle("Genes Per UMI")

metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

#' QC metrics plot some more examples
VlnPlot(merged_seurat, features = c("nGene", "nUMI", "mtcpRatio"), ncol = 3)+
  NoLegend()
FeatureScatter(merged_seurat, feature1 = "nUMI", feature2 = "nGene")

#'
# mean number of counts for each cell or gene with or without sample specificity
counts_per_cell <- Matrix::colSums(merged_seurat, slot = 'counts')
counts_per_gene <- Matrix::rowSums(merged_seurat, slot = 'counts', sample="kcl")

genes_per_cell <- Matrix::colSums(counts>0)
cells_per_gene <- Matrix::rowSums(counts>0)

hist(log10(counts_per_cell+1), main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('Counts vs Genes per Cell')

plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

# saveRDS(merged_seurat, file= "~/shruti/SNRIII/data/SeuratOut/output/mtCp/mrgdMtCpMnFt200.rds")

#' 2.3. Filtration
#' 2.3.1. Cell-level filtering
#' std is: nUMI > 200 & nFeature_RNA < 2500 & percent.mt < 5
#' This step is crucial: I filtered low quality reads at following thresholds
#' 
filtered_seurat <- subset(x = merged_seurat, subset= (nUMI >= 500) &
                            (nGene >= 200) & (log10GenesPerUMI > 0.9) &
                            (percent.mt < 5) & (percent.cp < 5))
table(filt0.8_seurat$sample)
# 28491 x 28217
# ctrl   kno 
# 15560 12931 

# For mt cp genome
# filt0.8_seurat <- subset(merged_seurat, subset= (nUMI >= 500) & (nGene >= 200)&
#                            (percent.mt < 5) & (percent.cp < 5))

filt0.8_seurat <- subset(merged_seurat, subset= (nUMI >= 1000) & (nGene >= 1000)&
                           (nGene < 7000) & (nUMI < 70000) &
                           (percent.mt < 5) & (percent.cp < 5))

# or you can filter at:
# filt0.8_seurat <- subset(x = merged_seurat, subset= (nUMI >= 500) & (nGene >= 200))
# ctrl   kno 
# 20539 14217

# If I set to filtration parameter in Chen et al., 2021
# subset= (nUMI >= 500) & (nUMI < 70000) & (nGene >= 200) & (nGene < 9000) & 
# (log10GenesPerUMI > 0.80))
# 34691 x 28217 
# ctrl   kno 
# 20537 14154 

# For mt cp genome,
# filt0.8_seurat <- subset(merged_seurat, subset= (nUMI >= 500) & (nGene >= 200))
# 37254x 34186
# kclmtcp knomtcp 
# 20581   13605

# 37254x 27855
# kclmtcp knomtcp 
# 15541   12314 

# saveRDS(filt0.9_seurat, file = "data/SeuratOut/output/mtCp/mtCpMnFilt0.9.rds")


#' 2.3.2. Gene-level filtering , a redundant step- decided not to do this 
#' Removed genes with zero expression in all cells 
filt_counts <- GetAssayData(object = filt0.8_seurat, slot = "counts")
nonzero <- filt_counts > 0
# keep only genes which are expressed in 3 or more cells.
keep_genes <- Matrix::rowSums(nonzero) >= 10
filt_counts <- filt_counts[keep_genes, ]
summary (keep_genes)
# Mode   FALSE    TRUE 
# logical     758   27459

#' 2.3.3. Re-assess QC metrics of the filtered Seurat object
#' 
filtered_seurat <- CreateSeuratObject(filt_counts,
                                      meta.data = filtered_seurat@meta.data)
metadata_clean <- filtered_seurat@meta.data
#' 
#' Visualize number of cell counts per sample and per cluster etc
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells after Filtration")
#'
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000) +
  ggtitle("UMI per Cell after Filtration")
#'
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500)+
  ggtitle("Genes per Cell after Filtration")
#'
#' via boxplot
metadata_clean %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("No. of Cells vs No. of Genes after Filtration")
#'
#' Overall complexity of the gene expression from genes detected per UMI
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.85) +
  ggtitle("Genes Per UMI after Filtration")
#'

dim(metadata_clean)

UMIGene <- filtered_seurat[[c("nUMI", "nGene")]]
VlnPlot(filtered_seurat, features = c("nGene", "nUMI"), ncol = 2)

# mean number of counts for each cell
counts_per_cell <- Matrix::colSums(filtered_seurat, slot = 'counts')
counts_per_gene <- Matrix::rowSums(filtered_seurat, slot = 'counts')

# mean number of counts for each cell you can run:
mean_counts_per_cell <- Matrix::colMeans(filtered_seurat, slot = 'counts')
# 
# save(filtered_seurat, file="/mnt/picea/home/schoudhary/shruti/SNRIII/data/SeuratOut/filt0.9_seurat.RData")
# save(filtered_seurat, file="~/shruti/SNR-u2023011/analysis/snrIII/filt0.9_mtcp.RData")
# Remove useless data

# Go to /mnt/picea/home/schoudhary/shruti/SNRIII/pipeline/doubletFinder.R to remove doublets
# and follow this script for merged and filtered samples (PART 1).
# There is no concrete evidence suggesting it is necessary to remove or include the doublets, 
# Based on discussion with Bert and Kedar, I will prefer to remove them
# 
#' 2.4. Normalize data
#' 
#'Count Normalization and Principal Component Analysis
#'
#' Next is clustering that separate different cell types into unique clusters.
#' For clustering, determine genes that have the most different expression between cells. 
#' These genes are then used to determine which correlated genes sets determine
# the largest expression differences between cells.

# Count normalization: essential to make accurate comparisons of gene expression
# between cells (or samples). The counts of mapped reads for each gene is 
# proportional to the expression of RNA ("interesting") in addition to other factors. 

# Normalization is the scaling of raw count values to account for "uninteresting" factors. 
# In this way the expression levels are more comparable between and/or within cells.
# main factors considered during normalization are:

# Sequencing depth: necessary for comparison of gene expression between cells. 
# Gene length: necessary for comparing expression between different genes within the same cell. 
# The number of reads mapped to a longer gene can appear to have equal 
# count/expression as a shorter gene that is more highly expressed.
#' 
#' 2.5. Cell Cycle scoring
#' 
#' Pick normalized pop_singlet data wherein doublets were removed (called seurat_phase)
#' 23720 x 27459 and continue from this point onwards
#' 
#' There are mixed opinions on removing the cell cycle associated variation and
#' since there is no solid explanation for removal or including them in the analysis,
#' I will anyways regress the variation due to cell cycle
#' 
#' Genes for cell cycle are picked from two arabidospis papers Yadav et al., 2009
#' and Gutierrez, 2009 as describe din Chen et al., 2021, Genome Biology
#' 
#' Find cell cycle genes (sphase and g2m pahse) list for poplar from CellTypeMarker.R
#' 
seurat_phase <- CellCycleScoring(seurat_phase, s.features = sphase, 
                                 g2m.features = c(g2phase,mphase), 
                                 set.ident = TRUE)
#' 
#' After scoring for cell cycle, determine if cell cycle could be a source of variation 
#' by looking at PCA. To do PCA, first choose the most variable features 
#' (default is 2000), then scale the data.

seurat_phase <- FindVariableFeatures(seurat_phase,
                                      selection.method = "vst", 
                                      nfeatures = 2000, mean.cutoff = c(0.01, 8), 
                                      dispersion.cutoff = c(0.5, Inf))

# ScaleData() function: Shifts expression of each gene, so that mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in object[["RNA"]]@scale.data

seurat_phase <- ScaleData(seurat_phase, features = rownames(seurat_phase))

#' Run and view PCA using expression of cell cycle genes. 
seurat_phase <- RunPCA(seurat_phase)
DimPlot(seurat_phase, reduction = "pca", 
        split.by = "Phase")

#' 3. SCTransform
#' 
#' The presented dataset of cells seem to group by cell cycle in the PCA, 
#' so I decided to regress out cell cycle variation
#' 
#' NOTE: (i) If cells are known to be differentiating and there is clear clustering 
#' differences between G2M and S phases, then maybe regress out by the difference
#' between the G2M and S phase scores, and still differentiating the cycling 
#' from the non-cycling cells.
#' 
#' (ii) ScaleData() function can remove unwanted sources of variation example, 
# ‘regress out’ heterogeneity associated with cell cycle stage
# However, for advanced users, it is recommended to use new normalization- SCTransform
# 
# Load mtcp 0.9 log singlets
Snglt0.9 <- readRDS("~/shruti/SNRIII/data/SeuratOut/output/mtCp/Snglt0.9.rds")
Snglt0.9 <- subset(Snglt0.9,  subset= (percent.mt < 5) & (percent.cp < 5))

# Adjust the memory
options(future.globals.maxSize = 4000 * 1024^2)

#' 
split_seurat <- SplitObject(pop.singlets, split.by = "sample")
split_seurat <- split_seurat[c("kclmtcp", "knomtcp")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], 
                                        g2m.features=c(g2phase,mphase), 
                                        s.features=sphase)
  split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], 
                                            selection.method = "vst", 
                                            nfeatures = 2000)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                   # vars.to.regress = c("S.Score", "G2M.Score"),
                                   variable.features.n = 3000
                                   )
}

# NOTE: 
# (i). By default, after normalizing, adjusting the variance, and regressing, 
# SCTransform will rank genes by residual variance and output 3000 most variant 
# genes. If the dataset has larger cell numbers, then adjust it to higher value
# using variable.features.n argument.

# (ii). can use seurat_phase object here instead of pop-singlets data

# (iii). can also use regress out protoplasting genes from the dataset 
# using Yadav_pp gene sets in CellTypeMarker.R file but I have not done this

#' 4. Integrate samples using shared highly variable genes

integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = 3000) 

split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")


# Run the standard workflow for visualization and clustering
# specify that downstream analysis will be done on the corrected data 
# original unmodified data still resides in the 'RNA' assay

DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50, reduction = "pca")
DimPlot(seurat_integrated)
#'
#' 4. Clustering of cells
#'
#' First explore the number of PC to use
ElbowPlot(object = seurat_integrated, ndims = 40)

# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/elbow_plot_metric.md
# Turns out to be 18 using the tutorial above

# Don't use too much PC
#' Find K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated,
                                   dims = 1:50)

# seurat_integrated <- FindNeighbors(object = seurat_integrated,
#                                    dims = 1:18)

#' Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

seurat_integrated@meta.data %>% View()

# Remove some of the objects: save only integ data
# save.image("~/shruti/SNRIII/data/SeuratOut/integ.RData")
save.image("~/shruti/SNRIII/data/SeuratOut/integafterdbltremoval.RData")

# seurat_integrated <- BuildClusterTree(seurat_integrated, dims = 1:50, 
#                                       reorder = F, reorder.numeric = F)
# 
# PlotClusterTree(seurat_integrated)

Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"

p3 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 2)
p4<- DimPlot(seurat_integrated, label = TRUE, split.by = "sample")  + NoLegend()
p4

# If cells group by cell cycle phase
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score")
p5 <- FeaturePlot(seurat_integrated, reduction = "umap", features = metrics,
            pt.size = 0.4, sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
p5

DimPlot(seurat_integrated, reduction = "umap", split.by = "Phase", label = TRUE,
        label.size = 6)

# Extract identity and sample info to determine no. of cells per cluster per sample
n_cells <- FetchData(integClCyc, 
                     vars = c("ident")) %>%
  dplyr::count(ident) %>%
  tidyr::spread(ident, n)

View(n_cells)

# Extract counts per gene
counts_per_gene <- Matrix::rowSums(filtered_seurat, slot = 'counts')

# To truly determine the identity of the clusters and whether the resolution is 
# appropriate, it is helpful to explore a handful of known gene markers for cell types expected.

#' 5. Finding DE markers for each cluster
#' 
#'# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize the RNA data 
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
# Note: raw and normalized counts are stored in counts and data slots of RNA assay. 

# By default, the functions for finding markers will use normalized data.
#'
markers <- FindAllMarkers(object = seurat_integrated,
                          only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25)

# Save markers as markersWilcox.RData

# 5.1. Goto enrichment.R for enrichment of markers clusterwise

write.table(markers, file = "output/afterDbltRemoval/allmarkercluster.txt", sep = "\t",
            row.names = T, col.names = T)

markers.sortedByPval = markers[order(markers$p_val),]

genes_uniquely_DE = markers.sortedByPval %>% 
  dplyr::filter(avg_log2FC >= 1) %>% group_by(gene) %>%  
  summarize(n=n()) %>%  dplyr::filter(n==1)

genes_uniquely_DE.markers.sortedByPval =
  markers.sortedByPval[markers.sortedByPval$gene %in% genes_uniquely_DE$gene & 
                         markers.sortedByPval$avg_log2FC >= 1,]

top_marker_each = genes_uniquely_DE.markers.sortedByPval %>% 
  dplyr::group_by(cluster) %>% do(head(., n=10))

write.table(top_marker_each, file = "output/afterDbltRemoval/top10markercluster.txt", sep = "\t",
            row.names = T, col.names = T)

png("plots/afterDbltRemoval/top10markercluster.png", width = 2480, height = 980,
    pointsize = 12, bg = "transparent")
for (i in length(top_marker_each$gene)){
  print(DotPlot(seurat_integrated, features = top_marker_each$gene)+RotatedAxis())
}
dev.off()
#' 
#' 5.2. or you can get conserved markers for each condition like as follows:
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE)
  }

conserved_markers <- map_dfr(c(0:20), get_conserved)

write.table(conserved_markers, file = "output/afterDbltRemoval/cons_marker.txt", sep = "\t",
            row.names = T, col.names = T)

# 6. decided not to regress cell cycle and find the markers
seurat_integrated <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(seurat_integrated) <- "RNA"
# seurat_integrated <- subset(seurat_integrated, integrated_snn_res.0.6 == '2', invert=T)
# seurat_integrated$oldclusters <- seurat_integrated$integrated_snn_res.0.6
# seurat_integrated$seurat_clusters <- seurat_integrated@active.ident
# seurat_integrated$orig.ident <- paste0("integ_",seurat_integrated$seurat_clusters)
# seurat_integrated$orig.samp <- "integ"
  
split_seurat <- SplitObject(seurat_integrated, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]

for (i in 1:length(split_seurat)) {
  # split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  # split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], 
  #                                       g2m.features=c(g2phase,mphase), 
  #                                       s.features=sphase)
  # split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], 
  #                                           selection.method = "vst", 
  #                                           nfeatures = 2000)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], variable.features.n = 3000
  )
}

integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = 3000) 

split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

integClCyc <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
# Remove unwanted objects
rm(split_seurat, i, integ_anchors, integ_features)

DefaultAssay(integClCyc) <- "integrated"
integClCyc <- ScaleData(integClCyc, verbose = FALSE)
integClCyc <- RunPCA(integClCyc)
integClCyc <- RunUMAP(integClCyc, dims = 1:50, reduction = "pca") #40
ElbowPlot(object = integClCyc, ndims = 40)
integClCyc <- FindNeighbors(integClCyc, dims = 1:50)
integClCyc <- FindClusters(integClCyc, resolution = c(0.1, 0.3, 0.4, 0.5, 0.6, 0.8))
DimPlot(integClCyc, label=T)

Idents(object = integClCyc) <- "integrated_snn_res.0.6"
integClCyc <- RunUMAP(integClCyc, dims = 1:30, reduction = "pca", seed.use = 5) #40
integClCyc <- FindNeighbors(integClCyc, dims = 1:30)
integClCyc <- FindClusters(integClCyc, resolution = c(0.1, 0.6), random.seed = 689213)
DimPlot(integClCyc, label=T)

# save.image("~/shruti/SNRIII/data/SeuratOut/integWthClCyc.RData")
Idents(object = integClCyc) <- "integrated_snn_res.0.6"

# Load in the data with mtcp genes
ctrl1.data <- Read10X_h5("~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/kcl2/outs/filtered_feature_bc_matrix.h5")
# 37254 x 20581 

kno1.data <- Read10X_h5("~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/kno2/outs/filtered_feature_bc_matrix.h5")
# 37254 x 13605

cdata1 <- CreateSeuratObject(counts = ctrl1.data, project = "kclmtcp", 
                             min.cells = 0, min.features = 200)
# 37254 x 20581 

# kclEmtyDrp <- readRDS("~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/kclEmtyDrp.rds")
# 37254 x 15234

tdata1 <- CreateSeuratObject(counts = kno1.data, project = "knomtcp", 
                             min.cells = 0, min.features = 200)

# 37254 x13605

merged_seurat <- merge(cdata1, y = tdata1, 
                       add.cell.ids = c("kclmtcp", "knomtcp"), project = "SNRIII")
# 37254 x 34186 

# merged_seurat <- merge(kclEmtyDrp, y = tdata1,
#                        add.cell.ids = c("kclmtcp", "knomtcp"), project = "SNRIII")
# 37254 x 28839

#' Now goto step 2.2 and then to 6 or follow ahead check how clustering varies with 
#' and without cell cycle

#' 7. Rename clusters 
# seurat_labelled <- RenameIdents(object = split_seurat$ctrl,
#                                 "4" = "Fibers1", "5" = "Rays1", 
#                                 "7" = "Rays2", "12" = "Fibers2", 
#                                 "15" = "Fibers3", "17" = "Vessels",
#                                 "19" = "Phloem", "20" = "Cambium")
# # 
# seurat_labelled$celltype.sample <- paste(Idents(seurat_labelled), seurat_labelled$sample,
#                                       sep = "_")

#' 8. See DE_analysis_scrnaseq.R for pseudobulk DE analysis and bulkDE analysis 
#' (from Part1 and 2, respectively)
#' 
#' 9. See the end of combined_cell_cycle.R for pseudo time trajectory: it works for only single sample
#' Working with pseudotime in trajectory.R
#' 
#' 10. See enrichment.R for enrichment
# 