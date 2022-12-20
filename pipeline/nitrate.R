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

#' This analysis is for single cell RNA-Seq with 10x standard 
#' where cDNA was amplified for 12 cycles at step2.2 while
#' sample index PCR for library prep for 11 (kno3) and 11 (kcl) cycles at step 3.5 
#' and sequenced on NovaSeq6000 with dual index kit TT set A well E5 and E6; 
#' 28nt(Read1)-10nt(Index1)-10nt(Index2) -90nt(Read2) setup using 'NovaSeqXp' workflow in 'SP' mode flowcell. 
#'
#' Cellranger is run by me with defaults and output is in ../data/CellrangerCount
#' See runCellRangerCount.sh in pipeline folder
#' Load packages
#' 
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
})

#' 1. Load cellranger output (.h5 matrix) for control (kcl) and treated (kno3)
#' 
setwd("~/shruti/SNRIII/data/SeuratOut/")
ctrl1.data <- Read10X_h5("../CellRangerCount/kcl/outs/filtered_feature_bc_matrix.h5")
kno1.data <- Read10X_h5("../CellRangerCount/kno/P27752_1002/outs/filtered_feature_bc_matrix.h5")
# ctrl.data <- Read10X_h5("../../../single_cell_analysis_poplar/data/nitrate_resp/cellranger_count/P26151_1005/outs/filtered_feature_bc_matrix.h5")
# kno.data <- Read10X_h5("../../../single_cell_analysis_poplar/data/nitrate_resp/cellranger_count/P26151_1006/outs/filtered_feature_bc_matrix.h5")
#'
#' Initialize Seurat object with the non-normalized data and filter:
#' 
cdata1 <- CreateSeuratObject(counts = ctrl1.data, project = "kcl1", min.cells = 3,
                            min.features = 200)
tdata1 <- CreateSeuratObject(counts = kno1.data, project = "kno1", min.cells = 3,
                            min.features = 200)
# cdata <- CreateSeuratObject(counts = ctrl.data, project = "kcl", min.cells = 3,
#                              min.features = 200)
# tdata <- CreateSeuratObject(counts = kno.data, project = "kno", min.cells = 3,
#                              min.features = 200)
# cdata$type = "ctrl"
cdata1$type = "ctrl"
# tdata$type = "kno"
tdata1$type = "kno"
#'
#' 2. Standard Pre-processing and plots
#' merge the two samples
# alldata <- merge(cdata, c(cdata1, tdata, tdata1), 
#                  add.cell.ids = c("kcl", "kcl1", "kno", "kno1"))

merged_seurat <- merge(x = cdata1, y = tdata1, 
                       add.cell.id = c("ctrl", "kno"))
rm(cdata, cdata1, tdata, tdata1, ctrl1.data, ctrl.data, kno1.data, kno.data)
gc()

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / 
  log10(merged_seurat$nCount_RNA)
#' 
#' Create a metadata frame by extracting metadata slot from the seurat object

metadata <- merged_seurat@meta.data
#'
#'Add cell IDs to metadata
metadata$cells <- rownames(metadata)
#'
#'Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
#'
#' Get sample names for each cell based on cell prefix:
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^kno"))] <- "kno"
#' Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
#'
#' Visualize number of cell counts per sample and per cluster etc.: not all plots are useful
#' Plot cell =barcode = orig.ident; nCount_RNA = nUMI; nFeature_RNA = nGene
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
#'
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  ggtitle("UMIperCell")
#'
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)+
  ggtitle("GeneperCell")
#'
#' via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
#'
#' Overall complexity of the gene expression from genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  ggtitle("Genes Per UMI")
#'
#' QC metrics plot example
VlnPlot(merged_seurat, features = c("nGene", "nUMI"), ncol = 2)
FeatureScatter(merged_seurat, feature1 = "nUMI", feature2 = "nGene")
as.data.frame(merged_seurat@assays$RNA@counts[1:10, 1:2])

#'
#' Percentage of Largest Gene
#'  
apply(merged_seurat@assays$RNA@counts,2, function(x)(100*max(x))/sum(x)) -> 
  merged_seurat$Percent.Largest.Gene
head(merged_seurat$Percent.Largest.Gene)

FeatureScatter(merged_seurat,feature1 = "nUMI", feature2 = "Percent.Largest.Gene")

as_tibble(merged_seurat[[c("nUMI","nGene","Percent.Largest.Gene")]],
          rownames="Cell.Barcode") -> qc.metrics
qc.metrics
qc.metrics %>%
  ggplot(aes(Percent.Largest.Gene)) + 
  geom_histogram(binwidth = 0.7, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Largest Gene") +
  geom_vline(xintercept = 20)
# 
# dim(merged_seurat)
# 28217 34756
#'
#' 2.1. Filtration
#' 2.1.1. Cell-level filtering
#' nUMI > 200 & nFeature_RNA < 2500 & percent.mt < 5
#' This step is crucial: I filtered out low quality reads using following thresholds
#' 
filtered_seurat <- subset(x = merged_seurat, subset= (nUMI >= 500) &
                            (nGene >= 200) & (log10GenesPerUMI > 0.90))

# selected_c <- WhichCells(alldata, expression = nFeature_RNA > 500)
# selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 10]
# data.filt <- subset(alldata, features = selected_f, cells = selected_c)

# Following filtration parameter in Chen et al., 2021
# subset= (nUMI >= 500) & (nUMI < 70000) & (nGene >= 200) & (nGene < 9000) & 
# (log10GenesPerUMI > 0.80))

#' 2.1.2. Gene-level filtering
#' Removed genes with zero expression in all cells and 
#' keep only those expressed in >=3 cells.
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
#'
#' Sums all TRUE values and returns TRUE if more than 3 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 3
filtered_counts <- counts[keep_genes, ]
#'
#' 2.2 Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, 
                                     meta.data = filtered_seurat@meta.data)
metadata_clean <- filtered_seurat@meta.data

#'
#' 2.3. Re-assess QC if you needed
#'
#' 2.4. Normalize data
seurat_phase <- NormalizeData(filtered_seurat, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)
#' 
#' 2.5. cell cycle scoring
#' verify if the genes are sorted into correct cell cycle phasing, goto marker.R
#'
seurat_phase <- CellCycleScoring(seurat_phase, s.features = sphase, 
                                 g2m.features = c(g2phase,mphase),
                                 set.ident = TRUE)

#'
#' After scoring for cell cycle, determine if cell cycle could be a source of variation 
#' by looking at PCA. 
#' To perform PCA, first choose the most variable features, then scale the data. 

seurat_phase <- FindVariableFeatures(seurat_phase,
                                      selection.method = "vst", 
                                      nfeatures = 2000, mean.cutoff = c(0.01, 8), 
                                      dispersion.cutoff = c(0.5, Inf))
VariableFeaturePlot(seurat_phase)
top10 <- head(VariableFeatures(seurat_phase), 10)

seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

DimPlot(seurat_phase, reduction = "pca", group.by= "Phase", split.by = "Phase")
#'
#' 3. SCTransform
options(future.globals.maxSize = 4000 * 1024^2)

#' Split seurat object by condition to (optional: perform cell cycle scoring)
#' and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]

# to regress out the cell cycle variation
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=c(g2phase,mphase), s.features=sphase)
  split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], selection.method = "vst", nfeatures = 2000)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score", "G2M.Score"))
}

#' 4. Integrate samples using shared highly variable genes

integ_features <- SelectIntegrationFeatures(object.list = split_seurat) 

split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

DefaultAssay(seurat_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(object = seurat_integrated, npcs = 50)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50, reduction = "pca")
DimPlot(seurat_integrated)
#'

#' 4. Cluster cells
#'
DimHeatmap(seurat_integrated, dims = 1:9, cells = 500, balanced = TRUE)
print(x = seurat_integrated[["pca"]], dims = 1:10, nfeatures = 5)
ElbowPlot(object = seurat_integrated, ndims = 50)

#' Find K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated,
                                   dims = 1:50)

#' Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

seurat_integrated@meta.data %>% View()

seurat_integrated <- BuildClusterTree(seurat_integrated, dims = 1:50, 
                                      reorder = F, reorder.numeric = F)

PlotClusterTree(seurat_integrated)

#' Assign identity of clusters and plot and explore resolutions
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 2)
#'
#' Explore QC if needed
#' Some useful data
DimPlot(seurat_integrated, label = TRUE, split.by = "sample")  + NoLegend()
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score")
FeaturePlot(seurat_integrated, reduction = "umap", features = metrics,
            pt.size = 0.4, sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

# Remove some of the objects: save only integ data
save.image("~/shruti/SNRIII/data/SeuratOut/integ.RData")

# Known markers (Tung et al., 2022)
FibOrgan <- c("Potra2n1c1345","Potra2n4c8430","Potra2n6c13614","Potra2n15c28766")
EarlyFib <- c("Potra2n568s35990","Potra2n7c16454")
InterFib <- c("Potra2n5c12126","Potra2n6c14930","Potra2n14c26481")
Fib <- c("Potra2n2c4078","Potra2n16c30073","Potra2n10c21170","Potra2n15c28689")
LateFib_Ves <- c("Potra2n5c10753","Potra2n8c16805","Potra2n17c31807","Potra2n19c33790")
RayPrecu <- c("Potra2n8c18231","Potra2n18c32244")
RayOrgan <- c("Potra2n5c11373","Potra2n8c16999")
RayParen <- c("Potra2n2c5467","Potra2n2c4942","Potra2n2c4144","Potra2n155s34744")
# 
FeaturePlot(seurat_integrated, reduction = "umap", features = RayParen,
            sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

#' 
#' 5. Finding DE markers for each cluster
#'
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- NormalizeData(seurat_integrated)
#'
markers <- FindAllMarkers(object = seurat_integrated)

markers.sortedByPval = markers[order(markers$p_val),]

genes_uniquely_DE = markers.sortedByPval %>% 
  dplyr::filter(avg_log2FC >= 1) %>% group_by(gene) %>%  
  summarize(n=n()) %>%  dplyr::filter(n==1)

genes_uniquely_DE.markers.sortedByPval =
  markers.sortedByPval[markers.sortedByPval$gene %in% genes_uniquely_DE$gene & 
                         markers.sortedByPval$avg_log2FC >= 1,]

top_marker_each = genes_uniquely_DE.markers.sortedByPval %>% 
  dplyr::group_by(cluster) %>% do(head(., n=10))

write.table(top_marker_each, file = "top10markercluster.txt", sep = "\t",
            row.names = F, col.names = T)

png("plots/top10markercluster.png", width = 2480, height = 980,
    pointsize = 12, bg = "transparent")
for (i in length(top_marker_each$gene)){
  print(DotPlot(seurat_integrated, features = top_marker_each$gene)+RotatedAxis())
}
dev.off()
#' 
#' or you can get conserved markers for each condition like as follows:
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE)
  }

conserved_markers <- map_dfr(c(0:21), get_conserved)

write.table(conserved_markers, file = "cons_marker.txt", sep = "\t",
            row.names = T, col.names = T)
#'
#' Rename clusters 
seurat_labelled <- RenameIdents(object = seurat_integrated,
                                "0" = "LateFib", 
                                "4" = "RayPrec", "5" = "Fib", "6" = "LateFib", 
                                "10" = "Fib", 
                                "12" = "Fib", "13" = "LateFib", "14" = "EarlyFib",
                                "15" = "FibOrgan", "16" = "EarlyFib")
# 
# seurat_labelled$celltype.sample <- paste(Idents(seurat_labelled), seurat_labelled$sample,
#                                       sep = "_")

#' 6. See DE_analysis_scrnaseq.R for pseudobulk DE analysis and bulkDE analysis (Part1 and 2, respectively)
#' 
#' 7. See the end of ./big/combined_cell_cycle.R for pseudo time trajectory: it works for only single sample
#' Working with pseudotime in trajectory.R
#' 
#' 8. See enrichment.R for enrichment