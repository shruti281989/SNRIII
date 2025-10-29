#' ---
#' title: "Single Cell data analysis for quick nitrate response in poplar SNRIII"
#' author: "Shruti"
#' #' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' 
#' Info:
#' 
#' Check if the paths to files are correct
#' 
#' This is Seurat analysis for single cell RNA-Seq data after short term (2h) nitrate response 
#' (SNRIII) in xylem of hybrid aspen (T89). Two pools of protoplast suspension - 
#' kcl (control) and kno (KNO3 treated) were generated from 5 plants per treatment.
#' 
#' The data was generated using 10x standard kit targeting 10000 cells 
#' where cDNA was amplified for 12 cycles at step2.2 while
#' sample index PCR for library prep for 11 (kno) and 11 (kcl) cycles at step 3.5 
#' and sequenced on NovaSeq6000 with dual index kit TT set A well E5 and E6; 
#' 28nt(Read1)-10nt(Index1)-10nt(Index2) -90nt(Read2) setup using 
#' 'NovaSeqXp' workflow in 'SP' mode flowcell. 
#'
#' Raw data is in: /mnt/ada/projects/aspseq/htuominen/SingleCellSeqSNRIII/raw
#' A link was created in the follwoing sample folders: 
#' /mnt/picea/home/schoudhary/shruti/SNRIII/data/raw/kcl for P27752_1001_S1
#' /mnt/picea/home/schoudhary/shruti/SNRIII/data/raw/kno for P27752_1002_S2
#'
#' For running cellranger see, /mnt/picea/home/schoudhary/shruti/SNRIII/pipeline/runCellRangerCount.sh
#' and move the output to the desired location
#' 
#' 
#' Load packages
#' 
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
  library(viridis)
  library(scCustomize)
  library(qs)
})
#'
#'
#' 1. Load cellranger output (.h5 matrix) for control (kcl) and treated (kno3)
ctrl1.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNucl/kcl2/outs/filtered_feature_bc_matrix.h5")
#' 20539 x 37184
#'
kno1.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNucl/kno2/outs/filtered_feature_bc_matrix.h5")
#'14217 x 37184 
#'
#' If you want to process the data from the previous runs as well (SNRII) and general LT
#' ctrl.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/single_cell_analysis_poplar/data/SNR1/cellranger_count/P26151_1005/outs/filtered_feature_bc_matrix.h5")
#' kno.data <- Read10X_h5("/mnt/picea/home/schoudhary/shruti/single_cell_analysis_poplar/data/SNR1/cellranger_count/P26151_1006/outs/filtered_feature_bc_matrix.h5")
#'
#'
#'Initialize Seurat object with the non-normalized or raw data and filter:
cdata1 <- CreateSeuratObject(counts = ctrl1.data, project = "kcl1", 
                             min.cells = 3, min.features = 200)
#' 20539 x 27687
#' 
tdata1 <- CreateSeuratObject(counts = kno1.data, project = "kno1", 
                             min.cells = 3, min.features = 200)
#' 14217 x 26963
#' 
#' 
#' 2. Standard Pre-processing and plots
#'
#' 2.1. Merge the two data sets together into a single Seurat object. 
#' This way it will be easier to run QC steps for all samples/ groups together
#' and to compare their data quality.
#'  
merged_seurat <- merge(x = cdata1, y = tdata1, add.cell.id = c("ctrl", "kno"))
#' 34756 x 28217
#' 
#' Retrieve specific values from the metadata
#' #' https://satijalab.org/seurat/articles/essential_commands.html
#' Remove the objects not needed after merging: 
rm(cdata1,tdata1,ctrl1.data,kno1.data)
#' garbage collect to free up memory
gc()
#'
#' 2.2. Quality metrics
#' Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / 
  log10(merged_seurat$nCount_RNA)
#'
#' Compute percent mito and  chloroplast ratio
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT")
merged_seurat[["percent.cp"]] <- PercentageFeatureSet(merged_seurat, pattern = "^CP")
#' 
#' To add more info to metadata slot in Seurat object for QC metrics, 
#' like cell IDs, condition etc., use the $ operator. But we first extract dataframe 
#' into a separate variable instead to avoid affecting the original merged_seurat object.
#' Create a metadata frame by extracting metadata slot from the seurat object
metadata <- merged_seurat@meta.data
#'
#'Add cell IDs to metadata
metadata$cells <- rownames(metadata)
#'
#'Rename columns
metadata <- metadata %>% dplyr::rename(seq_folder = orig.ident,
                                       nUMI = nCount_RNA, nGene = nFeature_RNA)
#'
#' Get sample names for each cell based on cell prefix:
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^kno1_"))] <- "kno"
#'
#' Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
#'
#' Assessing quality metrics
#' 2.2.1. Cell counts:
#' Visualize number of cell counts per sample and per cluster etc.: not all plots are useful
#' Plot cell =barcode = orig.ident; nCount_RNA = nUMI; nFeature_RNA = nGene
metadata %>% ggplot(aes(x=sample, fill=sample)) + geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells Per Sample")
#'
#' 2.2.2. Plot UMI (transcripts) per cell
metadata %>% ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  ylab("Cell density") + geom_vline(xintercept = 1000) + ggtitle("UMI Per Cell")
#'
#' 2.2.3. Plot Genes per cell
metadata %>% ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + theme_classic() + scale_x_log10() + 
  geom_vline(xintercept = 500)+ ggtitle("Genes Per Cell")
#'
#' 2.2.4. Plot genes per cell via boxplot
metadata %>% ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
#'
#' 2.2.5. Plot correlation between genes detected and number of UMIs and determine
#' whether strong presence of cells with low numbers of genes/UMIs
metadata %>% ggplot(aes(x=nUMI, y=nGene)) + geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) + scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(xintercept = 500) + geom_hline(yintercept = 250) +
  facet_wrap(~sample)
#'
#' 2.2.6. Overall complexity of the gene expression can be obtained from genes detected per UMI
metadata %>% ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) + theme_classic() + geom_vline(xintercept = 0.85) +
  ggtitle("Genes Per UMI")
#'
metadata %>% ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  geom_vline(xintercept = 0.2)
#'
#' 2.2.7. QC metrics plot more examples
VlnPlot(merged_seurat, features = c("nGene", "nUMI", "mtcpRatio"), ncol = 3)+
  NoLegend()
FeatureScatter(merged_seurat, feature1 = "nUMI", feature2 = "nGene")
#'
# 2.2.8. mean number of counts for each cell or gene with or without sample specificity
counts_per_cell <- Matrix::colSums(merged_seurat, slot = 'counts')
counts_per_gene <- Matrix::rowSums(merged_seurat, slot = 'counts', sample="kcl")
genes_per_cell <- Matrix::colSums(counts>0)
cells_per_gene <- Matrix::rowSums(counts>0)
hist(log10(counts_per_cell+1), main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('Counts vs Genes per Cell')
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
#'
#'
#' 2.3. Filtration
#' 2.3.1. Cell-level filtering
#' 
filtered_seurat <- subset(x = merged_seurat, subset= (nUMI >= 500) &
                            (nGene >= 200) & (log10GenesPerUMI > 0.9) &
                            (percent.mt < 5) & (percent.cp < 5))
#'
#'
#' 2.3.2. Gene-level filtering 
#' Removed genes with zero expression in all cells 
filt_counts <- GetAssayData(object = filt0.9_seurat, slot = "counts")
nonzero <- filt_counts > 0
#' keep only genes which are expressed in 3 or more cells.
keep_genes <- Matrix::rowSums(nonzero) >= 3
filt_counts <- filt_counts[keep_genes, ]
summary (keep_genes)
#' Mode   FALSE    TRUE 
#' logical     758   27459
#'
#' 2.3.3. Re-assess QC metrics of the filtered Seurat object
#' 
filtered_seurat <- CreateSeuratObject(filt_counts,
                                      meta.data = filtered_seurat@meta.data)
metadata_clean <- filtered_seurat@meta.data
#' 
#' Visualize number of cell counts per sample and per cluster etc
metadata_clean %>%  ggplot(aes(x=sample, fill=sample)) + geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells after Filtration")
#'
metadata_clean %>% ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  ylab("Cell density") + geom_vline(xintercept = 1000) +
  ggtitle("UMI per Cell after Filtration")
#'
metadata_clean %>% ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + theme_classic() + scale_x_log10() + 
  geom_vline(xintercept = 500)+ ggtitle("Genes per Cell after Filtration")
#'
#' via boxplot
metadata_clean %>% ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("No. of Cells vs No. of Genes after Filtration")
#'
#' Overall complexity of the gene expression from genes detected per UMI
metadata_clean %>% ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) + theme_classic() + geom_vline(xintercept = 0.85) +
  ggtitle("Genes Per UMI after Filtration")
#'
dim(metadata_clean)
#'
UMIGene <- filtered_seurat[[c("nUMI", "nGene")]]
VlnPlot(filtered_seurat, features = c("nGene", "nUMI"), ncol = 2)
#'
#' mean number of counts for each cell
counts_per_cell <- Matrix::colSums(filtered_seurat, slot = 'counts')
counts_per_gene <- Matrix::rowSums(filtered_seurat, slot = 'counts')
#'
#' mean number of counts for each cell you can run:
mean_counts_per_cell <- Matrix::colMeans(filtered_seurat, slot = 'counts')
#'
#'
save(filtered_seurat, file="data/SeuratOut/filt0.9_seurat.RData")
#' Remove objects not needed
#' 
#' 
#' Go to /mnt/picea/home/schoudhary/shruti/SNRIII/pipeline/doubletFinder.R 
#' to remove doublets for merged and filtered samples.
#' 
#' 2.4. Cell Cycle scoring
#' regress the variation due to cell cycle
#' 
#' Genes for cell cycle are picked from two arabidospis papers Yadav et al., 2009
#' and Gutierrez, 2009 as described in Chen et al., 2021, Genome Biology
#' 
#' Start from the normalized pop_singlet data wherein doublets were removed 
#' (called seurat_phase) 23720 x 27459 and continue from this point onwards
#' 
#' Cell cycle genes (sphase and g2m pahse) list for poplar in script CellTypeMarker.R
#' 
#' 
seurat_phase <- CellCycleScoring(seurat_phase, s.features = sphase, 
                                 g2m.features = c(g2phase,mphase), 
                                 set.ident = TRUE)
#' 
#' After scoring for cell cycle, check if cell cycle could be a source of variation 
#' by looking at PCA. To do PCA, first choose the most variable features 
#' (default is 2000), then scale the data.
#' 
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst", 
                                     nfeatures = 2000, mean.cutoff = c(0.01, 8), 
                                     dispersion.cutoff = c(0.5, Inf))
#' 
#' ScaleData() function: 
#' Shifts expression of each gene making mean expression across cells is 0
#' Scales the expression of each gene, so that the variance across cells is 1
#' This step gives equal weight in downstream analyses, 
#' so that highly-expressed genes do not dominate
#' The results of this are stored in object[["RNA"]]@scale.data
#' 
seurat_phase <- ScaleData(seurat_phase, features = rownames(seurat_phase))
#' 
#' Run and view PCA using expression of cell cycle genes
#' 
seurat_phase <- RunPCA(seurat_phase)
DimPlot(seurat_phase, reduction = "pca", split.by = "Phase")
#' 
#' 
#' 3. SCTransform
#' 
#' NOTE: (i) If cells are known to be differentiating and there is clear clustering 
#' differences between G2M and S phases, then regress out by the difference
#' between the G2M and S phase scores, and still differentiating the cycling 
#' from the non-cycling cells.
#' 
#' (ii) ScaleData() function can remove unwanted sources of variation example, 
#' ‘regress out’ heterogeneity associated with cell cycle stage
#' It is recommended to use new normalization- SCTransform
#' 
#' There seem no grouping by cell cycle in the PCA, still I regressed the 
#' variation due to cell cycle
#' 
#' Adjust the memory first
options(future.globals.maxSize = 4000 * 1024^2)
#'
split_seurat <- SplitObject(pop.singlets, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]],
                                        g2m.features=c(g2phase,mphase),
                                        s.features=sphase)
  split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], 
                                            selection.method = "vst", 
                                            nfeatures = 2000)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                   vars.to.regress = c("S.Score", "G2M.Score"),
                                   variable.features.n = 3000)
}
#' 
#' NOTE: 
#' 
#' (i). By default, after normalizing, adjusting the variance, and regressing, 
#' SCTransform will rank genes by residual variance and output 3000 most variant genes.
#' If the dataset has larger cell numbers, then adjust it to higher value
#' using variable.features.n argument.
#' 
#' (ii). We can also regress out protoplasting genes from the dataset (optional)
#' using Yadav_pp gene sets in CellTypeMarker.R file but I have not done this
#' 
#' 
#' 4. Integrate samples using shared highly variable genes
#' 
integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
#' 
#' Run the standard workflow for visualization and clustering
#' specify that downstream analysis will be done on the correct SCT integrated data
#' original unmodified data still resides in the 'RNA' assay
#' 
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50, reduction = "pca")
DimPlot(seurat_integrated)
#'
#'
#' 4.1. Clustering of cells
#' First explore the number of PC to use
#' 
ElbowPlot(object = seurat_integrated, ndims = 40)
#'
#' Find K-nearest neighbor graph
#' 
seurat_integrated <- FindNeighbors(object = seurat_integrated,
                                   dims = 1:50)
#'
#' Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
#'
seurat_integrated@meta.data %>% View()

saveRDS(seurat_integrated, file ="data/SeuratOut/integ.rds")
#'
#'
seurat_integrated <- BuildClusterTree(seurat_integrated, dims = 1:50,
                                      reorder = F, reorder.numeric = F)
PlotClusterTree(seurat_integrated)
#'
#'
#' Use the 0.6 resolution downstream
#' 
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"
#'
#'
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 2)
DimPlot(seurat_integrated, label = TRUE, split.by = "sample")  + NoLegend()
#'
#' Check If cells group by cell cycle phase
#' 
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score")
FeaturePlot(seurat_integrated, reduction = "umap", features = metrics,
            pt.size = 0.4, sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
DimPlot(seurat_integrated, reduction = "umap", split.by = "Phase", label = TRUE,
        label.size = 6)
#'
#' Extract identity and sample info to determine no. of cells per cluster per sample
integ <- readRDS("data/SeuratOut/integ.rds")
n_cells <- FetchData(integ, vars = c("ident")) %>% dplyr::count(ident) %>%
  tidyr::spread(ident, n)
View(n_cells)
#'
#'
#' Remove useless objects
#' 
#' 
#' 5. Markers
#' 5.1. Finding DE markers for each cluster
#' 
set.seed(42)
integ <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"
#' 
split_seurat <- SplitObject(integ, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]
ctrl <- split_seurat$ctrl
#' 
DefaultAssay(ctrl) <- "RNA"
ctrl <- NormalizeData(ctrl, verbose = FALSE)
#' 
#' find and report markers cluster wise for the control sample only in Table S2
#' 
marker <- FindAllMarkers( integ, logfc.threshold = -Inf, test.use = "wilcox", 
                          min.pct = 0.01, min.diff.pct = 0, only.pos = T, 
                          max.cells.per.ident = 20, assay = "RNA")
saveRDS(marker, file = "/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/ctrl_marker.rds")
#' 
background <- rownames(ctrl) # for GO background
saveRDS(background, file = "/mnt/ada/projects/aspseq/htuominen/SNR-results/ctrl_bg.rds")
#' 
markers.sortedByPval = marker[order(marker$p_val),]
genes_uniquely_DE = marker.sortedByPval %>% dplyr::filter(avg_log2FC >= 1) %>%
  group_by(gene) %>% summarize(n=n()) %>%  dplyr::filter(n==1)
genes_uniquely_DE.markers.sortedByPval =
  markers.sortedByPval[markers.sortedByPval$gene %in% genes_uniquely_DE$gene & 
                         markers.sortedByPval$avg_log2FC >= 1,]
top_marker_each = genes_uniquely_DE.markers.sortedByPval %>% 
  dplyr::group_by(cluster) %>% do(head(., n=20))
#'
#' See clusterGOnPlot.R for cluster wise GOs based on markers and their plots
#' 
#' 5.2. conserved markers for each condition can be found like:
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated, ident.1 = cluster,
                       grouping.var = "sample", only.pos = TRUE)}
conserved_markers <- map_dfr(c(0:20), get_conserved)
write.table(conserved_markers, file = "data/SeuratOut/output/afterDbltRemoval/cons_marker.txt",
            sep = "\t", row.names = T, col.names = T)
#'
#'  
#' 6. See degAnalysis.R for pseudobulk DE analysis between the samples
#' and degGOnPlot.R for GO enrichment in DEGs and plots
#' 
#' If you want to regress the protoplasting (optional)
#' 
# DefaultAssay(integ) <- "RNA"
# pp_genes <- readLines("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/pplast.txt")
#' 
#' Identify genes that are present in the Seurat object
# matching_pp_genes <- pp_genes %in% rownames(integ@assays$RNA@counts)
#' 
#' Subset the Seurat object to include only matching genes
#' 
# seurat_matching_pp_genes <- subset(integ, features = pp_genes[matching_pp_genes])
#' 
#' Calculate the percentage of matching genes in each cell
#' 
# percentage_matching_pp_genes <- Matrix::colMeans(
# seurat_matching_pp_genes@assays$RNA@counts) * 100
#' 
#' Create a new feature in the Seurat object to store the percentage information
# integ$percentage_pplast <- percentage_matching_pp_genes
#' 
#' 
#' Finally remove objects from integrated object for NCBI submission for the manuscript
#' 
integ <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"
integ@active.ident <- integ$integrated_snn_res.0.6
DefaultAssay(integ) <- "RNA"
integ$nUMI <- NULL
integ$nGene <- NULL
integ$type <- NULL
integ$pANN_0.25_0.26_4046 <- NULL
integ$pANN_0.25_0.3_3879 <- NULL
integ$doublet_finder <- NULL
integ$seurat_clusters <- NULL
integ$RNA_snn_res.0.1 <- NULL
# submit to NCBI
saveRDS(integ, file="data/SeuratOut/integDiffXyT89SNRIII.rds")
#'
#'
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
#' https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/elbow_plot_metric.md
#'