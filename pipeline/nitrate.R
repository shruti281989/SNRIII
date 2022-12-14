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
#' 
#' Firstly, the matrix has cellsxgenes: 1515x37158 (ctrl-kcl) and 2364x37158 (kno)
#' seurat object is created with filter of min. cells = 3 and min.features/genes = 200
#' So, we get for kcl- (1553x21681) and kno-(2402x22393) 
#' then the data is merged and filtered: 
#' 1. umi (<500), genes (<500) and log10GenesPerUMI (>0.8)
#' 2. Another filtration removes 0 count genes and those expressed in <3 cells 
#' 
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
  # library(SCINA)
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
  # library(SingleCellExperiment)
  # library(slingshot)
  # library(scater)
  # library(plotly)
  # library(scran)
  library(here)
})

suppressMessages({
  source(here("UPSCb-common/src/R/gopher.R"))
})
#'
#' 1. Load cellranger output (.h5 matrix) for control (kcl) and treated (kno3)
#' 
setwd("~/shruti/SNRIII/data/SeuratOut/")
ctrl.data <- Read10X_h5("../CellRangerCount/kcl/outs/filtered_feature_bc_matrix.h5")
kno.data <- Read10X_h5("../CellRangerCount/kno/P27752_1002/outs/filtered_feature_bc_matrix.h5")
#'
#' Initialize Seurat object with the non-normalized data and filter:
#' 
cdata <- CreateSeuratObject(counts = ctrl.data, project = "kcl", min.cells = 3,
                            min.features = 200)
tdata <- CreateSeuratObject(counts = kno.data, project = "kno", min.cells = 3,
                            min.features = 200)
#'
#' 2. Standard Pre-processing and plots
#' merge the two samples 
merged_seurat <- merge(x = cdata, y = tdata, 
                       add.cell.id = c("ctrl", "kno"))

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
#'
#' 2.1. Filtration
#' 2.1.1. Cell-level filtering
#' nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
#' This step is crucial: I filtered out low quality reads using following thresholds
#' 
filtered_seurat <- subset(x = merged_seurat, subset= (nUMI >= 500) &
                            (nGene >= 500) & (log10GenesPerUMI > 0.80))
# Following filtration parameter in Chen et al., 2021
# subset= (nUMI >= 500) & (nUMI < 70000) & (nGene >= 200) & (nGene < 9000) & 
# (log10GenesPerUMI > 0.80))

#' 2.1.2. Gene-level filtering
#' Removed genes with zero expression in all cells and keep only those expressed in >=3 cells.
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
#'
#' Sums all TRUE values and returns TRUE if more than 3 TRUE values per gene
# keep_genes <- Matrix::rowSums(nonzero) >= 3
keep_genes <- Matrix::rowSums(nonzero) >= 10
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
#' 
#' Split seurat object by condition to (optional: perform cell cycle scoring)
#' and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]

# for (i in 1:length(split_seurat)) {
#   split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
#   split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], selection.method = "vst", nfeatures = 3000)
#   split_seurat[[i]] <- SCTransform(split_seurat[[i]])
# }

# if you wnat to regress out the cell cycle variation
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=c(g2phase,mphase), s.features=sphase)
  split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], selection.method = "vst", nfeatures = 3000)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score", "G2M.Score"))
}
# 
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
# 
# I've also performed cell cycle regression and saved:
# saveRDS(seurat_integrated, "big/integ_seurat_cycle.rds")
#'
#' UMAP visualization
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50, reduction = "pca")
DimPlot(seurat_integrated)
#'
saveRDS(seurat_integrated, "integrated.rds")
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
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 4)
#'
#' Explore QC if needed
#' Some useful data
DimPlot(seurat_integrated, label = TRUE, split.by = "sample")  + NoLegend()
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score")
FeaturePlot(seurat_integrated, reduction = "umap", features = metrics,
            pt.size = 0.4, sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

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
# DefaultAssay(seurat_integrated) <- "RNA"
# FeaturePlot(seurat_integrated, reduction = "umap", features = RayParen, 
#             sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

#' Normalize RNA data for visualization
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
#' 
#' 5. Finding DE markers for each cluster
#'
markers <- FindAllMarkers(object = seurat_integrated, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 1, 
                          test.use = 'bimod', min.diff.pct = 0.25)

markWilcox <- FindAllMarkers(object = seurat_integrated, only.pos = TRUE,
                             logfc.threshold = 0.58, min.pct = 0.25)

markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC) 

markers %>% group_by(cluster) %>% slice(5) %>% pull(gene) -> best2.bimod.gene.per.cluster

markers.sortedByPval = markers[order(markers$p_val),]
genes_uniquely_DE = markers.sortedByPval %>% 
  dplyr::filter(avg_log2FC >= 1) %>% group_by(gene) %>%  
  summarize(n=n()) %>%  dplyr::filter(n==1)

genes_uniquely_DE.markers.sortedByPval =
  markers.sortedByPval[markers.sortedByPval$gene %in% genes_uniquely_DE$gene & 
                                      markers.sortedByPval$avg_log2FC >= 1,]

top_marker_each = genes_uniquely_DE.markers.sortedByPval %>% 
  dplyr::group_by(cluster) %>% do(head(., n=5))

write.table(top_marker_each, file = "top_marker.txt", sep = "\t",
            row.names = F, col.names = T)
png("plots/dot/top_marker_first.png", width = 2480, height = 980,
    pointsize = 12, bg = "transparent")
for (i in length(top_marker_each$gene)){
  print(DotPlot(seurat_phase, features = top_marker_each$gene)+RotatedAxis())
}
dev.off()
#' 
#' Finding conserved markers
Uniannot <- read.delim(file= "../../reference/Potra22_annotation.tsv",sep = "\t")
geneannot <- read.delim(file = "../../reference/gene_info.txt", sep = "\t")

get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) 
  # %>%
    # rownames_to_column(var = "gene") %>%
    # left_join(y = unique(annotations[, c("Potra22_id", "UniRef_annotation")]),
    #           by = c("gene" = "Potra22_id")) %>%
    # cbind(cluster_id = cluster, .)
}
#'
conserved_markers <- map_dfr(c(0:11), get_conserved)

top10_cons <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_log2FC + kno_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

write.table(top10_cons, file = "top10_cons_marker_annot.txt", sep = "\t",
            row.names = F, col.names = T)
#'
#' Rename clusters 
# Double check the following 
seurat_labelled <- RenameIdents(object = seurat_integrated, 
                                "0" = "Unknown", "1" = "Camb1", "2" = "RayPrec", "3" = "Camb2", "4" = "Camb3",
                                "5" = "Fibre", "6" = "Ray", "7" = "RayParenc","8" = "Vessel", "9" = "Ray2",
                                "10" = "Camb4","11" = "XylsideProlifCambium")
                                # "0" = "ILFib", "1" = "RFOrg", "2" = "FibORayPre",
                                # "3" = "LFibVes", "4" = "IFib", "5" = "Fibre",
                                # "6" = "RayParPre", "7" = "RayPar","8" = "Ves",
                                # "9" = "RayPar2","10" = "ILFibVes","11" = "XylProlCamb")

DimPlot(object = seurat_integrated, reduction = "umap", label = TRUE, split.by = "sample",
        label.size = 3,repel = TRUE)+ NoLegend()

seurat_labelled$celltype.sample <- paste(Idents(seurat_labelled), seurat_labelled$sample,
                                      sep = "_")


#' 7. See DE_analysis_scrnaseq.R for pseudobulk DE analysis and bulkDE analysis (Part1 and 2, respectively)
#' 
#' 8. See the end of ./big/combined_cell_cycle.R for pseudo time trajectory: it works for only single sample
#' Working with pseudotime in trajectory.R
#' 
#' 9. GO enrichment
deg.ls <- split(rownames(markWilcox), f = markWilcox$cluster)
# deg.ls is a list, still need to make a list for enrichment
gene.ls <- list(deg.ls)

#' Background to be used: Suggestions form Nico an Nat
#' use only the set of genes expressed in the tissue you are looking at
#' either devise that across all of your samples:
#' naïve approach: exclude genes that are always 0
#' more advanced, check the average expression independent filtering cutoff available in your DESeq results - if you use DESeq:
#' Cutoff off the independent filtering
#' cutoff <- metadata(res)$filterThreshold
#' We keep the rest of the data
# sel <- rowMeans(counts(dds_tx,normalized=TRUE)) > cutoff`
# independent filtering is how DESeq decides on the signal-to-noise ratio
# 2. your data is to sparse for 1. to work then use the set of genes expressed 
# in your tissue type (wood) from a bulk RNA resource

#' When I use rownames(filtered_seurat), the enrichment doesn't work
#' 
enr.list <- lapply(gene.ls,function(r){
  lapply(r,gopher,task=list("go","kegg","pfam"),background = iiit,url="potra2")
})

# Code from Aman about visualization
list <- enr.list[[1]]
df <- bind_rows(list$`0`$go, list$`1`$go, list$`2`$go, list$`3`$go, list$`4`$go, list$`5`$go, list$`6`$go, list$`7`$go, list$`8`$go, list$`9`$go, list$`10`$go, list$`11`$go, .id = "Cluster")
df %<>% rowwise() %>% mutate(Cluster = as.numeric(Cluster) - 1, .keep = "unused") 

BP <- df %>% filter(namespace %in% "BP")
MF <- df %>% filter(namespace %in% "MF")
CC <- df %>% filter(namespace %in% "CC")

#plotting function
pp <- function(df, go) {
  df %>% 
    ggplot(aes(x = factor(Cluster), y = name, color = padj, size = nt)) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    theme(text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) + 
    ylab(paste("GO Terms -", go)) + 
    xlab("Clusters") + 
    ggtitle("GO analysis")
}

#downloading function
dl <-  function(go) {
  pt %>% 
    download_this(output_name = go,
                  output_extension = ".pdf",
                  button_label = "Download plot",
                  button_type = "warning",
                  has_icon = TRUE,
                  icon = "fa fa-save")
} 
BP
MF
CC

pt <- pp(df = BP, go = "BP")
pt
# for i in length()
as.data.frame(enr.list[[1]]$`0`$go) %>% 
  ggplot(aes(x= id, y = name, color = padj, size = nt)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis")
