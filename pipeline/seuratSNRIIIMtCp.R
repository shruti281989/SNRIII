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
  library(clustree)
})

# 1. Basic filtering and merging
# Load in the data with mtcp genes
ctrl1.data <- Read10X_h5("~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/kcl2/outs/filtered_feature_bc_matrix.h5",
                         use.names = T)
# 37254 x 20581 

kno1.data <- Read10X_h5("~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/kno2/outs/filtered_feature_bc_matrix.h5",
                        use.names = T)
# 37254 x 13605

# Create seurat objects
cdata1 <- CreateSeuratObject(counts = ctrl1.data, project = "kclmtcp", 
                             min.cells = 3, min.features = 200)# min.cells = 0
# 37254 x 20581; 27695x 20581

tdata1 <- CreateSeuratObject(counts = kno1.data, project = "knomtcp", 
                             min.cells = 3, min.features = 200)# min.cells = 0
# 37254 x13605; 26931x13605

# Create a merged object from ctrl and treated (kno3)
merged_seurat <- merge(cdata1, y = tdata1, 
                       add.cell.ids = c("kclmtcp", "knomtcp"), project = "SNRIII")
# 37254 x 34186 ;28212 x 34186 

# Reomve unwanted objects
rm(cdata1,tdata1,ctrl1.data,kno1.data)
gc()

# Calculate QC
# log10 Genes per UMI
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / 
  log10(merged_seurat$nCount_RNA)

# Mito abundance 
merged_seurat[["percent_mito"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT")
 
# Chloroplast abundance
merged_seurat[["percent_cp"]] <- PercentageFeatureSet(merged_seurat, pattern = "^CP")

# 'ribosomal_genes' abundance
# ribosomal_genes <- readLines("analysis/markers/Ribo_gene_table.txt")
# matching_genes <- rownames(merged_seurat@assays$RNA@counts)[rownames(merged_seurat@assays$RNA@counts) %in% ribosomal_genes]
# seurat_ribo <- subset(merged_seurat,features = matching_genes)
# percent_ribo <- Matrix::colMeans(seurat_ribo@assays$RNA@counts) * 100
# merged_seurat$percent_ribo <- percent_ribo

# Rename columns
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)

# Add sample info
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^kclmtcp_"))] <- "kclmtcp"
metadata$sample[which(str_detect(metadata$cells, "^knomtcp_"))] <- "knomtcp"

merged_seurat@meta.data <- metadata
rm(metadata, ribosomal_genes, matching_genes, seurat_ribo, percent_ribo)
gc()

# PlotQC
feats <- c("nGene", "nUMI", "percent_mito", "percent_cp", "percent_ribo")
VlnPlot(merged_seurat, split.by = "sample", features = feats, pt.size = 0.1, 
        ncol = 3)

# Plot other correlations between the factors
# FeatureScatter(merged_seurat, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident",
#                pt.size = .5)
# # both percent.mt and cp are negatively correlated with both nCount and nFeature
# # both percent.mt and cp are very lowly correlated to each other

# Restart R

# 1. Cell and gene-level filtering
# Feb 27, 2024
filt0.9_seurat <- subset(merged_seurat, subset= (nCount_RNA >= 500) &(nFeature_RNA >= 300)&
                            (nFeature_RNA < 7000) & (nCount_RNA < 50000) &
                            (percent_mito < 2) & (percent_cp < 0.1) &
                            (log10GenesPerUMI > 0.9))

VlnPlot(filt0.9_seurat, split.by = "sample", features = feats, pt.size = 0.1, 
        ncol = 3)

filt_counts <- GetAssayData(object = filt0.9_seurat, slot = "counts")
nonzero <- filt_counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 3
summary (keep_genes)
filt_counts <- filt_counts[keep_genes, ]

filtered_seurat <- CreateSeuratObject(filt_counts,
                                      meta.data = filt0.9_seurat@meta.data)
rm(filt_counts, keep_genes, nonzero, filt0.9_seurat, merged_seurat)
gc()
VlnPlot(filtered_seurat, split.by = "sample", features = feats, pt.size = 0.1, 
        ncol = 3)

# 2. Remove doublets
# for kcl
pop.split <- SplitObject(filtered_seurat, split.by = "sample")
pop.sample <- NormalizeData(pop.split$kclmtcp)
rm(pop.split)
gc()

pop.sample <- FindVariableFeatures(pop.sample)
pop.sample <- ScaleData(pop.sample, 
                        vars.to.regress = c("nFeature_RNA", "percent_mito", "percent_cp", "log10GenesPerUMI"))
pop.sample <- RunPCA(pop.sample, nfeatures.print = 2, npcs = 50)

stdv <- pop.sample[["pca"]]@stdev
sum.stdv <- sum(pop.sample[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min.pc <- min(co1, co2) #20

pop.sample <- RunUMAP(pop.sample, dims = 1:min.pc)
# pop.sample <- FindNeighbors(object = pop.sample, dims = 1:min.pc)              
# pop.sample <- FindClusters(object = pop.sample, resolution = 0.1)

# pK identification (no ground-truth)
sweep.list <- paramSweep_v3(pop.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
nExp2 <- round(optimal.pk * nrow(pop.sample@meta.data))
nExp <- round(ncol(pop.sample) * 0.076) # expect 7.6% doublets for log8

pop.sample <- doubletFinder_v3(seu = pop.sample, PCs = 1:min.pc, 
                               pK = optimal.pk, nExp = nExp)

DF.name = colnames(pop.sample@meta.data)[grepl("DF.classification", colnames(pop.sample@meta.data))]
kclmtcp= pop.sample[, pop.sample@meta.data[, DF.name] == "Singlet"] #pk calculated manually
# 10472x 27142

kclmtcp2= pop.sample[, pop.sample@meta.data[, DF.name] == "Singlet"] 
# 8047x 27142 

VlnPlot(pop.sample, features = "nFeature_RNA", group.by = DF.name, pt.size = .1)

# remove unwanted objects
rm (stdv,sum.stdv, percent.stdv, cumulative, co1, co2, sweep.list, min.pc,
    sweep.stats, bcmvn, bcmvn.max, optimal.pk, nExp, pop.sample, DF.name)
gc()

saveRDS(kclmtcp, file = "~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/mtCpSeurat/kclmtcp.rds")
saveRDS(kclmtcp2, file = "~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/mtCpSeurat/kclmtcp2.rds")

# Run dblt finder for kno
pop.split <- SplitObject(filtered_seurat, split.by = "sample")
pop.sample <- NormalizeData(pop.split$knomtcp)
rm(pop.split)
gc()

pop.sample <- FindVariableFeatures(pop.sample)
pop.sample <- ScaleData(pop.sample, 
                        vars.to.regress = c("nFeature_RNA", "percent_mito", "percent_cp", "log10GenesPerUMI"))
pop.sample <- RunPCA(pop.sample, nfeatures.print = 2, npcs = 50)

stdv <- pop.sample[["pca"]]@stdev
sum.stdv <- sum(pop.sample[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min.pc <- min(co1, co2) #19

pop.sample <- RunUMAP(pop.sample, dims = 1:min.pc)
# pop.sample <- FindNeighbors(object = pop.sample, dims = 1:min.pc)              
# pop.sample <- FindClusters(object = pop.sample, resolution = 0.1)
# 
# pK identification (no ground-truth)
sweep.list <- paramSweep_v3(pop.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
nExp2 <- round(optimal.pk * nrow(pop.sample@meta.data))
# nExp <- round(ncol(pop.sample) * 0.076) # expect 7.6% doublets for log8

pop.sample <- doubletFinder_v3(seu = pop.sample, PCs = 1:min.pc, 
                               pK = optimal.pk, nExp = nExp2)
DF.name = colnames(pop.sample@meta.data)[grepl("DF.classification", colnames(pop.sample@meta.data))]
knomtcp= pop.sample[, pop.sample@meta.data[, DF.name] == "Singlet"]
# knomtcp -11190x27142; 8477x27142
saveRDS(knomtcp, file = "~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/mtCpSeurat/knomtcp.rds")
saveRDS(knomtcp, file = "~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/mtCpSeurat/knomtcp2.rds")
knomtcp
rm (stdv,sum.stdv, percent.stdv, cumulative, co1, co2, min.pc, sweep.list,
    sweep.stats, bcmvn, bcmvn.max, optimal.pk, DF.name, pop.split, min.pc,
    nExp, pop.sample)
gc()

# 4. run integration and clustering
# Restart R
# combine the singlets
pop.singlets <- merge(x = kclmtcp, y = knomtcp, project = "sngltSNRIII")
pop.singlets$pANN_0.25_0.29_861 <- NULL
pop.singlets$pANN_0.25_0.27_920 <- NULL
pop.singlets$DF.classifications_0.25_0.29_861 <- NULL
pop.singlets$DF.classifications_0.25_0.27_920 <- NULL

pop.singlets <- merge(x = kclmtcp2, y = knomtcp2, project = "sngltSNRIII")
pop.singlets$pANN_0.25_0.29_3287 <- NULL
pop.singlets$pANN_0.25_0.3_3633 <- NULL
pop.singlets$DF.classifications_0.25_0.29_3287 <- NULL
pop.singlets$DF.classifications_0.25_0.3_3633 <- NULL

rm(kclmtcp, knomtcp, kclmtcp2, knomtcp2)
gc()

# 3. check cell cycle phasing and mtio contribution
# Cell cycle genes (Zhang et al., 2021)
# g01phase <- c("Potra2n10c22352","Potra2n11c22967","Potra2n13c25403","Potra2n14c27234","Potra2n14c27322","Potra2n18c32466","Potra2n1c2151","Potra2n1c2633","Potra2n1c2837","Potra2n2c4043","Potra2n5c11567","Potra2n6c14588","Potra2n9c19242")
# g2phase <- c("Potra2n10c21070","Potra2n10c21695","Potra2n10c21833","Potra2n10c22093","Potra2n10c22145","Potra2n10c22346","Potra2n11c22815","Potra2n11c23463","Potra2n12c24927","Potra2n13c25321","Potra2n13c25579","Potra2n13c26121","Potra2n13c26208","Potra2n14c26403","Potra2n14c27462","Potra2n14c27475","Potra2n15c27972","Potra2n15c28096","Potra2n15c28656","Potra2n15c28758","Potra2n15c28945","Potra2n16c29301","Potra2n16c29689","Potra2n16c29960","Potra2n16c30306","Potra2n16c30563","Potra2n18c32638","Potra2n18c32677","Potra2n18c32912","Potra2n19c34004","Potra2n19c34102","Potra2n1c1559","Potra2n1c2018","Potra2n1c2288","Potra2n1c2573","Potra2n1c2758","Potra2n1c2942","Potra2n1c3087","Potra2n1c3141","Potra2n1c3180","Potra2n1c3295","Potra2n1c3611","Potra2n1c410","Potra2n2c4140","Potra2n2c4292","Potra2n2c4299","Potra2n2c4587","Potra2n2c4715","Potra2n2c4966","Potra2n2c5259","Potra2n2c5415","Potra2n2c5455","Potra2n2c5494","Potra2n2c5545","Potra2n2c5663","Potra2n2c6426","Potra2n372s35460","Potra2n3c7249","Potra2n3c7951","Potra2n3c7972","Potra2n3c8097","Potra2n3c8273","Potra2n4c8773","Potra2n4c9741","Potra2n5c10510","Potra2n5c10575","Potra2n5c10579","Potra2n5c11241","Potra2n5c11640","Potra2n5c12275","Potra2n680s36426","Potra2n6c14644","Potra2n6c14653","Potra2n6c14979","Potra2n6c15055","Potra2n6c15109","Potra2n7c15462","Potra2n7c16145","Potra2n8c17814","Potra2n8c17973","Potra2n8c18569","Potra2n8c18623","Potra2n9c18760","Potra2n9c19217","Potra2n9c19378","Potra2n9c20077")
# sphase <- c("Potra2n10c20260","Potra2n10c20508","Potra2n11c22684","Potra2n13c26037","Potra2n13c26116","Potra2n14c27103","Potra2n14c27106","Potra2n15c28275","Potra2n15c28844","Potra2n16c29833","Potra2n17c31198","Potra2n18c33119","Potra2n1c101","Potra2n1c2877","Potra2n1c3343","Potra2n1c3610","Potra2n2c4073","Potra2n2c4499","Potra2n2c5307","Potra2n2c6056","Potra2n3c7036","Potra2n3c7597","Potra2n4c9636","Potra2n5c10911","Potra2n5c12603","Potra2n6c13036","Potra2n6c14627","Potra2n7c15963","Potra2n7c16633","Potra2n8c16885","Potra2n8c17847","Potra2n9c19886")
# mphase <- c ("Potra2n10c20592","Potra2n10c20846","Potra2n10c21494","Potra2n10c21734","Potra2n10c22098","Potra2n11c22486","Potra2n11c22945","Potra2n11c23018","Potra2n12c24163","Potra2n13c25058","Potra2n13c25262","Potra2n13c25822","Potra2n13c25852","Potra2n14c27365","Potra2n14c27462","Potra2n14c27797","Potra2n15c29215","Potra2n16c29450","Potra2n16c30378","Potra2n16c30406","Potra2n16c30474","Potra2n17c30910","Potra2n17c31100","Potra2n17c31123","Potra2n18c32687","Potra2n19c34338","Potra2n1c2214","Potra2n1c2376","Potra2n1c2511","Potra2n1c2591","Potra2n1c3910","Potra2n1c4000","Potra2n1c671","Potra2n1c923","Potra2n2c4070","Potra2n2c4172","Potra2n3c7249","Potra2n3c7843","Potra2n4c8944","Potra2n4c8987","Potra2n5c10575","Potra2n5c10858","Potra2n5c10979","Potra2n5c11566","Potra2n5c11796","Potra2n5c11923","Potra2n5c12554","Potra2n6c12962","Potra2n6c13237","Potra2n6c13440","Potra2n6c13456","Potra2n6c14169","Potra2n6c15109","Potra2n7c15465","Potra2n7c16119","Potra2n8c17128","Potra2n8c17171","Potra2n8c17318","Potra2n8c18190","Potra2n8c18202","Potra2n8c18287","Potra2n9c19378","Potra2n9c19459","Potra2n9c19549","Potra2n6c13657", "Potra2n16c29642","Potra2n17c30791", "Potra2n6c14415")

split_seurat <- SplitObject(pop.singlets, split.by = "sample")
split_seurat <- split_seurat[c("kclmtcp", "knomtcp")]
rm(pop.singlets)
gc()

options(future.globals.maxSize = 4000 * 1024^2)

for (i in 1:length(split_seurat)) {
  # split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  # split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]],
  #                                       g2m.features=c(g2phase,mphase),
  #                                       s.features=sphase)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], variable.features.n = 3000,
                                   vst.flavor = "v2",
                                   vars.to.regress = c("percent_mito", "percent_cp", "log10GenesPerUMI"))
}

integ_features <- SelectIntegrationFeatures(split_seurat, nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

rm(split_seurat, i, integ_anchors, integ_features)
gc()

DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(seurat_integrated)
ElbowPlot(seurat_integrated, ndims = 50) #17  pc

# pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
# cumu <- cumsum(pct)
# co1 <- which(cumu > 90 & pct < 5)[1]
# co1
# co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# co2
# pcs <- min(co1, co2)
# pcs

integ <- RunUMAP(seurat_integrated, dims = 1:30, reduction = "pca")
integ <- FindNeighbors(integ, dims = 1:30)
integ <- FindClusters(integ, resolution = c(0.4, 0.6))

Idents(integ) <- integ$integrated_snn_res.0.6
integ <- RunUMAP(integ, dims = 1:30, reduction = "pca")
DimPlot(integ, reduction = "umap", label = TRUE, label.size = 2)

n_cells <- FetchData(integ, vars = c("ident")) %>%
  dplyr::count(ident) %>%
  tidyr::spread(ident, n)
clustree(integ, prefix = "integrated_snn_res.")

integ$seurat_clusters <- NULL

#' 5. Finding DE markers for each cluster
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
# Note: raw and normalized counts are stored in counts and data slots of RNA assay. 
markers <- FindAllMarkers(object = seurat_integrated,
                          only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25)
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

# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/elbow_plot_metric.md
