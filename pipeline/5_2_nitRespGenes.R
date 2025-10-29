setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(readxl)
  library(gplots)
})

set.seed(42)

# load seurat
integ <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"

# In case to plot for one sample, split the object
# split_seurat <- SplitObject(integ, split.by = "sample")
# split_seurat <- split_seurat[c("ctrl", "kno")]
# and plot using the object split_seurat$ctrl 

# Load markers
# selectedmarker <- read.table("/mnt/picea/home/schoudhary/SNR-u2023011/analysis/markers/selMarker.txt",
#                              sep = '\t',header=TRUE)
clustMarker <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/clustMarker.txt",
                          sep = '\t',header=TRUE)
# ploidy <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/endo_expPotra.txt",
#                      sep = '\t',header=TRUE)
# ribosomal <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/Ribo_gene_table.txt",
#                         sep = '\t',header=F)
# protoplasting <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/pplast.txt",
#                             sep = '\t',header=F)
degWilcox <- read.table("data/SeuratOut/output/afterDbltRemoval/markerWilcox_lfc1_fdr0.01_pct0.1.txt",
                        header = T)
nit <- read_excel("data/SeuratOut/bulkDegvsScDeg.xlsx", sheet = 2)

# degTime <- read.table("data/degTableS1E.txt",
#                        header = T, sep = '\t')

#Cluster wise heatmap: scale the data first: 
# integ <- ScaleData(integ, features = rownames(integ))
# DoHeatmap(integ, features = degAll$Potra, group.by="sample",
          # group.colors = viridis(100))

# # for label with own/AT annotation 
# nitRes <- nitRes[order(nitRes$AtName),]
# # nitRes <- nitRes[order(nitRes$Type),]
# # extract the genes first, since gene ids are not unique
# nrt <- nitRes[nitRes$Type == "Expansion",]
# markeranno <- paste(nrt$AtName)
# names(markeranno) <- nrt$GeneId

# DotPlot(subset(integ, subset = sample =="ctrl"), features = nrt$GeneId,
#         dot.scale = 10,idents=c("Fiber 2", "Fiber 1", "Early Fiber", 
#                                 "Fiber Precursor 2","Fiber Precursor 1",
#                                 "Late Vessel", "Early Vessel","Fusiform Initial",
#                                 "Ray 18","Ray 7","Ray 5","Ray/Fusiform Initial",
#                                 "Phloem like","Cambium"))+RotatedAxis()+ 
#   scale_x_discrete(labels = markeranno)+ylab(NULL)+xlab(NULL)+
#   scale_colour_gradient(low = "grey", high = "darkblue")
# 
# DotPlot(integ, selectedmarker[selectedmarker$CellType == "chen photosyn", ]$GeneId,
#         cols ="RdBu", dot.scale = 10) + RotatedAxis()+ coord_flip()
# # cols = c("#005AB5", "#DC3220"),
# DotPlot(integ, split.by = "sample", cols = c("#fde623ff", "#420051ff"),
#         features = "Potra2n15c29002",dot.scale = 10)+RotatedAxis()+ylab(NULL)+
#   xlab(NULL)+coord_flip()
#scale_colour_gradient(low = "grey", high = "#420051ff")
# subset(integ, subset = sample =="kno"
#c("#420051ff", "#fde623ff"),

# Expansion related: feronia, ralf, lrx, md 
# fer <- c("Potra2n16c30521", "Potra2n16c30523", "Potra2n6c14356")
# ralf <- c("Potra2n13c26153","Potra2n14c27428","Potra2n17c31433","Potra2n1c2781",
#           "Potra2n267s35162","Potra2n4c8790","Potra2n4c8791","Potra2n5c12610")
# lrx8 <- c("Potra2n1c2637", "Potra2n6c13964", "Potra2n6c14636", "Potra2n9c19254",
#           "Potra2n14c26497")
# md <- c("Potra2n14c27025","Potra2n6c14356","Potra2n8c17700","Potra2n16c30033",
#         "Potra2n10c21287")

# only metabolism genes in manuscript (Figure in the manuscript)
nitRes <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/fig3marker.txt",
                     sep = '\t',header=TRUE)
rownames(nitRes) <- nitRes$Potra

# DotPlot(integ, split.by = "sample", cols = c("#fde623ff", "#420051ff"),
#         features = nit$Potra, dot.scale = 8)+
#   RotatedAxis()+ylab(NULL)+xlab(NULL)+coord_flip()
# DotPlot(subset(integ, subset = sample =="kno"),
#         features = nit$Potra,dot.scale = 10)+
#         RotatedAxis()+ylab(NULL)+xlab(NULL)+
#           scale_colour_gradient(low = "grey", high = "#420051ff")+coord_flip()

integ <- ScaleData(integ, features = rownames(integ), assay = "RNA")

# use the following to draw heatmap for upregulated and
# downregulated genes in each cluster for each sample (supp Figure in the manuscript)

# features_to_plot <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/dnStrm/degSNRIII.txt",
#                                header = TRUE, fill = TRUE, sep = "\t", quote = "")
# features_to_plot <- features_to_plot %>% filter(status=="up in kno" & FDR <0.05) %>%
#   pull("GeneID")

# DoHeatmap(avgexp, features = features_to_plot,angle = 0) + guides(color="none") +
#   scale_fill_gradientn(colors = viridis(20))+
#   theme(axis.text.y = element_text(size = 1))

# use the following to draw heatmap for nitrogen metabolism genes
# (Figure in the manuscript)
features_to_plot <- nit$Potra
avgexp = AverageExpression(subset(integ, subset = sample =="ctrl"), assay="RNA",
                           features = nitRes$Potra,
                           return.seurat = T, group.by = 'integrated_snn_res.0.6')

# Calculate average expression for each sample's clusters
sample_list <- SplitObject(integ, split.by = "sample")
average_expression_list <- lapply(sample_list, function(x) {
  AverageExpression(x, assays = "RNA", return.seurat = TRUE)
})

# Combine data for each sample and filter for the list of genes
average_expression_combined <- lapply(names(average_expression_list), function(sample) {
  as.data.frame(GetAssayData(average_expression_list[[sample]], slot = "data")) %>%
    rownames_to_column("Gene") %>%
    filter(Gene %in% features_to_plot) %>%
    mutate(Sample = sample)
}) %>% bind_rows()

# Scale expression gene-wise for all genes
heatmap_data <- lapply(names(average_expression_list), function(sample) {
  as.data.frame(GetAssayData(average_expression_list[[sample]], slot = "data")) %>%
    rownames_to_column("Gene") %>%
    filter(Gene %in% features_to_plot) %>%
    pivot_longer(-Gene, names_to = "Cluster", values_to = "Expression") %>%
    group_by(Gene) %>%
    mutate(ScaledExpression = scale(Expression), Sample = sample) %>%
    ungroup()
}) %>% bind_rows()

# Aggregate duplicates by taking the mean
heatmap_data_clean <- heatmap_data %>%
  group_by(Gene, Cluster, Sample) %>%
  summarise(ScaledExpression = mean(ScaledExpression), .groups = 'drop')

# Extract and save gene order from clustering in kno
reference_sample <- unique(heatmap_data_clean$Sample)[2]  # kno as reference
reference_data <- heatmap_data_clean %>%
  filter(Sample == reference_sample) %>%
  select(Gene, Cluster, ScaledExpression) %>%
  pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Check for NA/NaN/Inf values in the reference_data
na_check <- apply(reference_data, 1, function(x) any(is.na(x) | is.infinite(x)))
clean_reference_data <- reference_data[!na_check, ]

gene_order <- hclust(dist(clean_reference_data))$order
ordered_genes <- rownames(clean_reference_data)[gene_order]

write.table(ordered_genes, paste0("gene_order_", reference_sample, ".txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# remove these genes from heatmap_data
heatmap_data_clean <- heatmap_data_clean %>% filter(Gene %in% ordered_genes)

# Generate heatmaps for two samples using the same gene order in manuscript
for (sample in unique(heatmap_data_clean$Sample)) {
  sample_data <- heatmap_data_clean %>%
    filter(Sample == sample) %>%
    select(Gene, Cluster, ScaledExpression) %>%
    pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # Reorder genes based on the reference sample
  sample_data <- sample_data[ordered_genes, ]
  
  svg(paste0("heatmap_", sample, ".svg"), width = 10, height = 8)
  heatmap.2(
    sample_data, scale = "none", col = viridis(20), trace = "none",
    margins = c(8, 8), dendrogram = "none", Rowv = F, Colv = F,
    key = TRUE, cexCol = 1, cexRow = 0.1
  )
  dev.off()
}

# for fiber Up and expansion (Figure in the manuscript)
fibUp <- deg %>% 
  filter(Level == "Upregulated", Cluster %in% c("1", "4", "10", "14", "15")) %>%
  pull("GeneId")
fib_clusters_of_interest <- c("1", "4", "10", "14", "15")
expsn <- readLines("data/Seurat_Out/expansion.txt")
fib_features_to_plot <- expsn[expsn %in% fibUp]

avg_expr_list <- lapply(sample_list, function(x) {
  x <- subset(x, idents = fib_clusters_of_interest)
  AverageExpression(x, assays = "RNA", return.seurat = TRUE)
})

heatmap_data <- lapply(names(avg_expr_list), function(sample) {
  expr <- GetAssayData(avg_expr_list[[sample]], slot = "data") %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    filter(Gene %in% fib_features_to_plot) %>%
    pivot_longer(-Gene, names_to = "Cluster", values_to = "Expression") %>%
    group_by(Gene) %>%
    mutate(ScaledExpression = scale(Expression), Sample = sample) %>%
    ungroup()
}) %>% bind_rows()

heatmap_data <- heatmap_data %>%
  group_by(Gene, Cluster, Sample) %>%
  summarise(ScaledExpression = mean(ScaledExpression), .groups = "drop")

ref_sample <- "kno"
ref_mat <- heatmap_data %>%
  filter(Sample == ref_sample) %>%
  pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

ref_mat <- ref_mat[complete.cases(ref_mat), ]
gene_order <- hclust(dist(ref_mat))$order
ordered_genes <- rownames(ref_mat)[gene_order]

heatmap_data <- heatmap_data %>% filter(Gene %in% ordered_genes)
for (sample in unique(heatmap_data$Sample)) {
  mat <- heatmap_data %>%
    filter(Sample == sample) %>%
    pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
    column_to_rownames("Gene")
  
  # fill missing ones with NA or 0
  missing_clusters <- setdiff(fib_clusters_of_interest, colnames(mat))
  if (length(missing_clusters) > 0) {
    mat[missing_clusters] <- NA  # or 0 if you prefer
  }
  
  # Reorder columns by cluster number
  mat <- mat[, sort(colnames(mat))]
  
  # Reorder rows by ordered genes
  mat <- mat[ordered_genes, , drop = FALSE]
  
  # handle NA
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0  
  svg(paste0("heatmap_", sample, ".svg"), width = 10, height = 8)
  heatmap.2(mat, scale = "none", col = viridis(20), trace = "none",
            margins = c(8, 8), dendrogram = "none", Rowv = FALSE, Colv = FALSE,
            key = TRUE, cexCol = 1, cexRow = 0.2)
  dev.off()
}

write.table(ordered_genes, paste0("gene_order_", ref_sample, ".txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# for vessel Up and expansion (Figure in the manuscript)
vesUp <- deg %>% 
  filter(Level == "Upregulated", Cluster %in% c("6", "14", "17")) %>%
  pull("GeneId")
ves_clusters_of_interest <- c("6", "14", "17")
expsn <- readLines("data/SeuratOut/expansion.txt")
ves_features_to_plot <- expsn[expsn %in% vesUp]

avg_expr_list <- lapply(sample_list, function(x) {
  x <- subset(x, idents = ves_clusters_of_interest)
  AverageExpression(x, assays = "RNA", return.seurat = TRUE)
})

heatmap_data <- lapply(names(avg_expr_list), function(sample) {
  expr <- GetAssayData(avg_expr_list[[sample]], slot = "data") %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    filter(Gene %in% ves_features_to_plot) %>%
    pivot_longer(-Gene, names_to = "Cluster", values_to = "Expression") %>%
    group_by(Gene) %>%
    mutate(ScaledExpression = scale(Expression), Sample = sample) %>%
    ungroup()
}) %>% bind_rows()

heatmap_data <- heatmap_data %>%
  group_by(Gene, Cluster, Sample) %>%
  summarise(ScaledExpression = mean(ScaledExpression), .groups = "drop")

ref_sample <- "kno"
ref_mat <- heatmap_data %>%
  filter(Sample == ref_sample) %>%
  pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

ref_mat <- ref_mat[complete.cases(ref_mat), ]
gene_order <- hclust(dist(ref_mat))$order
ordered_genes <- rownames(ref_mat)[gene_order]

heatmap_data <- heatmap_data %>% filter(Gene %in% ordered_genes)

# heatmap keeping the same gene order
for (sample in unique(heatmap_data$Sample)) {
  mat <- heatmap_data %>%
    filter(Sample == sample) %>%
    pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
    column_to_rownames("Gene")
  
  # fill missing ones with NA or 0
  missing_clusters <- setdiff(fib_clusters_of_interest, colnames(mat))
  if (length(missing_clusters) > 0) {
    mat[missing_clusters] <- NA  # or 0 if you prefer
  }
  
  # Reorder columns by cluster number
  mat <- mat[, sort(colnames(mat))]
  
  # Reorder rows by ordered genes
  mat <- mat[ordered_genes, , drop = FALSE]
  
  # handle NA
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0  
  svg(paste0("heatmap_", sample, ".svg"), width = 10, height = 8)
  heatmap.2(mat, scale = "none", col = viridis(20), trace = "none",
            margins = c(8, 8), dendrogram = "none", Rowv = FALSE, Colv = FALSE,
            key = TRUE, cexCol = 1, cexRow = 0.2)
  dev.off()
}

write.table(ordered_genes, paste0("gene_order_", ref_sample, ".txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
