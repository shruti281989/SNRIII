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

# load seurat onject
integ <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"

# figure 7D
DotPlot(integ, split.by = "sample", cols = c("#fde623ff", "#420051ff"),
        features = "Potra2n15c29002",dot.scale = 10)+RotatedAxis()+ylab(NULL)+
  xlab(NULL)+coord_flip()

# Figure 7B and 7C in the manuscript related to fib and vessel expansion
deg <-  read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/dnStrm/degSNRIII.txt",
                   header = T)

# for fiber Up and expansion related

fibUp <- deg %>% 
  filter(Level == "Upregulated", Cluster %in% c("1", "4", "10", "14", "15")) %>%
  pull("GeneId")
fib_clusters_of_interest <- c("1", "4", "10", "14", "15")

expsn <- readLines("data/SeuratOut/expansion.txt")
fib_features_to_plot <- expsn[expsn %in% fibUp]

# Calculate average expression for each sample's clusters
sample_list <- SplitObject(integ, split.by = "sample")
average_expression_list <- lapply(sample_list, function(x) {
  AverageExpression(x, assays = "RNA", return.seurat = TRUE)
})

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

#  Clear objects except integ seurat object
# similarly for vessel Up and expansion related
vesUp <- deg %>% 
  filter(Level == "Upregulated", Cluster %in% c("6", "14", "17")) %>%
  pull("GeneId")
ves_clusters_of_interest <- c("6", "14", "17")
expsn <- readLines("data/SeuratOut/expansion.txt")
ves_features_to_plot <- expsn[expsn %in% vesUp]

# Calculate average expression for each sample's clusters
sample_list <- SplitObject(integ, split.by = "sample")
average_expression_list <- lapply(sample_list, function(x) {
  AverageExpression(x, assays = "RNA", return.seurat = TRUE)
})

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


#' Restart R
library(dplyr)
library(tidyr)
library(readxl)
library(pheatmap)
library(viridis)
library(tibble)

# for nit metabolism related degs in scrna seq figure 5B
nit <- read.table(
  "/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/nitMetsc.txt",
  header=T)

deg <- read.table(
  "/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/dnStrm/degSNRIII.txt",
  header = TRUE)

df1 <- deg %>% filter(GeneId %in% nit$gene) %>% select(GeneId, Cluster, lfc)

gene_order1 <- unique(nit$gene)

df_sub1 <- df1 %>% group_by(GeneId, Cluster) %>% 
  summarise(lfc = mean(lfc), .groups = "drop")
df_complete1 <- df_sub1 %>%
  complete(GeneId = gene_order, Cluster, fill = list(lfc = 0))

mat_df1 <- df_complete1 %>%
  pivot_wider(names_from = Cluster, values_from = lfc)

rownames(mat_df1) <- mat_df1$GeneId
mat1 <- as.matrix(mat_df1[, -1])
mat1 <- apply(mat1, 2, as.numeric)
rownames(mat1) <- mat_df1$GeneId
mat1[is.na(mat1)] <- 0 # fill NAs

gene_order1 <- gene_order1[gene_order1 %in% rownames(mat1)]
mat1 <- mat1[gene_order1, , drop = FALSE]

pheatmap(mat1, cluster_rows = F, cluster_cols = F, color = viridis(100),
         display_numbers = T, fontsize_row = 7, fontsize_col = 2,
         border_color = NA)

# cell wall modifying genes in Fig 6B 
# file containing lfc of desired genes
df <- read_excel("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/figure6CellExpansionUpreg.xlsx") 

gene_order <- unique(df$gene)

df_sub <- df %>% group_by(gene, cluster) %>% 
  summarise(lfc = mean(lfc), .groups = "drop")
df_complete <- df_sub %>%
  complete(gene = gene_order, cluster, fill = list(lfc = 0))

mat_df <- df_complete %>%
  pivot_wider(names_from = cluster, values_from = lfc)

rownames(mat_df) <- mat_df$gene
mat <- as.matrix(mat_df[, -1])
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- mat_df$gene
mat[is.na(mat)] <- 0 # fill NAs

gene_order <- gene_order[gene_order %in% rownames(mat)]
mat <- mat[gene_order, , drop = FALSE]

pheatmap(mat, cluster_rows = F, cluster_cols = F, color = viridis(100),
         display_numbers = T, fontsize_row = 7, fontsize_col = 2,
         border_color = NA)

# supplementary deg figure in the manuscript S3
# Restart R.
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

lfc <- read.table(
  "/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/dnStrm/degSNRIII.txt",
  header = TRUE)

mat_df <- lfc %>% group_by(GeneId, Cluster) %>%
  summarise(lfc = mean(lfc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cluster, values_from = lfc, values_fill = 0)

mat <- as.matrix(mat_df[, -1])
storage.mode(mat) <- "numeric"
rownames(mat) <- mat_df$GeneId

stopifnot(!is.null(rownames(mat)))

mat <- mat[rowSums(mat != 0) > 0, , drop = FALSE]

lims   <- range(mat, na.rm = TRUE)
breaks <- seq(lims[1], lims[2], length.out = 101)
cols   <- viridis(100)

ph <- pheatmap(mat, cluster_rows = T, cluster_cols = F, show_rownames = F,
  color = cols, breaks = breaks, border_color = NA, fontsize_col = 7, silent = TRUE)

gene_order <- rownames(mat)[ph$tree_row$order]

write.table(gene_order, file = "deg_gene_order.txt", quote = F, row.names = F,
  col.names = F)
