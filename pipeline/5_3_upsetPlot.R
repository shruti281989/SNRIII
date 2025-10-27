suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyverse)
  library(tidyr)
  library(pheatmap)
  library(viridis)
  library(gplots)
  library(UpSetR)
  library(ComplexHeatmap)
  library(readxl)
})
set.seed(42)

# upset plots for deg set (Figure in the manuscript)
deg <- read.table("~/shruti/SNRIII/data/SeuratOut/output/afterDbltRemoval/markerWilcox_lfc1_fdr0.01_pct0.1.txt",
                             header = T)
deg <- deg %>% select("cluster", "avg_log2FC", "gene")

up = list(
  g1 = deg %>% filter(cluster=="Cluster_1" & avg_log2FC > 0),
  g3 = deg %>% filter(cluster=="Cluster_3" & avg_log2FC > 0),
  g4 = deg %>% filter(cluster=="Cluster_4" & avg_log2FC> 0),
  g5 = deg %>% filter(cluster=="Cluster_5" & avg_log2FC> 0),
  g6 = deg %>% filter(cluster=="Cluster_6" & avg_log2FC> 0),
  g7 = deg %>% filter(cluster=="Cluster_7" & avg_log2FC> 0),
  g8 = deg %>% filter(cluster=="Cluster_8" & avg_log2FC> 0),
  g10 = deg %>% filter(cluster=="Cluster_10" & avg_log2FC>0),
  g11 = deg %>% filter(cluster=="Cluster_11" & avg_log2FC>0),
  g12 = deg %>% filter(cluster=="Cluster_12" & avg_log2FC>0),
  g13 = deg %>% filter(cluster=="Cluster_13" & avg_log2FC>0),
  g14 = deg %>% filter(cluster=="Cluster_14" & avg_log2FC>0),
  g15 = deg %>% filter(cluster=="Cluster_15" & avg_log2FC>0),
  g16 = deg %>% filter(cluster=="Cluster_16" & avg_log2FC>0),
  g17 = deg %>% filter(cluster=="Cluster_17" & avg_log2FC>0),
  g18 = deg %>% filter(cluster=="Cluster_18" & avg_log2FC>0),
  g19 = deg %>% filter(cluster=="Cluster_19" & avg_log2FC>0),
  g20 = deg %>% filter(cluster=="Cluster_20" & avg_log2FC>0))

dn = list(
  g1 = deg %>% filter(cluster=="Cluster_1" & avg_log2FC< 0),
  g3 = deg %>% filter(cluster=="Cluster_3" & avg_log2FC < 0),
  g4 = deg %>% filter(cluster=="Cluster_4" & avg_log2FC< 0),
  g5 = deg %>% filter(cluster=="Cluster_5" & avg_log2FC< 0),
  g6 = deg %>% filter(cluster=="Cluster_6" & avg_log2FC< 0),
  g7 = deg %>% filter(cluster=="Cluster_7" & avg_log2FC< 0),
  g8 = deg %>% filter(cluster=="Cluster_8" & avg_log2FC< 0),
  g10 = deg %>% filter(cluster=="Cluster_10" & avg_log2FC< 0),
  g11 = deg %>% filter(cluster=="Cluster_11" & avg_log2FC< 0),
  g12 = deg %>% filter(cluster=="Cluster_12" & avg_log2FC< 0),
  g13 = deg %>% filter(cluster=="Cluster_13" & avg_log2FC< 0),
  g14 = deg %>% filter(cluster=="Cluster_14" & avg_log2FC< 0),
  g15 = deg %>% filter(cluster=="Cluster_15" & avg_log2FC< 0),
  g16 = deg %>% filter(cluster=="Cluster_16" & avg_log2FC< 0),
  g17 = deg %>% filter(cluster=="Cluster_17" & avg_log2FC< 0),
  g18 = deg %>% filter(cluster=="Cluster_18" & avg_log2FC< 0),
  g19 = deg %>% filter(cluster=="Cluster_19" & avg_log2FC< 0),
  g20 = deg %>% filter(cluster=="Cluster_20" & avg_log2FC< 0))

create_logFC_matrix <- function(gene_list) {
  all_genes <- unique(unlist(lapply(gene_list, function(df) df$gene)))
  all_clusters <- names(gene_list)
  logFC_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(all_clusters))
  rownames(logFC_matrix) <- all_genes
  colnames(logFC_matrix) <- all_clusters
  for (cluster in all_clusters) {
    cluster_data <- gene_list[[cluster]]
    logFC_matrix[cluster_data$gene, cluster] <- cluster_data$avg_log2FC
  }
  logFC_matrix[is.na(logFC_matrix)] <- 0
  return(logFC_matrix)
}

logFC_matrix_up <- create_logFC_matrix(up)
logFC_matrix_down <- create_logFC_matrix(dn)

hc_rows_up <- hclust(dist(logFC_matrix_up, method = "euclidean"), method = "complete")
hc_rows_down <- hclust(dist(logFC_matrix_down, method = "euclidean"), method = "complete")

align_matrix_with_dendrogram <- function(matrix, hc_rows) {
  dendro_order <- rownames(matrix)[hc_rows$order]
  aligned_matrix <- matrix[dendro_order, , drop = FALSE]
  return(aligned_matrix)
}

logFC_matrix_up_aligned <- align_matrix_with_dendrogram(logFC_matrix_up, hc_rows_up)
logFC_matrix_down_aligned <- align_matrix_with_dendrogram(logFC_matrix_down, hc_rows_down)

min_value <- min(logFC_matrix_up, logFC_matrix_down, na.rm = TRUE)
max_value <- max(logFC_matrix_up, logFC_matrix_down, na.rm = TRUE)
breaks <- seq(min_value, max_value, length.out = 51)

row_order <- rev(rownames(logFC_matrix_down_aligned)[hc_rows_down$order])
write.table(row_order,"dn.txt",quote = F,col.names = T,row.names = F)

row_order <- rev(rownames(logFC_matrix_up_aligned)[hc_rows_up$order])
write.table(row_order,"up.txt",quote = F,col.names = T,row.names = F)

dim(logFC_matrix_up_aligned_filtered)
dim(logFC_matrix_down_aligned)

filter_genes <- function(gene_list, genes_of_interest) {
  return(lapply(gene_list, function(df) df %>% filter(gene %in% genes_of_interest)))
}

filtered_up <- filter_genes(up, genes_of_interest)
filtered_dn <- filter_genes(dn, genes_of_interest)

logFC_matrix_up_filtered <- create_logFC_matrix(filtered_up)
logFC_matrix_down_filtered <- create_logFC_matrix(filtered_dn)
hc_rows_up_filtered <- hclust(dist(logFC_matrix_up_filtered, method = "euclidean"), method = "complete")
hc_rows_down_filtered <- hclust(dist(logFC_matrix_down_filtered, method = "euclidean"), method = "complete")

logFC_matrix_up_aligned_filtered <- align_matrix_with_dendrogram(logFC_matrix_up_filtered, hc_rows_up_filtered)
logFC_matrix_down_aligned_filtered <- align_matrix_with_dendrogram(logFC_matrix_down_filtered, hc_rows_down_filtered)

min_value_filtered <- min(logFC_matrix_up_filtered, logFC_matrix_down_filtered, na.rm = TRUE)
max_value_filtered <- max(logFC_matrix_up_filtered, logFC_matrix_down_filtered, na.rm = TRUE)
breaks_filtered <- seq(min_value_filtered, max_value_filtered, length.out = 51)

svg("heatmap_upregulated_fil.svg", width = 12, height = 10)
heatmap.2(logFC_matrix_up_aligned_filtered, col = viridis(50), breaks = breaks_filtered, 
          trace = "none", dendrogram = "both", Rowv = as.dendrogram(hc_rows_up_filtered),
          Colv = FALSE, cellnote = round(logFC_matrix_up_aligned_filtered, 3), notecol = "black",
          notecex = 0.1, scale = "none", cexCol = 1, cexRow = 0.1, margins = c(10, 10))
dev.off()

svg("heatmap_downregulated.svg", width = 12, height = 20)
heatmap.2(logFC_matrix_down_aligned, col = viridis(50), breaks = breaks,
          trace = "none", dendrogram = "both", Rowv = as.dendrogram(hc_rows_down),
          Colv = FALSE, cellnote = round(logFC_matrix_down_aligned, 3), notecol = "black",
          notecex = 0.01, scale = "none", cexCol = 1, cexRow = 0.1, margins = c(10, 10))
dev.off()

# upset plots
dnlt <- lapply(dn[c(1:12,14,15)], function(df) df$gene)
mdown = make_comb_mat(dnlt)
UpSet(mdown)

uplt <- lapply(up[c(1:18)], function(df) df$gene)
mup = make_comb_mat(uplt)
UpSet(mup)
