setwd("/mnt/picea/home/schoudhary/shruti/SNRIII")

library(svglite)
library(monocle3)
library(dplyr)
library(tibble)
library(pheatmap)
library(ggplot2)
library(ggplot2)
library(dplyr)
library(here)

# Module
deg_ids <- readLines("data/SeuratOut/fibTrajec/pseudotimeDeg0.01q_1.txt")
cds <- readRDS("data/SeuratOut/fibTrajec/cdsFib14/cds_object.rds")

# deg_ids <- readLines("data/SeuratOut/fusiformVes/pseudotimeDeg0.01q.txt")
# cds <- readRDS("data/SeuratOut/fusiformVes/cdsVes/cds_object.rds")

cds_subset <- cds[deg_ids, ]
cds_subset <- preprocess_cds(cds_subset, method = "PCA")
cds_subset <- reduce_dimension(cds_subset, reduction_method = "UMAP")
cds_subset <- cluster_cells(cds_subset, reduction_method = "UMAP")

cell_group_df <- tibble::tibble(cell = rownames(colData(cds_subset)),
                                cell_group = cds_subset@colData$seurat_clusters[colnames(cds_subset)])
gene_module_df <- find_gene_modules(cds_subset, resolution = 0.001, random_seed = 42)

trajectory_order <- c(14, 1, 10, 4, 15)
# trajectory_order <- c(14, 6, 17)
cluster_labels <- paste0("Cluster_", trajectory_order)

expr_matrix <- exprs(cds_subset)
cluster_info <- colData(cds_subset)$seurat_clusters
names(cluster_info) <- colnames(cds_subset)

avg_expr_by_cluster <- sapply(trajectory_order, function(cluster) {
  cells <- names(cluster_info)[cluster_info == cluster]
  rowMeans(expr_matrix[, cells, drop = FALSE])
})
colnames(avg_expr_by_cluster) <- cluster_labels

avg_expr_by_cluster <- as.data.frame(avg_expr_by_cluster)
avg_expr_by_cluster$gene <- rownames(avg_expr_by_cluster)

gene_module_df$gene <- gene_module_df$id
merged_df <- merge(gene_module_df, avg_expr_by_cluster, by = "gene")

dot_data_long <- tidyr::pivot_longer(merged_df, cols = starts_with("Cluster_"),
                                     names_to = "cluster", 
                                     values_to = "expression")

dot_data_long$cluster <- factor(dot_data_long$cluster, 
                                levels = paste0("Cluster_", c(14, 1, 10, 4, 15)))

module_expr <- dot_data_long %>% group_by(module, cluster) %>%
  summarise(avg_expr = mean(expression), .groups = "drop")

# line plot
ggplot(module_expr, aes(x = cluster, y = avg_expr, group = module, 
                        color = as.factor(module))) + geom_line(size = 1.2) +
  geom_point(size = 2) + theme_minimal() +
  labs(x = "Cluster", y = "Average Expression", color = "Module") +
  scale_x_discrete(limits = paste0("Cluster_", c(14, 1, 10, 4, 15)))

dir.create("data/SeuratOut/fibTrajec/clusterwise", showWarnings = FALSE)

# ptime gene plot according to module in each cluster
modules <- unique(merged_df$module)
for (mod in modules) {
  mod_genes <- merged_df$gene[merged_df$module == mod]
  if (length(mod_genes) < 2) next
  
  cluster_matrix <- avg_expr_by_cluster[avg_expr_by_cluster$gene %in% mod_genes,
                                        cluster_labels]
  rownames(cluster_matrix) <- avg_expr_by_cluster$gene[avg_expr_by_cluster$gene %in% mod_genes]
  cluster_matrix <- as.matrix(cluster_matrix)
  scaled_cluster_matrix <- t(scale(t(cluster_matrix)))
  gene_order <- rownames(scaled_cluster_matrix)[order(rowMeans(scaled_cluster_matrix),
                                                      decreasing = TRUE)]
  
  svg_filename_cluster <- paste0("data/SeuratOut/fibTrajec/clusterwise/module_", mod, "_clusterwise.svg")
  svg(svg_filename_cluster, width = 8, height = 10)
  pheatmap(scaled_cluster_matrix[gene_order, ], fontsize = 1, scale = "row", 
           clustering_method = "ward.D2", cluster_rows =T,
           cluster_cols = F, show_rownames = T, color = viridis(50), 
           main = paste("Module", mod, "- Cluster-wise"))
  dev.off()
  # writeLines(gene_order, paste0("data/SeuratOut/fibTrajec/orders/module_",
  #                               mod, "_gene_order.txt"))
}
