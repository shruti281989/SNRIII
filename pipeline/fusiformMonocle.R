set.seed(42)
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(monocle3)
  library(dplyr)
  library(SeuratWrappers)
  library(patchwork)
  library(ggplot2)
  library(ggridges)
  library(viridis)
})

# integrated object:
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
Idents(integ) <- "integrated_snn_res.0.6"
integ$seurat_clusters <- integ@active.ident
DefaultAssay(integ) <- "RNA"

# try these and check
#use this for paper if want to provide a common lineage
fusiform <- WhichCells(integ, ident= c("1","4","6","10","14","15","17")) 

fusiform <- WhichCells(integ, ident= c("1","4","10","14","15")) 
fusiform <- WhichCells(integ, ident= c("6","14","17")) 

integ1 = subset(integ, cells = fusiform)

cds <- as.cell_data_set(integ1)
cds <- cluster_cells(cds, resolution = 1e-3) 

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE)

cds <- order_cells(cds, reduction_method = "UMAP") #choose 14

plot_cells(cds, color_cells_by = "pseudotime",group_cells_by = "cluster",
           label_cell_groups = FALSE,label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE,
           label_roots = FALSE,trajectory_graph_color = "grey60")

integrated.sub <- as.Seurat(cds, assay = NULL)
FeaturePlot(integrated.sub, "monocle3_pseudotime")

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, 
                        reorder(seurat_clusters, monocle3_pseudotime), #reorder by ident
                        fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo, aes(x = monocle3_pseudotime, y = seurat_clusters, 
                        fill = seurat_clusters)) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_viridis_d(option = "C") + theme_classic() +
  labs(title = "Pseudotime distribution across clusters", 
       x = "Pseudotime", y = "Cluster") + theme(legend.position = "none")

cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph",
                                     cores = 8)
# saveRDS(cds_graph_test_results, "data/SeuratOut/fusiformTrajec/cdsfib_graph.rds")
saveRDS(cds_graph_test_results, "data/SeuratOut/fusiform-16Trajec/cdsfib_graph.rds")

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.01))

# write.table(deg_ids,"data/SeuratOut/fusiformTrajec/pseudotimeDeg0.01q.txt", 
#             sep="\t", row.names = F,col.names = F, quote = F)
write.table(deg_ids,"data/SeuratOut/fusiform-16Trajec/pseudotimeDeg0.01q.txt",
            sep="\t", row.names = F,col.names = F, quote = F)
# write.table(deg_ids,"data/SeuratOut/fusiformVes/pseudotimeDeg0.01q.txt", 
#             sep="\t", row.names = F,col.names = F, quote = F)

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, 
           cell_size = 1, x = 1, y = 2, label_branch_points=F, label_leaves=F, 
           label_groups_by_cluster=F, label_roots = F) + coord_fixed() + 
  theme_void() + labs(title = "Pseudotime") + 
  theme(plot.title = element_text(hjust = 0.5, size=15), 
        legend.position="bottom", legend.title=element_text (size=10), 
        legend.text=element_text(size=10)) + 
  guides(colour = guide_colourbar(barwidth = 10, barheight = 0.5, 
                                  title = "Pseudotime"))

cds@colData$sample <- integ@meta.data[colnames(cds), "sample"]
meta <- as.data.frame(colData(cds))
ggplot(meta, aes(x = monocle3_pseudotime, y = seurat_clusters, fill = seurat_clusters)) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_viridis_d(option = "C") +
  theme_classic() +
  facet_wrap(~sample, scales = "free_y") +
  labs(title = "Pseudotime distribution per cluster, split by sample", 
       x = "Pseudotime", y = "Cluster") +
  theme(legend.position = "none")

ggplot(meta, aes(x = monocle3_pseudotime, fill = seurat_clusters)) +
  geom_density(alpha = 0.6) +
  facet_grid(seurat_clusters ~ sample, scales = "free_y") +
  scale_fill_viridis_d(option = "C") +
  theme_classic() +
  labs(title = "Pseudotime distributions split by cluster and sample",
       x = "Pseudotime", y = "Density")

meta$sample_cluster <- paste0("C", meta$seurat_clusters, "_", meta$sample)

ggplot(meta, aes(x = monocle3_pseudotime, y = sample_cluster, fill = sample)) +
  geom_density_ridges(scale = 1.1, alpha = 0.8) +
  scale_fill_viridis_d(option = "C") +
  theme_classic() +
  labs(title = "Pseudotime Ridge Plot by Sample and Cluster",
       x = "Pseudotime", y = "Cluster_Sample") +
  theme(axis.text.y = element_text(size = 6))

ggplot(meta, aes(x = monocle3_pseudotime, color = sample, fill = sample)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~seurat_clusters, ncol = 1, scales = "free_y") +  # optional facet by cluster
  scale_fill_viridis_d(option = "C") +
  scale_color_viridis_d(option = "C") +
  theme_classic() +
  labs(title = "Pseudotime Density by Cluster and Sample",
       x = "Pseudotime", y = "Density")

ggplot(meta, aes(x = monocle3_pseudotime, y = seurat_clusters, fill = sample)) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01, alpha = 0.8, color = "black", size = 0.3) +
  scale_fill_viridis_d(option = "C") +
  theme_classic() +
  labs(title = "Pseudotime Ridge Plot by Cluster, Colored by Sample",
       x = "Pseudotime", y = "Cluster") +
  theme(legend.position = "right")

ggplot(meta, aes(x = monocle3_pseudotime, color = sample, fill = sample)) +
  geom_density(alpha = 0.4, size = 1) +
  scale_fill_viridis_d(option = "D") +
  scale_color_viridis_d(option = "D") +
  theme_classic() +
  labs(title = "Pseudotime Density Plot by Sample",
       x = "Pseudotime", y = "Density")

ggplot(meta, aes(x = seurat_clusters, y = monocle3_pseudotime, fill = sample)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) +
  theme_classic() +
  labs(title = "Pseudotime distribution by cluster and sample",
       x = "Cluster", y = "Pseudotime") +
  scale_fill_viridis_d(option = "C")

# deg_ids plot into modules
integ1 <- subset(integ, cells = fusiform)

cds <- as.cell_data_set(integ1)
cds <- cluster_cells(cds, resolution = 1e-3)
cds <- learn_graph(cds, use_partition = TRUE)
cds <- order_cells(cds, reduction_method = "UMAP")

cds$monocle3_pseudotime <- pseudotime(cds)
integ1$monocle3_pseudotime <- cds$monocle3_pseudotime
save_monocle_objects(cds, "data/SeuratOut/fusiform-16Trajec/cds/")

cds@colData$seurat_clusters <- integ1$seurat_clusters[colnames(cds)]
cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph", 
                                     cores = 8)

deg_ids <- rownames(subset(cds_graph_test_results, q_value < 0.01))
# deg_ids <- readLines("data/SeuratOut/fusiform-16Trajec/pseudotimeDeg0.01q_1.txt")
# cds <- readRDS("~/shruti/SNRIII/data/SeuratOut/fusiform-16Trajec/cds/cds_object.rds")

cds_subset <- cds[deg_ids, ]
cds_subset <- preprocess_cds(cds_subset, method = "PCA")
cds_subset <- reduce_dimension(cds_subset, reduction_method = "UMAP")
cds_subset <- cluster_cells(cds_subset, reduction_method = "UMAP")

cell_group_df <- tibble::tibble(cell = rownames(colData(cds_subset)),
                                cell_group = cds_subset@colData$seurat_clusters[colnames(cds_subset)])
gene_module_df <- find_gene_modules(cds_subset, resolution = 0.001, random_seed = 42)
write.table(gene_module_df, "data/SeuratOut/fusiform-16Trajec/gene_modules.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows = T, cluster_cols = T, scale = "column",
                   clustering_method = "ward.D2", fontsize = 6, color = viridis(50))

# heatmap of top genes
expr_matrix <- exprs(cds_subset)
cluster_info <- colData(cds_subset)$seurat_clusters

cluster_levels <- unique(cluster_info)
avg_expr_by_cluster <- sapply(cluster_levels, function(cluster) {
  cells_in_cluster <- names(cluster_info)[cluster_info == cluster]
  rowMeans(expr_matrix[, cells_in_cluster, drop = FALSE])
})

avg_expr_by_cluster <- as.data.frame(avg_expr_by_cluster)
colnames(avg_expr_by_cluster) <- paste0("Cluster_", cluster_levels)

top_genes <- gene_module_df %>% group_by(module) %>%
  group_map(~ {genes <- .x$id
  gene_avgs <- rowMeans(avg_expr_by_cluster[genes, , drop = FALSE])
  top_genes <- names(sort(gene_avgs, decreasing = TRUE))[1:min(10, length(gene_avgs))]
  return(top_genes)}) %>% unlist()

heatmap_matrix <- avg_expr_by_cluster[top_genes, ]
scaled_matrix <- t(scale(t(heatmap_matrix)))

desired_order <- c("Cluster_14", "Cluster_1", "Cluster_10", "Cluster_4",
                   "Cluster_15", "Cluster_6", "Cluster_17")
ordered_matrix <- scaled_matrix[, desired_order]

pheatmap(ordered_matrix, scale ="row", cluster_rows = T, cluster_cols = F,
         show_rownames = T, clustering_method = "ward.D2", fontsize = 2,
         color = viridis(100))
