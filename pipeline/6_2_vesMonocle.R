setwd(data/SeuratOut)
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

integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
Idents(integ) <- "integrated_snn_res.0.6"
integ$seurat_clusters <- integ@active.ident
DefaultAssay(integ) <- "RNA"

# trajectories for vessels
ves <- WhichCells(integ, ident= c("6","14","17"))
integ2 = subset(integ, cells = ves)

cds <- as.cell_data_set(integ2)
cds <- cluster_cells(cds, resolution = 1e-3) 

# p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
# p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
# wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
# plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE,
#            label_leaves=FALSE, label_branch_points=FALSE)

cds <- order_cells(cds, reduction_method = "UMAP") #choose 14 as root

plot_cells(cds, color_cells_by = "pseudotime",group_cells_by = "cluster",
           label_cell_groups = FALSE,label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE,
           label_roots = FALSE,trajectory_graph_color = "grey60")

integrated.sub <- as.Seurat(cds, assay = NULL)
# FeaturePlot(integrated.sub, "monocle3_pseudotime")

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, 
                        reorder(seurat_clusters, monocle3_pseudotime), #reorder by ident
                        fill = seurat_clusters)) + geom_boxplot()

cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph",
                                     cores = 8)
saveRDS(cds_graph_test_results, "vesTrajec/cdsVes_graph.rds")

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.01))
write.table(deg_ids,"vesTrajec/pseudotimeDeg0.01q.txt",
            sep="\t", row.names = F,col.names = F, quote = F)

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, 
           cell_size = 1, x = 1, y = 2, label_branch_points=F, label_leaves=F, 
           label_groups_by_cluster=F, label_roots = F) + coord_fixed() + 
  theme_void() + theme(plot.title = element_text(hjust = 0.5, size=15), 
                       legend.position="bottom", legend.title=element_text (size=10), 
                       legend.text=element_text(size=10)) + 
  guides(colour = guide_colourbar(barwidth = 10, barheight = 0.5))

cds@colData$sample <- integ@meta.data[colnames(cds), "sample"]
meta <- as.data.frame(colData(cds))

meta$sample_cluster <- paste0("C", meta$seurat_clusters, "_", meta$sample)

ggplot(meta, aes(x = monocle3_pseudotime, color = sample, fill = sample)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~seurat_clusters, ncol = 1, scales = "free_y") +
  scale_fill_viridis_d(option = "C") + scale_color_viridis_d(option = "C") +
  theme_classic() + labs(x = "Pseudotime", y = "Density")

# deg_ids plot into modules
cds <- as.cell_data_set(integ2)
cds <- cluster_cells(cds, resolution = 1e-3)
cds <- learn_graph(cds, use_partition = TRUE)
cds <- order_cells(cds, reduction_method = "UMAP")

cds$monocle3_pseudotime <- pseudotime(cds)
integ1$monocle3_pseudotime <- cds$monocle3_pseudotime

cds@colData$seurat_clusters <- integ1$seurat_clusters[colnames(cds)]
cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph", 
                                     cores = 8)

deg_ids <- rownames(subset(cds_graph_test_results, q_value < 0.01))
cds_subset <- cds[deg_ids, ]
cds_subset <- preprocess_cds(cds_subset, method = "PCA")
cds_subset <- reduce_dimension(cds_subset, reduction_method = "UMAP")
cds_subset <- cluster_cells(cds_subset, reduction_method = "UMAP")

cell_group_df <- tibble::tibble(cell = rownames(colData(cds_subset)),
                                cell_group = cds_subset@colData$seurat_clusters[colnames(cds_subset)])
gene_module_df <- find_gene_modules(cds_subset, resolution = 0.001, random_seed = 42)
write.table(gene_module_df, "data/SeuratOut/vesTrajec/gene_modules.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)