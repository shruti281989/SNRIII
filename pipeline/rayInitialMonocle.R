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
})

# integrated object:
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
Idents(integ) <- "integrated_snn_res.0.6"
integ$seurat_clusters <- integ@active.ident
DefaultAssay(integ) <- "RNA"

# try these and check
rayInitial <- WhichCells(integ, ident= c("16","18","8","5","7","3","9","19","20")) 
rayInitial <- WhichCells(integ, ident= c("16","18","8","5","7","3","9","11","12","13","19","20")) 
integ1 = subset(integ, cells = rayInitial)

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

cds <- order_cells(cds, reduction_method = "UMAP") #choose 19,20,18 for all

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
  labs(title = "Pseudotime distribution across fusiform clusters", 
       x = "Pseudotime", y = "Cluster") + theme(legend.position = "none")

cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph",
                                     cores = 8)
saveRDS(cds_graph_test_results, "data/SeuratOut/rayInitial1/cds_graph.rds")
rowData(cds)$gene_short_name <- row.names(rowData(cds))
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.01))
write.table(deg_ids,"data/SeuratOut/rayInitial1/pseudotimeDeg0.01q.txt", 
            sep="\t", row.names = F,col.names = F, quote = F)

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