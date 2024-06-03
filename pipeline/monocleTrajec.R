# https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
# https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# install.packages("devtools")
# # devtools::install_github('cole-trapnell-lab/monocle3')
# remotes::install_github('satijalab/seurat-wrappers')
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

# With mtcp genome and with cell cycle
# integMtCp <- readRDS("~/shruti/SNRIII/data/SeuratOut/integMtCp.rds")
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")

# Part A. trial with only control
split_seurat <- SplitObject(integ, split.by = "sample")
# kcl0 <- split_seurat[["kclmtcp"]]
kcl0 <- split_seurat[["ctrl"]]
DefaultAssay(kcl0) <- "SCT"
DimPlot(kcl0, reduction = "umap", pt.size = 0.01, label = TRUE)

# Take few clusters of the known identities
kcl1 <- subset (kcl0, subset= integrated_snn_res.0.6==c("16","5","14","18","13","15","17","4"))

DimPlot(kcl0, reduction = "umap", pt.size = 0.01, label = TRUE)
cds <- as.cell_data_set(kcl1)
# 1. Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# Assign cluster information
list.cluster <- kcl1@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <-
  kcl1@reductions$umap@cell.embeddings

plot_cells(cds, color_cells_by = "cluster", 
           label_groups_by_cluster = F, 
           group_label_size = 5,
           show_trajectory_graph = FALSE) +NoLegend()

# learn trajectory
cds <- learn_graph(cds, use_partition = F)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F,
           group_label_size = 5)

# Select 5 as the root node for kcl
# Select 1,10,14,0 as the root node for kcl1
cds <- order_cells(cds, reduction_method = "UMAP")
# cds <- order_cells(cds, reduction_method = "UMAP", 
      # root_cells = colnames(cds[, clusters(cds) == 20]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F)

# Cells ordered by Monocle3 Pseudotime
head(pseudotime(cds), 10)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, ident, 
                        fill = ident)) + geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident,
                                                     monocle3_pseudotime), 
                        fill = ident)) + geom_boxplot()
# save plots and restart R

# Find genes that change as a function of pseudotime
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
# FeaturePlot(integ, features = c("Rgs20", "Bend6", "Mcm3", "B3gat2"))

# PART B. trial with integrated object where:
# kcl (removed: empty drops first, high mt cells; and cluster 1,10,6 based on correlation)
# for kno (cluster 11,6,9 removed); 
# overall: (nUMI >= 1000) & (nGene >= 500)& (percent.mt < 5) & (percent.cp < 5) & (nGene < 9000))
# integ <- readRDS("~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/mtCp/filt0.8_seurat1integ.rds")

# for the integrated object with cell cycle regressed out and no mtcp genes
Idents(integ) <- "integrated_snn_res.0.6"
DefaultAssay(integ) <- "RNA"
DimPlot(integ, reduction = "umap", pt.size = 0.01, label = TRUE, split.by = "sample")

integ$seurat_clusters <- integ@active.ident
# Take few clusters of the known identities
integ1 <- WhichCells(integ, ident= c("1","10","6","14","17", "4","15"))
integ1 = subset(integ, cells = integ1)

# cds <- order_cells(cds, reduction_method = "UMAP", 
# root_cells = colnames(cds[, clusters(cds) == 20]))

# PART C. trial with integrated object where:
# cell cycle regressed out and no mtcp genes
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
Idents(integ) <- integ$integrated_snn_res.0.6
DefaultAssay(integ) <- "RNA"
DimPlot(integ, reduction = "umap", pt.size = 0.01, label = TRUE) # split.by = "sample"
integ$seurat_clusters <- integ@active.ident

# Take few clusters of the known identities

# fiber and vessel
integ1 <- WhichCells(integ, ident= c("1","10","6","14","17", "4","15"))

# fib
# integ1 <- WhichCells(integ, ident= c("1","10","4","20","12","15","8"))

# fib1
# integ1 <- WhichCells(integ, ident= c("1","10","4","12","15"))

# fib2
# integ1 <- WhichCells(integ, ident= c("1","10","4","15"))

# fib3 (seems ok to use)
# integ1 <- WhichCells(integ, ident= c("1","10","4","15","12","14"))

# fib4 (maybe good)
# integ1 <- WhichCells(integ, ident= c("1","10","4","15","20","14"))

# fibVes2 (maybe good)
# integ1 <- WhichCells(integ, ident= c("1","4","15","17","14","6","20"))
# remove 20, add 12 and 10

# fibVes3 (seems ok)
# integ1 <- WhichCells(integ, ident= c("1","4","15","17","14","10","12","6","20"))

# all- 16 as root
integ1 <- WhichCells(integ, ident= c("1","4","15","17","14","10","12","6","20","18","19","5","7","8","16"))

# fib only
# integ1 <- WhichCells(integ, ident= c("1","4","15","14"))

# fibVes4- 16 as root
# integ1 <- WhichCells(integ, ident= c("1","4","15","17","14","10","12","6","20","8","16"))

# integ1 <- WhichCells(integ, ident= c("17","14","6","16"))

# integ1 = subset(integ, cells = integ1)
cds <- as.cell_data_set(integ1)
# integ1$ident <- integ1@active.ident

# use partition
cds <- cluster_cells(cds, resolution=1e-3) #for fibVes, 2,; fib, 1,2,4,3; fibVes3,4, all

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE)

# cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 3]))
# for ves1

# cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 1]))
# for fib and ves; fib,1,2,

# cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 5])) 
#for fib4

# cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 10]))
# for fib3

# cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 6]))
# for fibVes2

# cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 9]))
# for fibVes2

# all (maybe)
# cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 6]))

# for fibVes4
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 15]))

plot_cells(cds,
           color_cells_by = "pseudotime",group_cells_by = "cluster",
           label_cell_groups = FALSE,label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE,
           label_roots = FALSE,trajectory_graph_color = "grey60")

integrated.sub <- as.Seurat(cds, assay = NULL)
FeaturePlot(integrated.sub, "monocle3_pseudotime")

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, 
                        reorder(seurat_clusters, monocle3_pseudotime), 
                        fill = seurat_clusters)) + geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, 
                        reorder(ident, monocle3_pseudotime), 
                        fill = ident)) + geom_boxplot()

# without partition
# Ves1 (ok if you try with part A) and 14 as the root
# vesOnly <- WhichCells(integ, ident= c("17","14","6"))
# vesOnly = subset(integ, cells = vesOnly)
# DimPlot(vesOnly, reduction = "umap", pt.size = 0.01, label = TRUE)

# fibVes
# fibVes <- WhichCells(integ, ident= c("1","10","6","14","17", "4","15", "12", "8"))
# fibVes = subset(integ, cells = fibVes)
# cds <- as.cell_data_set(fibVes)

# fib
# integ1 <- WhichCells(integ, ident= c("1","10","4","20","12","15","8", "14"))

# fib1
# integ1 <- WhichCells(integ, ident= c("1","10","4","12","15", "14"))

# fib2
# integ1 <- WhichCells(integ, ident= c("1","10","4","15", "14"))

# fib3 (maybe good)
# integ1 <- WhichCells(integ, ident= c("1","10","4","15","14"))

# fibVes2 (maybe good)
# integ1 <- WhichCells(integ, ident= c("1","10","4","15","14", "17","6"))

# all- 16 as root
# integ1 <- WhichCells(integ, ident= c("1","4","15","17","14","10","12","6","20","18","19","0","5","7","8","16"))

# fibVes3 - 16 as root
integ1 <- WhichCells(integ, ident= c("1","4","15","17","14","10","12","6","20","18","19","5","7","8","16"))
integ1 = subset(integ, cells = integ1)
DimPlot(integ1, reduction = "umap", pt.size = 0.01, label = TRUE)
cds <- as.cell_data_set(integ1)

# 1. Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- integ1@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <-
  integ1@reductions$umap@cell.embeddings
# plot_cells(cds, color_cells_by = "cluster", 
#            label_groups_by_cluster = F, group_label_size = 5,
#            show_trajectory_graph = FALSE) + NoLegend()

# learn trajectory
cds <- learn_graph(cds, use_partition = F)
# plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
#            label_branch_points = F, label_roots = F, label_leaves = F,
#            group_label_size = 5)

cds <- order_cells(cds, reduction_method = "UMAP")
# plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F,
#            label_branch_points = F, label_roots = F, label_leaves = F)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident,monocle3_pseudotime), 
                        fill = ident)) + geom_boxplot()
