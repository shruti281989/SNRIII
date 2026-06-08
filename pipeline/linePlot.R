setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
})
integ <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"
deg <- read.table("../single_cell_analysis_poplar/doc/lnePlotGene.txt",
                  header = T)
deg <- deg %>% filter(type =="expansionFig6B")
# deg <- deg %>% filter(type =="expansionFig6C")

# merge clusters 1 + 10
integ$seurat_clusters <- as.character(integ$integrated_snn_res.0.6)
integ$seurat_clusters[integ$seurat_clusters %in% c("1","10")] <- "1_10"
integ$seurat_clusters[integ$seurat_clusters %in% c("19","20")] <- "19_20"
clusters_keep <- c("19", "20", "14","1_10","4","15", "17", "12")
subset_obj <- subset(integ, subset = seurat_clusters %in% clusters_keep)
genes_use <- intersect(unique(deg$gene), rownames(integ))
avg_expr <- AverageExpression(subset_obj, features = genes_use,
  group.by = c("seurat_clusters","sample"), assays = "RNA", slot = "data")$RNA
summary_df <- avg_expr %>% as.data.frame() %>% 
  tibble::rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "cluster_sample", values_to = "expression") %>%
  separate(cluster_sample, into = c("cluster","sample"), sep = "_") %>%
  pivot_wider(names_from = sample, values_from = expression) %>%
  mutate(logFC = log2(kno + 1e-6) - log2(ctrl + 1e-6)) %>%
  group_by(cluster) %>% summarise(mean_logFC = mean(logFC, na.rm = TRUE),
    se = sd(logFC) / sqrt(n()), .groups = "drop")
summary_df$cluster <- factor(summary_df$cluster,
                             levels = c("g19-20","g14","g1-10","g4","g15","g17", "g12"))

ggplot(summary_df, aes(cluster, mean_logFC, group = 1)) + 
  geom_line(linewidth = 1.3) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean_logFC - se, ymax = mean_logFC + se), width = 0.1) +
  theme_classic() + labs(x = "Cluster", y = "Mean log2 Fold Change")
