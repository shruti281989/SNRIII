setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(tidyverse)
  library(pheatmap)
  library(viridis)
  library(gplots)
})

set.seed(42)

# load seurat onject
integ <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"
Idents(integ) <- "integrated_snn_res.0.6"
integ$seurat_clusters <- integ@active.ident

timeBinDeg <- read.table(here("doc/timeBinDeg.txt"), header = TRUE)
# time_bins <- c("2h","4h","8h","12h")
# levels <- c("Upregulated", "Downregulated") 
up2h <- timeBinDeg %>% filter(timeBin == "2h" & level =="Upregulated") 
avg <- AverageExpression(
  integ,
  assays = "RNA",
  # group.by = c("sample", "seurat_clusters"),
  features = up2h$gene,
  slot = "data"
)

mat <- avg$RNA
log_mat <- log2(mat+1)
# celltypes <- unique(integ$seurat_clusters)
# 
# mat_avg <- sapply(celltypes, function(ct) {
#   cols <- grep(paste0("_", ct, "$"), colnames(mat))
#   rowMeans(mat[, cols, drop = FALSE])
# })

# mat_avg <- t(mat_avg)

pheatmap(log_mat, scale = "none", cluster_rows = T, cluster_cols = F,
         display_numbers = T)

