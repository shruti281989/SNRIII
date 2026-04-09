setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/data/timeBinDegvsTungLCM/")
set.seed(42)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(here)

pal =rev(brewer.pal(9, "Spectral"))
hpal <- colorRampPalette(c("#3288BD","#FFFFBF","#D53E4F"))(50)
mar <- par("mar")

tungtpm <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/Tung/GSE180121_Cell_type_TPM_Ptr.txt", 
                      sep = "\t", header = TRUE)
tungtpm$Potri <- sub(".v4.1", "", tungtpm$Gene.ID)

potrapotri <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/potra_potri_BEST_DIAMOND_out.tsv.gz")
colnames(potrapotri) <- c("Potra", "Potri")

tungtpmPotra <- inner_join(tungtpm, potrapotri, by = "Potri")
tungtpmPotra$Potri <- NULL
tungtpmPotra$Gene.ID <- NULL

timeBinDeg <- read.table(here("doc/timeBinDeg.txt"), header = TRUE)
time_bins <- c("2h","4h","8h","12h")
levels <- c("Upregulated", "Downregulated") 

for (tb in time_bins) {
  
  for (lvl in levels) {
    subset_tb <- timeBinDeg %>% filter(timeBin == tb & level == lvl)
    if (nrow(subset_tb) == 0) next
    
    merged_data <- merge(tungtpmPotra, subset_tb, by.x = "Potra", by.y = "gene")
    if (nrow(merged_data) == 0) next
    
    merged_data <- merged_data[order(merged_data$timeBin), ]
    
    expression_data <- merged_data[, c("Ray.mean.TPM", "Fiber.mean.TPM", "Vessel.bio2.TPM")]
    rownames(expression_data) <- merged_data$Potra
    log_expression_data <- log2(expression_data + 1)
    
    # ray_expression <- log2(expression_data[ , "Ray.mean.TPM", drop = FALSE] + 1)
    # row_dist <- dist(ray_expression, method = "euclidean")
    # row_clust <- hclust(row_dist, method = "complete")
    
    svg_file <- paste0("heatmap_", tb, "_", lvl, ".svg")
    order_file <- paste0("order_", tb, "_", lvl, ".txt")
    
    svg(svg_file, width = 10, height = 12)
    p <- pheatmap(log_expression_data, cluster_rows = T, cluster_cols = F,
                  scale = "none", color = hpal, fontsize_row = 4,
                  fontsize_col = 8, fontsize = 10, display_numbers = T)

    # If you want genes to be clustered according to ray sample expression
    # p <- pheatmap(log_expression_data,cluster_rows = row_clust, cluster_cols = F, 
    #               scale = "none", display_numbers = T, color = hpal, 
    #               fontsize_row = 5, fontsize_col = 8, fontsize = 10)
    
    dev.off()
    
    ordered_genes <- rownames(log_expression_data)[p$tree_row$order]
    write.table(ordered_genes, file = order_file, col.names = F, 
                row.names = F, quote = F, sep = "\t")
    
    cat("Saved heatmap and gene order for", tb, lvl, "\n")
  }
}

