set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(ggplot2)
  library(here)
  library(hyperSpec)
  library(RColorBrewer)
  library(tidyverse)
  library(VennDiagram)
  library(readxl)
  library(pheatmap)
})

pal =rev(brewer.pal(9, "Spectral"))
hpal <- colorRampPalette(c("#3288BD","#FFFFBF","#D53E4F"))(100)
mar <- par("mar")

# Download the GSE180121_Cell_type_TPM_Ptr.txt from Tung et al., 2023
tungtpm <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/Tung/GSE180121_Cell_type_TPM_Ptr.txt", 
                      sep = "\t", header = TRUE)
tungtpm$Potri <- sub(".v4.1", "", tungtpm$Gene.ID)

potrapotri <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/potra_potri_BEST_DIAMOND_out.tsv.gz")
colnames(potrapotri) <- c("Potra", "Potri")

tungtpmPotra <- inner_join(tungtpm, potrapotri, by = "Potri")
tungtpmPotra$Potri <- NULL
tungtpmPotra$Gene.ID <- NULL

nitgene <- read.table("~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/correl/nitGene.txt")
colnames(nitgene) <- c("Potra", "Gene", "Family")

merged_data <- merge(tungtpmPotra, nitgene, by.x = "Potra", by.y = "Potra")
merged_data <- merged_data[order(merged_data$Family), ]  

expression_data <- merged_data[, c("Ray.mean.TPM", "Fiber.mean.TPM", "Vessel.bio2.TPM")]
rownames(expression_data) <- paste(merged_data$Potra, merged_data$Gene)
log_expression_data <- log2(expression_data + 1)

# first option :Order genes within each family by Ray.mean.TPM (or see second)
ordered_rows <- c()
for (fam in unique(merged_data$Family)) {
  family_data <- merged_data[merged_data$Family == fam, ]
  family_data <- family_data[order(-family_data$Ray.mean.TPM), ]
  ordered_rows <- c(ordered_rows, paste(family_data$Potra, family_data$Gene))
}

log_expression_data <- log_expression_data[ordered_rows, ]

annotation_row <- data.frame(Family = merged_data$Family[match(ordered_rows, rownames(expression_data))])
rownames(annotation_row) <- ordered_rows

svg(file.path("~/shruti/SNRIII/by_Ray_heatmap.svg"), width = 10, height = 12)
pheatmap(mat = log_expression_data, cluster_rows = FALSE,  cluster_cols = F,
  scale = "none", color = hpal, fontsize_row = 2, fontsize_col = 8, 
  fontsize = 10, annotation_row = annotation_row)
dev.off()

write.table(ordered_rows, file ="order.txt", col.names = F,row.names = F, 
            quote = F, sep="\t")

# second option: Perform clustering within each family
ordered_rows <- c()
for (fam in unique(merged_data$Family)) {
  family_data <- log_expression_data[merged_data$Family == fam, , drop = FALSE]
  
  if (nrow(family_data) > 1) {
    family_dist <- dist(family_data)
    family_hclust <- hclust(family_dist, method = "ward.D2")
    ordered_rows <- c(ordered_rows, rownames(family_data)[family_hclust$order])
  } else {
    ordered_rows <- c(ordered_rows, rownames(family_data))
  }
}

log_expression_data <- log_expression_data[ordered_rows, ]
annotation_row <- data.frame(Family = merged_data$Family[match(ordered_rows, rownames(expression_data))])
rownames(annotation_row) <- ordered_rows

svg(file.path("~/shruti/SNRIII/heatmap.svg"), width = 10, height = 12)
pheatmap(mat = log_expression_data, cluster_rows = F,  cluster_cols = F,
  scale = "none", color = hpal, fontsize_row = 2, fontsize_col = 8, 
  fontsize = 10, annotation_row = annotation_row)
dev.off()
