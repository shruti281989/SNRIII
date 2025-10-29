#' ---
#' title: "Single Cell pseudobulk analysis for quick nitrate repsonse in poplar"
#' author: "Shruti"
#' #' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#'
#' Info
#' Check if the paths to files are correct
#' 
setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/")
library(data.table)
library(here)
library(hyperSpec)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(reshape2)
library(tidyverse)
library(pheatmap)
library(viridis)
library(readxl)
#' 
#' 
#' Find markers to get DEGs
integ <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"
integ <- NormalizeData(integ, verbose = FALSE)
Idents(integ) <- "integrated_snn_res.0.6"
integ$seurat_clusters <- Idents(integ)
#' 
seurat_marker.list <- map(0:20,function(c){
  FindMarkers(subset(integ, subset = seurat_clusters == c), 
              ident.1 = "kno", group.by = "sample", min.pct = 0.1, only.pos=F,
              logfc.threshold = 1.0)})
names(seurat_marker.list) <- paste0("Cluster_",0:20)
filtered.marker <- map(names(seurat_marker.list),function(n){
  dplyr::filter(seurat_marker.list[[n]],p_val_adj < 0.01) %>%
    mutate(cluster = n, gene = rownames(.))
})
filtered_markers_df <- do.call(rbind, filtered.marker)
write.table(filtered_markers_df, "data/SeuratOut/degWilcoxpct0.1fdr0.01.txt", 
            row.names = F,col.names = T, quote = F, sep="\t")
#'
#' Figure in the manuscript
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")
#' 
#' DEGs in SNRIII
#' Modified data/SeuratOut/output/afterDbltRemoval/markerWilcox_lfc1_fdr0.01_pct0.1.txt
#' marker file in excel to have the desired coulmns as in the manuscript
degSnr <-read.delim(here("/data/SeuratOut/degWilcoxpct0.1fdr0.01.txt"), 
                    header = T,sep = "\t")
upKno <- degSnr %>% filter(Level.after.nitrate.tretament == "Upregulated") %>%
  pull(Gene) %>% unique()
#'
aspwood <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/publisheddatasets/AspWood_tpm.txt", 
                      header = TRUE)
aspwoodtpm <- dcast(aspwood, gene_id ~ sample_name)
orderaspwood <- c("T1-Phloem-01",
                  "T1-Phloem-02",
                  "T1-Phloem-03",
                  "T1-Phloem-04",
                  "T1-Phloem-05",
                  "T1-Cambium-06",
                  "T1-Cambium-07",
                  "T1-Cambium-08",
                  "T1-Cambium-09",
                  "T1-Cambium-10",
                  "T1-Cambium-11",
                  "T1-Cambium-12",
                  "T1-Expanding-xylem-13",
                  "T1-Expanding-xylem-14",
                  "T1-Expanding-xylem-15",
                  "T1-Expanding-xylem-16",
                  "T1-Expanding-xylem-17",
                  "T1-Expanding-xylem-18",
                  "T1-Expanding-xylem-19",
                  "T1-Lignified-xylem-20",
                  "T1-Lignified-xylem-21",
                  "T1-Lignified-xylem-22",
                  "T1-Lignified-xylem-23",
                  "T1-Lignified-xylem-24",
                  "T1-Lignified-xylem-25",
                  "T2-Phloem-01",
                  "T2-Phloem-02",
                  "T2-Phloem-03",
                  "T2-Phloem-04",
                  "T2-Phloem-05",
                  "T2-Cambium-06",
                  "T2-Cambium-07",
                  "T2-Cambium-08",
                  "T2-Cambium-09",
                  "T2-Cambium-10",
                  "T2-Cambium-11",
                  "T2-Expanding-xylem-12",
                  "T2-Expanding-xylem-13",
                  "T2-Expanding-xylem-14",
                  "T2-Expanding-xylem-15",
                  "T2-Expanding-xylem-16",
                  "T2-Expanding-xylem-17",
                  "T2-Expanding-xylem-18",
                  "T2-Expanding-xylem-19",
                  "T2-Lignified-xylem-20",
                  "T2-Lignified-xylem-21",
                  "T2-Lignified-xylem-22",
                  "T2-Lignified-xylem-23",
                  "T2-Lignified-xylem-24",
                  "T2-Lignified-xylem-25",
                  "T2-Lignified-xylem-26",
                  "T3-Phloem-01",
                  "T3-Phloem-02",
                  "T3-Phloem-03",
                  "T3-Phloem-04",
                  "T3-Phloem-05",
                  "T3-Cambium-06",
                  "T3-Cambium-07",
                  "T3-Cambium-08",
                  "T3-Cambium-09",
                  "T3-Cambium-10",
                  "T3-Cambium-11",
                  "T3-Cambium-12",
                  "T3-Cambium-13",
                  "T3-Cambium-14",
                  "T3-Expanding-xylem-15",
                  "T3-Expanding-xylem-16",
                  "T3-Expanding-xylem-17",
                  "T3-Expanding-xylem-18",
                  "T3-Expanding-xylem-19",
                  "T3-Expanding-xylem-20",
                  "T3-Expanding-xylem-21",
                  "T3-Lignified-xylem-22",
                  "T3-Lignified-xylem-23",
                  "T3-Lignified-xylem-24",
                  "T3-Lignified-xylem-25",
                  "T3-Lignified-xylem-26",
                  "T3-Lignified-xylem-27",
                  "T3-Lignified-xylem-28",
                  "T4-Phloem-01",
                  "T4-Phloem-02",
                  "T4-Phloem-03",
                  "T4-Phloem-04",
                  "T4-Phloem-05",
                  "T4-Cambium-06",
                  "T4-Cambium-07",
                  "T4-Cambium-08",
                  "T4-Cambium-09",
                  "T4-Cambium-10",
                  "T4-Cambium-11",
                  "T4-Cambium-12",
                  "T4-Expanding-xylem-13",
                  "T4-Expanding-xylem-14",
                  "T4-Expanding-xylem-15",
                  "T4-Expanding-xylem-16",
                  "T4-Expanding-xylem-17",
                  "T4-Expanding-xylem-18",
                  "T4-Expanding-xylem-19",
                  "T4-Expanding-xylem-20",
                  "T4-Lignified-xylem-21",
                  "T4-Lignified-xylem-22",
                  "T4-Lignified-xylem-23",
                  "T4-Lignified-xylem-24",
                  "T4-Lignified-xylem-25",
                  "T4-Lignified-xylem-26",
                  "T4-Lignified-xylem-27",
                  "T4-Lignified-xylem-28",
                  "gene_id")
aspwoodtpm <- aspwoodtpm[orderaspwood]
rownames(aspwoodtpm) <- aspwoodtpm$gene_id
# aspdata <- as.matrix(select(aspwoodtpm, c(1, 1:107)))
aspdata <- subset(aspwoodtpm, select = grep("T1-*", colnames(aspwoodtpm)))
#'
atnnotation <- read.delim(here("/mnt/picea/home/schoudhary/shruti/ERF85GeneExp/doc/potra_atgenes.txt"), 
                          header = FALSE, sep = "\t")
colnames(atnnotation) <- c("Potra_ID", "AT_Symbols")
degAnot <- atnnotation[match(rownames(aspwoodtpm), atnnotation$Potra_ID),]
#'
all(rownames(aspwoodtpm) == degAnot$Potra_ID)
#'
#'
#' Heatmap function:
#' 
hmap2 <- function(selGene, file_name) {
  tres1 <- aspdata[rownames(aspdata) %in% selGene, ]
  tres1 <- tres1[rowSums(tres1 != 0) > 0, ]
  if (nrow(tres1) == 0) stop("No genes found with non-zero expression.")
  
  svg(file.path("~/shruti/SNR-u2023011/analysis/plots/", paste0(file_name, ".svg")),
      width = 12, height = 12, pointsize = 8)
  
  heatmap_result <- heatmap.2(
    t(scale(t(tres1))), distfun = pearson.dist,
    hclustfun = function(X) hclust(X, method = "ward.D2"),
    trace = "none", col = hpal, margins = c(8, 8), cexCol = 0.1,
    cexRow = 0.1, key = TRUE, keysize = 1, main = file_name,
    Colv = FALSE, Rowv = TRUE, dendrogram = "row",
    labRow = paste(rownames(tres1), degAnot$AT_Symbols[match(rownames(tres1), degAnot$Potra_ID)])
  )
  
  dev.off()
  
  gene_order <- rownames(tres1)[rev(heatmap_result$rowInd)]
  write.table(gene_order, file = file.path("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/plots/", paste0(file_name, "_gene_order.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}
#'
hmap2(upKno, "upKno")
#'
#'
#' If the logTPM+1 is needed (optional)
#' 
# aspdatalog <- log2(aspdata + 1)
# tres <- aspdatalog[rownames(aspdatalog) %in% upKno, ]
# tres1 <- tres[rowSums(tres != 0) > 0, ]
# max(tres1)
# tres2 <- tres1[rowSums(tres1[])>1,]
# pheatmap(tres2, fontsize = 7, cluster_cols = FALSE, color = mako(100),
#          clustering_method = "ward.D2", border_color = NA)