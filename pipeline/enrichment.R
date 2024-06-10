#' 9. GO enrichment

set.seed(42)
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(RCurl)
  library(cowplot)
  library(viridis)
  library(qs)
  library(magrittr)
  library(Matrix)
  library(purrr)
  library(reshape2)
  library(S4Vectors)
  library(tibble)
  library(ComplexHeatmap)
  library(pheatmap)
  library(scales)
  library(ggplot2)
  library(here)
  library(RColorBrewer)
})

#' Load markers from step 5
integ <- readRDS("~/shruti/SNR-u2023011/analysis/snrIII/integNucl.rds")
DefaultAssay(integ) <- "RNA"

split_seurat <- SplitObject(integ, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]
ctrl <- split_seurat$ctrl
DefaultAssay(ctrl) <- "RNA"

ctrl <- NormalizeData(ctrl, verbose = FALSE)

# subset to 300 cells 
# sub <- subset(ctrl, cells = WhichCells(integ, downsample = 200))
# table(sub@active.ident)
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
# 300 300 300 300 300 300 300 300 300 300 287 261 289 159 206 241 226 
# 17  18  19  20 
# 191 119  87  24

# Compute differentiall expression
# markers_genes_sub <- FindAllMarkers(
#   sub, logfc.threshold = -Inf, test.use = "wilcox", min.pct = 0.05,
#   min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, assay = "RNA")

markers_genes <- FindAllMarkers(
  ctrl, logfc.threshold = -Inf, test.use = "wilcox", min.pct = 0.05,
  min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, assay = "RNA")

gene_rank <- setNames(markers_genes$avg_log2FC, 
                      casefold(rownames(markers_genes), upper = T))

# saveRDS(markers_genes_sub, file = "../SNR-u2023011/analysis/markers/crl_sub.rds")
# gene_rank <- setNames(crl_sub$avg_log2FC, 
#                       casefold(rownames(crl_sub), upper = T))

# gopher is down, use TopGO instead
suppressMessages({
  source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
  source(here("UPSCb-common/src/R/featureSelection.R"))
  source(here("UPSCb-common/src/R/volcanoPlot.R"))
  source(here("UPSCb-common/src/R/gopher.R"))
  source(here("UPSCb-common/src/R/topGoUtilities.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

deg.ls <- split(rownames(markers_genes), f = markers_genes$cluster)
# deg.ls <- split(rownames(crl_sub), f = crl_sub$cluster)
# deg.ls is a list, still need to make a list for enrichment
gene.ls <- list(deg.ls)

#' Background to be used: Suggestions from Nico and Nat
#' use only the set of genes expressed in the tissue you are looking at
#' either devise that across all of your samples:
#' naïve approach: exclude genes that are always 0
#' more advanced, check the average expression independent filtering cutoff available in your DESeq results - if you use DESeq:
#' Cutoff off the independent filtering
#' cutoff <- metadata(res)$filterThreshold
#' We keep the rest of the data
# sel <- rowMeans(counts(dds_tx,normalized=TRUE)) > cutoff`
# independent filtering is how DESeq decides on the signal-to-noise ratio
# 2. your data is to sparse for 1. to work then use the set of genes expressed 
# in your tissue type (wood) from a bulk RNA resource

# filt_counts <- GetAssayData(object = integ, slot = "counts")
#' nonzero <- filt_counts > 0
#' 
#' #' Sums all TRUE values and returns TRUE if more than 3 TRUE values per gene
#' keep_genes3 <- Matrix::rowSums(nonzero) >= 3
#' 
#' filt_seurat <- CreateSeuratObject(filt_counts,
#'                                       meta.data = integ@meta.data)

# bg <- list(rownames(integ))
#' When I use rownames(filtered_seurat), the enrichment doesn't work

#' gopher down
# enr.list <- lapply(gene.ls,function(r){
#   lapply(r,t,task=list("go","kegg","pfam"),background = bg,
#          url="potra2")
# })

background <- rownames(integ)
goannot <- prepAnnot(mapping = "/mnt/picea/storage/reference/Populus-tremula/v2.2/gopher/gene_to_go.tsv")

res.list <- list(deg.ls)
suppressMessages(enr.list <- lapply(res.list,function(r){
  lapply(r,topGO,background=background,annotation=goannot,alpha=0.1,p.adjust="none")
}))

# Code from Aman on visualization as in Chen et al., 2021
list <- enr.list[[1]]
df <- bind_rows(list$`0`$go, list$`1`$go, list$`2`$go, list$`3`$go, list$`4`$go, list$`5`$go, list$`6`$go, list$`7`$go, list$`8`$go, list$`9`$go, list$`10`$go, list$`11`$go, .id = "Cluster")
df %<>% rowwise() %>% mutate(Cluster = as.numeric(Cluster) - 1, .keep = "unused") 

BP <- df %>% filter(namespace %in% "BP")
MF <- df %>% filter(namespace %in% "MF")
CC <- df %>% filter(namespace %in% "CC")

#plotting function
pp <- function(df, go) {
  df %>% 
    ggplot(aes(x = factor(Cluster), y = name, color = padj, size = nt)) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    theme(text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5)) + 
    ylab(paste("GO Terms -", go)) + 
    xlab("Clusters") + 
    ggtitle("GO analysis")
}

#downloading function
dl <-  function(go) {
  pt %>% 
    download_this(output_name = go,
                  output_extension = ".pdf",
                  button_label = "Download plot",
                  button_type = "warning",
                  has_icon = TRUE,
                  icon = "fa fa-save")
} 
BP
MF
CC

pt <- pp(df = BP, go = "BP")
pt
# for i in length()
as.data.frame(enr.list[[1]]$`0`$go) %>% 
  ggplot(aes(x= id, y = name, color = padj, size = nt)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis")
