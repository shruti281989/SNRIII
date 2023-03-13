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
})

#' Load markers from step 5
load ("output/afterDbltRemoval/markersWilcox.RData")

suppressMessages({
  # source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
  # source(here("UPSCb-common/src/R/featureSelection.R"))
  # source(here("UPSCb-common/src/R/volcanoPlot.R"))
  source(here("UPSCb-common/src/R/gopher.R"))
})

DefaultAssay(seurat_integrated) <- "RNA"

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

deg.ls <- split(rownames(markers), f = markers$cluster)
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

filt_counts <- GetAssayData(object = seurat_integrated, slot = "counts")
nonzero <- filt_counts > 0

#' Sums all TRUE values and returns TRUE if more than 3 TRUE values per gene
keep_genes3 <- Matrix::rowSums(nonzero) >= 3

filt_seurat <- CreateSeuratObject(filt_counts,
                                      meta.data = seurat_integrated@meta.data)

bg <- list(rownames(seurat_integrated))
#' When I use rownames(filtered_seurat), the enrichment doesn't work
#' 
enr.list <- lapply(gene.ls,function(r){
  lapply(r,gopher,task=list("go","kegg","pfam"),background = bg,
         url="potra2")
})

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
