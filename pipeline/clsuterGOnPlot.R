setwd("data/SeuratOut/output/afterDbltRemoval/")
#' 
#' GO enrichment for clusters
#'
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
#' 
#' 
suppressMessages(source("~/shruti/SNRIII/UPSCb-common/src/R/topGoUtilities.R"))
goannot <- prepAnnot(mapping = "/mnt/reference/Populus-tremula/v2.2/gopher/gene_to_go.tsv")
ctrl_marker <- readRDS("~/shruti/SNR-u2023011/analysis/markers/ctrl_marker.rds")
background <- readRDS("/mnt/ada/projects/aspseq/htuominen/SNR-results/ctrl_bg.rds")
#' 
#' 
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")
deg.ls <- split(rownames(ctrl_marker), f = ctrl_marker$cluster)
res.list <- list(deg.ls)
suppressMessages(enr.list <- lapply(res.list,function(r){
  lapply(r,topGO,background=background,
         annotation=goannot,p.adjust="BH", alpha=0.05)
}))

suppressWarnings(extractEnrichmentResults(enr.list, count = 30))

extractEnrichmentResults <- function(enrichment,
                                     go.namespace=c("BP","CC","MF"),
                                     count=100,plot=TRUE){
  
  # sanity
  if(is.null(unlist(enrichment)) | length(unlist(enrichment)) == 0){
    message("No GO enrichment for",names(enrichment))
  } else {
    if(plot){
      gocatname <- c(BP="Biological Process",
                     CC="Cellular Component",
                     MF="Molecular Function")
      lapply(names(enrichment),function(n){
        lapply(names(enrichment[[n]]),function(gocat){
          dat <- enrichment[[n]][[gocat]]
          if(is.null(dat)){
            message("No GO enrichment for ",n," in category ",gocatname[gocat])
          } else {
            dat$GeneRatio <- dat$Significant/dat$Annotated
            dat$adjustedPvalue <- as.numeric(dat$FDR)
            dat$Count <- as.numeric(dat$Significant)
            dat <- dat[order(dat$GeneRatio),]
            if(nrow(dat) > count){ dat <- dat[1:count,] }
            dat$Term <- factor(dat$Term, levels = unique(dat$Term))
            ggplot(dat, aes(x =Term, y = GeneRatio, color = adjustedPvalue, size = Count)) + 
              geom_point() +
              scale_color_gradient(low = "red", high = "blue") +
              theme_bw() + 
              ylab("GeneRatio") + 
              xlab("") + 
              ggtitle(paste0("GO enrichment: ",n," ",gocatname[gocat])) +
              coord_flip()
          }
        })
      })
    }
  }
}