setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/data/")
#' 
#' GO enrichment for clusters
#'
set.seed(42)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(Seurat)
  library(purrr)
  library(tibble)
})
#' 
suppressMessages(source("../UPSCb-common/src/R/topGoUtilities.R"))
goannot <- prepAnnot(mapping = "/mnt/reference/Populus-tremula/v2.2/gopher/gene_to_go.tsv")
ctrl_marker <- readRDS("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/ctrl_marker.rds")
background <- readRDS("/mnt/ada/projects/aspseq/htuominen/SNR-results/ctrl_bg.rds")
#' 
#' 
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")
#'
#'
#' Gene list per cluster
deg.ls <- split(rownames(ctrl_marker), f = ctrl_marker$cluster)
#'
#'
#' Run GO enrichment
#' 
enr.list <- lapply(deg.ls, function(genes) {
  topGO(genes, background = background, annotation = goannot,
        p.adjust = "BH", alpha = 0.05)}) 
#'
#'
#' Function to extract and summarize GO results
#' 
extractEnrichmentResults <- function(enrichment, count = 20, plot = T, 
                                     save_csv = T) {
  all_results <- list()
  
  for (cluster in names(enrichment)) {
    for (gocat in names(enrichment[[cluster]])) {
      dat <- enrichment[[cluster]][[gocat]]
      if (!is.null(dat) && nrow(dat) > 0) {
        dat <- dat %>%mutate(Cluster = cluster, GO_Category = gocat,
            GeneRatio = Significant / Annotated, adjustedPvalue = as.numeric(FDR),
            Count = as.numeric(Significant)) %>%
          arrange(adjustedPvalue) %>% slice_head(n = count)
        all_results[[paste(cluster, gocat, sep = "_")]] <- dat
      }
    }
  }
  
  result_df <- bind_rows(all_results)
  
  if (nrow(result_df) == 0) {
    message("No enrichment found.")
    return(NULL)
  }
  
  if (save_csv) {
    write.csv(result_df, "GO_summary_table.csv", row.names = FALSE)
  }
  
  # Create a dot plot summary
  if (plot) {
    p <- ggplot(result_df, aes(x = Cluster, y = Term, color = adjustedPvalue,
      size = Count)) +
      geom_point() + 
      scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
      theme_bw() + ylab("GO Term") + xlab("Cluster") +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave("GO_summary_plot.pdf", p, width = 10, height = 8)
  }
  return(result_df)
}
#'
#'
#' Run summary extraction and plot
#' 
go_summary <- extractEnrichmentResults(enr.list, count = 10)
