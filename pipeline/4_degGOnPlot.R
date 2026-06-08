setwd("data/SeuratOut/output/afterDbltRemoval/")
#' 
#' GO enrichment for DEGs per sample
#'
setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/")
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
  library(RColorBrewer)
})
#' 
#' 
suppressMessages(source("UPSCb-common/src/R/topGoUtilities.R"))
goannot <- prepAnnot(mapping = "/mnt/reference/Populus-tremula/v2.2/gopher/gene_to_go.tsv")
background <- readRDS("/mnt/ada/projects/aspseq/htuominen/SNR-results/ctrl_bg.rds")

dir.create("~/shruti/SNRIII/data/SeuratOut/output/afterDbltRemoval/GO_tables_Wilcoxfdr0.01/", showWarnings = FALSE)
dir.create("~/shruti/SNRIII/data/SeuratOut/output/afterDbltRemoval/GO_plots_Wilcoxfdr0.01/", showWarnings = FALSE)

deg_data <- read.table("data/SeuratOut/degWilcoxpct0.1fdr0.01.txt", header = T)
deg_data <- deg_data %>% mutate(Direction = ifelse(lfc > 0, "up", "down"))

enrichment_results <- list()
summary_table <- data.frame()

#'
for (cl in unique(deg_data$cluster)) {
  for (dir in c("up", "down")) {
    genes <- deg_data %>% filter(cluster == cl, Direction == dir) %>% pull(gene)
    
    if (length(genes) < 10) next
    
    gores <- topGO(set = genes, background = background, annotation = goannot,
                   ontology = c("BP", "MF", "CC"), algorithm = "parentchild",
                   statistic = "fisher", p.adjust = "fdr", alpha = 0.05)
    
    enrichment_results[[cl]][[dir]] <- gores
    
    for (ont in names(gores)) {
      tab <- gores[[ont]]
      if (!is.null(tab) && nrow(tab) > 0) {
        tab$Cluster <- cl
        tab$Direction <- dir
        tab$Ontology <- ont
        summary_table <- bind_rows(summary_table, tab)
      }
    }
  }
}
#'
write.csv(summary_table, "~/shruti/SNRIII/data/SeuratOut/output/afterDbltRemoval/GO_tables_Wilcoxfdr0.01/GO_enrichment_summary.csv", 
          row.names = F)
#'
#' if you want gene list
enrichment_results <- list()
summary_table <- data.frame()

for (cl in unique(deg_data$cluster)) {
  
  dir.create(file.path("~/shruti/SNRIII/data/SeuratOut/output/afterDbltRemoval/GO_Wilcoxfdr0.01", paste0("Cluster_", cl)),
             showWarnings = FALSE)
  
  for (dir in c("up", "down")) {
    
    genes <- deg_data %>% filter(cluster == cl, Direction == dir) %>% pull(gene)
    
    if (length(genes) < 10) next
    
    gores <- topGO(set = genes, background = background, annotation = goannot,
                   ontology = c("BP", "MF", "CC"), algorithm = "parentchild",
                   statistic = "fisher", p.adjust = "fdr", alpha = 0.05)
    
    enrichment_results[[as.character(cl)]][[dir]] <- gores
    
    for (ont in names(gores)) {
      
      tab <- gores[[ont]]
      
      if (!is.null(tab) && nrow(tab) > 0) {
        
        ## ---- Extract genes per GO term ----
        term_genes <- lapply(tab$GO.ID, function(go) {
          intersect(goannot[[go]], genes)
        })
        
        tab$Genes <- sapply(term_genes, paste, collapse = ";")
        
        tab$Cluster <- cl
        tab$Direction <- dir
        tab$Ontology <- ont
        
        summary_table <- bind_rows(summary_table, tab)
        
        ## ---- Save cluster-specific enrichment table ----
        write.csv(
          tab,
          file = file.path("~/shruti/SNRIII/data/SeuratOut/output/afterDbltRemoval/GO_Wilcoxfdr0.01",
                           paste0("Cluster_", cl),
                           paste0("Cluster_", cl, "_", dir, "_", ont, "_GO.csv")
          ),
          row.names = FALSE
        )
      }
    }
  }
}

#'
#'
extractEnrichmentResults <- function(enrichment,
                                     go.namespace = c("BP", "CC", "MF"),
                                     save_plots = TRUE,
                                     output_dir = "GO_plots") {
  gocatname <- c(BP = "Biological Process", CC = "Cellular Component", 
                 MF = "Molecular Function")
  if (save_plots) dir.create(output_dir, showWarnings = FALSE)
  plot_list <- list()
  
  for (n in names(enrichment)) {
    for (direction in names(enrichment[[n]])) {
      for (gocat in intersect(names(enrichment[[n]][[direction]]), go.namespace)) {
        dat <- enrichment[[n]][[direction]][[gocat]]
        if (is.null(dat) || nrow(dat) == 0) next
        
        dat$GeneRatio <- dat$Significant / dat$Annotated
        dat$adjustedPvalue <- as.numeric(dat$FDR)
        dat$Count <- as.numeric(dat$Significant)
        dat <- dat[order(dat$GeneRatio), ]
        dat$Term <- factor(dat$Term, levels = unique(dat$Term))
        
        p <- ggplot(dat, aes(x = Term, y = GeneRatio, color = adjustedPvalue, 
                             size = Count)) + geom_point() +
          scale_color_viridis_c(direction = -1, end = 0.9) + theme_bw() +
          ylab("GeneRatio") + xlab("") +
          ggtitle(paste0("GO enrichment: ", n, " - ", direction, " - ", gocatname[gocat])) +
          coord_flip()
        
        plot_list[[paste(n, direction, gocat, sep = "_")]] <- p
        
        if (save_plots) {
          ggsave(filename = file.path(output_dir, paste0(n, "_", direction, "_", gocat, ".svg")),
                 plot = p, width = 8, height = 6)
        }
      }
    }
  }
  
  return(plot_list)
}
#'
#'
extractEnrichmentResults(enrichment_results)
#'
df <- read.csv("~/shruti/SNRIII/data/SeuratOut/output/afterDbltRemoval/GO_tables_Wilcoxfdr0.01/GO_enrichment_summary.csv")
df <- df %>% filter(FDR < 0.05)
ggplot(df, aes(x = Cluster, y = Term, size = Significant, color = FDR, 
               shape = Direction)) +
  geom_point(alpha = 0.8) + scale_color_viridis_c(direction = -1, end = 0.9) +
  scale_size(range = c(2, 10)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(size = "Gene Count", color = "FDR", shape = "Regulation")
#'
#'
#' GO plots in DEGs from scRNASeq data in figure 5A in manuscript
#' 
df <- read.table("GODeg.txt", sep = "\t", header = T)
df$Cluster.number <- factor(df$Cluster.number, levels = unique(df$Cluster.number))
ggplot(df, aes(x = Cluster.number, y = Description, size = Significant,
               color = q.value, shape = Status)) + geom_point(alpha = 0.8) +
  scale_color_viridis_c(direction = -1, end = 0.9) + 
  scale_size(range = c(2, 10)) + theme_bw() +
  theme(text=element_text(family="Arial"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(size = "Gene Count", color = "q-value", shape = "Regulation")