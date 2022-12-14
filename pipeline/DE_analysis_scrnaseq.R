#' ---
#' title: "Single Cell pseudobulk and bulk DE analysis for quick nitrate repsonse in poplar"
#' author: "Shruti"
#' #' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#'
#' Following parts of tutorials at:
#' https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
#' https://github.com/hbc/knowledgebase/blob/master/scrnaseq/pseudobulkDE_edgeR.md
#' https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#' https://bioconductor.riken.jp/packages/3.14/bioc/vignettes/Glimma/inst/doc/single_cell_edger.html
#' Seurat has been run per nitrate.R script
#'
# PART1: Single Cell pseudobulk DE analysis
suppressPackageStartupMessages({
  library(scater)
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(Matrix.utils)
  library(edgeR)
  library(Matrix)
  library(reshape2)
  library(S4Vectors)
  library(SingleCellExperiment)
  library(pheatmap)
  library(png)
  library(RColorBrewer)
  library(limma)
  library(magrittr)
  library(gridExtra)
  library(knitr)
  library(limma)
})

seurat <- readRDS("~/shruti/single_cell_analysis_poplar/data/nitrate_resp/integerated_seurat/seurat_integrated.rds")
setwd("~/shruti/single_cell_analysis_poplar/data/nitrate_resp/salmonQuant_out/DE_analysis_scrnaseq")

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = seurat_integrated, slot = "counts", assay="RNA")
metadata <- seurat_integrated@meta.data

# Set up metadata as desired for aggregation and DE analysis
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
metadata$cluster_id <- factor(seurat_integrated@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
#
# We have already filtered the cells
# Just Remove lowly expressed genes <10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
nk <- length(kids)

# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$sample)))
ns <- length(sids)

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sample)

## Turn class "table" into a named vector of cells per sample
n_cells <- table(sce$sample) %>%  as.vector()
names(n_cells) <- names(table(sce$sample))

## Match the named vector with metadata to combine it
m <- match(names(n_cells), sce$sample)

## Create the sample level metadata by selecting specific columns
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("sample", "n_cells")
kable(ei)

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample")]
groups$sample <- factor(groups$sample)

# Aggregate across cluster-sample groups
# Each row corresponds to aggregate counts for a cluster-sample combo
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",n = 2), `[`, 1)

# Split data and turn into a list
# Each component corresponds to a cluster; storing associated expression matrix (counts)
# Transform data i.e, so rows are genes and columns are samples 
pb <- split.data.frame(pb,factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), gsub(".*_", "", rownames(u))))

# Print out the table of cells in each cluster-sample group
options(width = 100)
kable(table(sce$cluster_id, sce$sample))

keepClusters <-as.character(c(0:11, 0:11))
(interestingClusters <- SingleCellExperiment(assays = pb[keepClusters]))

(design <- model.matrix(~ 0 + ei$sample) %>% 
    set_rownames(ei$sample) %>% 
    set_colnames(levels(factor(ei$sample))))

(contrast <- makeContrasts("ctrl-kno", levels = design))
bcv <- 0.2

# Run edgeR w/ default parameters- exact test as we have no reps
res <- lapply(keepClusters, function(k) {
  y <- assays(interestingClusters)[[k]]
  y <- DGEList(y, group=1:2, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  et <- exactTest(y, dispersion=bcv^2)
  # y <- estimateDisp(y, design)
  # fit <- glmQLFit(y, design)
  # fit <- glmQLFTest(fit, contrast = contrast)
  topTags(et, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)
})

# filter FDR < 0.05, |logFC| > 1 & sort by FDR
res_fil <- lapply(res, 
                  function(u)  u %>% 
                    dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% 
                    dplyr::arrange(p_adj))

## Count the number of differential gene findings by cluster.
# nb. & % of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
cbind(cluster=keepClusters, numDE_genes=n_de, 
      percentage = round(n_de / nrow(interestingClusters) * 100, digits =2)) %>%  kable()

for(cluster in 1:length(keepClusters)){
  # Full results
  filePath <- paste0("./all_genes/Cluster", keepClusters[cluster])
  out <- res[[cluster]][,c("gene", "logFC", "logCPM", "p_adj")]
  write.csv(out, file = paste0(filePath, "_ctrlkno.csv"), quote=F, row.names = F)
  
  # Sig genes
  filePath <- paste0("./sig0.05/Cluster", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "logCPM", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "ctrlkno.csv"), quote=F, row.names = F)
  
}
# Clear workspace and restart R
rm(ls())
# Part 2: bulk DE analysis
suppressPackageStartupMessages({
  library(tximport)
  library(GenomicFeatures)
  library(readr)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplot2)
  library(gplots)
  library(ggvenn)
  library(ggrepel)
  library(edgeR)
  library(data.table)
  library(here)
  library(hyperSpec)
  library(parallel)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(vsn)
  library(VennDiagram)
})

tx2gene <- suppressMessages(read_delim(here("~/shruti/single_cell_analysis_poplar/reference/annotation/tx2gene.tsv.gz"),delim="\t", col_names=c("TXID","GENE")))
setwd("~/shruti/single_cell_analysis_poplar/data/nitrate_resp/salmon/edger/")
files <- dir("~/shruti/single_cell_analysis_poplar/data/nitrate_resp/salmon/", recursive=TRUE, pattern="quant.sf", full.names =TRUE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi$counts)

cts <- txi$counts
colnames(cts) <- c("KCl","KNO")

#1. create a DGElist object:
y <- DGEList(counts=cts, group=1:2, genes = rownames(cts), remove.zeros = TRUE)
dim(y)
y$samples

#2.Filtering
y.full <- y
heatmap(cor(y.full$counts))
head(y$counts)

#retain genes represented at least 1cpm reads in at least 1 sample (cpm=counts per million)
countsPerMillion <- cpm(y)
summary(countsPerMillion)
head(cpm(y))

#find total gene counts per sample
# minimum of 1 count-per-million (CPM): ignores noisy background counts.
# keep genes that have minimum 1 CPM across 1 sample (since each of test groups consist of 2 samples).
apply(y$counts, 1, sum)
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep,]
dim(y)
heatmap(cor(y.full$counts))

# #reset the library sizes after filtering:
# y$samples$lib.size <- colSums(y$counts)
# y$samples

#3. Normalize
y <- calcNormFactors(y, method = "TMM")
y

eff.lib.size <- y$samples$lib.size*y$samples$norm.factors
normCounts <- cpm(y)
pseudoNormCounts <- log2(normCounts + 1)
boxplot(pseudoNormCounts, col="gray", las=3)
plotMDS(pseudoNormCounts)

# 4. DGE with exactTest() conducts tagwise tests using the exact negative binomial test. 
# default is to compare between first two groups:
bcv <- 0.2
y <- DGEList(y, group=1:2, remove.zeros = TRUE)
y <- calcNormFactors(y)
et <- exactTest(y, dispersion=bcv^2)

et$table$padj <- p.adjust(et$table$PValue, method="BH")
sum(et$table$padj < 0.01)

ctrl.up <- et$table$logFC > 0.5 & et$table$padj < 0.05
ctrl.dn <- et$table$logFC < 0.5 & et$table$padj < 0.05

#Save selected genes (FDR<0.05, lfc)
sel <- et$table$padj <= 0.05 & abs(et$table$logFC) >=0 & ! is.na(et$table$padj)

write.csv(et,file="results.csv")
write.csv(et[sel,],file="genes.csv")