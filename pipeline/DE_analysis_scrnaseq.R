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
#' Seurat has been run per seuratSNRIII.R script
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

seurat_integrated  <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
setwd("~/shruti/SNRIII/data/SeuratOut/output/afterDbltRemoval/pseudobulkDEClustWiseAftrDblt/")

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = seurat_integrated, slot = "counts", assay="RNA")
metadata <- seurat_integrated@meta.data

# Set up metadata as desired for aggregation and DE analysis
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6" 
#check it should be integ$integrated_snn_res.0.6

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

keepClusters <-as.character(c(0:20, 0:20))
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

# filter FDR < 0.01, |logFC| > 0.1 & sort by FDR
res_fil <- lapply(res, 
                  function(u)  u %>% 
                    dplyr::filter(p_adj < 0.01, abs(logFC) > 0.1) %>% 
                    dplyr::arrange(p_adj))

## Count the number of differential gene findings by cluster.
# nb. & % of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
cbind(cluster=keepClusters, numDE_genes=n_de, 
      percentage = round(n_de / nrow(interestingClusters) * 100, digits =2)) %>%  kable()

# mkdir clustWiseAllDE
for(cluster in 1:length(keepClusters)){
  # Full results
  # filePath <- paste0("data/SeuratOut/output/afterDbltRemoval/pseudobulkDEClustWiseAftrDblt/Cluster", keepClusters[cluster])
  # out <- res[[cluster]][,c("gene", "logFC", "logCPM", "p_adj")]
  # write.csv(out, file = paste0(filePath, "_ctrlkno.csv"), quote=F, row.names = F)
  # 
  # Sig genes
  # filePath <- paste0("data/SeuratOut/output/afterDbltRemoval/pseudobulkDEClustWiseAftrDblt/clustWiseSigDE/Cluster", keepClusters[cluster])
  filePath <- paste0("data/SeuratOut/output/afterDbltRemoval/", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "logCPM", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "ctrlkno.csv"), quote=F, row.names = F)
  
}
# Clear workspace and restart R

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

tx2gene <- suppressMessages(read_delim
                            (here("~/shruti/single_cell_analysis_poplar/reference/annotation/tx2gene.tsv.gz"),
                              delim="\t", col_names=c("TXID","GENE")))
setwd("~/shruti/SNRIII/data/edgeRscBulkDE/")
files <- dir(".", recursive=TRUE, 
             pattern="quant.sf", full.names =TRUE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi$countsFromAbundance)
txi

cts <- txi$counts
saveRDS(cts, file="cts_SNRIII.rds")
colnames(cts) <- c("KCl","KNO")
write.table(x, file, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

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

# Restart R
suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(hyperSpec)
  library(RColorBrewer)
  library(gplots)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(pheatmap)
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' DEGs in SNRIII
degSnr <-read.delim(here("~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNucl/degSNRIII.txt"), header = T,sep = "\t")
upKno <- degSnr %>% filter (degSnr$status == "up in kno")  %>%pull(gene)
dnkno <- degSnr %>% filter (degSnr$status == "down in kno")  %>%pull(gene)

aspwood <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/AspWood_tpm.txt", header = TRUE)
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
aspdata <- as.matrix(select(aspwoodtpm, c(1, 1:107)))

atnnotation <- read.delim(here("~/shruti/ERF85GeneExp/doc/potra_atgenes.txt"), header = FALSE, sep = "\t")
colnames(atnnotation) <- c("Potra_ID", "AT_Symbols")
degAnot <- atnnotation[match(rownames(aspwoodtpm), atnnotation$Potra_ID),]
all(rownames(aspwoodtpm) == degAnot$Potra_ID)

hmap2 <- function(selGene, file_name) {
  tres <- aspdata[rownames(aspdata) %in% selGene, ]
  tres1 <- tres[rowSums(tres != 0) > 0, ]
  # png(file.path(here("~/shruti/SNR-u2023011/analysis/plots/"),
  #               paste0(file_name,".png")), res= 250,height = 2000, width = 2000)
  # 
  svg(file.path(here("~/shruti/SNR-u2023011/analysis/plots/"),
                paste0(file_name,".svg")), pointsize = 8)

  heatmap.2(t(scale(t(tres1))),
            distfun = pearson.dist,
            hclustfun = function(X){hclust(X,method="ward.D2")},
            trace="none", col=hpal, margins =c(18,18), cexCol = 0.1,
            cexRow = 0.1, main = file_name, key = TRUE, keysize = 1,
            Colv = FALSE, Rowv = TRUE, dendrogram = "row",
            labRow = paste(rownames(tres1), degAnot$AT_Symbols[match(rownames(tres1), degAnot$Potra_ID)])
  )
  dev.off()
}

hmap2(dnkno,"dnKno")

# If the logTPM+1 is needed
aspdatalog <- log2(aspdata + 1)
tres <- aspdatalog[rownames(aspdatalog) %in% upKno, ]
tres1 <- tres[rowSums(tres != 0) > 0, ]
max(tres1)
tres2 <- tres1[rowSums(tres1[])>1,]

pheatmap(tres2, 
         fontsize = 7,
         cluster_cols = FALSE,
         color = mako(100),
         clustering_method = "ward.D2",
         border_color = NA)
