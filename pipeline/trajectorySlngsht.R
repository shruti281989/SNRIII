# https://rdrr.io/github/HectorRDB/condimentsPaper/f/vignettes/TGFB.Rmd
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
# http://barcwiki.wi.mit.edu/wiki/SOP/scRNA-seq/Slingshot
# https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html#using-slingshot
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#overview-1

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  # For analysis
  library(condiments)
  library(slingshot)
  # For data manipulation
  library(dplyr)
  library(tidyr)
  library(gam)
  library(destiny)
  # For visualization
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
  library(scran)
  library(scater)
  library(igraph)
  library(Polychrome)
  library(ggbeeswarm)
  library(ggthemes)
  library(grDevices)
  library(tradeSeq)
  library(uwot)
  library(mclust)
  library(magrittr)
  library(bioc2020trajectories)
  library(tradeSeq)
  library(UpSetR)
  library(pheatmap)
  library(gridExtra)
})
set.seed(42)
# theme_set(theme_classic())

integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
integ$seurat_clusters <- integ@active.ident

# convert to single cell experiment
sce <- as.SingleCellExperiment(integ, assay = "RNA")
rm(integ)

# go from factor to character
colData(sce)$active.ident <- as.character(sce$ident)
# saveRDS(sce1, file="~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/slngSht/sce1.rds")

# 1. first look pseudotime based on PCA
# load sce1
sce <- runPCA(sce, ncomponents = 50)
pca <- reducedDim(sce, "PCA")
sce$PC1 <- pca[, 1]
sce$PC2 <- pca[, 2]

colourCount = length(unique(sce$ident))
getPalette = colorRampPalette(brewer.pal(9, "21"))

sce$pseudotime_PC1 <- rank(sce$PC1)  # rank cells by their PC1 score
p2<- ggplot(as.data.frame(colData(sce1)), 
       aes(x = pseudotime_PC1, y = ident,colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_fill_manual(values = getPalette(colourCount), 
                    guide = guide_legend(nrow=2)) +
  # scale_color_tableau() + 
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

pdf("analysis/plots/pcaPstime.pdf")
# print(p1)     # Plot 1 --> in the first page of PDF
print(p2)     # Plot 2 ---> in the second page of the PDF
dev.off() 

# 2. pseudotime
sce <- readRDS("~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/slngSht/sce.rds")
sce1 <- slingshot(sce, clusterLabels = as.character(sce$ident), reducedDim = 'UMAP',
                 approx_points = 3000, allow.breaks = FALSE)

# Create a ggplot with many colors 
nb.cols <- 42
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
ggplot(as.data.frame(colData(sce)), 
       aes(x = sce1$slingPseudotime_1, y = ident, colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_minimal() +
  scale_fill_manual(values = mycolors)

# Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
# Identify the variable genes by ranking all genes by their variance.
Y <- log2(counts(sce) + 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:2000]
Y <- Y[var1K, ]  # only counts for variable genes

# Fit GAM for each gene using pseudotime as independent variable.
t <- sce$slingPseudotime_1
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]  

# # Filter and normalize
# # geneFilter <- apply(assays(sce)$counts,1,function(x){
# #   sum(x >= 3) >= 10
# # })
# # sce <- sce[geneFilter, ]
# 
# FQnorm <- function(counts){
#   rk <- apply(counts,2,rank,ties.method='min')
#   counts.sort <- apply(counts,2,sort)
#   refdist <- apply(counts.sort,1,median)
#   norm <- apply(rk,2,function(r){ refdist[r] })
#   rownames(norm) <- rownames(counts)
#   return(norm)
# }
# assays(sce)$norm <- FQnorm(assays(sce)$counts)

# 
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = sce$ident, #colData(sce)$sample
                 allow.breaks = FALSE, approx_points = 200)
# sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$sample,
#                  start.clus = c('0',"14","1","10","20", "16"), approx_points = 150)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

# Calculate pseudotime and lineages with Slingshot:
# use either PCA/UMAP/TSNE for reduction. Since Slingshot is designed where cells
# fall along a continuous trajectory, better to use PCA if clearly discrete 
# clusters presented in UMAP
lnes <- getLineages(reducedDim(sce,"UMAP"), sce$ident)
lnes

# this define the cluster color. You can change it with different color scheme.
my_color <- createPalette(length(levels(sce$ident)), c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(sce$ident))

slingshot_df <- data.frame(colData(sce)$ident)

# re-order y-axis for better figure: This should be tailored with your own cluster names
# slingshot_df$ident = factor(slingshot_df$ident, levels=c(0,20,14,1,10))

ggplot(slingshot_df, aes(x = sce$slingPseudotime_1, y = ident, 
                         colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)

plot(reducedDims(sce)$UMAP, col = my_color[as.character(sce$ident)], 
     pch=16, asp = 1)
legend("bottomleft",legend = names(my_color[levels(sce$ident)]),  
       fill = my_color[levels(sce$ident)])
lines(SlingshotDataSet(lnes), lwd=2, type = 'lineages', col = c("black"))

plotGeneCount(curve = sce, clusters = sce$ident)
topologyTest(SlingshotDataSet(sce), sce$sample, rep = 100,
            methods = "KS_mean", threshs = .01)

df <- bind_cols(
  as.data.frame(reducedDims(sce)$UMAP),
  as.data.frame(colData(sce)[,-3])
) %>%
  sample_frac(1)
curve <- slingCurves(sce)[[1]]
p4 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = slingPseudotime_1)) +
  geom_point(size = .7) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curve$s[curve$ord, ] %>% as.data.frame(),
            col = "black", size = 1.5)
p4
# Kolmogorov-Smirnov Test
ks.test(slingPseudotime(sce)[colData(sce)$sample == "ctrl", 1],
        slingPseudotime(sce)[colData(sce)$sample == "kno", 1])

icMat <- evaluateK(counts = sce,
                   conditions = factor(colData(sce)$sample),
                   nGenes = 300, k = 3:7)

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 2 # use 2 cores

sce <- fitGAM(counts = sce, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = TRUE, parallel=TRUE, BPPARAM = BPPARAM)


pseudotime <- slingPseudotime(sce, na=FALSE)
cellWeights <- slingCurveWeights(sce)
set.seed(42)
sce <- fitGAM(sce, conditions = factor(colData(sce)$sample),nknots = 5)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE)
mean(rowData(sce)$tradeSeq$converged)
rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce)$assocRes
ctrlGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionctrl, "fdr") <= 0.05)
]
knoGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionkno, "fdr") <= 0.05)
]

length(ctrlGenes)
length(knoGenes)

UpSetR::upset(fromList(list(ctrl = ctrlGenes, kno = knoGenes)))

# Visualization of DE genes
### based on mean smoother
yhatSmooth <- predictSmooth(sce, gene = ctrlGenes, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
cl <- sort(cutree(heatSmooth$tree_row, k = 6))
table(cl)

conditions <- colData(sce2)$sample
pt1 <- colData(sce2)$slingshot$pseudotime

### based on fitted values (plotting takes a while to run)
yhatCell <- predictCells(sce, gene=ctrlGenes)
yhatCellctrl <- yhatCell[,conditions == "Ctrl"]

# order according to pseudotime
# oo <- order(pt1[conditions == "ctrl"], decreasing=FALSE)
# yhatCellctrl <- yhatCellctrl[,Ctrl]
# pheatmap(t(scale(t(yhatCellctrl))), cluster_cols = FALSE,
#          show_rownames = FALSE, show_colnames=FALSE)

#  DE genes between conditions
plotSmoothers(sce2, assays(sce2)$counts, gene = "", alpha = 1, border = TRUE) + ggtitle("")

# conditionTest function tests the null hypothesis that genes have identical 
# expression patterns in each condition
condRes <- conditionTest(sce2, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)
conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]
# Visualize most and least significant gene
# plot genes
oo <- order(condRes$waldStat, decreasing = TRUE)

# most and least significant gene
plotSmoothers(sce2, assays(sce2)$counts,
              gene = rownames(assays(sce2)$counts)[oo[1]],
              alpha = 1, border = TRUE)
plotSmoothers(sce2, assays(sce2)$counts,
              gene = rownames(assays(sce2)$counts)[oo[nrow(sce2)]],
              alpha = 1, border = TRUE)

# Heatmaps of genes DE between conditions
### based on mean smoother
yhatSmooth <- predictSmooth(sce2, gene = conditionGenes, nPoints = 50, tidy = FALSE)
yhatSmoothScaled <- t(scale(t(yhatSmooth)))
heatSmooth_kno <- pheatmap(yhatSmoothScaled[, 51:100],
                           cluster_cols = FALSE,
                           show_rownames = FALSE, show_colnames = FALSE, main = "kno", legend = FALSE,
                           silent = TRUE)

matchingHeatmap_ctrl <- pheatmap(yhatSmoothScaled[heatSmooth_kno$tree_row$order, 1:50],
                                 cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = FALSE, show_colnames = FALSE, main = "ctrl",
                                 legend = FALSE, silent = TRUE
)

grid.arrange(heatSmooth_kno[[4]], matchingHeatmap_ctrl[[4]], ncol = 2)

# Gene set enrichment analysis
statsCond <- condRes$waldStat
names(statsCond) <- rownames(condRes)
# do the following
# eaRes <- fgsea(pathways = m_list, stats = statsCond, nperm = 5e4, minSize = 10)
# ooEA <- order(eaRes$pval, decreasing = FALSE)
# kable(head(eaRes[ooEA, 1:3], n = 20))

# 3.
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"
split_seurat <- SplitObject(integ, split.by = "sample")
kcl <- split_seurat[["ctrl"]]
rm(split_seurat)
dimred <- kcl@reductions$umap@cell.embeddings
clustering <- kcl$integrated_snn_res.0.6
counts <- as.matrix(kcl@assays$RNA@counts[kcl@assays$RNA@var.features, ])
lineages <- getLineages(data = dimred, clusterLabels = clustering,
                        end.clus = c("12","15","17")) # start.clus = "0") 
for(i in levels(clustering)){ 
  text( mean(dimred[clustering==i,1]),
        mean(dimred[clustering==i,2]), labels = i,font = 2) }

curves <- getCurves(lineages, approx_points = NULL, thresh = 0.01, stretch = 0.8,
                    allow.breaks = FALSE)
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]

dim(filt_counts)
# [1]     0 11514 why this is zero
sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)
# not working fitGAM
plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)
