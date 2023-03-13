# https://kkorthauer.org/fungeno2019/singlecell/vignettes/1.2-preprocess-droplet.html#4_quality_control_on_the_cells

# It should be run on unfiltered matrix always
setwd("data/SeuratOut/dropletUtilis/kcl/")

fname <- file.path("../../../CellRangerCount/kcl/outs/raw_feature_bc_matrix/")

library(DropletUtils)
library(scater)

sce <- read10xCounts(fname, col.names=TRUE)
sce
class(counts(sce))
head(rownames(sce))

#1. Calling cells from empty droplets
# 1.1. Testing for deviations from ambient expression
# https://kkorthauer.org/fungeno2019/singlecell/vignettes/1.2-preprocess-droplet.html#2_setting_up_the_data
bcrank <- barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# test whether expression profile for each cell barcode is significantly 
# different from ambient RNA pool. 
set.seed(100)
e.out <- emptyDrops(counts(sce))
is.cell <- e.out$FDR <= 0.001
# sum(is.cell, na.rm=TRUE)
sum(e.out$FDR <= 0.001, na.rm=TRUE)
table(Limited=e.out$Limited, Significant=is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability")

# Subset the sce object to retain only detected cells.
sce <- sce[,which(e.out$FDR <= 0.001)]

# 1.2. Examining cell-calling diagnostics
full.data <- read10xCounts(fname, col.names=TRUE)
set.seed(100)
limit <- 100   
all.out <- emptyDrops(counts(full.data), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
     xlab="P-value", main="", col="grey80") 

#3. Normalizing for cell-specific biases
library(scran)
library(BiocSingular)

set.seed(1000)
sce1 <- read10xCounts(fname, col.names=TRUE)
clusters <- quickCluster(sce1, BSPARAM=IrlbaParam())
table(clusters)

sce1 <- computeSumFactors(sce1, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce1))

plot(sce$total_counts, sizeFactors(sce), log="xy")







