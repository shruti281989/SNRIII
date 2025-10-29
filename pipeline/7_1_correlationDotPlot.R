setwd("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/dnStrm/correl/")

set.seed(42)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(tidyverse)
  library(pheatmap)
  library(scales)
  library(Matrix)
  library(viridis)
})

integ <- readRDS("/mnt/picea/home/schoudhary/shruti/SNRIII/data/SeuratOut/integ.rds")
load("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/SNRIII-FiltDbltRemIntegSplit.RData")
rc <- as.matrix(split_seurat$ctrl@assays$SCT@data)

merge.rownames <- function (x,y){
  dat <- merge(x = x, y = y, by = "row.names")
  rownames(dat) <- dat$Row.names
  dat <- dat[,-1]
  return(dat)
}

# 1. with Tung
tungtpm <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/publisheddatasets/Tung/GSE180121_Cell_type_TPM_Ptr.txt", sep = "\t", header = TRUE)
tungtpm$Potri <- sub(".v4.1","",tungtpm$Gene.ID)
potrapotri <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/publisheddatasets/potra_potri_BEST_DIAMOND_out.tsv.gz")
colnames(potrapotri) <- c("Potra","Potri")
tungtpmPotra <- inner_join(tungtpm,potrapotri)
tungtpmPotra$Potri <- NULL
tungtpmPotra$Gene.ID <- NULL
write.table(tungtpmPotra, 
            file = "dat/tungtpmPotraLcm.txt", 
            append = FALSE, sep = "\t", quote = F,row.names = F, col.names = TRUE)

TungAvg <- tungtpmPotra %>% select(Potra, Fiber.mean.TPM, Ray.mean.TPM, 
                                   Vessel.bio2.TPM) #vessel 2 probably right vessels
# TungAvg <- tungtpmPotra %>% select(Potra, Fiber.mean.TPM, Ray.mean.TPM, 
                                    # Vessel.bio3.TPM)
rownames(TungAvg) <- TungAvg$Potra
TungAvg$Potra <- NULL

exprAll <- Reduce(merge.rownames, list(TungAvg,rc))
exprAll_stat <- suppressWarnings(sapply(4:ncol(exprAll), function(i) sapply(1:3, function(j) cor.test(exprAll[,i],exprAll[,j],method = "pearson")[c(3,4)])))
exprAll_label=c("Fiber", "Ray", "Vessel")
exprAll_cor <- exprAll_stat[seq(2,nrow(exprAll_stat),2),]
exprAll_pvalue <- exprAll_stat[seq(1,nrow(exprAll_stat)-1,2),]
exprAll_max <- sapply(1:(ncol(exprAll)-3), function(i) max(as.numeric(exprAll_cor[,i])))
exprAll_ident <- sapply(1:(ncol(exprAll)-3), function(i) exprAll_label[which(as.numeric(exprAll_cor[,i])==max(as.numeric(exprAll_cor[,i])))])
exprAll_maxp <- sapply(1:(ncol(exprAll)-3), function(i) as.numeric(exprAll_pvalue[,i])[which(as.numeric(exprAll_cor[,i])==max(as.numeric(exprAll_cor[,i])))])
names(exprAll_max) <- exprAll_ident

split_seurat$ctrl@meta.data$exprAll.ID.P <- as.character(exprAll_ident)
split_seurat$ctrl@meta.data$exprAll.cor.P <- exprAll_max
split_seurat$ctrl@meta.data$exprAll.pvalue.P <- exprAll_maxp
split_seurat$ctrl@meta.data$exprAll.ID.P[which(split_seurat$ctrl@meta.data$exprAll.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order=c("Fiber", "Ray", "Vessel", "unknown")
split_seurat$ctrl$exprAll.ID.P <- factor(split_seurat$ctrl$exprAll.ID.P, levels = order[sort(match(unique(split_seurat$ctrl$exprAll.ID.P),order))]) 
palette <- c("seagreen","blue","pink","grey")
color <- palette[sort(match(unique(split_seurat$ctrl$exprAll.ID.P),order))]
DimPlot(split_seurat$ctrl, group.by="exprAll.ID.P", cols=color)+NoAxes()

# For fib
fibTungAvg <- tungtpmPotra %>% select(Potra, Fiber.mean.TPM)
rownames(fibTungAvg) <- fibTungAvg$Potra
fibTungAvg$Potra <- NULL

# For Ray
rayTungAvg <- tungtpmPotra %>% select(Potra, Ray.mean.TPM)
rownames(rayTungAvg) <- rayTungAvg$Potra
rayTungAvg$Potra <- NULL

# For Vessel2 and 3
# vesTungAvg <- tungtpmPotra %>% select(Potra, Vessel.bio2.TPM)
vesTungAvg <- tungtpmPotra %>% select(Potra, Vessel.bio3.TPM)
rownames(vesTungAvg) <- vesTungAvg$Potra
vesTungAvg$Potra <- NULL

fibData <- Reduce(merge.rownames, list(fibTungAvg,rc))
fib_stat <- sapply(2:ncol(fibData), function(i) sapply(1:1, function(j) cor.test(fibData[,i],fibData[,j],method = "pearson")[c(3,4)]))
fib_cor <- fib_stat[seq(2,nrow(fib_stat),2),]
fib_pvalue <- fib_stat[seq(1,nrow(fib_stat)-1,2),]
fib_cor_un <- unlist(fib_cor)
split_seurat$ctrl@meta.data$fib.cor.P <- fib_cor_un
integ@meta.data$fib.pvalue.P <- fib_pvalue
VlnPlot(integ, features = "fib.cor.P")

rayData <- Reduce(merge.rownames, list(rayTungAvg,rc))
ray_stat <- sapply(2:ncol(fibData), function(i) sapply(1:1, function(j) cor.test(fibData[,i],fibData[,j],method = "pearson")[c(3,4)]))
ray_cor <- ray_stat[seq(2,nrow(ray_stat),2),]
ray_pvalue <- ray_stat[seq(1,nrow(ray_stat)-1,2),]
ray_cor_un <- unlist(ray_cor)
integ@meta.data$ray.cor.P <- ray_cor_un
integ@meta.data$ray.pvalue.P <- ray_pvalue
VlnPlot(integ, features = "ray.cor.P")

vesData <- Reduce(merge.rownames, list(vesTungAvg,rc))
ves_stat <- sapply(2:ncol(fibData), function(i) sapply(1:1, function(j) cor.test(fibData[,i],fibData[,j],method = "pearson")[c(3,4)]))
ves_cor <- ves_stat[seq(2,nrow(ves_stat),2),]
ves_pvalue <- ves_stat[seq(1,nrow(ves_stat)-1,2),]
ves_cor_un <- unlist(ves_cor)
integ@meta.data$ves.cor.P <- ves_cor_un
integ@meta.data$ves.pvalue.P <- ves_pvalue
VlnPlot(integ, features = "ves.cor.P")+NoLegend()

# 2. Aspwood
aspwood <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/AspWood_tpm.txt", header = TRUE)
aspwoodtpm <- dcast(aspwood, gene_id ~ sample_name)
orderaspwood <- c("T1-Phloem-01","T1-Phloem-02","T1-Phloem-03","T1-Phloem-04",
                  "T1-Phloem-05","T1-Cambium-06","T1-Cambium-07","T1-Cambium-08",
                  "T1-Cambium-09","T1-Cambium-10","T1-Cambium-11","T1-Cambium-12",
                  "T1-Expanding-xylem-13","T1-Expanding-xylem-14",
                  "T1-Expanding-xylem-15","T1-Expanding-xylem-16",
                  "T1-Expanding-xylem-17","T1-Expanding-xylem-18",
                  "T1-Expanding-xylem-19","T1-Lignified-xylem-20",
                  "T1-Lignified-xylem-21","T1-Lignified-xylem-22",
                  "T1-Lignified-xylem-23","T1-Lignified-xylem-24",
                  "T1-Lignified-xylem-25","T2-Phloem-01","T2-Phloem-02",
                  "T2-Phloem-03","T2-Phloem-04","T2-Phloem-05","T2-Cambium-06",
                  "T2-Cambium-07","T2-Cambium-08","T2-Cambium-09","T2-Cambium-10",
                  "T2-Cambium-11","T2-Expanding-xylem-12","T2-Expanding-xylem-13",
                  "T2-Expanding-xylem-14","T2-Expanding-xylem-15",
                  "T2-Expanding-xylem-16","T2-Expanding-xylem-17",
                  "T2-Expanding-xylem-18","T2-Expanding-xylem-19",
                  "T2-Lignified-xylem-20","T2-Lignified-xylem-21",
                  "T2-Lignified-xylem-22","T2-Lignified-xylem-23",
                  "T2-Lignified-xylem-24","T2-Lignified-xylem-25",
                  "T2-Lignified-xylem-26","T3-Phloem-01","T3-Phloem-02",
                  "T3-Phloem-03","T3-Phloem-04","T3-Phloem-05","T3-Cambium-06",
                  "T3-Cambium-07","T3-Cambium-08","T3-Cambium-09","T3-Cambium-10",
                  "T3-Cambium-11","T3-Cambium-12","T3-Cambium-13","T3-Cambium-14",
                  "T3-Expanding-xylem-15","T3-Expanding-xylem-16",
                  "T3-Expanding-xylem-17","T3-Expanding-xylem-18",
                  "T3-Expanding-xylem-19","T3-Expanding-xylem-20",
                  "T3-Expanding-xylem-21","T3-Lignified-xylem-22",
                  "T3-Lignified-xylem-23","T3-Lignified-xylem-24",
                  "T3-Lignified-xylem-25","T3-Lignified-xylem-26",
                  "T3-Lignified-xylem-27","T3-Lignified-xylem-28",
                  "T4-Phloem-01","T4-Phloem-02","T4-Phloem-03","T4-Phloem-04",
                  "T4-Phloem-05","T4-Cambium-06","T4-Cambium-07","T4-Cambium-08",
                  "T4-Cambium-09","T4-Cambium-10","T4-Cambium-11","T4-Cambium-12",
                  "T4-Expanding-xylem-13","T4-Expanding-xylem-14",
                  "T4-Expanding-xylem-15","T4-Expanding-xylem-16",
                  "T4-Expanding-xylem-17","T4-Expanding-xylem-18",
                  "T4-Expanding-xylem-19","T4-Expanding-xylem-20",
                  "T4-Lignified-xylem-21","T4-Lignified-xylem-22",
                  "T4-Lignified-xylem-23","T4-Lignified-xylem-24",
                  "T4-Lignified-xylem-25","T4-Lignified-xylem-26",
                  "T4-Lignified-xylem-27","T4-Lignified-xylem-28","gene_id")

aspwoodtpm <- aspwoodtpm[orderaspwood]
rownames(aspwoodtpm) <- aspwoodtpm$gene_id

# For tree4: aspData <- aspwoodtpm %>% select(80:108)
# For tree1:
aspData <- aspwoodtpm %>% select(1:25,108)
colnames(aspData)
aspData$gene_id <- NULL

# For aspwood with Tung
exprAll <- Reduce(merge.rownames, list(aspData,TungAvg))
exprAll <- Reduce(merge.rownames, list(exprAll,rc))
exprAll_stat <- suppressWarnings(sapply(32:ncol(exprAll), function(i) sapply(1:31, function(j) cor.test(exprAll[,i],exprAll[,j],method = "pearson")[c(3,4)])))
exprAll_label=c("P1", "P2", "P3", "P4", "P5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "Exp13", "Exp14", "Exp15", "Exp16", "Exp17","Exp18", "Exp19", "Exp20", "Lig21", "Lig22","Lig23", "Lig24","Lig25", "Lig26","Lig27", "Lig28", "TungFib", "TungRay", "TungVes")
exprAll_cor <- exprAll_stat[seq(2,nrow(exprAll_stat),2),]
exprAll_pvalue <- exprAll_stat[seq(1,nrow(exprAll_stat)-1,2),]
exprAll_max <- sapply(1:(ncol(exprAll)-31), function(i) max(as.numeric(exprAll_cor[,i])))
exprAll_ident <- sapply(1:(ncol(exprAll)-31), function(i) exprAll_label[which(as.numeric(exprAll_cor[,i])==max(as.numeric(exprAll_cor[,i])))])
exprAll_maxp <- sapply(1:(ncol(exprAll)-31), function(i) as.numeric(exprAll_pvalue[,i])[which(as.numeric(exprAll_cor[,i])==max(as.numeric(exprAll_cor[,i])))])
names(exprAll_max) <- exprAll_ident

integ@meta.data$exprAll.ID.P <- as.character(exprAll_ident)
integ@meta.data$exprAll.cor.P <- exprAll_max
integ@meta.data$exprAll.pvalue.P <- exprAll_maxp
integ@meta.data$exprAll.ID.P[which(integ@meta.data$exprAll.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order=c("P1", "P2", "P3", "P4", "P5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "Exp13", "Exp14", "Exp15", "Exp16", "Exp17","Exp18", "Exp19", "Exp20", "Lig21", "Lig22","Lig23", "Lig24","Lig25", "Lig26","Lig27", "Lig28", "TungFib", "TungRay", "TungVes")
integ$exprAll.ID.P <- factor(integ$exprAll.ID.P, levels = order[sort(match(unique(integ$exprAll.ID.P),order))]) 
palette <- c("darkviolet","orange","yellow","pink","lightblue", "maroon", "green",
             "thistle2", "grey", "violet","burlywood4","blue", "seagreen", "wheat", 
             "lightgreen", "magenta2", "peachpuff", "thistle4","aquamarine",
             "slateblue","purple4","rosybrown","sandybrown","peru","palevioletred",
             "salmon","turquoise","red","olivedrab1","black","skyblue")
color <- palette[sort(match(unique(integ$exprAll.ID.P),order))]
DimPlot(integ, group.by="exprAll.ID.P", cols=color)

# For expanding Xylem
expXAsp <- aspwoodtpm %>% select(1,88:103)
rownames(expXAsp) <- expXAsp$gene_id
expXAsp$gene_id <- NULL
expX <- Reduce(merge.rownames, list(expXAsp,rc))
expX_label=c("Exp13", "Exp14", "Exp15", "Exp16", "Exp17", "Exp18", "Exp19","Exp20")
             #"Lig21", "Lig22", "Lig23", "Lig24", "Lig25", "Lig26", "Lig27", "Lig28"
expX_stat <- suppressWarnings(sapply(9:ncol(expX), function(i) sapply(1:8, function(j) cor.test(expX[,i],expX[,j],method = "pearson")[c(3,4)])))
expX_cor <- expX_stat[seq(2,nrow(expX_stat),2),]
expX_pvalue <- expX_stat[seq(1,nrow(expX_stat)-1,2),]
expX_max <- sapply(1:(ncol(expX)-8), function(i) max(as.numeric(expX_cor[,i])))
expX_ident <- sapply(1:(ncol(expX)-8), function(i) expX_label[which(as.numeric(expX_cor[,i])==max(as.numeric(expX_cor[,i])))])
expX_maxp <- sapply(1:(ncol(expX)-8), function(i) as.numeric(expX_pvalue[,i])[which(as.numeric(expX_cor[,i])==max(as.numeric(expX_cor[,i])))])
names(expX_max) <- expX_ident
integ@meta.data$expX.ID.P <- as.character(expX_ident)
integ@meta.data$expX.cor.P <- expX_max
integ@meta.data$expX.pvalue.P <- expX_maxp
integ@meta.data$expX.ID.P[which(integ@meta.data$expX.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order <- c("Exp13", "Exp14", "Exp15", "Exp16", "Exp17", "Exp18", "Exp19", "Exp20",
           "Lig21", "Lig22", "Lig23", "Lig24", "Lig25", "Lig26", "Lig27", "Lig28","unknown")
palette <- c("lightgreen","lightblue","royalblue","lightpink","grey", "turquoise", "seagreen", "red", "lightgrey")
palette <- c("darkorange","lightgrey","maroon","yellow", "lightpink", "lightgreen",
             "red", "blue", "thistle2","burlywood4","lightgrey", "seagreen", "wheat", 
             "lightgreen", "magenta2", "olivedrab1","black","skyblue")
integ$expX.ID.P <- factor(integ$expX.ID.P, levels = order[sort(match(unique(integ$expX.ID.P),order))]) 
color <- palette[sort(match(unique(integ$expX.ID.P),order))]
DimPlot(integ, group.by="expX.ID.P", cols=color)+NoAxes()

# For lignified Xylem
ligXAsp <- aspwoodtpm %>% select(1,102:103)
rownames(ligXAsp) <- ligXAsp$gene_id
ligXAsp$gene_id <- NULL
ligX <- Reduce(merge.rownames, list(ligXAsp,rc))
ligX_label=c("T27", "T28")
ligX_stat <- suppressWarnings(sapply(3:ncol(ligX), function(i) sapply(1:2, function(j) cor.test(ligX[,i],ligX[,j],method = "pearson")[c(3,4)])))
ligX_cor <- ligX_stat[seq(2,nrow(ligX_stat),2),]
ligX_pvalue <- ligX_stat[seq(1,nrow(ligX_stat)-1,2),]
ligX_max <- sapply(1:(ncol(ligX)-2), function(i) max(as.numeric(ligX_cor[,i])))
ligX_ident <- sapply(1:(ncol(ligX)-2), function(i) ligX_label[which(as.numeric(ligX_cor[,i])==max(as.numeric(ligX_cor[,i])))])
ligX_maxp <- sapply(1:(ncol(ligX)-2), function(i) as.numeric(ligX_pvalue[,i])[which(as.numeric(ligX_cor[,i])==max(as.numeric(ligX_cor[,i])))])
names(ligX_max) <- ligX_ident
integ@meta.data$ligX.ID.P <- as.character(ligX_ident)
integ@meta.data$ligX.cor.P <- ligX_max
integ@meta.data$ligX.pvalue.P <- ligX_maxp
integ@meta.data$ligX.ID.P[which(integ@meta.data$ligX.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order <- c("T27", "T28","unknown")
palette <- c("skyblue", "black", "grey")
integ$ligX.ID.P <- factor(integ$ligX.ID.P, levels = order[sort(match(unique(integ$ligX.ID.P),order))]) 
color <- palette[sort(match(unique(integ$ligX.ID.P),order))]
DimPlot(integ, group.by="ligX.ID.P", cols=color)
VlnPlot(integ, features="ligX.cor.P")+NoLegend()

# for section 27 and 28 separately
ligXAsp <- aspwoodtpm %>% select(1,102)
# ligXAsp <- aspwoodtpm %>% select(1,103)
rownames(ligXAsp) <- ligXAsp$gene_id
ligXAsp$gene_id <- NULL
ligX <- Reduce(merge.rownames, list(ligXAsp,rc))
# ligX_stat <- sapply(2:ncol(ligX), function(i) sapply(1:1, function(j) cor.test(ligX[,i],ligX[,j],method = "pearson")[c(3,4)]))
ligX_stat <- readRDS("~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/correl/ligX27_stat.rds")
ligX_cor <- ligX_stat[seq(2,nrow(ligX_stat),2),]
ligX_pvalue <- ligX_stat[seq(1,nrow(ligX_stat)-1,2),]
ligX_cor_un <- unlist(ligX_cor)
integ@meta.data$ligX.cor.P <- ligX_cor_un
integ@meta.data$ligX.pvalue.P <- ligX_pvalue
VlnPlot(integ, features = "ligX.cor.P")+NoLegend()
FeaturePlot(integ, features = "ligX.cor.P")

# For Phloem
phlmAsp <- aspwoodtpm %>% select (1,105:108)
rownames(phlmAsp) <- phlmAsp$gene_id
phlmAsp$gene_id <- NULL
phlm <- Reduce(merge.rownames, list(phlmAsp,rc))
phlm_stat <- suppressWarnings(sapply(5:ncol(phlm), function(i) sapply(1:4, function(j) cor.test(phlm[,i],phlm[,j],method = "pearson")[c(3,4)])))
phlm_cor <- phlm_stat[seq(2,nrow(phlm_stat),2),]
phlm_pvalue <- phlm_stat[seq(1,nrow(phlm_stat)-1,2),]
phlm_max <- sapply(1:(ncol(phlm)-4), function(i) max(as.numeric(phlm_cor[,i])))
phlm_label=c("P2", "P3", "P4", "P5")
phlm_ident <- sapply(1:(ncol(phlm)-4), function(i) phlm_label[which(as.numeric(phlm_cor[,i])==max(as.numeric(phlm_cor[,i])))])
phlm_maxp <- sapply(1:(ncol(phlm)-4), function(i) as.numeric(phlm_pvalue[,i])[which(as.numeric(phlm_cor[,i])==max(as.numeric(phlm_cor[,i])))])
names(phlm_max) <- phlm_ident
integ@meta.data$phlm.ID.P <- as.character(phlm_ident)
integ@meta.data$phlm.cor.P <- phlm_max
integ@meta.data$phlm.pvalue.P <- phlm_maxp
integ@meta.data$phlm.ID.P[which(integ@meta.data$phlm.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order <- c("P2", "P3", "P4", "P5","unknown")
palette <- c("seagreen","red","blue", "grey", "green")
integ$phlm.ID.P <- factor(integ$phlm.ID.P, levels = order[sort(match(unique(integ$phlm.ID.P),order))]) 
color <- palette[sort(match(unique(integ$phlm.ID.P),order))]
DimPlot(integ, group.by="phlm.ID.P", cols=color)
VlnPlot(integ, features="phlm.cor.P")+NoLegend()
DotPlot(integ, features="phlm.cor.P", cols="RdBu")+NoLegend()

# For Cambium
cambAsp <- aspwoodtpm %>% select(1,81:87)
rownames(cambAsp) <- cambAsp$gene_id
cambAsp$gene_id <- NULL
camb <- Reduce(merge.rownames, list(cambAsp,rc))
camb_stat <- suppressWarnings(sapply(8:ncol(camb), function(i) sapply(1:7, function(j) cor.test(camb[,i],camb[,j],method = "pearson")[c(3,4)])))
camb_cor <- camb_stat[seq(2,nrow(camb_stat),2),]
camb_pvalue <- camb_stat[seq(1,nrow(camb_stat)-1,2),]
camb_max <- sapply(1:(ncol(camb)-7), function(i) max(as.numeric(camb_cor[,i])))
camb_label=c("C6", "C7", "C8", "C9", "C9", "C10", "C11", "C12")
camb_ident <- sapply(1:(ncol(camb)-7), function(i) camb_label[which(as.numeric(camb_cor[,i])==max(as.numeric(camb_cor[,i])))])
camb_maxp <- sapply(1:(ncol(camb)-7), function(i) as.numeric(camb_pvalue[,i])[which(as.numeric(camb_cor[,i])==max(as.numeric(camb_cor[,i])))])
names(camb_max) <- camb_ident
integ@meta.data$camb.ID.P <- as.character(camb_ident)
integ@meta.data$camb.cor.P <- camb_max
integ@meta.data$camb.pvalue.P <- camb_maxp
integ@meta.data$camb.ID.P[which(integ@meta.data$camb.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order <- c("C6", "C7", "C8", "C9", "C9", "C10", "C11", "C12","unknown")
palette <- c("seagreen","orange","blue","pink","skyblue", "green", "red", "grey")
integ$camb.ID.P <- factor(integ$camb.ID.P, levels = order[sort(match(unique(integ$camb.ID.P),order))]) 
color <- palette[sort(match(unique(integ$camb.ID.P),order))]
DimPlot(integ, group.by="camb.ID.P", cols=color)
VlnPlot(integ, features="camb.cor.P")+NoLegend()

# For aspwood 
exprAll <- Reduce(merge.rownames, list(aspData,rc))
exprAll_stat <- suppressWarnings(sapply(26:ncol(exprAll), function(i) sapply(1:25, function(j) cor.test(exprAll[,i],exprAll[,j],method = "pearson")[c(3,4)])))
exprAll_label=c("Phloem 1", "Phloem 2", "Cambium 3", "Cambium 4", "Cambium 5",
                 "Cambium 6", "Cambium 7", "Expanding xylem 8", 
                 "Expanding xylem 9", "Expanding xylem 10", "Expanding xylem 11",
                 "Lignified xylem 12", "Lignified xylem 13", "Lignified xylem 14",
                 "Lignified xylem 15", "Lignified xylem 16", "Lignified xylem 17",
                 "Lignified xylem 18", "Lignified xylem 19", "Lignified xylem 20",
                 "Lignified xylem 21", "Lignified xylem 22","Lignified xylem 23",
                 "Lignified xylem 24","Lignified xylem 25")
exprAll_cor <- exprAll_stat[seq(2,nrow(exprAll_stat),2),]
exprAll_pvalue <- exprAll_stat[seq(1,nrow(exprAll_stat)-1,2),]
exprAll_max <- sapply(1:(ncol(exprAll)-25), function(i) max(as.numeric(exprAll_cor[,i])))
exprAll_ident <- sapply(1:(ncol(exprAll)-25), function(i) exprAll_label[which(as.numeric(exprAll_cor[,i])==max(as.numeric(exprAll_cor[,i])))])
exprAll_maxp <- sapply(1:(ncol(exprAll)-25), function(i) as.numeric(exprAll_pvalue[,i])[which(as.numeric(exprAll_cor[,i])==max(as.numeric(exprAll_cor[,i])))])
names(exprAll_max) <- exprAll_ident

split_seurat$ctrl@meta.data$exprAll.ID.P <- as.character(exprAll_ident)
split_seurat$ctrl@meta.data$exprAll.cor.P <- exprAll_max
split_seurat$ctrl@meta.data$exprAll.pvalue.P <- exprAll_maxp
split_seurat$ctrl@meta.data$exprAll.ID.P[which(split_seurat$ctrl@meta.data$exprAll.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order=c("Phloem 1", "Phloem 2", "Cambium 3", "Cambium 4", "Cambium 5",
        "Cambium 6", "Cambium 7", "Expanding xylem 8", 
        "Expanding xylem 9", "Expanding xylem 10", "Expanding xylem 11",
        "Lignified xylem 12", "Lignified xylem 13", "Lignified xylem 14",
        "Lignified xylem 15", "Lignified xylem 16", "Lignified xylem 17",
        "Lignified xylem 18", "Lignified xylem 19", "Lignified xylem 20",
        "Lignified xylem 21", "Lignified xylem 22","Lignified xylem 23",
        "Lignified xylem 24","Lignified xylem 25")
split_seurat$ctrl$exprAll.ID.P <- factor(split_seurat$ctrl$exprAll.ID.P, levels = order[sort(match(unique(split_seurat$ctrl$exprAll.ID.P),order))]) 
palette <- c("violet","orange","wheat","green","grey","blue","pink", "aquamarine",
             "grey", "turquoise","darkorange","purple", "turquoise","pink", 
             "magenta2","slateblue","olivedrab1", "purple4","yellow","rosybrown",
             "sandybrown","peru","palevioletred", "peru","red")
color <- palette[sort(match(unique(split_seurat$ctrl$exprAll.ID.P),order))]
DimPlot(split_seurat$ctrl, group.by="exprAll.ID.P", cols=color, label = TRUE)+NoAxes()