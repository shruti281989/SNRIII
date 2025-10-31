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

# The big integrated object
# integ <- readRDS("/mnt/picea/home/schoudhary/shruti/SNRIII/data/SeuratOut/integ.rds")
# load("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/SNRIII-FiltDbltRemIntegSplit.RData")

integ <- readRDS("/mnt/picea/home/schoudhary/shruti/SNRIII/data/SeuratOut/integDiffXyT89SNRIII.rds")
# Remove useless object 
integ$integrated_snn_res.0.4 <- NULL
integ$integrated_snn_res.0.8 <- NULL
integ$integrated_snn_res.1 <- NULL
integ$integrated_snn_res.1.4 <- NULL
integ$Phase <- NULL
integ$S.Score <- NULL
integ$G2M.Score <- NULL

split_seurat <- SplitObject(integ, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]

merge.rownames <- function (x,y){
  dat <- merge(x = x, y = y, by = "row.names")
  rownames(dat) <- dat$Row.names
  dat <- dat[,-1]
  return(dat)
}

# 1. with Tung et al., 2022 dataset for LCM cells in aspen
# For control sample
rc <- as.matrix(split_seurat$ctrl@assays$SCT@data)
# rc <- as.matrix(integ@assays$SCT@data)
tungtpm <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/publisheddatasets/Tung/GSE180121_Cell_type_TPM_Ptr.txt", sep = "\t", header = TRUE)
tungtpm$Potri <- sub(".v4.1","",tungtpm$Gene.ID)
potrapotri <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/publisheddatasets/potra_potri_BEST_DIAMOND_out.tsv.gz")
colnames(potrapotri) <- c("Potra","Potri")
tungtpmPotra <- inner_join(tungtpm,potrapotri)
tungtpmPotra$Potri <- NULL
tungtpmPotra$Gene.ID <- NULL
write.table(tungtpmPotra, file = "tungtpmPotraLcm.txt", 
            append = F, sep = "\t", quote = F,row.names = F, col.names = T)

TungAvg <- tungtpmPotra %>% select(Potra, Fiber.mean.TPM, Ray.mean.TPM, 
                                   Vessel.bio2.TPM) #vessel 2 might be purely vessels
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

split_seurat$ctrl@meta.data$exprTungAll.ID.P <- as.character(exprAll_ident)
split_seurat$ctrl@meta.data$exprTungAll.cor.P <- exprAll_max
split_seurat$ctrl@meta.data$exprTungAll.pvalue.P <- exprAll_maxp
split_seurat$ctrl@meta.data$exprTungAll.ID.P[which(split_seurat$ctrl@meta.data$exprTungAll.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order=c("Fiber", "Ray", "Vessel", "unknown")
split_seurat$ctrl$exprTungAll.ID.P <- factor(split_seurat$ctrl$exprTungAll.ID.P, levels = order[sort(match(unique(split_seurat$ctrl$exprTungAll.ID.P),order))]) 
palette <- c("seagreen","blue","pink","grey")
color <- palette[sort(match(unique(split_seurat$ctrl$exprTungAll.ID.P),order))]
DimPlot(split_seurat$ctrl, group.by="exprTungAll.ID.P", cols=color)+NoAxes()
saveRDS(split_seurat, "split_seurat.rds")

# integ@meta.data$exprTungAll.ID.P <- as.character(exprAll_ident)
# integ@meta.data$exprTungAll.cor.P <- exprAll_max
# integ@meta.data$exprTungAll.pvalue.P <- exprAll_maxp
# integ@meta.data$exprTungAll.ID.P[which(integ@meta.data$exprTungAll.ID.P=='character(0)')]="unknown"
# options(repr.plot.width=30, repr.plot.height=24)
# order=c("Fiber", "Ray", "Vessel", "unknown")
# integ$exprTungAll.ID.P <- factor(integ$exprTungAll.ID.P, levels = order[sort(match(unique(integ$exprTungAll.ID.P),order))]) 
# palette <- c("seagreen","blue","pink","grey")
# color <- palette[sort(match(unique(integ$exprTungAll.ID.P),order))]
# DimPlot(integ, group.by="exprTungAll.ID.P", cols=color)+NoAxes()
# saveRDS(integ, "integCorrel.rds")

# Clear useless objects

# For nitrate sample:
rn <- as.matrix(split_seurat$kno@assays$SCT@data)
exprAll <- Reduce(merge.rownames, list(TungAvg,rn))
exprAll_stat <- suppressWarnings(sapply(4:ncol(exprAll), function(i) sapply(1:3, function(j) cor.test(exprAll[,i],exprAll[,j],method = "pearson")[c(3,4)])))
exprAll_label=c("Fiber", "Ray", "Vessel")
exprAll_cor <- exprAll_stat[seq(2,nrow(exprAll_stat),2),]
exprAll_pvalue <- exprAll_stat[seq(1,nrow(exprAll_stat)-1,2),]
exprAll_max <- sapply(1:(ncol(exprAll)-3), function(i) max(as.numeric(exprAll_cor[,i])))
exprAll_ident <- sapply(1:(ncol(exprAll)-3), function(i) exprAll_label[which(as.numeric(exprAll_cor[,i])==max(as.numeric(exprAll_cor[,i])))])
exprAll_maxp <- sapply(1:(ncol(exprAll)-3), function(i) as.numeric(exprAll_pvalue[,i])[which(as.numeric(exprAll_cor[,i])==max(as.numeric(exprAll_cor[,i])))])
names(exprAll_max) <- exprAll_ident

split_seurat$kno@meta.data$exprTungAll.ID.P <- as.character(exprAll_ident)
split_seurat$kno@meta.data$exprTungAll.cor.P <- exprAll_max
split_seurat$kno@meta.data$exprTungAll.pvalue.P <- exprAll_maxp
split_seurat$kno@meta.data$exprTungAll.ID.P[which(split_seurat$kno@meta.data$exprTungAll.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order=c("Fiber", "Ray", "Vessel", "unknown")
split_seurat$kno$exprTungAll.ID.P <- factor(split_seurat$kno$exprTungAll.ID.P, levels = order[sort(match(unique(split_seurat$kno$exprTungAll.ID.P),order))]) 
palette <- c("seagreen","blue","pink","grey")
color <- palette[sort(match(unique(split_seurat$kno$exprTungAll.ID.P),order))]
DimPlot(split_seurat$kno, group.by="exprTungAll.ID.P", cols=color)+NoAxes()

saveRDS(split_seurat, "split_seurat.rds")

# 2. with Aspwood
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

# For correlation of conntrl sample with aspwood
exprAll <- Reduce(merge.rownames, list(aspData,rc))
saveRDS(exprAll, "exprAllAspwdinteg.rds")

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

split_seurat$ctrl@meta.data$tree1.ID.P <- as.character(exprAll_ident)
split_seurat$ctrl@meta.data$tree1.cor.P <- exprAll_max
split_seurat$ctrl@meta.data$tree1.pvalue.P <- exprAll_maxp
split_seurat$ctrl@meta.data$tree1.ID.P[which(split_seurat$ctrl@meta.data$tree1.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order=c("Phloem 1", "Phloem 2", "Cambium 3", "Cambium 4", "Cambium 5",
        "Cambium 6", "Cambium 7", "Expanding xylem 8", 
        "Expanding xylem 9", "Expanding xylem 10", "Expanding xylem 11",
        "Lignified xylem 12", "Lignified xylem 13", "Lignified xylem 14",
        "Lignified xylem 15", "Lignified xylem 16", "Lignified xylem 17",
        "Lignified xylem 18", "Lignified xylem 19", "Lignified xylem 20",
        "Lignified xylem 21", "Lignified xylem 22","Lignified xylem 23",
        "Lignified xylem 24","Lignified xylem 25")
split_seurat$ctrl$tree1.ID.P <- factor(split_seurat$ctrl$tree1.ID.P, levels = order[sort(match(unique(split_seurat$ctrl$tree1.ID.P),order))]) 
palette <- c("violet","orange","wheat","green","grey","blue","pink", "aquamarine",
             "grey", "turquoise","darkorange","purple", "turquoise","pink", 
             "magenta2","slateblue","olivedrab1", "purple4","yellow","rosybrown",
             "sandybrown","peru","palevioletred", "peru","red")
color <- palette[sort(match(unique(split_seurat$ctrl$tree1.ID.P),order))]
DimPlot(split_seurat$ctrl, group.by="tree1.ID.P", cols=color, label = TRUE)+NoAxes()
saveRDS(split_seurat, "split_seurat.rds")

# integ@meta.data$tree1.ID.P <- as.character(exprAll_ident)
# integ@meta.data$tree1.cor.P <- exprAll_max
# integ@meta.data$tree1.pvalue.P <- exprAll_maxp
# integ@meta.data$tree1.ID.P[which(integ@meta.data$tree1.ID.P=='character(0)')]="unknown"
# options(repr.plot.width=30, repr.plot.height=24)
# order=c("Phloem 1", "Phloem 2", "Cambium 3", "Cambium 4", "Cambium 5",
#         "Cambium 6", "Cambium 7", "Expanding xylem 8", 
#         "Expanding xylem 9", "Expanding xylem 10", "Expanding xylem 11",
#         "Lignified xylem 12", "Lignified xylem 13", "Lignified xylem 14",
#         "Lignified xylem 15", "Lignified xylem 16", "Lignified xylem 17",
#         "Lignified xylem 18", "Lignified xylem 19", "Lignified xylem 20",
#         "Lignified xylem 21", "Lignified xylem 22","Lignified xylem 23",
#         "Lignified xylem 24","Lignified xylem 25")
# integ$tree1.ID.P <- factor(integ$tree1.ID.P, levels = order[sort(match(unique(integ$tree1.ID.P),order))]) 
# palette <- c("violet","orange","wheat","green","grey","blue","pink", "aquamarine",
#              "grey", "turquoise","darkorange","purple", "turquoise","pink", 
#              "magenta2","slateblue","olivedrab1", "purple4","yellow","rosybrown",
#              "sandybrown","peru","palevioletred", "peru","red")
# color <- palette[sort(match(unique(integ$tree1.ID.P),order))]
# DimPlot(integ, group.by="tree1.ID.P", cols=color, label = TRUE)+NoAxes()
# saveRDS(integ, "integ.rds")
# For nitrate
exprAll <- Reduce(merge.rownames, list(aspData,rn))
# saveRDS(exprAll, "exprAllAspwdkno.rds")
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

split_seurat$kno@meta.data$tree1.ID.P <- as.character(exprAll_ident)
split_seurat$kno@meta.data$tree1.cor.P <- exprAll_max
split_seurat$kno@meta.data$tree1.pvalue.P <- exprAll_maxp
split_seurat$kno@meta.data$tree1.ID.P[which(split_seurat$kno@meta.data$tree1.ID.P=='character(0)')]="unknown"
options(repr.plot.width=30, repr.plot.height=24)
order=c("Phloem 1", "Phloem 2", "Cambium 3", "Cambium 4", "Cambium 5",
        "Cambium 6", "Cambium 7", "Expanding xylem 8", 
        "Expanding xylem 9", "Expanding xylem 10", "Expanding xylem 11",
        "Lignified xylem 12", "Lignified xylem 13", "Lignified xylem 14",
        "Lignified xylem 15", "Lignified xylem 16", "Lignified xylem 17",
        "Lignified xylem 18", "Lignified xylem 19", "Lignified xylem 20",
        "Lignified xylem 21", "Lignified xylem 22","Lignified xylem 23",
        "Lignified xylem 24","Lignified xylem 25")
split_seurat$kno$tree1.ID.P <- factor(split_seurat$kno$tree1.ID.P, levels = order[sort(match(unique(split_seurat$kno$tree1.ID.P),order))]) 
palette <- c("violet","orange","wheat","green","grey","blue","pink", "aquamarine",
             "grey", "turquoise","darkorange","purple", "turquoise","pink", 
             "magenta2","slateblue","olivedrab1", "purple4","yellow","rosybrown",
             "sandybrown","peru","palevioletred", "peru","red")
color <- palette[sort(match(unique(split_seurat$kno$tree1.ID.P),order))]
DimPlot(split_seurat$kno, group.by="tree1.ID.P", cols=color, label = TRUE)+NoAxes()
saveRDS(split_seurat, "split_seurat.rds")