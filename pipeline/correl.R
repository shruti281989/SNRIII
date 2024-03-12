#Comparison of kcl0 (control sample from previous experiment) with publicly available datasets
#Here are pearson pairwise correlations between cluster-based pseudobulk of kcl0 and public bulk data

set.seed(42)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(tidyverse)
  library(pheatmap)
  library(DESeq2)
  library(tximport)
  library(viridis)
})

options(bitmapType = "cairo")

# load data Without mtcp genome and without cell cycle
integ<- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
split_seurat <- SplitObject(integ, split.by = "sample")
kcl0 <- split_seurat[["ctrl"]]
rm(split_seurat, integ)

# With mtcp genome and with cell cycle
integ8 <- readRDS("~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/mtCp/integ8.rds")
Idents(integ8) <- integ8$integrated_snn_res.0.6

# decided to remove cluster 1 from kcl
integ8 <- subset(integ8, subset = seurat_clusters == '1', invert=T)

# Split the data you loaded and extract only the control
split_seurat <- SplitObject(integ8, split.by = "sample")
kcl0 <- split_seurat[["kclmtcp"]]
# kno0 <- split_seurat[["knomtcp"]]
# kcl0 <- split_seurat[["ctrl"]]
rm(split_seurat)

#Setting default assay back to RNA
DefaultAssay(kcl0) <- "RNA"
kcl0 <- NormalizeData(kcl0)
kcl0 <- FindVariableFeatures(kcl0, selection.method = "vst", nfeatures = 2000)
kcl0 <- ScaleData(kcl0)

#Renaming and adding metadata
kcl0$oldclusters <- kcl0$seurat_clusters
kcl0$seurat_clusters <- kcl0@active.ident
kcl0$orig.ident <- paste0("kcl0_",kcl0$seurat_clusters)
kcl0$orig.samp <- "kcl0"

#Check umap once again
# DimPlot(kcl0, reduction = "umap", pt.size = 0.01, label = TRUE)
#FeaturePlot(kcl0, features = c("nCount_RNA","nFeature_RNA"))

#pseudobulk
kcl0bulkscale <- AggregateExpression(kcl0, group.by = "orig.ident") #slot = "data"
#kcl0bulkcount <- AggregateExpression(kcl0, group.by = "orig.ident", slot = "counts")
# 'scale' is normalized count which I'll use. If raw count is needed in any case, one can add slot = "counts"

# saveRDS(kcl0bulkscale, file = "~/shruti/SNR-u2023011/analysis/kcl0bulkscale.rds")
# saveRDS(kcl0bulkscale, file = "~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/kclbulkMtCp5final.rds")
# saveRDS(kcl0bulkscale, file = "~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNucl/kclbulkClCyc.rds")
# saveRDS(kcl0bulkscale, file = "~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNucl/kclbulkClCycFinal.rds")
# saveRDS(kcl0bulkscale, file = "~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNucl/kclbulkFinal.rds")
# saveRDS(kcl0bulkscale, file = "~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/kclbulkMtCp.rds")
# saveRDS(kcl0bulkscale, file = "~/shruti/SNR-u2023011/analysis/snrIII/clRngrCntNuclMtCp/kclbulkMtCpFinal.rds")

# Clear workspace and restart R
##1. Aspwood
#Load AspWood TPM
#Ref: https://plantgenie.org/FTP?dir=Data%2FPlantGenIE%2FPopulus_tremula%2Fv2.2%2FExpression
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

#Join tables and correlation
pseudobulk <- as.data.frame(kcl0bulkscale$RNA)
pseudobulk$gene_id <- rownames(pseudobulk)
exp <- inner_join(aspwoodtpm,pseudobulk)
rownames(exp) <- exp$gene_id
exp$gene_id <- NULL
pcor <- cor(exp, method = "pearson")
pheatmap(pcor, 
         clustering_method = "ward.D2",
         fontsize = 8) #Heatmap with clustering
pheatmap(pcor,
         fontsize = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE) #Heatmap with aspwood order preserved

#Checking pseudobulk in other objects
#SCT
pseudobulk <- as.data.frame(kcl0bulkscale$SCT)
pseudobulk$gene_id <- rownames(pseudobulk)
exp <- inner_join(aspwoodtpm,pseudobulk)
rownames(exp) <- exp$gene_id
exp$gene_id <- NULL
pcor <- cor(exp, method = "pearson")
pheatmap(pcor,
         fontsize = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
#similar to lognormalized

##2. Shi
#Load Shi's bulk LCM dataset
#Ref: https://link.springer.com/article/10.1007/s00425-016-2640-1
shilcm <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/Tung/GSE81077_Tissue_LCM_rawcount.txt", header = TRUE)
#Normalize
rownames(shilcm) <- shilcm$gene
shilcm$gene <- NULL
shilcmmeta <- data.frame(colnames(shilcm),sub(".[123]","",colnames(shilcm)))
colnames(shilcmmeta) <- c("sampleID","celltype")
dds <- DESeqDataSetFromMatrix(countData = shilcm, colData = shilcmmeta, design = ~ celltype)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
boxplot(split(sizes,dds$celltype),las=2,
        main="Sequencing libraries size factor by cell type")
shilcmnorm <- as.data.frame(counts(dds, normalized=TRUE))

#Load ortholog list
potrapotri <- read.table(gzfile("~/shruti/SNR-u2023011/analysis/publisheddatasets/potra_potri_BEST_DIAMOND_out.tsv.gz"))
colnames(potrapotri) <- c("Potra","Potri")
#length(Reduce(intersect,list(potrapotri$Potri,rownames(shilcmnorm)))) #26577

shilcmnorm$Potri <- rownames(shilcmnorm)
shilcmnormPotra <- inner_join(shilcmnorm,potrapotri)
shilcmnormPotra$Potri <- NULL
pseudobulk <- as.data.frame(kcl0bulkscale$RNA)
pseudobulk$Potra <- rownames(pseudobulk)
exp <- inner_join(shilcmnormPotra,pseudobulk)
rownames(exp) <- exp$Potra
exp$Potra <- NULL
nrow(exp) #25772
pcor <- cor(exp, method = "pearson")
pheatmap(pcor,
         fontsize = 8)

##3. Tung
#Tung bulk data
#Ref: https://link.springer.com/article/10.1186/s13059-022-02845-1 
tungtpm <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/Tung/GSE180121_Cell_type_TPM_Ptr.txt", sep = "\t", header = TRUE)
tungtpm$Potri <- sub(".v4.1","",tungtpm$Gene.ID)

#Load potra to potri diamond hits
potrapotri <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/potra_potri_BEST_DIAMOND_out.tsv.gz")
colnames(potrapotri) <- c("Potra","Potri")

tungtpmPotra <- inner_join(tungtpm,potrapotri)
tungtpmPotra$Potri <- NULL
tungtpmPotra$Gene.ID <- NULL
pseudobulk <- as.data.frame(kcl0bulkscale$RNA)
pseudobulk$Potra <- rownames(pseudobulk)
exp <- inner_join(tungtpmPotra,pseudobulk)
rownames(exp) <- exp$Potra
exp$Potra <- NULL
nrow(exp) #24578
pcor <- cor(exp, method = "pearson")
pheatmap(pcor,
         fontsize = 8)

#Compare aspwood with Potri datasets
aspwoodtpm$Potra <- aspwoodtpm$gene_id

#Reorder the kcl clusters
#And.. combine them all!!
# kclorder <- c("kcl0_8")
# pseudobulk <- as.data.frame(kcl0bulkscale$RNA)
# pseudobulk <- pseudobulk[kclorder]
# pseudobulk$gene_id <- rownames(pseudobulk)

#4. Load Du's clusterwise enriched genes from: https://pgx.zju.edu.cn/stRNAPal/
duSTmain <- read.csv("../SNR-u2023011/analysis/publisheddatasets/Du_main_rnamean.csv", header = TRUE)
duSTmain <- subset(duSTmain, select = -geneID)

# in the Du's subclusters: C1.0: Procambium-a	; C1.1: Protoxylem ; C1.2: Protophloem-a;
# C1.3: Procambium-b;	C1.4: Metacambium;	C1.5: Protophloem-b;	C11.0: Cambium zone
# C11.1: Differentiating xylem;	C14.0: Procambium-like;	C14.1: Differentiating phloem
duSTsub<- read.csv("../SNR-u2023011/analysis/publisheddatasets/Du_sub_rnamean.csv", header = TRUE)
duSTsub <- subset(duSTsub, select = -c(geneID,description))

#Load potra to potri diamond hits
potrapotri <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/potra_potri_BEST_DIAMOND_out.tsv.gz")
colnames(potrapotri) <- c("Potra","Potri")

duSTmain$Potri <- sub(".v4.1","",duSTmain$gene)
duSTmainPotra <- inner_join(duSTmain,potrapotri)
duSTmainPotra$Potri <- NULL
duSTmainPotra$gene <- NULL

exp <- inner_join(duSTmainPotra,pseudobulk)
rownames(exp) <- exp$Potra
exp$Potra <- NULL
nrow(exp)
pcor <- cor(exp, method = "pearson")
pheatmap(pcor,
         fontsize = 8)

duSTsub$Potri <- sub(".v4.1","",duSTsub$gene)
duSTsubPotra <- inner_join(duSTsub,potrapotri)
duSTsubPotra$Potri <- NULL
duSTsubPotra$gene <- NULL

exp <- inner_join(duSTsubPotra,pseudobulk)
rownames(exp) <- exp$Potra
exp$Potra <- NULL
nrow(exp)
pcor <- cor(exp, method = "pearson")
pheatmap(pcor,
         fontsize = 8)

# make one big correl heatmap
# aspwood
aspwood <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/AspWood_tpm.txt", header = TRUE)
aspwoodtpm <- dcast(aspwood, gene_id ~ sample_name)

##Potri datasets
potrapotri <- read.table(gzfile("~/shruti/SNR-u2023011/analysis/publisheddatasets/potra_potri_BEST_DIAMOND_out.tsv.gz"))
colnames(potrapotri) <- c("Potra","Potri")

#Tung bulk data
tungtpm <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/Tung/GSE180121_Cell_type_TPM_Ptr.txt", sep = "\t", header = TRUE)
tungtpm$Potri <- sub(".v4.1","",tungtpm$Gene.ID)
tungtpmPotra <- inner_join(tungtpm,potrapotri)
tungtpmPotra$Potri <- NULL
tungtpmPotra$Gene.ID <- NULL

#Load Du's clusterwise enriched genes from:
duSTmain <- read.csv("~/shruti/SNR-u2023011/analysis/publisheddatasets/Du_main_rnamean.csv", header = TRUE)
duSTmain <- subset(duSTmain, select = -geneID)
duSTmain$Potri <- sub(".v4.1","",duSTmain$gene)
duSTmainPotra <- inner_join(duSTmain,potrapotri)
duSTmainPotra$Potri <- NULL
duSTmainPotra$gene <- NULL

duSTsub<- read.csv("~/shruti/SNR-u2023011/analysis/publisheddatasets/Du_sub_rnamean.csv", header = TRUE)
duSTsub <- subset(duSTsub, select = -c(geneID,description))

duSTsub$Potri <- sub(".v4.1","",duSTsub$gene)
duSTsubPotra <- inner_join(duSTsub,potrapotri)
duSTsubPotra$Potri <- NULL
duSTsubPotra$gene <- NULL
# C1.0"  "C1.1"  "C1.2"  "C1.3"  "C1.4"  "C1.5"  "C11.0", "C11.1" "C14.0" "C14.1"
# "C1.0 Procambium-a", "C1.1 Primary Xylem", "C1.2 Primary Phloem-a", "C1.3 Procambium-b",
# "C1.4 Metacambium", "C1.5 Primary Phloem-b", "C11.0 Cambium Zone", "C11.1 Differentiating Xylem"
# "C14.0 Procamium-like", "C14.1 Differentiating Phloem" 

colnames(duSTsubPotra) <- c("C1.0 Procambium apex", "C1.1 Primary Xylem apex", 
                       "C1.2 Primary Phloem apex", "C1.3 Procambium-b",
                       "C1.4 Metacambium", "C1.5 Primary Phloem-b",
                       "C11.0 Cambium Zone", "C11.1 Differentiating Xylem",
                       "C14.0 Procamium-like", "C14.1 Differentiating Phloem",
                       "Potra")

#Shi's bulk LCM dataset
shilcm <- read.table("~/shruti/SNR-u2023011/analysis/publisheddatasets/Tung/GSE81077_Tissue_LCM_rawcount.txt", header = TRUE)
rownames(shilcm) <- shilcm$gene
shilcm$gene <- NULL
shilcmmeta <- data.frame(colnames(shilcm),sub(".[123]","",colnames(shilcm)))
colnames(shilcmmeta) <- c("sampleID","celltype")
dds <- DESeqDataSetFromMatrix(countData = shilcm, colData = shilcmmeta, design = ~ celltype)
dds <- estimateSizeFactors(dds)
shilcmnorm <- as.data.frame(counts(dds, normalized=TRUE))
remove(dds,shilcm,shilcmmeta)
shilcmnorm$Potri <- rownames(shilcmnorm)
shilcmnormPotra <- inner_join(shilcmnorm,potrapotri)
shilcmnormPotra$Potri <- NULL

#Drop unused columns
aspwooddrop <- subset(aspwoodtpm, select = grep("T4-*", colnames(aspwoodtpm)))
aspwooddrop$Potra <- aspwoodtpm$gene_id
tungdrop <- subset(tungtpmPotra, select = -grep("*.mean.TPM", colnames(tungtpmPotra)))
shidrop <- subset(shilcmnormPotra, select = -grep("Three*", colnames(shilcmnormPotra)))
shidrop <- subset(shidrop, select = -grep("Leaf*", colnames(shidrop)))
shidrop <- subset(shidrop, select = -grep("Shoot*", colnames(shidrop)))
shidrop <- subset(shidrop, select = -grep("Root*", colnames(shidrop)))
duSTsubpick <- subset(duSTsubPotra, select = c(4:11))
duSTmainpick <- subset(duSTmainPotra, select = c(8,9,10,12,13,14,15,16,18))

#Join tables and correlation
pseudobulk <- as.data.frame(kcl0bulkscale$RNA)
pseudobulk$Potra <- rownames(pseudobulk)
exp <- inner_join(pseudobulk,aspwooddrop)
exp2 <- inner_join(exp,tungdrop)
exp3 <- inner_join(exp2,shidrop)
exp4 <- inner_join(exp3,duSTsubpick)
exp5 <- inner_join(exp4,duSTmainpick)
exp5$Potra <- NULL
pcor <- cor(exp5, method = "pearson")
#pheatmap(pcor, 
#         fontsize = 7,
#         cluster_rows = FALSE,
#         cluster_cols = FALSE) #Order preserved

#Small heatmap
pcorsubset <- pcor[-grep("kcl*", rownames(pcor)),grep("kcl*", rownames(pcor))]
pheatmap(pcorsubset, 
         fontsize = 7,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = round(pcorsubset, digits = 2))

#Set min clolor range at 0.3/ 0.2
mat_breaks <- seq(0.3, max(pcorsubset), length.out = 100)
# mat_breaks <- seq(0.2, max(pcorsubset), length.out = 100)

pheatmap(pcorsubset, 
         fontsize = 16,
         cluster_rows = FALSE,
         # cluster_cols = FALSE,
         #display_numbers = round(pcorsubset, digits = 2),
         # color = inferno(100), 
         color = magma(100),
         # color = viridis(100), 
         # color = mako(100),
         border_color = NA,
         breaks = mat_breaks,
         gaps_row = c(28,37, 49))
