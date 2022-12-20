# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/
# seurat_01_qc.html#Predict_doublets

suppressMessages(require(DoubletFinder))

# before doublet detection run scaling, variable gene selection and pca, 
# UMAP for visualization

data.filt = FindVariableFeatures(data.filt, verbose = F)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA"),
                      verbose = F)
data.filt = RunPCA(data.filt, verbose = F, npcs = 50)
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)

# Then run doubletFinder, selecting first 10 PCs and a pK value of 0.9. 
# To optimize the parameters, run paramSweep function first

sweep.res <- paramSweep_v3(data.filt) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# define the expected number of doublet cells
nExp <- round(ncol(data.filt) * 0.076)  # expect 4% doublets
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, 
                              PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", 
                                              colnames(data.filt@meta.data))]


cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + 
                     NoAxes(),
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())

# two cells can have more detected genes than a single cell, so check if the
# predicted doublets also have more detected genes in general.
VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

# Now, remove all predicted doublets from data
data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
dim(data.filt)

saveRDS(data.filt, "seurat_filterd.rds")