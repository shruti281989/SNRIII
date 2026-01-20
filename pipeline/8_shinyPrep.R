setwd("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/dnStrm/correl/")

set.seed(42)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

#  Prepare the object for shiny app
integ <- readRDS("/mnt/picea/home/schoudhary/shruti/SNRIII/data/SeuratOut/integDiffXyT89SNRIII.rds")

# Remove useless metadata
integ$integrated_snn_res.0.4 <- NULL
integ$integrated_snn_res.0.8 <- NULL
integ$integrated_snn_res.1 <- NULL
integ$integrated_snn_res.1.4 <- NULL
integ$Phase <- NULL
integ$S.Score <- NULL
integ$G2M.Score <- NULL
integ$seq_folder <- NULL
integ$nCount_SCT <- NULL
integ$nFeature_SCT <- NULL
integ$log10GenesPerUMI <- NULL

# add metadata for public data correlations of each sample to the seurat object
split_seurat <- readRDS("/mnt/ada/projects/aspseq/htuominen/SNR-results/snrIII/dnStrm/correl/split_seurat.rds")
samples <- c("ctrl", "kno")
correlCol <- c("tree1.ID.P", "tree1.cor.P", "exprTungAll.ID.P","exprTungAll.cor.P")
for (s in samples) {
  obj <- split_seurat[[s]]
  inter <- intersect(colnames(integ), colnames(obj))
  meta <- bind_rows(split_seurat[["ctrl"]]@meta.data[, correlCol, drop = FALSE], 
                    split_seurat[["kno"]]@meta.data[, correlCol, drop = FALSE])
  integ <- AddMetaData(integ, meta)
}

# add cell types ID to a different column
celltypeID <- c("0" = "Unknown Cluster 0", "1" = "Fiber precursor Cluster 1",
                "2" = "Unknown Cluster 2", "3" = "Putative Early Ray Cluster 3",
                "4" = "Fiber Cluster 4", "5" = "Ray Cluster 5", "6" = "Vessel precursor Cluster 6",
                "7" = "Ray Cluster 7", "8" = "Putative Early Ray Cluster 8",
                "9" = "Putative Paratracheal Parenchyma Cluster 9", "10" = "Fiber precursor Cluster 10",
                "11" = "Putative Early Ray Cluster 11", "12" = "Fiber Cluster 12",
                "13" = "Putative Early Ray Cluster 13", "14" = "Fusiform initials Cluster 14",
                "15" = "Fiber Cluster 15", "16" = "Fusiform initials Cluster 16",
                "17" = "Vessel Cluster 17", "18" = "Ray Cluster 18",
                "19" = "Cambium Cluster 19", "20" = "Cambium Cluster 20")

integ$CellTypeId <- dplyr::recode(as.character(Idents(integ)),!!!celltypeID)

celltype <- c("0" = "Unknown", "1" = "Fiber precursor", "2" = "Unknown", 
              "3" = "Putative Early Ray", "4" = "Fiber", "5" = "Ray", 
              "6" = "Vessel precursor", "7" = "Ray", "8" = "Putative Early Ray",
              "9" = "Putative Paratracheal Parenchyma", "10" = "Fiber precursor",
              "11" = "Putative Early Ray", "12" = "Fiber", 
              "13" = "Putative Early Ray", "14" = "Fusiform initials",
              "15" = "Fiber", "16" = "Fusiform initials",
              "17" = "Vessel", "18" = "Ray", "19" = "Cambium", "20" = "Cambium")

integ$CellType <- dplyr::recode(as.character(Idents(integ)),!!!celltype)

# change the sample names 
integ$sample <- gsub("^ctrl$", "Control", integ$sample)
integ$sample <- gsub("^kno$", "Nitrate-treated", integ$sample)

# change the names in correlation metadata with Tung LCM cells
integ$exprTungAll.ID.P <- recode(integ$exprTungAll.ID.P, "Fiber" = "Fiber mean",
                                 "Ray" = "Ray mean", "Vessel" = "Vessel replicate 2")

# change name of integrated_snn_res.0.6 metadata column to Clusters
# names(integ@meta.data)[names(integ@meta.data) == "integrated_snn_res.0.6"] <- "Clusters"
integ$Clusters <- integ$integrated_snn_res.0.6
integ$integrated_snn_res.0.6 <- NULL

saveRDS(integ, "integAllLayers.rds")

# remove layers if needed:
integ[["SCT"]] <- NULL
integ[["integrated"]] <- NULL
saveRDS(integ, "integTiny.rds")

# shiny files prep
library(ShinyCell2)

seu <- readRDS("integAllLayers.rds")

scConf <- createConfig(seu)

palette <- c("green","grey","blue","pink", 
             "darkorange","purple","turquoise","magenta2",
             "olivedrab1","yellow","red")

scConf = modColours(scConf, meta.to.mod = "tree1.ID.P", 
                    new.colours= palette)

scConf = modMetaName(scConf, 
                     meta.to.mod = c("exprTungAll.ID.P", "exprTungAll.cor.P"), 
                     new.name = c("Cell types from LCM in Tung et al., 2023", 
                                  "Highest Pearson Correl Value with cell types from LCM in Tung et al., 2023"))

scConf = modMetaName(scConf, 
                     meta.to.mod = c("tree1.ID.P", "tree1.cor.P"), 
                     new.name = c("Tissue section from Aspwood", 
                                  "Highest Pearson Correl Value with tissue section from Aspwood"))

scConf = modColours(scConf, meta.to.mod = "exprTungAll.ID.P", 
                    new.colours= c("seagreen","blue","pink"))

showLegend(scConf)


checkConfig(scConf, seu)

makeShinyFiles(seu, scConf = scConf, dimred.to.use = "umap",
               shiny.prefix = "scPop", shiny.dir = "shinyApp/")

makeShinyCodes(shiny.prefix = "scPop", 
               shiny.dir = "shinyApp/", 
               shiny.title="Populus wood single cell RNA-seq")

shiny::runApp("shinyApp")
