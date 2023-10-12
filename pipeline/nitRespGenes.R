library(Seurat)
library(dplyr)
library(ggplot2)
library(here)

integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")

# In case to plot for one sample, split the object
split_seurat <- SplitObject(seurat_integrated, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]
# and plot using the object split_seurat$ctrl

#' Plot using different methods
#' DefaultAssay(GnFilt0.8Integ) <- "integrated"
# integ0.8 <- subset(GnFilt0.8Integ, subset = seurat_clusters == '21',invert=T)
DefaultAssay(integ) <- "RNA"

DotPlot(integ8Flt, features = 
          selectedmarker[selectedmarker$CellType == "chen photosyn", ]$GeneId,
        cols ="RdBu") + RotatedAxis()+ coord_flip()

DotPlot(integ, features = nitRes[nitRes$Type == "LGO", ]$GeneId,
        cols ="RdBu") + RotatedAxis()+ coord_flip()

DefaultAssay(integ) <- "RNA"
DimPlot(integ8Flt, reduction = "umap", label = TRUE, label.size = 2)

DotPlot(integ, features = resTFVarala, cols = c("#005AB5", "#DC3220"),
        dot.scale = 8, split.by = "sample") + RotatedAxis() + coord_flip()+ 
  ylab(NULL) +xlab(NULL)

DotPlot(integ, features = nitRes[nitRes$Type == "005min Varala; up kno SNRIII", ]$GeneId,
        cols = c("#005AB5", "#DC3220"), dot.scale = 8, split.by = "sample") + RotatedAxis() + coord_flip()+ 
  ylab(NULL) +xlab(NULL)

DotPlot(integ, features = "Potra2n16c29671", cols ="RdBu") +
  RotatedAxis()+ ylab("Cluster Number") +xlab("Genes")+coord_flip()

# Split Violin 
plots <- VlnPlot(integ0.8, features = "Potra2n1c723", cols = c("blue", "red"),
                 split.by = "sample", pt.size = 0, combine = FALSE, 
                 split.plot = TRUE)
wrap_plots(plots = plots, ncol = 1)

# DoHeatmap(subset(seurat_phase, downsample = 100),
#           features = ces_synthase, size = 3)

# Stacked violin
a <- VlnPlot(integ, features =nlp7fam, stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("nlp7fam")
plot_grid(a)

# png("lac_fam.png")
# for (i in length(lac_fam)) {
#   print(DotPlot(seurat_phase, features = lac_fam)+RotatedAxis())}
# dev.off()
# 
#' Find markers based on literature and plot and discuss if they 
#' could be good markers for a cell type
#' 
#' selected cell type markers for figure
selectedmarker <- read.table("data/SeuratOut/selMarker.txt",sep = '\t',
                             header=TRUE)

generateDotPlot <- function(data, markerfile, cellType) {
  markers <- markerfile[markerfile$CellType == cellType, ]$GeneId
  
  # png(file.path(here("data/SeuratOut/plots/withClCyc/"),
  png(file.path(here("data/SeuratOut/plots/mtCp/"),
                paste0(cellType,".png")), res= 250,height = 4000, width = 3000)
  
  p <- DotPlot(data, features = markers, cols ="RdBu") +
    RotatedAxis()+  ylab(NULL) +xlab(NULL)+coord_flip() +
    labs(title = cellType)
  
  print(p)
  dev.off()
}

unique(selectedmarker$CellType)

#' 
#' selected nitrate repsonsive genes
nitRes <- read.table("~/shruti/SNRIII/doc/NitInt.txt",sep = '\t',header=TRUE)

# "005min Varala" "005min Varala; up kno SNRIII" "010min Varala" 
# "010min Varala; up kno SNRIII" "015min Varala" "015min Varala; up kno SNRIII"           
# [7] "020min Varala" "020min Varala; up kno SNRIII" "030min Varala"
# [10] "030min Varala; up kno SNRIII" "045min Varala"
# "045min Varala; up kno SNRIII" "Cooke" "ACC synthase Li review"
# [13] "060min Varala" "060min Varala; up kno SNRIII" "090min Varala"
# [16] "090min Varala; up kno SNRIII" "120min Varala"
# "120min Varala; down kno SNRIII" "120min Varala; up kno SNRIII" 
# [22] "DEG Lng Trm" "DEG Lng Trm; down kno SNRIII" "DEG Lng Trm; up kno SNRIII"             
# [25] "down at high N Lng Trm" "down at high N Lng Trm; down kno SNRIII"
# [27] "down at high N Lng Trm; up kno SNRIII" "down kno SNRIII"                        
# [29] "Early nitrate responsive TF" "ERF" "Expansin Potra Li review"               
# [32] "k Trnsprt Li review" "LGO" "Nitrate inducible genes" "NLP" "NRT"
# "Other signifcant ERFs" "PME Li review" "Reductase" "Responsive TFs in Varala"               
# [41] "up at high N Lng Trm" "up at high N Lng Trm; down kno SNRIII"  
# [43] "up at high N Lng Trm; up kno SNRIII" "up kno SNRIII" "Varala"
# [46] "Varala TF target" "Varala; down kno SNRIII" "Varala; up kno SNRIII"
# [49] "XET Li review" "Nit Related upregGO in SNRIII"

# lgo interactors in plantgenie
# krp2, RNAPolII, SNZ,CDF1, PHL1, LBD38/37/39, NAC4, HRS1, HHO5/6)
krp2ni <- c("Potra2n7c16280","Potra2n5c11608","Potra2n4c9630","Potra2n16c30563",
            "Potra2n10c21425","Potra2n4c9273","Potra2n4c9276", "Potra2n4c9277",
            "Potra2n8c17012","Potra2n10c20485","Potra2n17c31186","Potra2n4c9546",
            "Potra2n1c535","Potra2n10c21811","Potra2n8c17705","Potra2n14c26447",
            "Potra2n15c28472","Potra2n9c19330","Potra2n7c16165","Potra2n5c11527",
            "Potra2n8c17705","Potra2n5c11630","Potra2n9c19181")