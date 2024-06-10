library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(here)
library(viridis)
# library(scCustomize)
library(qs)

# load seurat
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")

# In case to plot for one sample, split the object
split_seurat <- SplitObject(integ, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]
# and plot using the object split_seurat$ctrl

#' Plot using different methods
DefaultAssay(integ) <- "RNA"

# In the paper, make figure with control:
DefaultAssay(split_seurat$ctrl) <- "RNA"
rm(integ)

#' Find markers based on literature and plot and discuss if they 
#' could be good markers for a cell type
#' selected cell type markers for figure
selectedmarker <- read.table("~/shruti/SNR-u2023011/analysis/selMarker.txt",sep = '\t',
                             header=TRUE)

#' selected nitrate repsonsive genes
nitRes <- read.table("~/shruti/SNR-u2023011/analysis/nitRes.txt",sep = '\t',header=TRUE)
# for label with own/AT annotation 
nitRes <- nitRes[order(nitRes$AtName),]
# nitRes <- nitRes[order(nitRes$Type),]

# extract the genes first, since gene ids are not unique
nrt <- nitRes[nitRes$Type == "Expansion",]

# extract the atnames
markeranno <- paste(nrt$AtName)
names(markeranno) <- nrt$GeneId

# reorder the idents alphabetically
# Idents(integ) <- factor(integ@active.ident, sort(levels(integ@active.ident)))

# rename idents
integ <- RenameIdents(object = integ,
                      "0" = "Unknown 0","1" = "Fiber Precursor 1",
                      "2" = "Unknown 2","3" = "Unknown 3",
                      "4" = "Early Fiber 4", "5" = "Ray 5",
                      "6" = "Early Vessel 6","7" = "Ray 7",
                      "8" = "Unknown 8", "9" = "Unknown 9",
                      "10" = "Fiber Precursor 19","11" = "Unknown 11",
                      "12" = "Fiber 12","13" = "Unknown 13",
                      "14" = "Fusiform Initial 14", "15" = "Fiber 15",
                      "16" = "Ray/Fusiform Initial 16","17" = "Late Vessel 17",
                      "18" = "Ray 18","19" = "Phloem-like 19",
                      "20" = "Cambium 20")


# order idents
integ@active.ident <- factor(integ@active.ident,
                             levels=c("Phloem like 19","Cambium 20","Ray/ Fusiform Initial 16",
                                      "Ray 5","Ray 7","Ray 18", "Fusiform Initial 14",
                                      "Early Vessel 6","Late Vessel 17",
                                      "Fiber Precursor 10","Fiber Precursor 1",
                                      "Early Fiber 4", "Fiber 12", "Fiber 15", 
                                      "Unknown 0", "Unknown 2","Unknown 3",
                                      "Unknown 8", "Unknown 9","Unknown 11", "Unknown 13"))

DotPlot(subset(integ, subset = sample =="ctrl"), features = nrt$GeneId,dot.scale = 10,
        idents=c("Fiber 2", "Fiber 1", "Early Fiber", "Fiber Precursor 2",
                 "Fiber Precursor 1","Late Vessel", "Early Vessel","Fusiform Initial",
                 "Ray 18","Ray 7","Ray 5","Ray/Fusiform Initial", "Phloem like",
                 "Cambium"))+RotatedAxis()+ scale_x_discrete(labels = markeranno)+
  ylab(NULL)+xlab(NULL)+scale_colour_gradient(low = "grey", high = "darkblue")

DotPlot(subset(integ, subset = sample =="kno"), features = nrt$GeneId,dot.scale = 10,
        idents=c("Fiber 2", "Fiber 1", "Early Fiber", "Fiber Precursor 2",
                 "Fiber Precursor 1","Late Vessel", "Early Vessel","Fusiform Initial",
                 "Ray 18","Ray 7","Ray 5","Ray/Fusiform Initial", "Phloem like",
                 "Cambium"))+ RotatedAxis()+ scale_x_discrete(labels = markeranno)+
  ylab(NULL)+xlab(NULL)+scale_colour_gradient(low = "grey", high = "seagreen")

DotPlot(integ, features = nrt$GeneId,dot.scale = 10,
        idents=c("Fiber 2", "Fiber 1", "Early Fiber", "Fiber Precursor 2",
                 "Fiber Precursor 1","Late Vessel", "Early Vessel","Fusiform Initial",
                 "Ray 18","Ray 7","Ray 5","Ray/Fusiform Initial", "Phloem like",
                 "Cambium"),cols = c("darkblue", "seagreen"), split.by = "sample")+
  RotatedAxis()+scale_x_discrete(labels = markeranno)+ylab(NULL)+xlab(NULL)
  
# plotting function 2
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


# basic plotting commands
DotPlot_scCustom(integ, features = nitRes[nitRes$Type == "NRT", ]$GeneId,
                 colors_use = viridis_plasma_dark_high)+RotatedAxis()+coord_flip()

DotPlot(integ8, features = 
          selectedmarker[selectedmarker$CellType == "chen photosyn", ]$GeneId,
        cols ="RdBu") + RotatedAxis()+ coord_flip()

DotPlot(integ, features = md, cols = c("#005AB5", "#DC3220"),
        dot.scale = 8, split.by = "sample") + RotatedAxis() + coord_flip()+ 
  ylab(NULL) +xlab(NULL)

DotPlot(integ, features = nitRes[nitRes$Type == "NRT", ]$GeneId,
        cols = c("#005AB5", "#DC3220"), dot.scale = 8, split.by = "sample") + RotatedAxis() + coord_flip()+ 
  ylab(NULL) +xlab(NULL)

DotPlot(integ, features = "Potra2n16c29671", cols ="RdBu") +
  RotatedAxis()+ ylab("Cluster Number") +xlab("Genes")+coord_flip()

# Split Violin 
plots <- VlnPlot(integ, features = "Potra2n1c723", cols = c("blue", "red"),
                 split.by = "sample", pt.size = 0, combine = FALSE, 
                 split.plot = TRUE)
wrap_plots(plots = plots, ncol = 1)

# Stacked violin
a <- VlnPlot(integ, features =c("Potra2n1c723","Potra2n16c29671"), stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("nlp7fam")
plot_grid(a)

# png("lac_fam.png")
# for (i in length(lac_fam)) {
#   print(DotPlot(seurat_phase, features = lac_fam)+RotatedAxis())}
# dev.off()
# 

# Potential expansion related
erf<- c("Potra2n6c14920","Potra2n6c15202","Potra2n741s36648","Potra2n7c15468","Potra2n7c15915","Potra2n7c15962","Potra2n7c16238","Potra2n7c16270","Potra2n8c17255","Potra2n8c18165","Potra2n8c18167","Potra2n8c18558","Potra2n8c18593","Potra2n900s36884","Potra2n9c19227, Potra2n10c20778","Potra2n10c20805","Potra2n10c21813","Potra2n10c21814","Potra2n10c21815","Potra2n10c22291","Potra2n11c22537","Potra2n11c22843","Potra2n11c23270","Potra2n11c23301","Potra2n11c23302","Potra2n12c23983","Potra2n12c24633","Potra2n12c24822","Potra2n13c24952","Potra2n13c25101","Potra2n13c25435","Potra2n13c25436","Potra2n13c25445","Potra2n13c25901","Potra2n14c26504","Potra2n14c26681","Potra2n14c26682","Potra2n14c26683","Potra2n14c26684","Potra2n14c26758","Potra2n14c27085","Potra2n14c27132","Potra2n14c27376","Potra2n15c28051","Potra2n15c29002","Potra2n16c29403","Potra2n16c29714","Potra2n16c30086","Potra2n17c31501","Potra2n17c31814","Potra2n18c32534","Potra2n18c32848","Potra2n18c32875","Potra2n18c32911","Potra2n19c33278","Potra2n19c33433","Potra2n19c33580","Potra2n19c33705","Potra2n19c33706","Potra2n19c33749","Potra2n19c33808","Potra2n19c34263","Potra2n1c1345","Potra2n1c1346","Potra2n1c1362","Potra2n1c1377","Potra2n1c1598","Potra2n1c1648","Potra2n1c2725","Potra2n1c3138","Potra2n1c3461","Potra2n1c3816","Potra2n1c541","Potra2n1c558","Potra2n1c643","Potra2n1c644","Potra2n1c646","Potra2n1c766","Potra2n1c788","Potra2n1c934","Potra2n2c4582","Potra2n2c4847","Potra2n2c4848","Potra2n2c4898","Potra2n2c5129","Potra2n2c5283","Potra2n2c5554","Potra2n2c5640","Potra2n2c6054","Potra2n2c6088","Potra2n2c6089","Potra2n2c6090","Potra2n2c6091","Potra2n2c6178","Potra2n3c7035","Potra2n3c7050","Potra2n3c7135","Potra2n3c7136","Potra2n3c7238","Potra2n3c7265","Potra2n3c7404","Potra2n3c7800","Potra2n3c7807","Potra2n3c7820","Potra2n3c7834","Potra2n3c7987","Potra2n3c8021","Potra2n4c8816","Potra2n4c8818","Potra2n4c8857","Potra2n4c8859","Potra2n4c9699","Potra2n5c10756","Potra2n5c10848","Potra2n5c10850","Potra2n5c10852","Potra2n5c10853","Potra2n5c10884","Potra2n5c11289","Potra2n5c11291","Potra2n5c11364","Potra2n5c11580","Potra2n5c12062","Potra2n5c12161","Potra2n6c13204","Potra2n6c13400","Potra2n6c13894","Potra2n6c14084","Potra2n6c14085","Potra2n6c14086")
# feronia
fer <- c("Potra2n16c30521", "Potra2n16c30523", "Potra2n6c14356")
ralf <- c("Potra2n13c26153","Potra2n14c27428","Potra2n17c31433","Potra2n1c2781",
          "Potra2n267s35162","Potra2n4c8790","Potra2n4c8791","Potra2n5c12610")
lrx8 <- c("Potra2n1c2637", "Potra2n6c13964", "Potra2n6c14636", "Potra2n9c19254", "Potra2n14c26497")
md <- c("Potra2n14c27025","Potra2n6c14356","Potra2n8c17700","Potra2n16c30033","Potra2n10c21287")