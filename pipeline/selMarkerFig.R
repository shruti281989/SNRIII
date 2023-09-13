library(Seurat)
library(dplyr)
library(ggplot2)
library(here)

integMtCp <- readRDS("~/shruti/SNRIII/data/SeuratOut/integMtCp.rds")
load("~/shruti/SNRIII/data/SeuratOut/integafterdbltremoval.RData")
load("~/shruti/SNRIII/data/SeuratOut/integWthClCyc.RData")

Idents(object = integMtCp) <- "integrated_snn_res.0.6"
DefaultAssay(integMtCp) <- "RNA"

DefaultAssay(seurat_integrated) <- "RNA"

DefaultAssay(integClCyc) <- "RNA"
DimPlot(integClCyc, reduction = "umap", label = TRUE)

DotPlot(integClCyc, features = selectedmarker[selectedmarker$CellType == "vessel", ]$GeneId, 
        cols ="RdBu")+
  RotatedAxis()+ ylab("Cluster Number") +xlab("Genes")+coord_flip()
DotPlot(seurat_integrated, features = c("Potra2n13c25437","Potra2n13c25438","Potra2n158s34776","Potra2n484s35780"), 
        cols ="RdBu")+RotatedAxis()+ ylab("Cluster Number") +xlab("Genes")+coord_flip()

phloem <- c("Potra2n10c20892","Potra2n10c21280","Potra2n10c21469","Potra2n10c21994","Potra2n11c22484","Potra2n11c22540","Potra2n11c22649","Potra2n11c22821","Potra2n12c23967","Potra2n12c24479","Potra2n12c24713","Potra2n12c24714","Potra2n13c25317","Potra2n13c26179","Potra2n16c29487","Potra2n17c31828","Potra2n18c33107","Potra2n19c33601","Potra2n19c34419","Potra2n1c1686","Potra2n1c2159","Potra2n1c267","Potra2n1c3642","Potra2n1c700","Potra2n2c4235","Potra2n2c5371","Potra2n2c5643","Potra2n2c5644","Potra2n391s35520","Potra2n3c6746","Potra2n4c9160","Potra2n4c9835","Potra2n5c11135","Potra2n5c12223","Potra2n5c12320","Potra2n6c13736","Potra2n6c14623","Potra2n6c15124","Potra2n7c15682","Potra2n7c15685","Potra2n7c15689","Potra2n81s34525","Potra2n81s34526","Potra2n8c17353","Potra2n8c18612","Potra2n9c18825","Potra2n9c19910","Potra2n4c10378","Potra2n16c30506","Potra2n6c14369","Potra2n11c23289","Potra2n13c25854","Potra2n17c31798","Potra2n18c32257","Potra2n18c32258","Potra2n18c32505","Potra2n18c32506","Potra2n18c32791","Potra2n19c33677","Potra2n19c34194","Potra2n1c1091")
cambium <- c("Potra2n10c20235","Potra2n10c20736","Potra2n10c20760","Potra2n10c20846","Potra2n10c21451","Potra2n10c21673","Potra2n10c22031","Potra2n11c22555","Potra2n11c22567","Potra2n11c22811","Potra2n11c22828","Potra2n11c23322","Potra2n12c23826","Potra2n12c23983","Potra2n12c24164","Potra2n12c24650","Potra2n12c24764","Potra2n13c24952","Potra2n13c25185","Potra2n13c25262","Potra2n13c25541","Potra2n13c25822","Potra2n13c25852","Potra2n14c26391","Potra2n14c26506","Potra2n14c26580","Potra2n14c27085","Potra2n14c27215","Potra2n14c27365","Potra2n14c27708","Potra2n15c27987","Potra2n15c28802","Potra2n15c28861","Potra2n15c29002","Potra2n16c29899","Potra2n16c30563","Potra2n16c30576","Potra2n17c30987","Potra2n17c31100","Potra2n17c31267","Potra2n18c32073","Potra2n18c32470","Potra2n18c32483","Potra2n18c32702","Potra2n18c33107","Potra2n19c33278","Potra2n19c33517","Potra2n19c33790","Potra2n1c1635","Potra2n1c17","Potra2n1c2087","Potra2n1c2195","Potra2n1c2732","Potra2n1c2900","Potra2n1c3147","Potra2n1c357","Potra2n1c3842","Potra2n1c605","Potra2n1c788","Potra2n1c998","Potra2n2c4070","Potra2n2c4139","Potra2n2c4218","Potra2n2c4898","Potra2n2c5282","Potra2n2c5745","Potra2n2c5977","Potra2n2c6181","Potra2n2c6375","Potra2n2c6439","Potra2n3c6864","Potra2n3c7095","Potra2n3c7265","Potra2n3c7460","Potra2n3c7528","Potra2n4c10039","Potra2n4c8814","Potra2n4c9561","Potra2n4c9798","Potra2n5c10575","Potra2n5c10753","Potra2n5c10979","Potra2n5c11241","Potra2n5c11446","Potra2n5c11602","Potra2n5c11916","Potra2n5c11923","Potra2n5c12554","Potra2n6c12962","Potra2n6c13323","Potra2n6c13849","Potra2n6c14336","Potra2n6c14711","Potra2n6c15109","Potra2n7c15470","Potra2n7c15818","Potra2n7c15935","Potra2n7c16113","Potra2n7c16264","Potra2n7c16566","Potra2n8c16773","Potra2n8c17318","Potra2n8c17825","Potra2n8c17900","Potra2n8c18011","Potra2n8c18019","Potra2n9c19378","Potra2n9c19384","Potra2n9c19454","Potra2n9c19718","Potra2n9c19735", "Potra2n6c14967","Potra2n9c19651","Potra2n1c2250","Potra2n6c13361","Potra2n6c15120")
ray5_7 <- c("Potra2n10c20444","Potra2n11c22600","Potra2n14c27801","Potra2n14c27834","Potra2n15c28669","Potra2n16c30308","Potra2n18c33011","Potra2n1c3480","Potra2n2c4118","Potra2n2c4509","Potra2n3c6926","Potra2n3c6944","Potra2n3c8005","Potra2n4c8449","Potra2n5c11003","Potra2n6c13067","Potra2n6c15285","Potra2n7c16335","Potra2n9c19477","Potra2n13c25349","Potra2n14c27607","Potra2n16c29243","Potra2n16c30464","Potra2n18c32429","Potra2n18c32564","Potra2n1c1049","Potra2n1c2010","Potra2n1c232","Potra2n1c90","Potra2n2c4144","Potra2n2c4942","Potra2n2c5467","Potra2n2c6010","Potra2n3c8020","Potra2n4c10139","Potra2n4c9969")
ray2<- c("Potra2n13c25215","Potra2n18c32373","Potra2n1c394","Potra2n2c4494","Potra2n2c5697","Potra2n2c5875","Potra2n4c9970","Potra2n5c10874","Potra2n5c11107","Potra2n5c11849","Potra2n8c17115","Potra2n9c18880")
xpc_ray <- c("Potra2n15c28663","Potra2n16c29896","Potra2n1c3834","Potra2n1c749","Potra2n4c9930","Potra2n6c13057","Potra2n6c13538","Potra2n6c15220","Potra2n14c26690","Potra2n14c26986","Potra2n14c27082","Potra2n14c27101","Potra2n16c30018","Potra2n16c30450","Potra2n17c30920","Potra2n18c32844","Potra2n1c1897","Potra2n1c2429","Potra2n1c394","Potra2n2c4928","Potra2n3c6993","Potra2n4c9205","Potra2n5c11091","Potra2n5c11250","Potra2n5c11849","Potra2n6c14448","Potra2n8c17012","Potra2n8c17375","Potra2n9c18721")
fiber <- c("Potra2n11c22673","Potra2n11c23199","Potra2n15c28689","Potra2n16c30073","Potra2n16c30273","Potra2n1c546","Potra2n2c4078","Potra2n4c9297","Potra2n6c14160","Potra2n7c15752","Potra2n7c16228","Potra2n19c33344","Potra2n10c20411","Potra2n5c11573","Potra2n13c24959","Potra2n5c11765","Potra2n7c16495","Potra2n1c3619","Potra2n7c16234","Potra2n5c11571","Potra2n3c7305","Potra2n4c8952","Potra2n5c12319","Potra2n12c24762","Potra2n15c28411","Potra2n1c307","Potra2n1c953","Potra2n2c6235","Potra2n16c30051","Potra2n9c19651","Potra2n1c4016","Potra2n11c22498","Potra2n6c14201","Potra2n8c16778","Potra2n8c16946","Potra2n8c17885","Potra2n9c20072")
celldeath <- c("Potra2n14c26495","Potra2n14c27047","Potra2n1c3263","Potra2n5c11905", "Potra2n4c10301","Potra2n6c15163","Potra2n16c29448","Potra2n10c21346","Potra2n1c3148","Potra2n11c23290","Potra2n1c2692")

selectedmarker <- read.table("data/SeuratOut/selMarker.txt",sep = '\t',header=TRUE)

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
# [1] "aspwood"               "chen photosyn"        
# [3] "etc"                   "wilcox"               
# [5] "gus"                               
# [7] "S"                     "fiber"                
# [9] "cambium"               "ray"                  
# [11] "lignin"                "du cell exp"          
# [13] "chen lit"              "du dna rep cell cyc"  
# [15] "conde"                 "M"                    
# [17] "chen fv"               "du sec wall"          
# [19] "laccase"               "cluster 6"            
# [21] "peroxidase"            "du phloem diff"       
# [23] "phloem"                "G2"                   
# [25] "Tung7"                 "vessel"               
# [27] "du meristem identity"  "du aux"               
# [29] "gutirrez"              "cluster 14"           
# [31] "xylan"                 "g10"                  
# [33] "pectin"                "cluster 5"            
# [35] "cluster 18"            "chen ce"              
# [37] "cluster 11"            "fiber late"           
# [39] "erf"                   "chen sieve"           
# [41] "chen companion"        "chen cork"            
# [43] "cluster 20"            "ligx"                 
# [45] "cluster 7"             "cellulose"            
# [47] "chen cambium"          "Tung2"                
# [49] "kucukoglu"             "du xylem diff"        
# [51] "cluster 9"             "Tung8"                
# [53] "cluster 8"             "Tung6"                
# [55] "chen ce "              "Tung1"                
# [57] "chen phloem mother"    "Tung5"                
# [59] "cluster 2"             "chen xylem mother"    
# [61] "cluster 16"            "not exp"              
# [63] "seyfferth"             "chen e"               
# [65] "Tung4"                 "Tung3"                
# [67] "chen xylem parenchyma"