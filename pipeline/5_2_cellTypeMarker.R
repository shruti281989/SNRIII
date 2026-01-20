setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/")
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(viridis)
library(scCustomize)

set.seed(42)

integ <- readRDS("data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"

# To plot for one sample, split the object
split_seurat <- SplitObject(integ, split.by = "sample")
split_seurat <- split_seurat[c("ctrl", "kno")]
# and plot using the object split_seurat$ctrl

DimPlot_scCustom(split_seurat$ctrl, figure_plot = T, pt.size = 0.3, label=F,
                 colors_use = DiscretePalette_scCustomize(num_colors = 30,
                                                          palette = "varibow", 
                                                          shuffle_pal = T))+ 
  NoLegend()

# figure 4A
DimPlot_scCustom(integ, pt.size = 0.3, label=F,
                 colors_use = DiscretePalette_scCustomize(num_colors = 30,
                                                          palette = "varibow", 
                                                          shuffle_pal = T)) 

# figure 4B
fig4 <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/fig4.txt",
                             sep = '\t',header=TRUE)

# for label with Populus/Arabidospsis annotation 
fig4 <- fig4[order(fig4$cluster), ]
# Combine gene and name for labels
markeranno <- paste(fig4$gene, fig4$name, sep = " ")
names(markeranno) <- fig4$gene

DotPlot(subset(integ, subset = sample == "ctrl"), features = fig4$gene,
  dot.scale = 10) + RotatedAxis() + coord_flip()+ xlab(NULL) + ylab(NULL) +
  scale_x_discrete(labels = markeranno) + scale_colour_viridis()

# Load markers
# ploidy <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/endo_expPotra.txt",
#                      sep = '\t',header=TRUE)
# ribosomal <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/Ribo_gene_table.txt",
#                         sep = '\t',header=F)
# protoplasting <- read.table("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/markers/pplast.txt",
#                             sep = '\t',header=F)
#' Plot using different methods
# Cluster wise heatmap: scale the data first: 
# integ <- ScaleData(integ, features = rownames(integ))
# DoHeatmap(integ, features=clustMarker$GeneId, group.by="sample",
#           group.colors = viridis(100))

# DotPlot(subset(integ, subset = sample =="ctrl"), features = clustMarker$GeneId,
#         dot.scale = 10)+RotatedAxis()+ 
#   scale_x_discrete(labels = markeranno)+ylab(NULL)+xlab(NULL)+
#   scale_colour_gradient(low = "grey", high = "darkblue")
# 
# DotPlot(integ, selectedmarker[selectedmarker$CellType == "chen photosyn", ]$GeneId,
#         cols ="RdBu", dot.scale = 10) + RotatedAxis()+ coord_flip()
# # cols = c("#005AB5", "#DC3220"),

# grep( "^mt-", rownames(seurat.object), value = T)
# DoHeatmap(subset(integ, downsample = 100), features = ces_synthase, size = 3)

#' test other markers based on literature for a cell type
#' 
#' 5.1 from literature
#' phloem <- c("Potra2n10c20892","Potra2n10c21280","Potra2n10c21469","Potra2n10c21994",
#'             "Potra2n11c22484","Potra2n11c22540","Potra2n11c22649","Potra2n11c22821",
#'             "Potra2n12c23967","Potra2n12c24479","Potra2n12c24713","Potra2n12c24714",
#'             "Potra2n13c25317","Potra2n13c26179","Potra2n16c29487","Potra2n17c31828",
#'             "Potra2n18c33107","Potra2n19c33601","Potra2n19c34419","Potra2n1c1686",
#'             "Potra2n1c2159","Potra2n1c267","Potra2n1c3642","Potra2n1c700","Potra2n2c4235",
#'             "Potra2n2c5371","Potra2n2c5643","Potra2n2c5644","Potra2n391s35520",
#'             "Potra2n3c6746","Potra2n4c9160","Potra2n4c9835","Potra2n5c11135",
#'             "Potra2n5c12223","Potra2n5c12320","Potra2n6c13736","Potra2n6c14623",
#'             "Potra2n6c15124","Potra2n7c15682","Potra2n7c15685","Potra2n7c15689",
#'             "Potra2n81s34525","Potra2n81s34526","Potra2n8c17353","Potra2n8c18612",
#'             "Potra2n9c18825","Potra2n9c19910","Potra2n4c10378","Potra2n16c30506",
#'             "Potra2n6c14369","Potra2n11c23289","Potra2n13c25854","Potra2n17c31798",
#'             "Potra2n18c32257","Potra2n18c32258","Potra2n18c32505","Potra2n18c32506",
#'             "Potra2n18c32791","Potra2n19c33677","Potra2n19c34194","Potra2n1c1091")
#' cambium <- c("Potra2n10c20235","Potra2n10c20736","Potra2n10c20760","Potra2n10c20846",
#'              "Potra2n10c21451","Potra2n10c21673","Potra2n10c22031","Potra2n11c22555",
#'              "Potra2n11c22567","Potra2n11c22811","Potra2n11c22828","Potra2n11c23322",
#'              "Potra2n12c23826","Potra2n12c23983","Potra2n12c24164","Potra2n12c24650",
#'              "Potra2n12c24764","Potra2n13c24952","Potra2n13c25185","Potra2n13c25262",
#'              "Potra2n13c25541","Potra2n13c25822","Potra2n13c25852","Potra2n14c26391",
#'              "Potra2n14c26506","Potra2n14c26580","Potra2n14c27085","Potra2n14c27215",
#'              "Potra2n14c27365","Potra2n14c27708","Potra2n15c27987","Potra2n15c28802",
#'              "Potra2n15c28861","Potra2n15c29002","Potra2n16c29899","Potra2n16c30563",
#'              "Potra2n16c30576","Potra2n17c30987","Potra2n17c31100","Potra2n17c31267",
#'              "Potra2n18c32073","Potra2n18c32470","Potra2n18c32483","Potra2n18c32702",
#'              "Potra2n18c33107","Potra2n19c33278","Potra2n19c33517","Potra2n19c33790",
#'              "Potra2n1c1635","Potra2n1c17","Potra2n1c2087","Potra2n1c2195","Potra2n1c2732",
#'              "Potra2n1c2900","Potra2n1c3147","Potra2n1c357","Potra2n1c3842","Potra2n1c605",
#'              "Potra2n1c788","Potra2n1c998","Potra2n2c4070","Potra2n2c4139","Potra2n2c4218",
#'              "Potra2n2c4898","Potra2n2c5282","Potra2n2c5745","Potra2n2c5977","Potra2n2c6181",
#'              "Potra2n2c6375","Potra2n2c6439","Potra2n3c6864","Potra2n3c7095","Potra2n3c7265",
#'              "Potra2n3c7460","Potra2n3c7528","Potra2n4c10039","Potra2n4c8814","Potra2n4c9561",
#'              "Potra2n4c9798","Potra2n5c10575","Potra2n5c10753","Potra2n5c10979","Potra2n5c11241",
#'              "Potra2n5c11446","Potra2n5c11602","Potra2n5c11916","Potra2n5c11923","Potra2n5c12554",
#'              "Potra2n6c12962","Potra2n6c13323","Potra2n6c13849","Potra2n6c14336","Potra2n6c14711",
#'              "Potra2n6c15109","Potra2n7c15470","Potra2n7c15818","Potra2n7c15935","Potra2n7c16113",
#'              "Potra2n7c16264","Potra2n7c16566","Potra2n8c16773","Potra2n8c17318","Potra2n8c17825",
#'              "Potra2n8c17900","Potra2n8c18011","Potra2n8c18019","Potra2n9c19378","Potra2n9c19384",
#'              "Potra2n9c19454","Potra2n9c19718","Potra2n9c19735", "Potra2n6c14967","Potra2n9c19651",
#'              "Potra2n1c2250","Potra2n6c13361","Potra2n6c15120")
#' ray5_7 <- c("Potra2n10c20444","Potra2n11c22600","Potra2n14c27801","Potra2n14c27834",
#'             "Potra2n15c28669","Potra2n16c30308","Potra2n18c33011","Potra2n1c3480",
#'             "Potra2n2c4118","Potra2n2c4509","Potra2n3c6926","Potra2n3c6944",
#'             "Potra2n3c8005","Potra2n4c8449","Potra2n5c11003","Potra2n6c13067",
#'             "Potra2n6c15285","Potra2n7c16335","Potra2n9c19477","Potra2n13c25349",
#'             "Potra2n14c27607","Potra2n16c29243","Potra2n16c30464","Potra2n18c32429",
#'             "Potra2n18c32564","Potra2n1c1049","Potra2n1c2010","Potra2n1c232","Potra2n1c90",
#'             "Potra2n2c4144","Potra2n2c4942","Potra2n2c5467","Potra2n2c6010","Potra2n3c8020",
#'             "Potra2n4c10139","Potra2n4c9969","Potra2n13c25215","Potra2n18c32373",
#'             "Potra2n1c394","Potra2n2c4494", "Potra2n2c5697","Potra2n2c5875",
#'             "Potra2n4c9970","Potra2n5c10874", "Potra2n5c11107","Potra2n5c11849","Potra2n8c17115","Potra2n9c18880")
#' xpc_ray <- c("Potra2n15c28663","Potra2n16c29896","Potra2n1c3834","Potra2n1c749",
#'              "Potra2n4c9930","Potra2n6c13057","Potra2n6c13538","Potra2n6c15220",
#'              "Potra2n14c26690","Potra2n14c26986","Potra2n14c27082","Potra2n14c27101",
#'              "Potra2n16c30018","Potra2n16c30450","Potra2n17c30920","Potra2n18c32844",
#'              "Potra2n1c1897","Potra2n1c2429","Potra2n1c394","Potra2n2c4928",
#'              "Potra2n3c6993","Potra2n4c9205","Potra2n5c11091","Potra2n5c11250",
#'              "Potra2n5c11849","Potra2n6c14448","Potra2n8c17012","Potra2n8c17375","Potra2n9c18721")
#' fiber <- c("Potra2n11c22673","Potra2n11c23199","Potra2n15c28689","Potra2n16c30073",
#'            "Potra2n16c30273","Potra2n1c546","Potra2n2c4078","Potra2n4c9297",
#'            "Potra2n6c14160","Potra2n7c15752","Potra2n7c16228","Potra2n19c33344",
#'            "Potra2n10c20411","Potra2n5c11573","Potra2n13c24959","Potra2n5c11765",
#'            "Potra2n7c16495","Potra2n1c3619","Potra2n7c16234","Potra2n5c11571",
#'            "Potra2n3c7305","Potra2n4c8952","Potra2n5c12319","Potra2n12c24762",
#'            "Potra2n15c28411","Potra2n1c307","Potra2n1c953","Potra2n2c6235",
#'            "Potra2n16c30051","Potra2n9c19651","Potra2n1c4016","Potra2n11c22498",
#'            "Potra2n6c14201","Potra2n8c16778","Potra2n8c16946","Potra2n8c17885","Potra2n9c20072")
#' celldeath <- c("Potra2n14c26495","Potra2n14c27047","Potra2n1c3263","Potra2n5c11905", 
#'                "Potra2n4c10301","Potra2n6c15163","Potra2n16c29448","Potra2n10c21346",
#'                "Potra2n1c3148","Potra2n11c23290","Potra2n1c2692")
#' 
#' #' 5.2. From network gene aspwood paper
#' # Expansins family
#' expansin_fam <- c("Potra2n10c20623","Potra2n10c20953","Potra2n13c24992",
#'                   "Potra2n13c25752","Potra2n16c30156","Potra2n16c30169",
#'                   "Potra2n16c30491","Potra2n17c30722","Potra2n17c31117",
#'                   "Potra2n17c31179","Potra2n19c33898","Potra2n1c11",
#'                   "Potra2n1c2087","Potra2n1c3505","Potra2n1c960","Potra2n2c4734",
#'                   "Potra2n2c4737","Potra2n2c6293","Potra2n4c9142","Potra2n4c9551",
#'                   "Potra2n5c10644","Potra2n6c14375","Potra2n6c14588","Potra2n8c17125",
#'                   "Potra2n8c17409","Potra2n9c18642","Potra2n9c19851","Potra2n1c1267",
#'                   "Potra2n1c1318","Potra2n3c7743","Potra2n3c7783","Potra2n4c10039",
#'                   "Potra2n9c18880")
#' 
#' #' cellulase synthase genes: cesA4, A7-A, A7-B, A8-A, A8-B
#' ces <- c("Potra2n2c4078","Potra2n6c13749", "Potra2n18c32336", "Potra2n11c23199",
#'          "Potra2n4c8952")
#' 
#' #' Phloem: shouldn't exist in present dataset
#' phloem <- c("Potra2n17c30730","Potra2n12c24024","Potra2n15c28954","Potra2n4c9149",
#'             "Potra2n8c17353","Potra2n10c20892","Potra2n17c30726","Potra2n4c9146",
#'             "Potra2n10c20934","Potra2n8c17387","Potra2n1c1172","Potra2n17c30744",
#'             "Potra2n4c9160","Potra2n8c18492")
#' bsp <- c("Potra2n13c25437", "Potra2n13c25438", "Potra2n13c25439", "Potra2n13c25440",
#'          "Potra2n13c25589", "Potra2n158s34776", "Potra2n484s35780")
#' 
#' #' Cambium/phloem: PtOPS-A, PtOPS-B, PtBRX-A, PtBRX-B, CLE41, MYB83/PtMYB20/21
#' camb_or_phloem <- c("Potra2n6c14490","Potra2n16c30270","Potra2n1c1136",
#'                     "Potra2n3c7606","Potra2n2c4235","Potra2n9c19455")
#' 
#' #' Cambium/expanding xylem:
#' #' PtCDC2, PtAINT, PtAINT, PtAINT, PtAINT, LOG4/PtLOG6, PtPM2, PtPL1-27
#' cambium <- c("Potra2n16c30563","Potra2n2c5371","Potra2n6c14967","Potra2n227s35019",
#'              "Potra2n1c1493","Potra2n16c30088", "Potra2n5c10720","Potra2n3c6899")
#' 
#' #' SNARE-like, MYB83/PtMYB20/21, MYB46/PtMYB21, PtTMO5-A, PtTMO5-B, PtEXPA1
#' cambium_expandingxylem <- c ("Potra2n1c1686","Potra2n1c2305", "Potra2n9c19651",
#'                              "Potra2n10c21280","Potra2n8c17692", "Potra2n1c2087")
#' 
#' #' Expanding xylem: PtTMO5-A, PtTMO5-B, snd, vnd, vnd, snd, AtLAC11, LAC, LAC12
#' expx <- c ("Potra2n10c21280","Potra2n17c31797", "Potra2n12c24750","Potra2n3c7486",
#'            "Potra2n7c15489","Potra2n1c2159","Potra2n317s35299","Potra2n11c22797")
#' 
#' #' Lignified xylem or everywhere except expanding xylem
#' #' PtPM1, PtBFN1, PtBAM3-A, PtBAM3-B, PtBAM3-C, PtLHW-A, PtLHW-B
#' lignified_xylem <- c ("Potra2n2c6214","Potra2n10c21346", "Potra2n1c593", 
#'                       "Potra2n3c7083","Potra2n4c10378","Potra2n280s35193","Potra2n1c1878")
#' 
#' #' Xyloglucan and pectin biosynthesis genes
#' xylo <- c("Potra2n267s35161","Potra2n10c22184","Potra2n10c22183",
#'           "Potra2n12c24288","Potra2n15c28674")
#' 
#' pectin <- c("Potra2n14c26890","Potra2n2c5051","Potra2n4c10217","Potra2n16c29244",
#'             "Potra2n6c15382","Potra2n16c30575","Potra2n252s35126","Potra2n5c10503",
#'             "Potra2n3c7089","Potra2n1c601","Potra2n2c4911","Potra2n14c27069",
#'             "Potra2n12c24734","Potra2n114s34573")
#' 
#' #' 5.3. Laccases family (as per aspwood) specific to cluster or cell type could exist:
#' lac_fam <- c("Potra2n10c20714","Potra2n10c20806","Potra2n11c22797","Potra2n11c23187",
#'              "Potra2n12c24132","Potra2n13c25011","Potra2n14c27140","Potra2n8c17196",
#'              "Potra2n15c28842","Potra2n16c30017","Potra2n16c30239","Potra2n16c30241",
#'              "Potra2n16c30273","Potra2n16c30274","Potra2n6c13254",
#'              "Potra2n19c33334","Potra2n19c33575","Potra2n1c1775","Potra2n1c2108",
#'              "Potra2n1c2159","Potra2n1c2996","Potra2n1c3502","Potra2n1c447",
#'              "Potra2n317s35299","Potra2n4c9830","Potra2n5c11061","Potra2n5c11062",
#'              "Potra2n5c11063","Potra2n6c14486","Potra2n6c14515","Potra2n6c14571",
#'              "Potra2n6c14577","Potra2n7c16431","Potra2n8c17197","Potra2n8c17280",
#'              "Potra2n9c18755","Potra2n9c19221","Potra2n9c19763","Potra2n9c19826")
#' #' 
#' #' 5.4. For rays
#' #' PIRIN2 could be a marker of lignified xylem or rays
#' prn2 <- c("Potra2n2c4342","Potra2n14c27610", "Potra2n9c19440")
#' 
#' #' 5.5. For vessels 
#' #' right homolog:  "Potra2n6c13361",
#' acl5 <- c("Potra2n10c21673", "Potra2n6c13361", "Potra2n8c18019")
#' wat1 <- c("Potra2n5c10753","Potra2n2c6181")
#' #' XCP: xylem cysteine protease inhibitor: xcp2, xcp2, xcp1, xbcp3
#' #' XCPs were used as markers in Li et al., 2021; Chen et al, 2021 as well
#' xcp <-c("Potra2n5c10536","Potra2n2c6410", "Potra2n4c10232","Potra2n5c10760")
#' 
#' #' vnd
#' vnd <- c("Potra2n114s34556","Potra2n12c24750","Potra2n13c25337","Potra2n19c33624",
#'          "Potra2n1c1031","Potra2n3c7486","Potra2n5c11773","Potra2n7c16511")
#' 
#' rnaseiii <- c("Potra2n8c17393")
#' 
#' #' Kubo et al, 2005- highly exp vessel genes
#' # 
#' #' Metacaspase5
#' mc5 <- c("Potra2n8c17345","Potra2n10c20883")
#' #' 
#' #' 5.6. Core cell cycle TAIR
#' coreCyc <- c("Potra2n2c4390","Potra2n1c2176","Potra2n1c2151","Potra2n13c25822",
#'              "Potra2n1c1559","Potra2n8c17361","Potra2n12c24142","Potra2n12c24163",
#'              "Potra2n5c10575","Potra2n2c6426","Potra2n735s36629","Potra2n5c11241",
#'              "Potra2n2c5307","Potra2n2c6231","Potra2n9c19212","Potra2n8c17316",
#'              "Potra2n10c21618","Potra2n8c17791","Potra2n5c11923","Potra2n5c11433",
#'              "Potra2n5c11441","Potra2n7c16280","Potra2n6c13237","Potra2n1c2376",
#'              "Potra2n9c20077","Potra2n5c10858","Potra2n1c2727","Potra2n6c13528",
#'              "Potra2n3c7661","Potra2n16c30563","Potra2n6c14644","Potra2n12c24354",
#'              "Potra2n6c15109","Potra2n1c257","Potra2n2c4226","Potra2n4c9630",
#'              "Potra2n5c11567","Potra2n13c25044","Potra2n3c7422","Potra2n9c19647",
#'              "Potra2n18c32557","Potra2n14c27002","Potra2n1c2633","Potra2n7c16154",
#'              "Potra2n6c14563","Potra2n16c30120","Potra2n12c23916","Potra2n6c13129",
#'              "Potra2n8c17128","Potra2n13c25612","Potra2n669s36351","Potra2n4c10173",
#'              "Potra2n18c32954","Potra2n5c12590","Potra2n2c4525","Potra2n1c724","Potra2n7c15742")
#' # Cell cycle genes (Zhang et al., 2021)
#' g01phase <- c("Potra2n10c22352","Potra2n11c22967","Potra2n13c25403","Potra2n14c27234",
#'               "Potra2n14c27322","Potra2n18c32466","Potra2n1c2151","Potra2n1c2633",
#'               "Potra2n1c2837","Potra2n2c4043","Potra2n5c11567","Potra2n6c14588","
#'               Potra2n9c19242")
#' g2phase <- c("Potra2n10c21070","Potra2n10c21695","Potra2n10c21833","Potra2n10c22093","Potra2n10c22145","Potra2n10c22346","Potra2n11c22815","Potra2n11c23463","Potra2n12c24927","Potra2n13c25321","Potra2n13c25579","Potra2n13c26121","Potra2n13c26208","Potra2n14c26403","Potra2n14c27462","Potra2n14c27475","Potra2n15c27972","Potra2n15c28096","Potra2n15c28656","Potra2n15c28758","Potra2n15c28945","Potra2n16c29301","Potra2n16c29689","Potra2n16c29960","Potra2n16c30306","Potra2n16c30563","Potra2n18c32638","Potra2n18c32677","Potra2n18c32912","Potra2n19c34004","Potra2n19c34102","Potra2n1c1559","Potra2n1c2018","Potra2n1c2288","Potra2n1c2573","Potra2n1c2758","Potra2n1c2942","Potra2n1c3087",
#'              "Potra2n1c3141","Potra2n1c3180","Potra2n1c3295","Potra2n1c3611","Potra2n1c410","Potra2n2c4140","Potra2n2c4292","Potra2n2c4299",
#'              "Potra2n2c4587","Potra2n2c4715","Potra2n2c4966","Potra2n2c5259","Potra2n2c5415","Potra2n2c5455","Potra2n2c5494","Potra2n2c5545",
#'              "Potra2n2c5663","Potra2n2c6426","Potra2n372s35460","Potra2n3c7249","Potra2n3c7951","Potra2n3c7972","Potra2n3c8097","Potra2n3c8273","Potra2n4c8773","Potra2n4c9741","Potra2n5c10510","Potra2n5c10575","Potra2n5c10579","Potra2n5c11241","Potra2n5c11640","Potra2n5c12275","Potra2n680s36426","Potra2n6c14644","Potra2n6c14653","Potra2n6c14979","Potra2n6c15055","Potra2n6c15109","Potra2n7c15462","Potra2n7c16145","Potra2n8c17814","Potra2n8c17973","Potra2n8c18569","Potra2n8c18623","Potra2n9c18760","Potra2n9c19217","Potra2n9c19378","Potra2n9c20077")
#' sphase <- c("Potra2n10c20260","Potra2n10c20508","Potra2n11c22684","Potra2n13c26037","Potra2n13c26116","Potra2n14c27103","Potra2n14c27106","Potra2n15c28275","Potra2n15c28844","Potra2n16c29833","Potra2n17c31198","Potra2n18c33119","Potra2n1c101","Potra2n1c2877","Potra2n1c3343","Potra2n1c3610","Potra2n2c4073","Potra2n2c4499","Potra2n2c5307","Potra2n2c6056","Potra2n3c7036","Potra2n3c7597","Potra2n4c9636","Potra2n5c10911","Potra2n5c12603","Potra2n6c13036","Potra2n6c14627","Potra2n7c15963","Potra2n7c16633","Potra2n8c16885","Potra2n8c17847","Potra2n9c19886")
#' mphase <- c ("Potra2n10c20592","Potra2n10c20846","Potra2n10c21494","Potra2n10c21734","Potra2n10c22098","Potra2n11c22486","Potra2n11c22945","Potra2n11c23018","Potra2n12c24163","Potra2n13c25058","Potra2n13c25262","Potra2n13c25822","Potra2n13c25852","Potra2n14c27365","Potra2n14c27462","Potra2n14c27797","Potra2n15c29215","Potra2n16c29450","Potra2n16c30378","Potra2n16c30406","Potra2n16c30474","Potra2n17c30910","Potra2n17c31100","Potra2n17c31123","Potra2n18c32687","Potra2n19c34338","Potra2n1c2214","Potra2n1c2376","Potra2n1c2511","Potra2n1c2591","Potra2n1c3910","Potra2n1c4000","Potra2n1c671","Potra2n1c923","Potra2n2c4070","Potra2n2c4172","Potra2n3c7249","Potra2n3c7843","Potra2n4c8944","Potra2n4c8987","Potra2n5c10575","Potra2n5c10858","Potra2n5c10979","Potra2n5c11566","Potra2n5c11796","Potra2n5c11923","Potra2n5c12554","Potra2n6c12962","Potra2n6c13237","Potra2n6c13440","Potra2n6c13456","Potra2n6c14169","Potra2n6c15109","Potra2n7c15465","Potra2n7c16119","Potra2n8c17128","Potra2n8c17171","Potra2n8c17318","Potra2n8c18190","Potra2n8c18202","Potra2n8c18287","Potra2n9c19378","Potra2n9c19459","Potra2n9c19549","Potra2n6c13657", "Potra2n16c29642","Potra2n17c30791", "Potra2n6c14415")
#' # s.genes <- c("Potra2n10c20508", "Potra2n1c2877")
#' # g2m.genes <- c("Potra2n8c18287", "Potra2n6c12962","Potra2n5c11241","Potra2n2c6426")
#' # g1.genes <- c("Potra2n14c27322", "Potra2n1c2151")
#' 
#' #' 5.7. Genes in suppl table of Shi et al., 2017, laser microdissection paper
#' Shifib <- c("Potra2n1c546","Potra2n8c16946","Potra2n10c20318","Potra2n13c24959","Potra2n1c3619","Potra2n2c4078","Potra2n5c12319","Potra2n6c14201","Potra2n7c16495","Potra2n7c15752","Potra2n8c17885","Potra2n11c23199","Potra2n11c22673","Potra2n16c30051","Potra2n16c30091","Potra2n1c452","Potra2n1c829","Potra2n1c1020","Potra2n1c2246","Potra2n1c3530","Potra2n2c4212","Potra2n3c7978","Potra2n3c7705","Potra2n3c7477","Potra2n4c8837","Potra2n4c9512","Potra2n5c11659","Potra2n5c11186","Potra2n7c15765","Potra2n7c15466","Potra2n8c16778","Potra2n9c20072","Potra2n10c21057","Potra2n13c24976","Potra2n15c28920","Potra2n16c29749","Potra2n18c32675","Potra2n1c1059","Potra2n1c2112","Potra2n1c2242","Potra2n584s36051","Potra2n5c12633","Potra2n14c27757","Potra2n17c31154","Potra2n18c32560","Potra2n758s36699")
#' ShiCelWal <- c("Potra2n1c307","Potra2n3c6817","Potra2n1c1541","Potra2n2c5206","Potra2n3c7945","Potra2n3c6847","Potra2n3c6783","Potra2n4c8952","Potra2n5c11571","Potra2n6c15127","Potra2n6c14160","Potra2n6c13749","Potra2n5c11670","Potra2n9c20060","Potra2n9c19275","Potra2n9c19246","Potra2n12c23766","Potra2n12c24640","Potra2n14c26547","Potra2n14c26628")
#' ShiXyl <- c("Potra2n1c953","Potra2n1c1657","Potra2n1c2250","Potra2n1c2305","Potra2n1c3191","Potra2n1c3797","Potra2n1c4016","Potra2n2c6375","Potra2n2c6070","Potra2n2c5746","Potra2n2c5385","Potra2n2c5151","Potra2n2c5075","Potra2n14c27102","Potra2n2c4785","Potra001157g10045","Potra2n2c4080","Potra2n3c8027","Potra2n3c7305","Potra2n3c6563","Potra2n4c9127","Potra2n4c9297","Potra2n4c9868","Potra2n5c12307","Potra2n5c11773","Potra2n5c10901","Potra2n5c10869","Potra2n5c10648","Potra2n6c13965","Potra2n7c16575","Potra2n7c16511","Potra2n203s34957","Potra2n7c15498","Potra2n7c15489","Potra2n8c18444","Potra2n369s35436","Potra2n9c19727","Potra2n9c19651","Potra2n10c22314","Potra2n10c21389","Potra2n11c23287","Potra2n4c8839","Potra2n11c22968","Potra2n12c23697","Potra2n12c24046","Potra2n12c24750","Potra2n12c24760","Potra2n13c25576","Potra2n14c26711","Potra2n14c26818","Potra2n14c26918","Potra2n14c27178","Potra2n15c28915","Potra2n15c28466","Potra2n15c28118","Potra2n16c29403","Potra2n6c14928","Potra2n16c30501","Potra2n17c31797","Potra2n17c31126","Potra2n17c30794","Potra2n19c33603")
#' ShiVes <- c("Potra2n2c5125","Potra2n6c13849","Potra2n14c26495","Potra2n12c24354","Potra2n17c31793","Potra2n3c6715","Potra2n4c9874","Potra2n4c10104","Potra2n6c13368","Potra2n7c16311","Potra2n9c19893","Potra2n9c18811","Potra2n13c26225","Potra2n14c26600","Potra2n1c2330","Potra2n1c3263","Potra2n11c22702","Potra2n2c5282","Potra2n2c5031","Potra2n2c4928","Potra2n3c7393","Potra2n3c7342","Potra2n5c11905","Potra2n6c15090","Potra2n6c14723","Potra2n6c12972","Potra2n8c17258","Potra2n10c21280","Potra2n10c21180","Potra2n11c23465","Potra2n443s35699","Potra2n12c24560","Potra2n14c26506","Potra2n14c26919","Potra2n14c27047","Potra2n14c27085","Potra2n16c29522","Potra2n16c29709","Potra2n18c32911","Potra2n18c32086")
#' Shietc <- c("Potra2n8c16843","Potra2n1542s37129","Potra2n246s35111","Potra2n10c20411","Potra2n19c33285","Potra2n2c5125","Potra2n2c4880","Potra2n6c13849","Potra2n11c22498","Potra2n14c26495","Potra2n17c31793","Potra2n2c4256","Potra2n17c30701","Potra2n3c6715","Potra2n11c23285","Potra2n14c26600","Potra2n7c16311","Potra2n2c4256",
#'             "Potra2n10c21180","Potra2n11c22702","Potra2n11c23465","Potra2n12c24560","Potra2n14c26919","Potra2n14c27047","Potra2n14c27085","Potra2n14c27828","Potra2n16c29522","Potra2n16c29709","Potra2n18c32086","Potra2n18c32911","Potra2n18c33178","Potra2n1c2330","Potra2n1c3597","Potra2n1c923","Potra2n2c5031","Potra2n3c7040","Potra2n3c7342","Potra2n4c10104","Potra2n6c15090","Potra2n7c16348","Potra2n8c17258","Potra2n9c18811")
#' 
#' #' 5.8. Markers in SuppData (Li et al., 2021)
#' #' SuppData6 (Li et al., 2021) cluster 0 to 11: 
#' #' Of these, Cluster 0,2,3,4,5,6,8,9,10,11 genes were selected, others do not exist
#' li1ST6 <- c("Potra2n1c596","Potra2n19c33855", "Potra2n1c1565", "Potra2n1c885","Potra2n13c26019", "Potra2n1c1049","Potra2n13c26037","Potra2n14c27483","Potra2n18c33199","Potra2n9c19467","Potra2n16c29966")
#' 
#' #' Genes in SuppData7 (Li et al., 2021):
#' li2ST7 <- c("Potra2n10c20736","Potra2n10c21184","Potra2n10c22101","Potra2n10c22334","Potra2n11c22555","Potra2n11c22739","Potra2n12c24388","Potra2n14c26962","Potra2n16c30018","Potra2n16c30222","Potra2n16c30520","Potra2n17c30922","Potra2n17c30987","Potra2n17c31220","Potra2n1c3290","Potra2n1c3599", "Potra2n1c777","Potra2n1c892","Potra2n2c4139","Potra2n2c4470","Potra2n2c4509","Potra2n2c5397","Potra2n2c5623","Potra2n2c5742","Potra2n2c5999","Potra2n2c6010","Potra2n3c6548","Potra2n3c6709","Potra2n3c6863","Potra2n3c7085","Potra2n3c8387","Potra2n5c10517","Potra2n5c11113","Potra2n5c11775","Potra2n5c12062","Potra2n6c15294","Potra2n7c16270","Potra2n7c16346","Potra2n8c17128","Potra2n8c18019","Potra2n9c18964","Potra2n9c19384","Potra2n9c19859","Potra2n5c10911","Potra2n14c27106")
#' 
#' liSelMarkST7<- c("Potra2n2c4994","Potra2n14c27682","Potra2n8c18380","Potra2n7c15772","Potra2n4c9797","Potra2n2c6426","Potra2n5c11306","Potra2n16c30290","Potra2n14c27721","Potra2n1c1438","Potra2n11c23043","Potra2n5c11668","Potra2n5c10940","Potra2n5c10926","Potra2n1c3835","Potra2n6c14358","Potra2n6c14530","Potra2n16c29448","Potra2n9c19538","Potra2n4c9403","Potra2n3c7365","Potra2n9c19373","Potra2n2c5831","Potra2n5c12253","Potra2n16c29767","Potra2n4c10402","Potra2n7c15743","Potra2n9c19117","Potra2n820s36802","Potra2n4c9402")
#' #'
#' #' 5.9. Genes for clusters and cell types (Chen et al., 2021)
#' mark <- c("Potra2n6c14918","Potra2n4c9802","Potra2n15c29012","Potra2n1c986","Potra2n5c12689",
#'           "Potra2n2c5262","Potra2n18c32485","Potra2n1c2071","Potra2n3c7164","Potra2n5c11923",
#'           "Potra2n4c8852","Potra2n19c34474","Potra2n1c1024","Potra2n7c16246","Potra2n14c27521",
#'           "Potra2n369s35437","Potra2n12c24741","Potra2n5c10890","Potra2n1c3322","Potra2n17c31073",
#'           "Potra2n12c24778","Potra2n1c769","Potra2n1c957","Potra2n369s35445","Potra2n6c13032",
#'           "Potra2n3c7257","Potra2n4c8671")
#' 
#' #' 5.10. For fibers
#' #' Late fiber gene list sent by Hannele
#' lateFibHanle <- c ("Potra2n6c14170","Potra2n1c1849","Potra2n6c15100","Potra2n6c15101","Potra2n19c33746","Potra2n5c12551","Potra2n1c94","Potra2n4c9157","Potra2n12c24599","Potra2n4c9039","Potra2n8c17428","Potra2n19c34385","Potra2n15c28085","Potra2n6c14603","Potra2n15c28689","Potra2n12c24486","Potra2n633s36227","Potra2n12c24762","Potra2n16c29514","Potra2n19c33745","Potra2n125s34627","Potra2n13c25019","Potra2n13c25443","Potra2n15c28411","Potra2n16c30073","Potra2n19c33329","Potra2n19c33344","Potra2n19c33747","Potra2n3c6760","Potra2n4c9502")
#' lowexp <- c("Potra2n6c14170","Potra2n1c1849","Potra2n6c15100","Potra2n6c15101","Potra2n19c33746","Potra2n12c24599","Potra2n4c9039","Potra2n6c14603","Potra2n12c24762","Potra2n19c33745","Potra2n15c28411","Potra2n19c33344")
#' 
#' #' 5.11. Tung et al., 2022
#' lignin <- c("Potra2n3c6783","Potra2n6c14201","Ptra2n8c16946","Potra2n14c26628","Potra2n16c30091","Potra2n5c10696","Potra2n8c16891","Potra2n13c25045","Potra2n17c30922")
#' figs13 <- c("Potra2n19c33452","Potra2n1c856","Potra2n1c1639","Potra2n1c1897","Potra2n2c5835","Potra2n2c4490","Potra2n3c8296","Potra2n4c9930","Potra2n4c10139","Potra2n8c17951","Potra2n11c22655","Potra2n15c28183","Potra2n15c28184")
#' figs19 <- c("Potra2n17c31807","Potra2n8c16805","Potra2n5c12126","Potra2n5c11373","Potra2n14c26481", "Potra2n7c16454")
#' s14fib <- c("Potra2n10c20156","Potra2n10c22093","Potra2n12c23702","Potra2n12c23837","Potra2n12c24890","Potra2n13c25019","Potra2n13c25704","Potra2n13c26198","Potra2n14c27637","Potra2n15c28014","Potra2n15c29105","Potra2n15c29113","Potra2n15clength29164","Potra2n16c29484","Potra2n16c30218","Potra2n16c30267","Potra2n17c31215","Potra2n17c31695","Potra2n18c32051","Potra2n18c32369","Potra2n19c33344","Potra2n1c1134","Potra2n1c1591","Potra2n1c1997","Potra2n1c2382","Potra2n1c3857","Potra2n2c4696","Potra2n2c6000","Potra2n2c6306","Potra2n3c6819","Potra2n3c7291","Potra2n4c8742","Potra2n4c9502","Potra2n4c9523","Potra2n4c9719","Potra2n5c11070","Potra2n5c11913","Potra2n5c12161","Potra2n6c14379","Potra2n6c14530","Potra2n6c14780","Potra2n6c14781","Potra2n7c16243","Potra2n8c17402","Potra2n8c17650","Potra2n8c18023","Potra2n8c18128","Potra2n8c18294","Potra2n9c18710","Potra2n9c19633","Potra2n9c19741","Potra2n9c20007","Potra2n6c15285","Potra2n6c13790")
#' s14ves <- c("Potra001029g08632","Potra2n10c20716","Potra2n10c20846","Potra2n10c21715","Potra2n11c22924","Potra2n12c23924","Potra2n12c24430","Potra2n12c24758","Potra2n13c25053","Potra2n13c25938","Potra2n13c26020","Potra2n14c26398","Potra2n14c26658","Potra2n14c27455","Potra2n14c27675","Potra2n15c28703","Potra2n16c29368","Potra2n16c29956","Potra2n18c32909","Potra2n1c2383","Potra2n1c3919","Potra2n1c483","Potra2n1c631","Potra2n203s34942","Potra2n2c4081","Potra2n2c5047","Potra2n2c5927","Potra2n2c6217","Potra2n2c6318","Potra2n3c7065","Potra2n3c7465","Potra2n4c8586","Potra2n4c8788","Potra2n4c9127","Potra2n4c9402","Potra2n514s35859","Potra2n5c10700","Potra2n5c11603","Potra2n5c11639","Potra2n5c12176","Potra2n5c12351","Potra2n5c12441","Potra2n6c14479","Potra2n6c14597","Potra2n7c15960","Potra2n8c17079","Potra2n8c17589","Potra2n8c17613","Potra2n8c18257","Potra2n9c19793","Potra2n9c19878")
#' s14ray <- c("Potra001886g15019","Potra2n10c20181","Potra2n10c20234","Potra2n10c20683","Potra2n10c20932","Potra2n10c21140","Potra2n10c21228","Potra2n10c21351","Potra2n10c21638","Potra2n10c21713","Potra2n10c21725","Potra2n10c21912","Potra2n10c22232","Potra2n11c22873","Potra2n11c23444","Potra2n12c23733","Potra2n12c23990","Potra2n12c24140","Potra2n13c24999","Potra2n13c25225","Potra2n13c25636","Potra2n13c25697","Potra2n13c25934","Potra2n13c26148","Potra2n14c26461","Potra2n14c26462","Potra2n14c26464","Potra2n14c26497","Potra2n14c26529","Potra2n14c27128","Potra2n14c27196","Potra2n14c27357","Potra2n14c27545","Potra2n14c27691","Potra2n14c27805","Potra2n14c27850","Potra2n154s34741","Potra2n15c28231","Potra2n15c29127","Potra2n16c29491","Potra2n16c29513","Potra2n16c30395","Potra2n17c30597","Potra2n19c33348","Potra2n19c33407","Potra2n19c33924","Potra2n19c34230","Potra2n1c1749","Potra2n1c2034","Potra2n1c2037","Potra2n1c2350","Potra2n1c3008","Potra2n1c3440","Potra2n1c3617","Potra2n1c3711","Potra2n1c3908","Potra2n1c531","Potra2n1c677","Potra2n1c693","Potra2n1c714","Potra2n258s35141","Potra2n2c4261","Potra2n2c4724")
#' s14ray1 <- c("Potra2n2c5150","Potra2n2c5325","Potra2n2c5691","Potra2n2c5770","Potra2n2c5997","Potra2n2c6123","Potra2n2c6195","Potra2n3c6552","Potra2n3c6865","Potra2n3c7178","Potra2n3c7315","Potra2n3c7363","Potra2n3c7903","Potra2n4c8432","Potra2n4c9135","Potra2n4c9752","Potra2n4c9884","Potra2n5c10715","Potra2n5c10780","Potra2n5c11137","Potra2n5c11557","Potra2n5c12296","Potra2n5c12429","Potra2n5c12447","Potra2n5c12522","Potra2n5c12541","Potra2n5c12687","Potra2n5c12804","Potra2n620s36177","Potra2n6c12828","Potra2n6c13140","Potra2n6c13191","Potra2n6c13220","Potra2n6c13318","Potra2n6c13467","Potra2n6c13491","Potra2n6c13586","Potra2n6c13868","Potra2n6c13967","Potra2n6c14775","Potra2n6c14862","Potra2n6c15097","Potra2n6c15141","Potra2n6c15233","Potra2n6c15336","Potra2n747s36678","Potra2n7c15481","Potra2n7c15868","Potra2n7c15955","Potra2n7c16004","Potra2n7c16369","Potra2n7c16410","Potra2n7c16446","Potra2n8c16988","Potra2n8c17145","Potra2n8c17251","Potra2n8c17293","Potra2n8c17331","Potra2n9c18794","Potra2n9c18942","Potra2n9c18949","Potra2n9c19022","Potra2n9c19045","Potra2n9c19321","Potra2n9c19718","Potra2n9c19805")
#' s18 <- c("Potra2n1c1345","Potra2n2c5467","Potra2n2c4942","Potra2n2c4144","Potra2n2c4078","Potra2n4c8430","Potra2n568s35990","Potra2n5c12126","Potra2n5c11373","Potra2n5c10753","Potra2n6c14930","Potra2n155s34744","Potra2n16c30073","Potra2n6c13614","Potra2n7c16454","Potra2n8c16805","Potra2n8c16999","Potra2n8c18231","Potra2n10c21170","Potra2n14c26481","Potra2n15c28766","Potra2n15c28689","Potra2n17c31807","Potra2n18c32244","Potra2n19c33790")
#' 
#' # Selected markers (Tung et al., 2022)
#' LateFib_Ves <- c("Potra2n5c10753","Potra2n8c16805","Potra2n17c31807","Potra2n19c33790")
#' InterFib <- c("Potra2n5c12126","Potra2n6c14930","Potra2n14c26481")
#' RayOrgan <- c("Potra2n5c11373","Potra2n8c16999")
#' EarlyFib <- c("Potra2n568s35990","Potra2n7c16454")
#' RayPrecu <- c("Potra2n8c18231","Potra2n18c32244")
#' FibOrgan <- c("Potra2n1c1345","Potra2n4c8430","Potra2n6c13614","Potra2n15c28766")
#' Fib <- c("Potra2n2c4078","Potra2n16c30073","Potra2n10c21170","Potra2n15c28689")
#' RayParen <- c("Potra2n2c5467","Potra2n2c4942","Potra2n2c4144","Potra2n155s34744")
#' 
#' #' 5.12. Lignin genes from Mikko's paper
#' ligMiko <- c("Potra2n19c33285","Potra2n1c351","Potra2n10c20411","Potra2n9c19275","Potra2n9c19246","Potra2n8c17885","Potra2n8c16946","Potra2n7c16495","Potra2n6c15127","Potra2n6c14201","Potra2n5c11765","Potra2n3c7945","Potra2n3c6847","Potra2n3c6817","Potra2n3c6783","Potra2n1c307","Potra2n1c2649","Potra2n1c1541","Potra2n16c30091","Potra2n16c29966","Potra2n13c24959","Potra2n12c23766")
#' 
#' # 5.13. O-Carillo et al cell death genes
#' death <- c("Potra2n4c10301","Potra2n16c29448","Potra2n10c21346","Potra2n1c3148","Potra2n4c8470")
#' 
#' #5.14. ERF 81-88
#' erf <- c("Potra2n15c29002","Potra2n1c788","Potra2n2c4898","Potra2n3c7265",
#'          "Potra2n13c24952","Potra2n14c27085","Potra2n19c33278","Potra2n12c23983")
#' 
#' #5.15.
#' GT43Mello<- c("Potra2n2c5445","Potra2n5c11571","Potra2n6c14160","Potra2n6c13189",
#'               "Potra2n7c16228","Potra2n16c30051","Potra2n18c32894")
#' 
#' # 5.16
#' # GUS lines
#' gus <- c("Potra2n3c7095","Potra2n4c9818","Potra2n5c11792","Potra2n9c19178",
#'          "Potra2n18c32992","Potra2n19c33746","Potra2n10c20181","Potra2n13c26298",
#'          "Potra2n13c25538","Potra2n15c28411","Potra2n2c5389","Potra2n4c9202",
#'          "Potra2n9c18829","Potra2n10c22224","Potra2n17c31277","Potra2n10c20842",
#'          "Potra2n16c30091","Potra2n8c17315","Potra2n2c5294")
#' # issa Validation
#' issa <- c("Potra2n4c8952","Potra2n5c10536","Potra2n2c6410","Potra2n14c27598", "Potra2n4c9696", "Potra2n2c5606", "Potra2n7c16495", "Potra2n5c11765", "Potra2n7c16234", "Potra2n5c11573", "Potra2n3c7305", "Potra2n2c6235")
#' 
#' # 5.17 Larisch et al., 2012
#' larish <- c("Potra2n13c26124","Potra2n15c28479","Potra2n13c26252","Potra2n3c7481","Potra2n10c21640","Potra2n5c11320","Potra2n10c20915","Potra2n8c17983","Potra2n2c5606","Potra2n4c9696","Potra2n9c19920","Potra2n18c32485","Potra2n8c18284","Potra2n10c22245","Potra2n17c31716","Potra2n19c33927","Potra2n16c29957","Potra2n2c4577","Potra2n10c21132","Potra2n1c1290","Potra2n17c30652","Potra2n18c32561","Potra2n4c9927","Potra2n15c28361","Potra2n10c20361","Potra2n5c12183","Potra2n8c18507","Potra2n18c32183","Potra2n13c26115","Potra2n8c18472","Potra2n8c17720","Potra2n15c28125","Potra2n7c15917","Potra2n5c12157","Potra2n8c17568","Potra2n7c15596","Potra2n10c21949","Potra2n14c27777","Potra2n14c27633")
#' 
#' # 5.19. Expansion related: feronia, ralf, lrx, md 
#' fer <- c("Potra2n16c30521", "Potra2n16c30523", "Potra2n6c14356")
#' ralf <- c("Potra2n13c26153","Potra2n14c27428","Potra2n17c31433","Potra2n1c2781",
#'           "Potra2n267s35162","Potra2n4c8790","Potra2n4c8791","Potra2n5c12610")
#' lrx8 <- c("Potra2n1c2637", "Potra2n6c13964", "Potra2n6c14636", "Potra2n9c19254",
#'           "Potra2n14c26497")
#' md <- c("Potra2n14c27025","Potra2n6c14356","Potra2n8c17700","Potra2n16c30033",
#'         "Potra2n10c21287")