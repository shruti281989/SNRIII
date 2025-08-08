suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(readxl)
  # library(qs)
  library(gplots)
})

set.seed(42)
# load seurat
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"

# In case to plot for one sample, split the object
# split_seurat <- SplitObject(integ, split.by = "sample")
# split_seurat <- split_seurat[c("ctrl", "kno")]
# and plot using the object split_seurat$ctrl 

# Load markers
selectedmarker <- read.table("~/shruti/SNR-u2023011/analysis/markers/selMarker.txt",
                             sep = '\t',header=TRUE)
clustMarker <- read.table("~/shruti/SNR-u2023011/analysis/markers/clustMarker.txt",
                          sep = '\t',header=TRUE)
ploidy <- read.table("~/shruti/SNR-u2023011/analysis/markers/endo_expPotra.txt",
                     sep = '\t',header=TRUE)
ribosomal <- read.table("~/shruti/SNR-u2023011/analysis/markers/Ribo_gene_table.txt",
                        sep = '\t',header=F)
protoplasting <- read.table("~/shruti/SNR-u2023011/analysis/markers/pplast.txt",
                            sep = '\t',header=F)
degWilcox <- read.table("~/shruti/SNRIII/data/SeuratOut/output/markerWilcox_lfc1_fdr0.01_pct0.5.txt",
                        header = T)
nit <- read_excel("~/shruti/SNRIII/data/SeuratOut/bulkDegvsScDeg.xlsx", sheet = 2)

degTime <- read.table("~/shruti/SNRIII/data/degTableS1E.txt",
                       header = T, sep = '\t')

#Cluster wise heatmap: scale the data first: 
integ <- ScaleData(integ, features = rownames(integ))
DoHeatmap(integ, features = degAll$Potra, group.by="sample",
          group.colors = viridis(100))

# # for label with own/AT annotation 
# nitRes <- nitRes[order(nitRes$AtName),]
# # nitRes <- nitRes[order(nitRes$Type),]
# # extract the genes first, since gene ids are not unique
# nrt <- nitRes[nitRes$Type == "Expansion",]
# markeranno <- paste(nrt$AtName)
# names(markeranno) <- nrt$GeneId

DotPlot(subset(integ, subset = sample =="ctrl"), features = nrt$GeneId,
        dot.scale = 10,idents=c("Fiber 2", "Fiber 1", "Early Fiber", 
                                "Fiber Precursor 2","Fiber Precursor 1",
                                "Late Vessel", "Early Vessel","Fusiform Initial",
                                "Ray 18","Ray 7","Ray 5","Ray/Fusiform Initial",
                                "Phloem like","Cambium"))+RotatedAxis()+ 
  scale_x_discrete(labels = markeranno)+ylab(NULL)+xlab(NULL)+
  scale_colour_gradient(low = "grey", high = "darkblue")


# plotting function
generateDotPlot <- function(data, markerfile, cellType) {
  markers <- markerfile[markerfile$CellType == cellType, ]$GeneId
  
  # png(file.path(here("data/SeuratOut/plots/withClCyc/"),
  png(file.path(here("data/SeuratOut/plots/mtCp/"),
                paste0(cellType,".png")), res= 250,height = 4000, width = 3000)
  
  p <- DotPlot(data, features = markers, cols ="RdBu", dot.scale = 10) +
    RotatedAxis()+  ylab(NULL) +xlab(NULL)+coord_flip() +labs(title = cellType)
  
  dev.off()
}

DotPlot(integ, selectedmarker[selectedmarker$CellType == "chen photosyn", ]$GeneId,
        cols ="RdBu", dot.scale = 10) + RotatedAxis()+ coord_flip()
# cols = c("#005AB5", "#DC3220"),
DotPlot(integ, split.by = "sample", cols = c("#fde623ff", "#420051ff"),
        features = "Potra2n15c29002",dot.scale = 10)+RotatedAxis()+ylab(NULL)+
  xlab(NULL)+coord_flip()
#scale_colour_gradient(low = "grey", high = "#420051ff")
# subset(integ, subset = sample =="kno"
#c("#420051ff", "#fde623ff"),

# Split Violin 
plots <- VlnPlot(integ, features = "Potra2n1c723", cols = c("blue", "red"),
                 split.by = "sample", pt.size = 0, combine = FALSE, split.plot = T)
wrap_plots(plots = plots, ncol = 1)

# Stacked violin
a <- VlnPlot(integ, features =c("Potra2n1c723","Potra2n16c29671"), stack = T,
             sort = T) + theme(legend.position = "none") + ggtitle("nlp7fam")
plot_grid(a)

# png("lac_fam.png")
# for (i in length(lac_fam)) {
#   print(DotPlot(seurat_phase, features = lac_fam)+RotatedAxis())}
# dev.off()
# 

# Potential expansion related: feronia, ralf, lrx, md 
fer <- c("Potra2n16c30521", "Potra2n16c30523", "Potra2n6c14356")
ralf <- c("Potra2n13c26153","Potra2n14c27428","Potra2n17c31433","Potra2n1c2781",
          "Potra2n267s35162","Potra2n4c8790","Potra2n4c8791","Potra2n5c12610")
lrx8 <- c("Potra2n1c2637", "Potra2n6c13964", "Potra2n6c14636", "Potra2n9c19254",
          "Potra2n14c26497")
md <- c("Potra2n14c27025","Potra2n6c14356","Potra2n8c17700","Potra2n16c30033",
        "Potra2n10c21287")

# only metabolism related in Figure 3 of scrnaseq paper
nitRes <- read.table("~/shruti/SNR-u2023011/analysis/markers/fig3marker.txt",
                     sep = '\t',header=TRUE)
rownames(nitRes) <- nitRes$Potra

DotPlot(integ, split.by = "sample", cols = c("#fde623ff", "#420051ff"),
        features = nit$Potra, dot.scale = 8)+
  RotatedAxis()+ylab(NULL)+xlab(NULL)+coord_flip()
DotPlot(subset(integ, subset = sample =="kno"), 
        features = nit$Potra,dot.scale = 10)+
        RotatedAxis()+ylab(NULL)+xlab(NULL)+
          scale_colour_gradient(low = "grey", high = "#420051ff")+coord_flip()

##E66100; #5D3A9B; #fde623ff; #420051ff; #1A85FF; #D41159; #1AFF1A; #5D3A9; 
##994F00; #006CD1

# 2h time series genes heatmap
degTimeSeries <- read.csv("~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/degSNRTime.csv",)
features_to_plot <- degTimeSeries %>% 
  filter(`level.in.KNO3` == "Upregulated" & timepoint == "KNO3_2h_vs_KCL_2h") %>%
  pull(X)
features_to_plot <- readLines("deg2hTF.txt")

# dn2h <- degTimeSeries %>% 
#   filter(`level.in.KNO3` == "Downregulated" & timepoint == "KNO3_2h_vs_KCL_2h") %>% 
#   pull(X)

integ <- ScaleData(integ, features = rownames(integ), assay = "RNA")
avgexp = AverageExpression(subset(integ, subset = sample =="ctrl"), assay="RNA",
                           features = features_to_plot,
                           return.seurat = T, group.by = 'integrated_snn_res.0.6')

DoHeatmap(avgexp, features = features_to_plot,angle = 0) + guides(color="none") +
  scale_fill_gradientn(colors = viridis(20))+
  theme(axis.text.y = element_text(size = 1))

# PRN2 upeg in PRN2::GFP positive vs all cells
# prnUp<- read_xlsx("~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/SupplementaryTable_PRN2CellsDeseq.xlsx",
#                   sheet = "S1B", skip = 1, col_names = T)
# features_to_plot <- prnUp %>% filter(`Upregulated in` == "PRN2::GFP posiitve cells") %>%
#   pull(`P. tremula Gene IDs`)

# heatmap of the average of upregulated genes in each cluster of SNRIII-scrnaseq
# separately for each sample
sample_list <- SplitObject(integ, split.by = "sample")

features_to_plot <- read.table("~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/degSNRIII.txt",
                               header = TRUE, fill = TRUE, sep = "\t", quote = "")
features_to_plot <- features_to_plot %>% filter(status=="up in kno" & FDR <0.05) %>%
  pull("GeneID")

features_to_plot <- nit$Potra

# Calculate average expression for each sample's clusters
average_expression_list <- lapply(sample_list, function(x) {
  AverageExpression(x, assays = "RNA", return.seurat = TRUE)
})

# Combine data for each sample and filter for the list of genes
average_expression_combined <- lapply(names(average_expression_list), function(sample) {
  as.data.frame(GetAssayData(average_expression_list[[sample]], slot = "data")) %>%
    rownames_to_column("Gene") %>%
    filter(Gene %in% features_to_plot) %>%
    mutate(Sample = sample)
}) %>% bind_rows()

# Scale expression gene-wise for all genes
heatmap_data <- lapply(names(average_expression_list), function(sample) {
  as.data.frame(GetAssayData(average_expression_list[[sample]], slot = "data")) %>%
    rownames_to_column("Gene") %>%
    filter(Gene %in% features_to_plot) %>%
    pivot_longer(-Gene, names_to = "Cluster", values_to = "Expression") %>%
    group_by(Gene) %>%
    mutate(ScaledExpression = scale(Expression), Sample = sample) %>%
    ungroup()
}) %>% bind_rows()

# Aggregate duplicates by taking the mean
heatmap_data_clean <- heatmap_data %>%
  group_by(Gene, Cluster, Sample) %>%
  summarise(ScaledExpression = mean(ScaledExpression), .groups = 'drop')

# Extract and save gene order from clustering in kno
reference_sample <- unique(heatmap_data_clean$Sample)[2]  # kno as reference
reference_data <- heatmap_data_clean %>%
  filter(Sample == reference_sample) %>%
  select(Gene, Cluster, ScaledExpression) %>%
  pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Check for NA/NaN/Inf values in the reference_data
na_check <- apply(reference_data, 1, function(x) any(is.na(x) | is.infinite(x)))
clean_reference_data <- reference_data[!na_check, ]

# hierarchical clustering on the cleaned data
gene_order <- hclust(dist(clean_reference_data))$order
ordered_genes <- rownames(clean_reference_data)[gene_order]

# Save the gene order
write.table(ordered_genes, paste0("gene_order_", reference_sample, ".txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Ensure to remove these genes from your heatmap_data as well
heatmap_data_clean <- heatmap_data_clean %>%
  filter(Gene %in% ordered_genes)

# Generate heatmaps for all samples using the same gene order in manuscript
for (sample in unique(heatmap_data_clean$Sample)) {
  sample_data <- heatmap_data_clean %>%
    filter(Sample == sample) %>%
    select(Gene, Cluster, ScaledExpression) %>%
    pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # Reorder genes based on the reference sample
  sample_data <- sample_data[ordered_genes, ]
  
  # Plot the heatmap
  svg(paste0("heatmap_", sample, ".svg"), width = 10, height = 8)
  heatmap.2(
    sample_data, scale = "none", col = viridis(20), trace = "none",
    margins = c(8, 8), dendrogram = "none", Rowv = F, Colv = F,
    key = TRUE, key.title = "Scaled Expression", key.xlab = "Expression",
    main = paste("Clustered Heatmap for", sample), cexCol = 1, cexRow = 0.1
  )
  dev.off()
}

# for special clusters deg 
set.seed(42)
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"

sample_list <- SplitObject(integ, split.by = "sample")
rm(integ)

features_to_plot <- read.table("~/shruti/SNRIII/data/SeuratOut/degTableS2F.txt",
                               header = TRUE, fill = TRUE, sep = "\t", quote = "")
features_to_plot <- features_to_plot %>% 
  filter(Level == "Upregulated", Cluster %in% c("4", "5", "7", "12", "15", "17", "18")) %>%
  pull("GeneId")

clusters_of_interest <- c("4", "5", "7", "12", "15", "17", "18")

# Subset samples to desired clusters and calculate average expression
avg_expr_list <- lapply(sample_list, function(x) {
  x <- subset(x, idents = clusters_of_interest)
  AverageExpression(x, assays = "RNA", return.seurat = TRUE)
})

heatmap_data <- lapply(names(avg_expr_list), function(sample) {
  expr <- GetAssayData(avg_expr_list[[sample]], slot = "data") %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    filter(Gene %in% features_to_plot) %>%
    pivot_longer(-Gene, names_to = "Cluster", values_to = "Expression") %>%
    group_by(Gene) %>%
    mutate(ScaledExpression = scale(Expression), Sample = sample) %>%
    ungroup()
}) %>% bind_rows()

heatmap_data <- heatmap_data %>%
  group_by(Gene, Cluster, Sample) %>%
  summarise(ScaledExpression = mean(ScaledExpression), .groups = "drop")

ref_sample <- "kno"
ref_mat <- heatmap_data %>%
  filter(Sample == ref_sample) %>%
  pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

ref_mat <- ref_mat[complete.cases(ref_mat), ]
gene_order <- hclust(dist(ref_mat))$order
ordered_genes <- rownames(ref_mat)[gene_order]

heatmap_data <- heatmap_data %>% filter(Gene %in% ordered_genes)

# Plot heatmaps using the same gene order
for (sample in unique(heatmap_data$Sample)) {
  mat <- heatmap_data %>%
    filter(Sample == sample) %>%
    pivot_wider(names_from = Cluster, values_from = ScaledExpression) %>%
    column_to_rownames("Gene")
  
  # Ensure matrix has all clusters (fill missing ones with NA or 0)
  missing_clusters <- setdiff(clusters_of_interest, colnames(mat))
  if (length(missing_clusters) > 0) {
    mat[missing_clusters] <- NA  # or 0 if you prefer
  }
  
  # Reorder columns by cluster number
  mat <- mat[, sort(colnames(mat))]
  
  # Reorder rows by ordered genes
  mat <- mat[ordered_genes, , drop = FALSE]
  
  # Convert to numeric matrix, handle NA
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0  # Optional: impute missing with 0
  
  # Plot
  svg(paste0("heatmap_", sample, ".svg"), width = 10, height = 8)
  heatmap.2(mat, scale = "none", col = viridis(20), trace = "none",
            margins = c(8, 8), dendrogram = "none", Rowv = FALSE, Colv = FALSE,
            key = TRUE, key.title = "Scaled Expression", key.xlab = "Expression",
            main = paste("Heatmap:", sample), cexCol = 1, cexRow = 0.5)
  dev.off()
}

write.table(ordered_genes, paste0("gene_order_", ref_sample, ".txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# groups <- list(ray = c("5", "7", "18"), fiber = c("4", "12", "15"),
#                vessel = "17", fusiformInitial = c("14", "16"), vesselPrecursor = "6",
#                fiberPrecursor = c("1", "10"), cambium = c("19", "20"),
#                earlyRay = c("3", "8", "9", "11", "13"),
#                unknown = c("0", "2"))
# 
# samples <- sapply(strsplit(colnames(integ), "_"), `[`, 1)
# integ$seurat_clusters <- integ$integrated_snn_res.0.6
# clusters <- as.character(integ$seurat_clusters)
# sample_cluster <- paste(samples, clusters, sep = "_")
# integ <- AddMetaData(integ, metadata = sample_cluster, col.name = "group")
# Idents(integ) <- "group"
# table(Idents(integ))

# specific cluster groups
set.seed(42)
integ <- readRDS("~/shruti/SNRIII/data/SeuratOut/integ.rds")
DefaultAssay(integ) <- "RNA"
Idents(integ) <- "integrated_snn_res.0.6"

sample_list <- SplitObject(integ, split.by = "sample")

features_to_plot <- read.table("~/shruti/SNRIII/data/SeuratOut/degTableS2F.txt",
                               header = TRUE, fill = TRUE, sep = "\t", quote = "")
features_to_plot <- features_to_plot %>% 
  filter(Level == "Upregulated", Cluster %in% c("4", "5", "7", "12", "15", "17", "18")) %>%
  pull("GeneId")

features_to_plot <- up_shared_2_only$fiber_vessel_not_ray

clusters_of_interest <- c("4", "5", "7", "12", "15", "17", "18")
cluster_groups <- list(`5,7,18` = c("5", "7", "18"), `4,12,15` = c("4", "12", "15"),
  `17` = c("17"))

avg_expr_list <- lapply(sample_list, function(sample_obj) {
  expr_data <- GetAssayData(sample_obj, assay = "RNA", slot = "data")
  meta <- sample_obj@meta.data
  meta$cluster <- Idents(sample_obj)
  
  avg_expr_by_group <- lapply(names(cluster_groups), function(group_name) {
    clusters <- cluster_groups[[group_name]]
    cells_in_group <- rownames(meta)[meta$cluster %in% clusters]
    if (length(cells_in_group) > 0) {
      rowMeans(expr_data[, cells_in_group, drop = FALSE])
    } else {
      rep(NA, nrow(expr_data))  # handle missing clusters
    }
  })
  
  avg_expr_mat <- do.call(cbind, avg_expr_by_group)
  colnames(avg_expr_mat) <- names(cluster_groups)
  rownames(avg_expr_mat) <- rownames(expr_data)
  as.data.frame(avg_expr_mat) %>%
    rownames_to_column("Gene")
})

heatmap_data <- lapply(names(avg_expr_list), function(sample) {
  expr <- avg_expr_list[[sample]] %>%
    filter(Gene %in% features_to_plot) %>%
    pivot_longer(-Gene, names_to = "ClusterGroup", values_to = "Expression") %>%
    group_by(Gene) %>%
    mutate(ScaledExpression = scale(Expression), Sample = sample) %>%
    ungroup()
}) %>% bind_rows()

heatmap_data <- heatmap_data %>%
  group_by(Gene, ClusterGroup, Sample) %>%
  summarise(ScaledExpression = mean(ScaledExpression), .groups = "drop")

ref_sample <- "kno"
ref_mat <- heatmap_data %>% filter(Sample == ref_sample) %>%
  pivot_wider(names_from = ClusterGroup, values_from = ScaledExpression) %>%
  column_to_rownames("Gene") %>% as.matrix()

ref_mat <- ref_mat[complete.cases(ref_mat), ]
gene_order <- hclust(dist(ref_mat))$order
ordered_genes <- rownames(ref_mat)[gene_order]

heatmap_data <- heatmap_data %>% filter(Gene %in% ordered_genes)

# Plot heatmaps using the same gene order
for (sample in unique(heatmap_data$Sample)) {
  mat <- heatmap_data %>%
    filter(Sample == sample) %>%
    pivot_wider(names_from = ClusterGroup, values_from = ScaledExpression) %>%
    column_to_rownames("Gene")
  
  # Ensure matrix has all clusters (fill missing ones with NA or 0)
  missing_clusters <- setdiff(clusters_of_interest, colnames(mat))
  if (length(missing_clusters) > 0) {
    mat[missing_clusters] <- NA  # or 0 if you prefer
  }
  
  # Reorder columns by cluster number
  mat <- mat[, sort(colnames(mat))]
  
  # Reorder rows by ordered genes
  mat <- mat[ordered_genes, , drop = FALSE]
  
  # Convert to numeric matrix, handle NA
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0  # Optional: impute missing with 0
  
  # Plot
  svg(paste0("heatmap_", sample, ".svg"), width = 10, height = 8)
  heatmap.2(mat, scale = "none", col = viridis(20), trace = "none",
            margins = c(8, 8), dendrogram = "none", Rowv = FALSE, Colv = FALSE,
            key = TRUE, key.title = "Scaled Expression", key.xlab = "Expression",
            main = paste("Heatmap:", sample), cexCol = 1, cexRow = 0.5)
  dev.off()
}

write.table(ordered_genes, paste0("gene_order_", ref_sample, ".txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# for group wise deg
library(readxl)
library(dplyr)
sheets <- list(ray = 1, vessel = 2, fiber = 3)
deg <- lapply(sheets, function(sheet) read_xlsx("~/shruti/SNRIII/data/SeuratOut/degGroupWiseWilcox.xlsx", sheet = sheet))

get_genes <- function(df, direction = "up") {
  if (direction == "up") {
    return(df %>% filter(avg_log2FC >= 1) %>% pull(gene))
  } else {
    return(df %>% filter(avg_log2FC < 1) %>% pull(gene))
  }
}

gene_lists_up <- lapply(deg, get_genes, direction = "up")
gene_lists_dn <- lapply(deg, get_genes, direction = "dn")

get_common_unique <- function(gene_lists) {
  common <- Reduce(intersect, gene_lists)
  unique <- lapply(names(gene_lists), function(name) {
    setdiff(gene_lists[[name]], unlist(gene_lists[names(gene_lists) != name]))
  })
  names(unique) <- paste0(names(gene_lists), "_only")
  return(list(common = common, unique = unique))
}

up_genes   <- get_common_unique(gene_lists_up)
dn_genes   <- get_common_unique(gene_lists_dn)

get_pairwise_shared <- function(gene_lists) {
  list(
    ray_fiber_not_vessel   = intersect(gene_lists$ray, gene_lists$fiber) %>% setdiff(gene_lists$vessel),
    ray_vessel_not_fiber   = intersect(gene_lists$ray, gene_lists$vessel) %>% setdiff(gene_lists$fiber),
    fiber_vessel_not_ray   = intersect(gene_lists$fiber, gene_lists$vessel) %>% setdiff(gene_lists$ray)
  )
}

up_shared_2_only <- get_pairwise_shared(gene_lists_up)
dn_shared_2_only <- get_pairwise_shared(gene_lists_dn)

# GOTO LINE 356 to plot