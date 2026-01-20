setwd("/mnt/picea/home/schoudhary/shruti/SNRIII/")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyverse)
  library(tidyr)
  library(pheatmap)
  library(viridis)
  library(gplots)
  library(UpSetR)
  library(ComplexHeatmap)
  library(readxl)
})
set.seed(42)

# upset plots for deg set (Figure in the manuscript)
deg <- read.table("data/SeuratOut/output/afterDbltRemoval/degWilcox_lfc1_fdr0.01_pct0.1.txt", header = T)
deg <- deg %>% select("cluster", "avg_log2FC", "gene")

up = list(
  g1 = deg %>% filter(cluster=="Cluster_1" & avg_log2FC > 0),
  g3 = deg %>% filter(cluster=="Cluster_3" & avg_log2FC > 0),
  g4 = deg %>% filter(cluster=="Cluster_4" & avg_log2FC> 0),
  g5 = deg %>% filter(cluster=="Cluster_5" & avg_log2FC> 0),
  g6 = deg %>% filter(cluster=="Cluster_6" & avg_log2FC> 0),
  g7 = deg %>% filter(cluster=="Cluster_7" & avg_log2FC> 0),
  g8 = deg %>% filter(cluster=="Cluster_8" & avg_log2FC> 0),
  g8 = deg %>% filter(cluster=="Cluster_9" & avg_log2FC> 0),
  g10 = deg %>% filter(cluster=="Cluster_10" & avg_log2FC>0),
  g11 = deg %>% filter(cluster=="Cluster_11" & avg_log2FC>0),
  g12 = deg %>% filter(cluster=="Cluster_12" & avg_log2FC>0),
  g13 = deg %>% filter(cluster=="Cluster_13" & avg_log2FC>0),
  g14 = deg %>% filter(cluster=="Cluster_14" & avg_log2FC>0),
  g15 = deg %>% filter(cluster=="Cluster_15" & avg_log2FC>0),
  g16 = deg %>% filter(cluster=="Cluster_16" & avg_log2FC>0),
  g17 = deg %>% filter(cluster=="Cluster_17" & avg_log2FC>0),
  g18 = deg %>% filter(cluster=="Cluster_18" & avg_log2FC>0),
  g19 = deg %>% filter(cluster=="Cluster_19" & avg_log2FC>0),
  g20 = deg %>% filter(cluster=="Cluster_20" & avg_log2FC>0))

dn = list(
  g1 = deg %>% filter(cluster=="Cluster_1" & avg_log2FC< 0),
  g3 = deg %>% filter(cluster=="Cluster_3" & avg_log2FC < 0),
  g4 = deg %>% filter(cluster=="Cluster_4" & avg_log2FC< 0),
  g5 = deg %>% filter(cluster=="Cluster_5" & avg_log2FC< 0),
  g6 = deg %>% filter(cluster=="Cluster_6" & avg_log2FC< 0),
  g7 = deg %>% filter(cluster=="Cluster_7" & avg_log2FC< 0),
  g8 = deg %>% filter(cluster=="Cluster_8" & avg_log2FC< 0),
  g8 = deg %>% filter(cluster=="Cluster_9" & avg_log2FC> 0),
  g10 = deg %>% filter(cluster=="Cluster_10" & avg_log2FC< 0),
  g11 = deg %>% filter(cluster=="Cluster_11" & avg_log2FC< 0),
  g12 = deg %>% filter(cluster=="Cluster_12" & avg_log2FC< 0),
  g13 = deg %>% filter(cluster=="Cluster_13" & avg_log2FC< 0),
  g14 = deg %>% filter(cluster=="Cluster_14" & avg_log2FC< 0),
  g15 = deg %>% filter(cluster=="Cluster_15" & avg_log2FC< 0),
  g16 = deg %>% filter(cluster=="Cluster_16" & avg_log2FC< 0),
  g17 = deg %>% filter(cluster=="Cluster_17" & avg_log2FC< 0),
  g18 = deg %>% filter(cluster=="Cluster_18" & avg_log2FC< 0),
  g19 = deg %>% filter(cluster=="Cluster_19" & avg_log2FC< 0),
  g20 = deg %>% filter(cluster=="Cluster_20" & avg_log2FC< 0))

# upset plots in figure 5C and 5D
dnlt <- lapply(dn[c(1:12,14,15)], function(df) df$gene)
mdown = make_comb_mat(dnlt)
UpSet(mdown)

uplt <- lapply(up[c(1:18)], function(df) df$gene)
mup = make_comb_mat(uplt)
UpSet(mup)
