setwd("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/snrIII/dnStrm/correl/")

suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

exprAll <- readRDS("/mnt/picea/home/schoudhary/shruti/SNR-u2023011/shruti/SNR-u2023011/analysis/snrIII/dnStrm/correl/exprAllAspwdctrl.rds")
exprAll_stat <- suppressWarnings(sapply(32:ncol(exprAll), function(i) sapply(1:31, function(j) cor.test(exprAll[,i],exprAll[,j],method = "pearson")[c(3,4)])))
saveRDS(exprAll_stat, file ="exprAllctrl_stat.rds")