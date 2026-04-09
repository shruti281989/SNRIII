# goto doc/
# tail -n +2 timeBinDeg.txt | cut -f1 | sort -u | uniq > genesTime.txt
# tail -n +2 ../data/SeuratOut/output/afterDbltRemoval/degWilcox_lfc1_fdr0.01_pct0.1.txt | cut -f7 | sort -u | uniq > geneSc.txt
# cat geneSc.txt genesTime.txt | sort -u | uniq > allDeg.txt
# rm geneSc.txt genesTime.txt
# grep -F -f allDeg.txt ../data/Potra02_NuclMtCpGenes.gtf > allDeg.gtf

# Run from within R
setwd("~/shruti/SNRIII")
use_python("/usr/bin/python3", required = TRUE)
library(reticulate)

# without promoter
py_run_string("
import pandas as pd

gtf_file = 'doc/allDeg.gtf'
fimo_file = 'data/nreVariant/fimo.tsv'
deg_file = 'doc/allDeg.txt'

# Load DEG list
deg_genes = set()
with open(deg_file) as f:
    for line in f:
        gene = line.strip().split('.')[0]
        deg_genes.add(gene)

# Parse GTF
gene_coords = {}
with open(gtf_file) as f:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue
        fields = line.strip().split('\\t')
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attr = fields[8]
        gene_id = None
        for x in attr.split(';'):
            if 'gene_id' in x:
                gene_id = x.replace('gene_id','').replace('\"','').strip()
        if gene_id and gene_id in deg_genes:
            if gene_id not in gene_coords:
                gene_coords[gene_id] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}
            else:
                gene_coords[gene_id]['start'] = min(gene_coords[gene_id]['start'], start)
                gene_coords[gene_id]['end'] = max(gene_coords[gene_id]['end'], end)

# Load FIMO output
fimo = pd.read_csv(fimo_file, sep='\t', comment='#')
fimo['Start'] = fimo['start'].astype(int)
fimo['End'] = fimo['stop'].astype(int)
fimo['Strand'] = fimo['strand']
fimo['Motif ID'] = fimo['motif_id']
fimo['Sequence Name'] = fimo['sequence_name']
fimo['Matched Sequence'] = fimo['matched_sequence']

# Map motifs to genes
motif_to_gene = []
for idx, row in fimo.iterrows():
    chrom = row['Sequence Name']
    start = row['Start']
    end = row['End']
    strand = row['Strand']
    motif = row['Motif ID']
    matched = row['Matched Sequence']
    for gene_id, info in gene_coords.items():
        if chrom == info['chrom'] and start >= info['start'] and end <= info['end']:
            motif_to_gene.append([gene_id, motif, start, end, strand, matched])
            break

# Save results
df = pd.DataFrame(motif_to_gene, columns=['GeneID','MotifID','Start','End','Strand','Sequence'])
df.to_csv('data/nre/nreVariantInDeg.tsv', sep='\\t', index=False)
print('Mapping complete. Output saved to motif_gene_table.tsv')
")

# with promoter

py_run_string("
import pandas as pd

gtf_file = 'doc/allDeg.gtf'
fimo_file = 'data/nreVariant/fimo.tsv'         
deg_file = 'doc/allDeg.txt' # DEG list

deg_genes = set()
with open(deg_file) as f:
    for line in f:
        gene = line.strip().split('.')[0]
        deg_genes.add(gene)
print(f'Total DEGs loaded: {len(deg_genes)}')

gene_coords = {}
with open(gtf_file) as f:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue
        fields = line.strip().split('\\t')
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attr = fields[8]
        gene_id = None
        for x in attr.split(';'):
            if 'gene_id' in x:
                gene_id = x.replace('gene_id','').replace('\"','').strip()
        if gene_id and gene_id in deg_genes:
            if gene_id not in gene_coords:
                gene_coords[gene_id] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}
            else:
                gene_coords[gene_id]['start'] = min(gene_coords[gene_id]['start'], start)
                gene_coords[gene_id]['end'] = max(gene_coords[gene_id]['end'], end)
print(f'Total DEG genes found in GTF: {len(gene_coords)}')

promoter_coords = {}
promoter_length = 2000

for gene_id, info in gene_coords.items():
    chrom = info['chrom']
    strand = info['strand']
    tss = info['start'] if strand == '+' else info['end']
    if strand == '+':
        start = max(0, tss - promoter_length)
        end = tss - 1
    else:
        start = tss + 1
        end = tss + promoter_length
    promoter_coords[gene_id] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}

fimo = pd.read_csv(fimo_file, sep='\\t', comment='#')
fimo['Start'] = fimo['start'].astype(int)
fimo['End'] = fimo['stop'].astype(int)
fimo['Strand'] = fimo['strand']
fimo['MotifID'] = fimo['motif_id']
fimo['SequenceName'] = fimo['sequence_name']
fimo['MatchedSequence'] = fimo['matched_sequence']

motif_in_promoter = []

for idx, row in fimo.iterrows():
    chrom = row['SequenceName']
    start = row['Start']
    end = row['End']
    strand = row['Strand']
    motif = row['MotifID']
    matched = row['MatchedSequence']
    
    for gene_id, info in promoter_coords.items():
        if chrom == info['chrom'] and start >= info['start'] and end <= info['end']:
            motif_in_promoter.append([gene_id, motif, start, end, strand, matched])
            break

df = pd.DataFrame(motif_in_promoter,
                  columns=['GeneID','MotifID','Start','End','Strand','Sequence'])
df.to_csv('data/nre/nreVariantInDeg_promoter.tsv', sep='\\t', index=False)
print('Mapping complete. Output saved to DEG_with_NRE_promoter_hits.tsv')
")

# 
library(txdbmaker)
library(ChIPseeker)
txdb <- makeTxDbFromGFF("/mnt/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff",format="gff3")

# 1. if using the mac2 outputs directly done in different ways
input_dir <- "data/analysis/additionalAnalysis/macs2_nucChr"
output_dir <- "data/annotated_results_macs2_noNuChr"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

peak_files <- list.files(input_dir, pattern = "\\.narrowPeak$", 
                         full.names = TRUE)

annotate_and_save <- function(file_path, output_dir, txdb) {
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  annotated_peak <- annotatePeak(file_path, tssRegion = c(-3000, 3000), 
                                 TxDb = txdb, flankDistance = 5000, 
                                 sameStrand = FALSE)
  
  annotated_df <- as.data.frame(annotated_peak)
  output_file <- file.path(output_dir, paste0(sample_name, "_annotated.tsv"))
  write.table(annotated_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  return(annotated_df)
}

annotated_dfs <- lapply(peak_files, function(file) {
  annotate_and_save(file, output_dir, txdb)})

# 2. Look at one sample at a time.
library(tidyverse)
library(dplyr)
library(ChIPseeker)
library(ggplot2)
library(viridis)
library(txdbmaker)
setwd("/mnt/picea/projects/aspseq/htuominen/ChipSeqERF85OE/analysis/additionalAnalysis/macs2_nucChr/idr")

# original macs narrow peaks
# X173macs <- read_tsv("Xylem-173_X1_peaks.narrowPeak", show_col_types = F, col_names = c("chr173","start173","stop173","name173","score173","strand173","signalValue173","pValue173","qValue173","summit173"))
# V8 has the p-value and V9 has the q-value

txdb <- makeTxDbFromGFF("/mnt/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff",format="gff3")

samples <- c("173", "206", "349")
input_paths <- paste0("Xylem-", samples, "_X1_peaks.narrowPeak")
output_paths <- paste0("/mnt/picea/home/schoudhary/shruti/chipERF/data/", samples, ".tsv")

input_paths <- paste0("Phloem-", samples, "_peaks.narrowPeak")
output_paths <- paste0("/mnt/picea/home/schoudhary/shruti/chipERF/data/", samples, ".tsv")

for (i in seq_along(samples)) {
  annot <- annotatePeak(input_paths[i], tssRegion = c(-3000, 3000), 
                        TxDb = txdb, flankDistance = 5000, sameStrand = FALSE)
  annot_df <- as.data.frame(annot)
  annot_df$real_p_value <- 10^(-annot_df$V8)
  annot_df$real_q_value <- 10^(-annot_df$V9)
  # annot_df <- annot_df[, c("annotation", "geneId", "real_p_value", "real_q_value")]
  # annot_df <- annot_df[annot_df$annotation != "Distal Intergenic" & annot_df$real_p_value < 0.05, ]
  write.table(annot_df, output_paths[i], sep = "\t", row.names = FALSE, quote = FALSE)
}

# 3. for consensusIDR:
annot <- annotatePeak("Phloem-206-349-idrP_005.tsv", 
                      tssRegion=c(-3000, 3000), TxDb=txdb,
                      flankDistance = 5000, sameStrand = FALSE)
annot <- annot %>% as.data.frame()
write.table(annot,"~/shruti/chipERF/data/idrPh206-349_p0.05.tsv", quote = F,
            row.names = F, sep ="\t", col.names = F)
annot <- annotatePeak("consensusPeaks_idr.tsv", 
                      tssRegion=c(-3000, 3000), TxDb=txdb,
                      flankDistance = 5000, sameStrand = FALSE)
annot <- annot %>% as.data.frame()
write.table(annot,"~/shruti/chipERF/data/consensusPeaks_XyIdr.tsv", quote = F,
            row.names = F, sep ="\t", col.names = T)

annotate_all_tsv <- function(input_dir, output_dir, txdb) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  tsv_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  annotate_and_save <- function(file_path, output_dir, txdb) {
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    annot <- annotatePeak(file_path, tssRegion = c(-3000, 3000), 
                          TxDb = txdb, flankDistance = 5000, sameStrand = FALSE)
    annot_df <- as.data.frame(annot)
    output_file <- file.path(output_dir, paste0(sample_name, "_annotated.tsv"))
    write.table(annot_df, output_file, quote = FALSE, row.names = FALSE, sep = "\t", col.names = TRUE)
    return(annot_df)
  }
  annotated_dfs <- lapply(tsv_files, function(file) {
    annotate_and_save(file, output_dir, txdb)
  })
}

# Example usage
input_directory <- "."
output_directory <- "~/shruti/chipERF/data/annotated/"
annotate_all_tsv(input_directory, output_directory, txdb)
