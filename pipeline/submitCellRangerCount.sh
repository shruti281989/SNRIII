#!/bin/bash

## be verbose and print
set -eux

# singularity
export SINGULARITY_BINDPATH="/mnt:/mnt"
singc=/mnt/picea/home/schoudhary/shruti/single_cell_analysis_poplar/singularity/cellranger_7.0.0.sif

# proj=u2022027
out=$(realpath ../data)/CellRangerCount

if [ ! -d $out ]; then
    mkdir -p $out
fi

# Run the script one by one for each sample
# kclfolder=kcl
knofolder=kno

# kclreads=/mnt/picea/home/schoudhary/shruti/SNRIII/data/raw/kcl
knoreads=/mnt/picea/home/schoudhary/shruti/SNRIII/data/raw/kno

# kclsample=P27752_1001
knosample=P27752_1002

# Use potra genome as T89 genome is not well annotated
ref_transc=/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/reference/Potra02_NuclMtCp
# ref_transc=/mnt/ada/projects/aspseq/htuominen/SNR-results/potrxref/Potrx01

# sbatch -A $proj -p node -w picea --mem 120G \
# -o $out/kcl.out -e $out/kcl.err \

# ~/shruti/SNRIII/pipeline/runCellRangerCount.sh \
# $singc $kclfolder $kclreads $kclsample $ref_transc

~/shruti/SNRIII/pipeline/runCellRangerCount.sh \
$singc $knofolder $knoreads $knosample $ref_transc

# The output is generated in the folder from wgere the script is run
# Move the output files to the desired location