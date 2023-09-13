#!/bin/bash

## be verbose and print
set -eux

# singularity
export SINGULARITY_BINDPATH="/mnt:/mnt"
# singc=$(realpath ../../single_cell_analysis_poplar/singularity/cellranger_7.0.0.sif)

singc=/mnt/picea/home/schoudhary/shruti/single_cell_analysis_poplar/singularity/cellranger_7.0.0.sif

# proj=u2022027
# mail=huge_ashes@yahoo.com
# 
# out=$(realpath ../data)/CellRangerCount
# 
# if [ ! -d $out ]; then
#     mkdir -p $out
# fi

# kclfolder=P27752_1001
# kclfolder=kcl0

# knofolder=P27752_1002
knofolder=kno2

# kclreads=$(realpath ../data/raw/kcl)
# kclreads=/mnt/picea/home/schoudhary/shruti/SNRIII/data/raw/kcl

# knoreads=$(realpath ../data/raw/kno)
knoreads=/mnt/picea/home/schoudhary/shruti/SNRIII/data/raw/kno

# kclsample=P27752_1001
knosample=P27752_1002

# ref_transc=$(realpath ../../single_cell_analysis_poplar/data/reference/cellrangeref/Potra02_genome)
# ref_transc=/mnt/picea/home/schoudhary/shruti/SNRIV/Potra02_genome
ref_transc=/mnt/picea/home/schoudhary/shruti/SNR-u2023011/analysis/reference/Potra02_NuclMtCp

# sbatch -A $proj -p node -w picea --mem 120G --mail-user=$mail \
# -o $out/kcl.out -e $out/kcl.err \

# ~/shruti/SNRIII/pipeline/runCellRangerCount.sh \
# $singc $kclfolder $kclreads $kclsample $ref_transc

~/shruti/SNRIII/pipeline/runCellRangerCount.sh \
$singc $knofolder $knoreads $knosample $ref_transc