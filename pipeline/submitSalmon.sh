#!/bin/bash
set -eu 

# singularity
export SINGULARITY_BINDPATH="/mnt:/mnt"
sings=$(realpath ../../single_cell_analysis_poplar/singularity/salmon-1.6.0.simg)

# in=$(realpath ../data/raw)
# in=$(realpath ../data/raw/kcl.fq.gz)
in=$(realpath ../data/raw/kno.fq.gz)
inx=$(realpath ../../single_cell_analysis_poplar/reference/indices/salmon1.6.0/)
# out=$(realpath ../data)/salmonQuant_out
out=$(realpath ../data)/salmon_cat/kno
# out=$(realpath ../data)/salmon_cat/kcl

[[ ! -d $out ]] && mkdir -p $out

sbatch -A u2022027 -o $out/kno.out -e $out/kno.err \
--mail-user huge_ashes@yahoo.com \
$(realpath ../../single_cell_analysis_poplar/UPSCb-common/pipeline/runSalmonSE.sh) \
$inx $in $out 250 25

# sbatch -A u2022003 -o $out/kcl.out -e $out/kcl.err \
# --mail-user huge_ashes@yahoo.com \
# $(realpath ../../single_cell_analysis_poplar/UPSCb-common/pipeline/runSalmonSE.sh) \
# $inx $in $out 250 25

# for f in $(find $in -name "*_R2_001.fastq.gz"); do
# 
# sname=$(basename ${f/_R2_001.fastq.gz/})
# sbatch -A u2022003 -o $out/$sname.out -e $out/$sname.err \
# --mail-user huge_ashes@yahoo.com \
# $(realpath ../UPSCb-common/pipeline/runSalmonSE.sh) \
# $inx $f $out 250 25
#   
# done
