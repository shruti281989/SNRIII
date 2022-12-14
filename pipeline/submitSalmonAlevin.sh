#!/bin/bash
set -eu 

in=$(realpath ../data/raw)
#fwd=$(realpath ../data/raw/P24151_1001_S1_L002_R1_001.fastq.gz)
#rev=$(realpath ../data/raw/P24151_1001_S1_L002_R2_001.fastq.gz)
inx=$(realpath ../data/reference/salmon_GFP_indices)
tgmap=$(realpath ../data/reference/tx2gene.tsv)
out=$(realpath ../data)/alevin_output/P24151

[[ ! -d $out ]] && mkdir -p $out
for f in $(find $in -name "*_R1_001.fastq.gz"); do 
  sname=$(basename ${f/_R1_001.fastq.gz/})
  sbatch -A u2022003 -o $out/$sname.out -e $out/$sname.err \
  --mail-user huge_ashes@yahoo.com $(realpath ../UPSCb-common/pipeline/runAlevin.sh) \
  $inx $f ${f/_R1_001.fastq.gz/_R2_001.fastq.gz} $out $tgmap

done 