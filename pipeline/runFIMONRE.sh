#!/bin/bash -l
#SBATCH -p main
#SBATCH -t 12:00:00
#SBATCH -A u2022003
#SBATCH -n 1

## stop on error and be verbose in the output
set -eux

export SINGULARITY_BINDPATH="/mnt:/mnt"

cd /mnt/ada/projects/aspseq/htuominen/SingleCellSeqSNRIII/miniEx/meme

# mkdir -p /mnt/picea/home/schoudhary/shruti/SNRIII/data/nre
mkdir -p /mnt/picea/home/schoudhary/shruti/SNRIII/data/nreVariant

# apptainer exec -B "/mnt:/mnt" meme.sif fimo --max-stored-scores 10000 \
# --oc /mnt/picea/home/schoudhary/shruti/SNRIII/data/nre \
# /mnt/picea/home/schoudhary/shruti/SNRIII/doc/nre.meme Potra02_genome.fasta

apptainer exec -B "/mnt:/mnt" meme.sif fimo --max-stored-scores 10000 \
--oc /mnt/picea/home/schoudhary/shruti/SNRIII/data/nreVariant --thresh 1e-3 \
/mnt/picea/home/schoudhary/shruti/SNRIII/doc/nreVariant.meme Potra02_genome.fasta