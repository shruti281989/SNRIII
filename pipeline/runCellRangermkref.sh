#!/bin/bash -l
#SBATCH -p node
#SBATCH -t 05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem 120G

## stop on error and be verbose in the output
set -eux

singi=$1
ref_dir=$2
assembly=$3
gtf=$4

# run the command
singularity exec $singi cellranger mkref \
--genome=$ref_dir \
--fasta=$assembly \
--genes=$gtf


