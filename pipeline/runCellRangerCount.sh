#!/bin/bash -l
#SBATCH -p node
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem 120G

## stop on error and be verbose in the output
set -eux

singi=$1
foldname=$2
reads=$3
sample_id=$4
ref_transc=$5

# vars
# OPTIONS="--chemistry=SC3Pv3LT"

# run the command
singularity exec $singi cellranger count \
--id=$foldname --fastqs=$reads \
--sample=$sample_id --transcriptome=$ref_transc \

# $OPTIONS

