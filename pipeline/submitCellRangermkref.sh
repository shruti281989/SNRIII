#!/bin/bash
set -eu 

proj=u2022027
mail=huge_ashes@yahoo.com

# singularity
export SINGULARITY_BINDPATH="/mnt:/mnt"
singc=$(realpath ../singularity/cellranger_7.0.0.sif)

gname=Potra02_genome

gnome=$(realpath ../data/reference/Potra02_genome.fasta)
gtf=$(realpath ../data/reference/Potra02_genes.gtf)

out=$(realpath ../data/reference/cellrangeref)

if [ ! -d $out ]; then
    mkdir -p $out
fi

sbatch -A $proj -e $out/index.err -o /$out/index.out --mail-user $mail \
./runCellRangermkref.sh $singc $gname $gnome $gtf
 