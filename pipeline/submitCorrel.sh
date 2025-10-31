#!/bin/bash
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=10GB
#SBATCH --account=u2023011

cd ~/shruti/SNR-u2023011/analysis/snrIII/dnStrm/correl/

apptainer exec -B /mnt:/mnt -B /mnt/picea/Modules/apps/compilers/R/4.3.1/lib/R/library:/usr/local/lib/R/library /mnt/picea/storage/singularity/R-4.3.1.sif Rscript --vanilla correl.R 