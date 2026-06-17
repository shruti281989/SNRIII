#!/bin/bash
#SBATCH -A hpc2n2026-155
#SBATCH -n 4
#SBATCH -t 12:00:00

ml GCC/14.2.0 SAMtools/1.22 

samtools coverage $1 > $2
samtools flagstat -@ 3 $1 >> $2
