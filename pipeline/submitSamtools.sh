#!/bin/bash

sbatch -J kcl ./runSamtools.sh /pfs/proj/nobackup/fs/projnb10/hpc2nstor2025-047/shruti/temp_scRNA/SNR-results/snrIII/clRngrCntNuclMtCp/kcl2/outs/possorted_genome_bam.bam kcl.txt
sbatch -J kno ./runSamtools.sh /pfs/proj/nobackup/fs/projnb10/hpc2nstor2025-047/shruti/temp_scRNA/SNR-results/snrIII/clRngrCntNuclMtCp/kno2/outs/possorted_genome_bam.bam kno.txt
