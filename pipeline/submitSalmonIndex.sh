#!/bin/bash
set -eu 
in=$(realpath ../data/reference/Potra02_transcripts_GFP.fasta)
out=$(realpath ../data/reference)/salmon_GFP_indices
decoy=$(realpath ../reference/fasta/Potra02_genome.fasta)

[[ ! -d $out ]] && mkdir -p $out

sbatch -A u2022003 -o $out/index.out -e $out/index.err \
--mail-user huge_ashes@yahoo.com $(realpath ../UPSCb-common/pipeline/runSalmonIndex.sh) \
-d $decoy $in $out


 