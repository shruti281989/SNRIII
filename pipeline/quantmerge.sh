module load bioinfo-tools Salmon

srun -A u2022027 -w picea -t 2:00:00 salmon quantmerge --quants {P27752_1001_S1_L001_R2_001,P27752_1001_S1_L002_R2_001} --names {kcl1,kcl2} --output kcl_quant.sf
srun -A u2022027 -w picea -t 2:00:00 salmon quantmerge --quants {P27752_1002_S2_L001_R2_001,P27752_1002_S2_L002_R2_001} --names {kno1,kno2} --output kno_quant.sf

srun -A u2022027 -w picea -t 2:00:00 cat kcl/P27752_1001_S1_L001_R2_001.fastq.gz kcl/P27752_1001_S1_L002_R2_001.fastq.gz > kcl.fq.gz
srun -A u2022027 -w picea -t 2:00:00 cat kno/P27752_1002_S2_L001_R2_001.fastq.gz kno/P27752_1002_S2_L002_R2_001.fastq.gz > kno.fq.gz
