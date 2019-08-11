#!/bin/bash
salmon index -t data/Mus_musculus.GRCm38.cdna.all.fa -i data/chr_index
for i in *_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
salmon quant -i data/chr_index -l A -1 ${prefix}_R1.fastq.gz -2 
${prefix}_R2.fastq.gz -o quant/${prefix};
done
