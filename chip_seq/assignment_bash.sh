#!/bin/bash
if [ $# -eq 0 ]
  then
    echo "Script needs a Chr(fa) file, Chip(fastq) and control(fastq) file"
    exit
fi
chr=$1
chip=$2
control=$3
time fastqc $chip
read
time fastqc $control
read
module load bowtie
module load samtools
bowtie2-build
head $control
read
head $chip
read
head $chr
read
bowtie2-build $chr hg19
read
time bowtie2 -x hg19 -U $control -S input.sam
read
time bowtie2 -x hg19 -U $chip -S chip.sam
read
diff -y input.fastq input.sam | head
read
samtools
samtools view
time samtools view -Sb input.sam > input.bam
time samtools view -Sb chip.sam > chip.bam
samtools rmdup input.bam input.nodup.bam
samtools rmdup chip.bam chip.nodup.bam
samtools sort 
time samtools sort chip.nodup.bam chip.nodup.sorted
time samtools sort input.nodup.bam input.nodup.sorted
time samtools index input.nodup.sorted.bam
time samtools index chip.nodup.sorted.bam
time samtools flagstat chip.nodup.sorted.bam > chip.flagstat.txt
time samtools flagstat input.nodup.sorted.bam > input.flagstat.txt
diff -y chip.flagstat.txt input.flagstat.txt
read
macs2
macs2 callpeak
time macs2 callpeak -t chip.nodup.sorted.bam -c input.nodup.sorted.bam -f BAM -g hs -n macs_out --call-summits -B
read
scp chip.nodup.sorted.* input.nodup.sorted.* chip.bam input.bam macs_out* skumar@smgate.nuigalway.ie:/home/skumar/assignment/
