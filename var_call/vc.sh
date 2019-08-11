#!/bin/sh
#$ -N SAMTOOLS
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -v PATH
#$ -v LD_LIBRARY_PATH

#module load samtools
#module load bcftools

##FUNCTIONS:
parallel_samtools() ##Runs variant calling in parallel
{
##No need to index reference seqs using samtools faidx because they're already indexed
##fixmate checks the two mates from a paired-end bam (name sorted) and then updates the flags and insert sizes
samtools fixmate $1.bam $1.fixed.bam 
samtools sort $1.fixed.bam $1.sorted ##Sort bam file
samtools rmdup $1.sorted.bam $1.sorted.rmdup.bam ##Remove duplicates
samtools index $1.sorted.rmdup.bam ##Index bam files
##multi-way 'pile' seqs based on read groups
bcftools mpileup -O u --threads 16 -f /data4/pilib/genomes/hs/hg19/hg19.fa $1.sorted.rmdup.bam | bcftools call -m  > $1.raw.vcf
#coverage depth formula is : (read count*Avg read length)/total bases
bcftools view $1.raw.vcf | vcfutils.pl varFilter -d 10 > $1.flt.vcf ##Filter variants
}

##SUBSCRIPT
##ENTRYPOINT
##Runs variant calling
##WAIT for bwa.sh output

proc_ids=() ##Used to store process IDs

##Alternate way of saying "nohup" because nohup doesn't directly support calling functions in a script
( trap "true" HUP ; parallel_samtools HN51 ) &> ./HN51.sam.out &
proc_ids+=("$!")
( trap "true" HUP ; parallel_samtools HN60 ) &> ./HN60.sam.out &
proc_ids+=("$!")

##Wait until process completion
for i in ${proc_ids[@]}; do
        wait $i
done

