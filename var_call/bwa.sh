#!/bin/bash

#module load bwa

##SUBSCRIPT
##ENTRYPOINT
##Runs bwa on samples
##WAIT until user has validated fast_qc.sh output

proc_ids=() ##Used to store process IDs
##Runs bwa mem in background on HN51 sequences
nohup bwa mem -t 16 -R '@RG\tID:HN51\tSM:BD1LYPACXX.lane_3_I12\tLB:Lib1'  /data4/pilib/genomes/hs/hg19/hg19.fa /data4/nextgen2015/pilib/MA5112_Assign3/HN51_S2_normal.BD1LYPACXX.lane_3_P1_I12.hg19.sequence.fastq /data4/nextgen2015/pilib/MA5112_Assign3/HN51_S2_normal.BD1LYPACXX.lane_3_P2_I12.hg19.sequence.fastq | samtools view -Sb - > HN51.bam&
proc_ids+=("$!") ##Store process ID
##Runs bwa mem in bg on HN60 sequences
nohup bwa mem -t 16 -R '@RG\tID:HN60\tSM:BD1LYPACXX.lane_7_I16\tLB:Lib2' /data4/pilib/genomes/hs/hg19/hg19.fa /data4/nextgen2015/pilib/MA5112_Assign3/HN60_s2_normal.BD1LYPACXX.lane_7_P1_I16.hg19.sequence.fastq /data4/nextgen2015/pilib/MA5112_Assign3/HN60_s2_normal.BD1LYPACXX.lane_7_P2_I16.hg19.sequence.fastq | samtools view -Sb - > HN60.bam&
proc_ids+=("$!") ##Store process ID

##Wait until all processes are complete
for id in ${proc_ids[@]}; do
  wait $id
done

