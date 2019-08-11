#!/bin/sh
#$ -N VAR_CALL
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -v PATH
#$ -v LD_LIBRARY_PATH

#[P] {script runs commands in parallel}

##MAIN SCRIPT
##ENTRYPOINT
##Runs the variant calling pipeline [both {samtools,bcftools},{gatk}]

proc_ids=() ##Used to store process IDs

##SUBSCRIPT-fast_qc.sh[P] [Runs fastqc on samples]
##REQUIRES user validation of output
nohup fast_qc.sh /data4/nextgen2015/pilib/MA5112_Assign3/
wait $!

##Following lines are run after fast_qc.sh output validation
##no need to trim adapters because adapter contamination is < 0.1% according to fastqc

##SUBSCRIPT-bwa.sh[P] [Runs BWA for the samples]
##WAIT until user validation of fast_qc.sh output
nohup bwa.sh
wait $!

##SUBSCRIPT-vc.sh[P] [Variant calling-{samtools,bcftools}]
##DEPENDS on bwa.sh output
nohup vc.sh
proc_ids+=("$!") 

##SUBSCRIPT-vc.gatk.sh[P] [Variant calling[P]-{gatk}]
##DEPENDS on bwa.sh output
nohup vc.gatk.sh HN51 HN60 #sample file prefixes
proc_ids+=("$!")

for i in ${proc_ids[@]}; do
	wait $i
done

echo "Done!"
