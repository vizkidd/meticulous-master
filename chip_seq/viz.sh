##!/bin/bash

##Functions
##Gimmick
show_spinner()
{
  local -r pid="${1}"
  local -r delay='0.75'
  local spinstr='\|/-'
  local temp
  while ps a | awk '{print $1}' | grep -q "${pid}"; do
    temp="${spinstr#?}"
    printf " [%c]  " "${spinstr}"
    spinstr=${temp}${spinstr%"${temp}"}
    sleep "${delay}"
    printf "\b\b\b\b\b\b"
  done
  printf "    \b\b\b\b"
}

file_checks()
{

##Check file extensions
##PATTERN MATCHING
if [[ ${files[0]} != *.fa ]]; then
        echo "Give proper Reference File"
        exit 111
fi

##COMMAND SUBSTITUTION
for i in $(seq 1 1 2);
do
if [[ ${files[$i]} != *.fq && ${files[$i]} != *.fastq ]]; then
        echo "Give proper fastq file for ${files[$i]}"
        exit 112
fi
done

##Check if all files exist and are readable
filechk=()
for i in ${files[@]} ;
do
##   echo $i
   [[ -f $i && -r $i ]] ##EXPRESSION TO CHECK IF FILE EXISTS AND IS READABLE
   filechk+=($?)
done

##filechk array contains the exit codes, if it has 1 it means some file 
##failed to read. We perform this check here by doing a simple
##'if filechk contains 1'

if [[ " ${filechk[@]} " =~ "1" ]]; then
   echo "Some files are not readable or do not exist. Exiting..."
   exit 121 ##Use different exit codes for different errors
else
   echo "All files exist and are readable"
fi
}

##CODE ENTRYPOINT
##Check parameter count
if [ $# != 3 ]
  then
    echo "Script needs a Chr(fa) file, Chip(fastq) and control(fastq) file"
    exit 111
fi

chr=$1
chip=$2
control=$3

files=($chr $chip $control)

file_checks ##${files[@]}

##Put .fastq files in an array
fastq=($chip $control)

module load fastqc
[ -x "$(command -v fastqc)" ]
if [ $? != 0 ]; then
	echo "fastqc not found!"
	exit 211
else
	echo "Starting QC on data"
fi

for item in ${fastq[@]};
do
  nohup fastqc $item &> fastqc.out&
  if [ $? != 0 ];
  then
	echo "Error in fastqc step. Exiting with status $?"
	exit $?
  fi
  show_spinner $!
  wait $!
done
echo "Writing output to fastqc.out"

##Making sure bowtie & samtools exist
module load bowtie
[ -x "$(command -v bowtie2)" ]
if [ $? != 0 ]; then
        echo "bowtie2 not found!"
        exit 311
fi
module load samtools
[ -x "$(command -v samtools)" ]
if [ $? != 0 ]; then
        echo "samtools not found!"
        exit 311
fi

##View data, commented to remove clutter
##less $control
##read -p "Press Enter to continue"
##less $chip
##read -p "Press Enter to continue"
##less $chr
##read -p "Press Enter to continue"

#Parallelizing bowtie because it's cool
echo "Building index"
nohup time bowtie2-build $chr hg19 &> bowtie-build.out& 
nbb=$!
echo "Starting sub-process [$nbb]"

##Fancy gimmicks
show_spinner "$nbb"

wait $nbb 

echo "Starting alignment"
nohup time bowtie2 -x hg19 -U $control -S input.sam &> bowtie-control.out& 
nbco=$!
nohup time bowtie2 -x hg19 -U $chip -S chip.sam &> bowtie-chip.out&
nbch=$!
##diff -y input.fastq input.sam | head
echo "Writing output to bowtie-chip.out & bowtie-control.out"

##Compressing sam files to bam files & post-processing
echo "Starting sub-process [$nbco]"
echo "Starting sub-process [$nbch]"

show_spinner "$nbco"

wait $nbco ##wait until process is complete because other functions 
##require the index to be built and aligned
echo "Compressing .sam files to .bam files"
samtools view -Sb input.sam > input.bam

show_spinner "$nbch"

wait $nbch
samtools view -Sb chip.sam > chip.bam
echo "Removing duplicates"
samtools rmdup input.bam input.nodup.bam
samtools rmdup chip.bam chip.nodup.bam 
echo "Sorting .bam files"
samtools sort chip.nodup.bam chip.nodup.sorted
samtools sort input.nodup.bam input.nodup.sorted
echo "Indexing"
samtools index input.nodup.sorted.bam
samtools index chip.nodup.sorted.bam
echo "Estimating stats of sorted .bam files"
samtools flagstat chip.nodup.sorted.bam > chip.flagstat.txt
samtools flagstat input.nodup.sorted.bam > input.flagstat.txt
##diff -y chip.flagstat.txt input.flagstat.txt
echo "Peak calling using macs2"
nohup time macs2 callpeak -t chip.nodup.sorted.bam -c input.nodup.sorted.bam -f BAM -g hs -n macs_out --call-summits -B &> macs2.out&
macsproc=$!
show_spinner $!
wait $macsproc
echo "Writing output to macs2.out"
echo "Finished..."
read -p "ENTER smgate username:" usr
scp chip.nodup.sorted.* input.nodup.sorted.* chip.bam input.bam macs_out* *.out "$usr"@smgate.nuigalway.ie:/home/"$usr"/assignment/
