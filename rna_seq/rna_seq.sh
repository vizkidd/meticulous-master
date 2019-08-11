#!/bin/sh
#$ -N RNA_SEQ_HISAT_STIE_BGOWN
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -v PATH
#$ -v LD_LIBRARY_PATH

if [ $# -eq 0 ]
  then
    echo "Script needs a file with assenscion numbers"
    exit
fi
module load samtools
asc_file=$1
line_count=$(wc -l < $1)
echo $line_count 
#hisat2_extract_splice_sites.py data/Mus_musculus.GRCm38.95.gtf > mus.ss
#hisat2_extract_exons.py data/Mus_musculus.GRCm38.95.gtf > mus.exon
hisat2-build --ss mus.ss --exon mus.exon -p 8 -f data/Mus_musculus.GRCm38.dna.primary_assembly.fa data/indexes/chr_trans
echo "Genome Indexed"
for ASN in $(< $asc_file); do
    echo "Processing $ASN"
    #hisat2 -p 8 --dta -x data/indexes/chr_tran -1 data/fastq/"$ASN"_1.fastq -2 data/samples/"$ASN"_chrX_2.fastq.gz -S "$ASN"_chrX.sam
    hisat2 -p 8 --dta -x data/indexes/chr_trans --sra-acc "$ASN" -S "$ASN"_chr.sam
    echo "Completed Hisat2"
    samtools view -bS "$ASN"_chr.sam >"$ASN"_chr_unsorted.bam
    samtools sort "$ASN"_chr_unsorted.bam "$ASN"_chr
    echo "Completed Samtools"
    stringtie -p 8 -l "$ASN" -G data/Mus_musculus.GRCm38.95.gtf -o "$ASN"_chr.gtf "$ASN"_chr.bam
    echo "Completed stringtie"
    echo "$ASN"_chr.gtf>>mergelist.txt
    wait
done
stringtie --merge -p 8 -G data/Mus_musculus.GRCm38.95.gtf -o stringtie_merged.gtf mergelist.txt
gffcompare -r data/Mus_musculus.GRCm38.95.gtf –G -o merged stringtie_merged.gtf
for ASN in $(< $asc_file); do
	stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/"$ASN"/"$ASN"_chr.gtf "$ASN"_chr.bam
done
rm *.sam *unsorted.bam 
