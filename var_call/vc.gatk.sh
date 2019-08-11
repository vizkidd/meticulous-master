#!/bin/bash
#$ -N GATK
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -v PATH
#$ -v LD_LIBRARY_PATH

#module load jdk
#module load GATK or use custom GATK

##SUBSCRIPT
##ENTRYPOINT
##Runs gatk variant calling
##WAIT for bwa.sh output

proc_ids=()

##Runs gatk Haplotype Calling
##nohup gatk HaplotypeCaller --native-pair-hmm-threads 16 -R /data4/pilib/genomes/hs/hg19/hg19.fa -I $1.sorted.rmdup.bam -O 
##$1.gatk.vcf -ERC GVCF -L $1.flt.vcf &> $1.gatk.out&
nohup gatk HaplotypeCaller --native-pair-hmm-threads 16 -R /data4/pilib/genomes/hs/hg19/hg19.fa -I $1.sorted.rmdup.bam -O "$1"_gatk.vcf -ERC GVCF &> "$1"_gatk.out&
proc_ids+=("$!")
##nohup gatk HaplotypeCaller --native-pair-hmm-threads 16 -R /data4/pilib/genomes/hs/hg19/hg19.fa -I $2.sorted.rmdup.bam -O 
##$2.gatk.vcf -ERC GVCF -L $2.flt.vcf &> $2.gatk.out&
nohup gatk HaplotypeCaller --native-pair-hmm-threads 16 -R /data4/pilib/genomes/hs/hg19/hg19.fa -I $2.sorted.rmdup.bam -O "$2"_gatk.vcf  -ERC GVCF &> "$2"_gatk.out&
proc_ids+=("$!")

##Wait until process completion
for id in ${proc_ids[@]}; do
  wait $id
done

##MERGE VCF files
##CANT use DBImport needs intervals
##gatk GenomicsDBImport --reader-threads 16 -V "$1"_gatk.vcf   -V "$2"_gatk.vcf --genomicsdb-workspace-path . -L $1.flt.vcf -L 
##$2.flt.vcf
#gatk CombineGVCFs -R /data4/pilib/genomes/hs/hg19/hg19.fa -V 
#"$1"_gatk.vcf   -V "$2"_gatk.vcf -O HN_merged.vcf
##gatk GenomicsDBImport --reader-threads 16 -V $1.gatk.vcf   -V $2.gatk.vcf --genomicsdb-workspace-path . -L $1.flt.vcf -L $2.flt.vcf

##Genotype the called haplotypes for further processing

nohup gatk GenotypeGVCFs -R /data4/pilib/genomes/hs/hg19/hg19.fa -V "$1"_gatk.vcf -O "$1"_gatk.g.vcf  &> geno_"$1"_gatk.out
nohup gatk GenotypeGVCFs -R /data4/pilib/genomes/hs/hg19/hg19.fa -V "$2"_gatk.vcf -O "$2"_gatk.g.vcf  &> geno_"$2"_gatk.out

##Extract data from vcf files- Chromosome,Location,genotype
#bcftools query HN51.gatk.gz -f '%CHROM\t%POS[\t%GT]\n' > HN51_matrix.txt
##Extract sample names
#bcftools query -l HN51.gatk.gz > HN51_names.txt

##COMPARISON-(Move to a different script)
##Compress outputs to save space, -c param prints the compressed output to STDOUT
##Sample 1 from GATK & samtools
bgzip -c "$1".flt.vcf > "$1".sam.gz
bgzip -c "$1"_gatk.g.vcf > "$1".gatk.gz
##Sample 2 from GATK ad samtools
bgzip -c "$2".flt.vcf > "$2".sam.gz
bgzip -c "$2"_gatk.g.vcf > "$2".gatk.gz

##Index the SNP location regions in the format "chr:beginPos-endPos" 
##using tabix. -p param is the file type 
#Sample 1
tabix -p vcf "$1".sam.gz
tabix -p vcf "$1".gatk.gz
#Sample 2
tabix -p vcf "$2".sam.gz
tabix -p vcf "$2".gatk.gz

##Compare the sample files from different pipelines {GATK,samtools}
##Comparison tells which pipeline is more effective in finding SNPs
#Sample 1
bcftools stats "$1".sam.gz "$1".gatk.gz  > "$1".stats.txt
#Sample 2
bcftools stats "$2".sam.gz "$2".gatk.gz  > "$2".stats.txt

##Plot the stats using plot-vcfstats
#module load python
mkdir plots
cd plots
mkdir "$1"
mkdir "$2"

plot-vcfstats -p "$1"/ -r "$1".stats.txt
plot-vcfstats -p "$2"/ -r "$2".stats.txt

echo "Done!"

