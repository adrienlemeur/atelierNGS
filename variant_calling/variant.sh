#!/bin/bash

#Trimming des reads
mkdir -p trimmed_fastq
for file_name in fastq/*"_r1F.fastq.gz"
do
	base=$(basename $file_name "_r1F.fastq.gz")

	trimmomatic PE "fastq/$base""_r1F.fastq.gz" "fastq/$base""_r2F.fastq.gz" -baseout "trimmed_fastq/$base.fastq" \
	LEADING:20 TRAILING:20 MINLEN:50
done

mkdir -p sam
#Mapping des reads
for file_name in trimmed_fastq/*"_1P.fastq"
do
	base=$(basename $file_name "_1P.fastq")
	bwa mem -M -t 2 -A 2 -E 1 human_genome/chr16.fa.gz trimmed_fastq/$base"_1P.fastq" trimmed_fastq/$base"_2P.fastq" > \
	sam/$base.sam
done
mkdir -p pileup
mkdir -p bam

#Processing SAM files
for file_name in sam/*.sam
do
	base=$(basename $file_name ".sam")
	samtools view -S -b $file_name  > bam/$base.bam
	samtools sort bam/$base.bam -o bam/$base.sorted.bam
	samtools index bam/$base.sorted.bam
	samtools mpileup -B -A -f human_genome/chr16.fa bam/$base.sorted.bam > pileup/$base.msf
done

#Variant Calling
mkdir -p vcf
varscan somatic pileup/TCRBOA7-N-WEX-chr16.msf pileup/TCRBOA7-T-WEX-chr16.msf vcf/TCRBOA7-WEX-chr16 --variants --p-value 0.001 --min-avg-qual 15 --output-vcf 1

#Téléchargement et extraction du fichier d'annotation
mkdir -p annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz -O annotation/gencode.v24lift37.basic.annotation.gtf.gz
gunzip -f annotation/gencode.v24lift37.basic.annotation.gtf.gz

#VCF annotation
mkdir -p output

grep vcf/* -e ';SOMATIC;' | cut -d ':' -f2- > output/somatic_SNP_indel.table
awk '{OFS="\t"; if (!/^#/){print "chr16",$2-1,$2,$4"/"$5,"+"}}' output/somatic_SNP_indel.table > output/somatic_SNP_indel.bed

#On ne garde garde que les annotations du chromosome 16
rm -f annotation/chr16.gtf
grep "chr16" annotation/*.gtf > annotation/chr16.gtf

#Bedtools ne marchait pas donc j'ai fait la jointure à la main
cut output/somatic_SNP_indel.bed -f 2,3 | grep - annotation/chr16.gtf | grep "\sgene\s" | awk '{print " " $1 " " $4 " " $5 " " $16}' > output/gene_variant_list.txt
