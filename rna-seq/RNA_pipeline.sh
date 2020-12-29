#!/bin/bash


#Téléchargement des fastq
mkdir -p fastq
wget -nc -O fastq/TPrnaseq.tar.gz http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz
tar -zxvf fastq/TPrnaseq.tar.gz -C fastq
rm -f fastq/TPrnaseq.tar.gz


#Téléchargement du chromosome 18 et de ses annotations
mkdir -p human_genome
if [ ! -f human_genome/chr18.fa ]
then
	wget -nc -O human_genome/chr18.fa.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz
	gunzip human_genome/chr18.fa.gz
fi

if [ ! -f human_genome/gencode.v24lift37.basic.annotation.gtf ]
then
	wget -nc -O human_genome/gencode.v24lift37.basic.annotation.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz
	gunzip human_genome/gencode.v24lift37.basic.annotation.gtf.gz
fi

#Indexation du génome
#--genomeSAindexNbases 12 > option demandé par STAR quand je l'ai lancé
STAR --runMode genomeGenerate --runThreadN 4 \
	--genomeDir human_genome \
	--genomeFastaFiles human_genome/chr18.fa \
	--sjdbGTFfile human_genome/gencode.v24lift37.basic.annotation.gtf \
	--genomeSAindexNbases 12  



#Trimming des reads
mkdir -p trimmed_fastq
for file_name in fastq/*".R1.fastq"
do
	base=$(basename $file_name ".R1.fastq")
	trimmomatic PE fastq/$base.R1.fastq fastq/$base.R2.fastq -baseout trimmed_fastq/$base.fastq LEADING:10 TRAILING:10 MINLEN:50
	#Note : j'ai changé les valeurs de LEADING et TRAILING parce que l'alignement ne fonctionnait pas pour les valeurs proposées dans le TP
	#STAR affichait l'alignement de 90% des reads comme "too short", c'est à dire qu'il n'arrivait à les aligner suffisamment nulle part en paired-end
done

#Alignement des reads
mkdir -p bam
for file_name in trimmed_fastq/*_1P.fastq
do
	base=$(basename $file_name "_1P.fastq")

	STAR --runThreadN 4 --outFilterMultimapNmax 1\
		--genomeDir human_genome \
		--outSAMattributes All --outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix bam/$base.bam \
		--readFilesIn trimmed_fastq/$base"_1P.fastq" trimmed_fastq/$base"_2P.fastq"
done

#Comptage des reads avec feature count
mkdir -p count
featureCounts -p -t exon -g gene_id -a human_genome/gencode.v24lift37.basic.annotation.gtf -o count/resultat_feature_counts.txt bam/*.bam
perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' human_genome/gencode.v24lift37.basic.annotation.gtf | sort | uniq > equivalence_table.txt


#Reformatage des résultats et ajout des noms de gène
mkdir -p results
#grep -v ";" permet de se débarasser des reads multimappées -> probablement à cause des modifications des paramètres que j'ai effectué mais c'est le seul moyen que j'ai trouvé pour obtenir des résultats
grep -v ";" count/resultat_feature_counts.txt | grep "^[^#;]" | grep "chr18" | sort > results/temp1
sort equivalence_table.txt > results/temp2
join results/temp1 results/temp2 | grep "chr18" > results/temp3
awk '{print $13,$S6,$7,$8,$9,$10,$11,$12}' results/temp3 > results/table_comptage_read_for_R.table
rm -f temp1, temp2, temp3, temp4
