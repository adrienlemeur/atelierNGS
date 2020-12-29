# Scripts pour l'atelier NGS



## RNA-Seq :
sh RNA_pipeline.sh

- Le script télécharge les différents fichiers et indexe le chromosome (ne prend pas beaucoup de temps)
- Résultats dans results/table_comptage_read_for_R.table


## Analyse Variants :

sh variant.sh
- ATTENTION : nécessite un dossier human_genome/ contenant le chromosome 18 indexé et un dossier fastq contenant les fastq
- Le script télécharge l'annotation
