# Scripts pour l'atelier NGS



## RNA-Seq :
Commande : sh RNA_pipeline.sh

- Le script télécharge les différents fichiers et index le chromosome (ne prend pas beaucoup de temps)
- Résultats dans results/table_comptage_read_for_R.table


## Analyse Variants :

Commande : sh variant.sh
- ATTENTION : nécessite un dossier human_genome/ contenant le chromosome 18 indexé et un dossier fastq contenant les fastq
- Le script télécharge le fichier d'annotation
- Résultats dans output/gene_variant_list.txt
