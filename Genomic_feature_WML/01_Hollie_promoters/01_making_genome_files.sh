#!/bin/bash

#SBATCH --job-name=gene
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=01_making_genome_files.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

grep "gene" /mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.augustus.hints.gff3 > genes.txt
cut -f1,4,5,7,9 genes.txt > genes_next.txt
sed 's/ID=//g' genes_next.txt > new.txt
echo -e "chr\tstart\tend\tstrand\tgene_id" | cat - new.txt > genes_with_start_and_end.txt
