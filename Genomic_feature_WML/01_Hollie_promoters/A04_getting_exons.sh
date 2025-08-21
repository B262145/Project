#!/bin/bash

#SBATCH --job-name=exon
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=A04_getting_exons.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

( echo -e "chr\tstart\tend\tstrand\tgene_id" && \
grep "exon" /mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.augustus.hints.gff3 | \
cut -f1,4,5,7,9 | \
sed 's/ID=//g' ) > exons.txt
