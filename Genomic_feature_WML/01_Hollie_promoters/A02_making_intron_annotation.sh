#!/bin/bash

#SBATCH --job-name=intron
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=A02_making_intron_annotation.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

#GFF file doesn't contain intron regions so need to define them and add to a bed file with all annotation so can plot methylation one features.

#Found link which explains how to define introns from a gff:https://davetang.org/muse/2013/01/18/defining-genomic-regions/

#In brief:
#Use subtractBed command (part of bedtools) to subtract the exon regions from the overall gene region, leaving introns behind.

#First define all exon without overlap then take from gene:
cat /mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.augustus.hints.gff3 |
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' |
sortBed |
mergeBed -i - | gzip > exon_merged.bed.gz

gunzip exon_merged.bed.gz

cat /mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.augustus.hints.gff3 |
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
sortBed |
subtractBed -a stdin -b exon_merged.bed > intron.bed

#Then can put this file into R along with the main annotation file and use SQL to add in intron info with gene name.