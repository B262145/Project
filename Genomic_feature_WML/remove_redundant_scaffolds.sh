#!/bin/bash


#SBATCH --job-name=remove_redundant_scaffolds
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=remove_redundant_scaffolds.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

awk '
    BEGIN {keep=0}
    /^>/ {
        # If current header is one of the excluded scaffolds, set keep=0
        if ($0 ~ />OX465514.1/ || $0 ~ />CATLOI010000001.1/ || $0 ~ />CATLOI010000004.1/) {
            keep=0
        } else {
            keep=1
        }
    }
    keep==1 {print}
' /mnt/loki/ross/assemblies/scales/pseudococcidae/Planococcus_citri/GCA_950023065.1_ihPlaCitr1.1_genomic.fna > genome_filtered.fna

# grep -v ">" removes all header lines (those starting with '>'). 

# wc outputs line count, word count, and character count. 

# awk '{print $3-$1}' subtracts the number of lines from the character count to remove the newline characters, 

# resulting in the total number of nucleotide bases in the file. 

grep -v ">" genome_filtered.fna | wc | awk '{print $3-$1}' 