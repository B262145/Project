#!/bin/bash

# This script strips away parts of a .gff to leave only the TE annotations we need 

#SBATCH --job-name=get_clean_gffAnnot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=get_clean_gffAnnot.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

# Step 1: Remove comment lines starting with '#'
# Step 2: Remove lines where the first column equals 'OX465514.1'
# Step 3: Extract columns 1 (chromosome), 3 (feature type), 4 (start), and 5 (end)

grep -v "^#" /mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.filteredRepeats.renamed.Wes.trimmed.gff \
| awk '$1 != "OX465514.1"' \
| cut -f 1,3,4,5 > TEAnnot.gff
