#!/bin/bash


#SBATCH --job-name=remove_repeats
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=remove_repeats.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address


# This script removes the Low_complexity, Satellite, Simple_repeat and Unknown features from the TE annotation 
# and only retains TEs.
# Additionally, Remove comment lines starting with '#' and Remove lines where the first column equals 'OX465514.1'

awk '$3 != "Low_complexity" && $3 != "Satellite" && $3 != "Simple_repeat" && $3 != "Unknown"' /mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.filteredRepeats.renamed.Wes.trimmed.gff \
 > filtered_TEannot.gff

# This script filters a TE annotation GFF file with the following steps:
# 1. Remove comment lines starting with '#'.
# 2. Remove lines where the first column equals 'OX465514.1' (mitochondria).
# 3. Exclude features annotated as 'Low_complexity', 'Satellite', 'Simple_repeat', or 'Unknown' in the third column.

grep -v "^#" /mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.filteredRepeats.renamed.Wes.trimmed.gff \
| awk '$1 != "OX465514.1" && $3 != "Low_complexity" && $3 != "Satellite" && $3 != "Simple_repeat" && $3 != "Unknown"' \
> filtered_TEannot.gff