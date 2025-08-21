#!/bin/bash


#SBATCH --job-name=remove_scaffolds
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=remove_unplaced_scaffolds.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address


# This script removes the unplaced scaffolds from the filtered TE annotation
awk '$1 != "CATLOI010000001.1" && $1 != "CATLOI010000003.1"' /mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Genomic_feature_WML/filtered_TEannot.gff > filtered_no_unplaced_TE.gff