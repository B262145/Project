#!/bin/bash

#SBATCH --job-name=01_cov_column_extractor          # Job name
#SBATCH --ntasks=1                               # Number of tasks (Rscript runs as a single process)
#SBATCH --mem=8G                                 # Total memory required
#SBATCH --time=2:00:00                           # Maximum run time (HH:MM:SS)
#SBATCH --output=01_cov_column_extractor.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                     # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk  

# This is a modified version of a script originally written by Hollie Marshall,
# available at: https://github.com/MooHoll/physalia-DNAm-EcoEvo/tree/main/5_Friday
input_dir="/home/s2343706/tamsin/diff_exp_meth/WGBS_newassembly/output/"

for i in $(seq -w 1 24); do
    cut -f1,2,4,5 "${input_dir}/DNA${i}/DNA${i}.CpG_report.txt" > "DNA${i}_coverage.txt"
done

# Remove mitochondrial chromosome (OX465514.1) and CATLOI010000001.1 to CATLOI010000004.1
#for i in $(seq -w 1 24); do
    #cut -f1,2,4,5 "${input_dir}/DNA${i}/DNA${i}.CpG_report.txt" | \
    #awk '$1 != "OX465514.1" && 
         #$1 != "CATLOI010000001.1" && 
         #$1 != "CATLOI010000002.1" && 
         #$1 != "CATLOI010000003.1" && 
         #$1 != "CATLOI010000004.1"' \
    #> "DNA${i}_coverage.txt"
#done
