#!/bin/bash
#SBATCH --job-name=step4   # Job name
#SBATCH --ntasks=1                                 # Number of tasks 
#SBATCH --cpus-per-task=1                          # CPUs required
#SBATCH --mem=4G                                  # Total memory required
#SBATCH --output=step4_merge_weighted_methylation_results.%j.out  # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

#Run the script
Rscript step4_merge_weighted_methylation_results.R
