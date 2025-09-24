#!/bin/bash
#SBATCH --job-name=pooled_group_level_genomewide_WML    # Job name
#SBATCH --ntasks=1                                 # Number of tasks 
#SBATCH --cpus-per-task=12                         # CPUs required
#SBATCH --mem=16G                                  # Total memory required
#SBATCH --time=24:00:00                            # Maximum run time (HH:MM:SS)
#SBATCH --output=04_pooled_group_level_genomewide_WML%j.out  # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

#Run the script
Rscript 04_pooled_group_level_genomewide_WML.R
