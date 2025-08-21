#!/bin/bash
#SBATCH --job-name=test    # Job name
#SBATCH --ntasks=1                                 # Number of tasks 
#SBATCH --cpus-per-task=2                         # CPUs required
#SBATCH --mem=100G                                  # Total memory required
#SBATCH --output=test.%j.out  # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

#Run the script
Rscript test.R
