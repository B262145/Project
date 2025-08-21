#!/bin/bash
#SBATCH --job-name=New_Split_By_Feature    # Job name
#SBATCH --ntasks=1                                 # Number of tasks 
#SBATCH --cpus-per-task=1                          # CPUs required
#SBATCH --mem=4G                                  # Total memory required
#SBATCH --output=Split.%j.out  # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

#Run the script
Rscript Split.R
