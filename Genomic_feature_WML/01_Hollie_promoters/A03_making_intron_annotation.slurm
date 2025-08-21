#!/bin/bash
#SBATCH --job-name=A03_making_intron_annotation    # Job name
#SBATCH --ntasks=1                                 # Number of tasks 
#SBATCH --cpus-per-task=1                          # CPUs required
#SBATCH --mem=4G                                  # Total memory required
#SBATCH --time=2:00:00                            # Maximum run time (HH:MM:SS)
#SBATCH --output=A03_making_intron_annotation.%j.out  # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

#Run the script
Rscript A03_making_intron_annotation.R
