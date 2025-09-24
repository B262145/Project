#!/bin/bash
#SBATCH --job-name=format_methylation_coverage_by_sample       # Job name
#SBATCH --ntasks=1                                             # Number of tasks (Rscript runs as a single process)
#SBATCH --mem=8G                                               # Total memory required
#SBATCH --time=2:00:00                                         # Maximum run time (HH:MM:SS)
#SBATCH --output=02_format_methylation_coverage_by_sample.%j.out  # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                                   # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk                          # Email address

#Run the script
Rscript 02_format_methylation_coverage_by_sample.R
