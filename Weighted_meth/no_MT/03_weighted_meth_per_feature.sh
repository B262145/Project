#!/bin/bash
#SBATCH --job-name=New_WML_per_TE       # Job name
#SBATCH --ntasks=1                                 # Number of tasks 
#SBATCH --cpus-per-task=24                         # CPUs required
#SBATCH --mem=100G                                  # Total memory required
#SBATCH --output=03_weighted_meth_per_feature.%j.out  # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

#Run the script
Rscript 03_weighted_meth_per_feature.R
