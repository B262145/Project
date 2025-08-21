#!/bin/bash

#SBATCH --job-name=3rd_female_14
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G

#SBATCH --output=3rd_female_14.%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk

Rscript 3rd_female_14.R
