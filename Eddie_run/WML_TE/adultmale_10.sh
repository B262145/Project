#!/bin/bash

#SBATCH --job-name=adultmale_10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G

#SBATCH --output=adultmale_10.%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk

Rscript adultmale_10.R