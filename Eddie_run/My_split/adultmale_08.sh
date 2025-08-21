#!/bin/bash

#SBATCH --job-name=adultmale_08
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G

#SBATCH --output=adultmale_08.%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk

Rscript adultmale_08.R
