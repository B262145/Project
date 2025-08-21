#!/bin/bash

#SBATCH --job-name=B02_getting_intergenic_regions
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=B02_getting_intergenic_regions.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address


python B02_getting_intergenic_regions.py annotations_for_getting_intergenic_regions.gtf > intergenic.txt
