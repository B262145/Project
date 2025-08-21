#!/bin/bash

#SBATCH --job-name=get_scaffold_lengths
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=get_scaffold_lengths.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

# Define the fasta file path
FASTA="/mnt/loki/ross/assemblies/scales/pseudococcidae/Planococcus_citri/GCA_950023065.1_ihPlaCitr1.1_genomic.fna"

# Run AWK to calculate lengths of specific scaffolds
awk '
BEGIN {
  seq = ""; name = ""
}
/^>/ {
  if (name ~ /^OX465509\.1$|^OX465510\.1$|^OX465511\.1$|^OX465512\.1$|^OX465513\.1$|^OX465514\.1$|^CATLOI010000001\.1$|^CATLOI010000002\.1$|^CATLOI010000003\.1$|^CATLOI010000004\.1$/) {
    printf "%s\t%.3f Mb\n", name, length(seq)/1e6
  }
  name = substr($0, 2)
  split(name, a, " "); name = a[1]
  seq = ""
  next
}
{
  seq = seq $0
}
END {
  if (name ~ /^OX465509\.1$|^OX465510\.1$|^OX465511\.1$|^OX465512\.1$|^OX465513\.1$|^OX465514\.1$|^CATLOI010000001\.1$|^CATLOI010000002\.1$|^CATLOI010000003\.1$|^CATLOI010000004\.1$/) {
    printf "%s\t%.3f Mb\n", name, length(seq)/1e6
  }
}
' "$FASTA"
