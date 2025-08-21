#!/bin/bash

#SBATCH --job-name=CpG
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=process_CpG_reports.%j.out     # Standard output and error log (%j = job ID)
#SBATCH --mail-type=END,FAIL                       # Send email on job END or FAIL
#SBATCH --mail-user=s2673561@ed.ac.uk              # Email address

input_base="/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/diff_exp_meth/WGBS_newassembly/output"
output_base="/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Cov_files"

mkdir -p "$output_base"

for i in $(seq -w 1 24); do
    sample="DNA$i"
    input_file="$input_base/$sample/$sample.CpG_report.txt"
    output_file="$output_base/$sample.cov.txt"

    if [[ ! -f "$input_file" ]]; then
        echo "$input_file not found, skipping..."
        continue
    fi

    awk '
    BEGIN {
        OFS = "\t"
    }
    $1 != "OX465514.1" && 
    $1 != "CATLOI010000001.1" &&
    $1 != "CATLOI010000002.1" &&
    $1 != "CATLOI010000003.1" &&
    $1 != "CATLOI010000004.1" &&
    ($4 + $5) > 0 {
        meth_perc = ($4 / ($4 + $5)) * 100
        print $1, $2, $2, meth_perc, $4, $5
    }' "$input_file" > "$output_file"

    echo "Processed $sample → $output_file"
done
