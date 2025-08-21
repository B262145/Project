import csv
import os

# Base directory containing DNA01–DNA24 folders
base_dir = "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/diff_exp_meth/WGBS_newassembly/output"

# Output file path
output_file = "methylation_summary.tsv"

# Define methylation level function
def calc_level(m, u):
    return m / (m + u) * 100 if (m + u) > 0 else 0

# Open output file for writing
with open(output_file, "w") as out_f:
    out_f.write("Sample\tCpG_methylation(%)\tCHG_methylation(%)\tCHH_methylation(%)\n")

    # Loop over DNA01 to DNA24
    for i in range(1, 25):
        sample_id = f"DNA{i:02d}"
        input_file = os.path.join(base_dir, sample_id, f"{sample_id}.cytosine_context_summary.txt")

        # Initialize counters
        cpg_m, cpg_u = 0, 0
        chg_m, chg_u = 0, 0
        chh_m, chh_u = 0, 0

        try:
            with open(input_file) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    context = row['C-context'].upper()
                    try:
                        m = int(row['count methylated'])
                        u = int(row['count unmethylated'])
                    except ValueError:
                        continue

                    if context.startswith('CG'):
                        cpg_m += m
                        cpg_u += u
                    elif context[0] == 'C' and context[2] == 'G' and context[1] != 'G':
                        chg_m += m
                        chg_u += u
                    elif context[0] == 'C' and context[1] != 'G' and context[2] != 'G':
                        chh_m += m
                        chh_u += u

            # Calculate methylation percentages
            cpg_pct = calc_level(cpg_m, cpg_u)
            chg_pct = calc_level(chg_m, chg_u)
            chh_pct = calc_level(chh_m, chh_u)

            # Print to screen
            print(f"{sample_id}")
            print(f"  CpG methylation: {cpg_pct:.2f}%")
            print(f"  CHG methylation: {chg_pct:.2f}%")
            print(f"  CHH methylation: {chh_pct:.2f}%\n")

            # Write to file
            out_f.write(f"{sample_id}\t{cpg_pct:.2f}\t{chg_pct:.2f}\t{chh_pct:.2f}\n")

        except FileNotFoundError:
            print(f"{sample_id} - File not found: {input_file}")
            out_f.write(f"{sample_id}\tNA\tNA\tNA\n")
