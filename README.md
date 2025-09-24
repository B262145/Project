# *Planococcus citri* DNA Methylation and TE Analysis

This project investigates DNA methylation patterns in *Planococcus citri* with a focus on transposable elements (TEs). It integrates WGBS data, genomic annotations, and statistical analysis to quantify methylation levels, identify differentially methylated TEs, and explore sex- and developmental stage-specific patterns.

---

## Project Structure

This project is organized into six repositories, each with a specific role in the workflow:

1. **Cov_files** – Generate and organize CpG coverage files from Bismark-outputted CpG reports.  
2. **Eddie_run** – Run weighted_meth_per_feature.R separately for each sample. This was done on Eddie and could avoid the memory and runtime issues caused while processing all 24 samples at once on AC3.
3. **Data_preprocessing** – Preprocess methylation and genomic data for downstream analysis. Scripts include calculation of CpG/CHG/CHH methylation levels per sample and extraction of scaffold lengths from the reference genome.
5. **Genomic_feature_WML** – Calculate weighted methylation levels for genomic features.
6. **TE_WML** – Calculate weighted methylation levels for TEs.    
7. **TE_differential_methylation** – Perform differential methylation analysis at TE copy level.  

---
