#This is a modified version of a script originally written by Hollie Marshall,
#available at: https://github.com/MooHoll/physalia-DNAm-EcoEvo/blob/main/5_Friday/weighted_meth_per_feature.R

## -------------------------------------------------------------------------
# Weighted methylation per genomic feature
## -------------------------------------------------------------------------

# This is quite an intensive script and may need running on a cluster for your real samples
library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(foreach)
library(doParallel)

## -------------------------------------------------------------------------

# STEP 1: Read in sample methylation count files
message("STEP 1: Reading in sample methylation count files...")

# Define the directory path of the methylation data
meth_dir <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/TE_WML/no_MT/"

file.list <- list.files(meth_dir, pattern = "*_final_coverage.txt", full.names = TRUE)

sample_names <- gsub("_final_coverage.txt", "", basename(file.list))
ordering <- order(sample_names)
file.list <- file.list[ordering]
sample_names <- sample_names[ordering]

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
names(samples) <- sample_names

# Check if sample_names matches file.list
cat("File check:\n")
for (i in seq_along(sample_names)) {
  cat(sample_names[i], " --> ", basename(file.list[i]), "\n")
}

# STEP 2: Read in annotations
message("STEP 2: Reading in annotations...")

# Read in TE annotations with start/end
annotation <- read.delim("merged_annotations.txt", header=TRUE, sep="\t")
##Guess I don't need width and strand columns
annotation <- annotation[, c("chr", "start", "end", "feature")]

# STEP 3: Calculate weighted methylation per feature per sample
message("STEP 3: Calculating weighted methylation per feature for each sample...")

# Tell R how many cores it has available to use
registerDoParallel(cores = 24)

# Define the output directory
output_dir <- "24_cores"

# Calculate weighted meth for each feature for each sample (Windows)
foreach(i = seq_along(samples), .packages=c("sqldf","doBy","dplyr"),
        .export = ls(globalenv())) %dopar% {
  start_time <- Sys.time()         
  df <- samples[[i]]
  df <- subset(df, total_coverage > 10) ###Changed the minimal coverage to 10
  output <- sqldf("SELECT sample.chr AS chr,
                    sample.cpg,
                    sample.count_c,
                    sample.total_coverage,
                    annot.chr,
                    annot.feature,
                    annot.start,
                    annot.end
                    FROM df AS sample
                    LEFT JOIN annotation AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
  output <- output[!is.na(output$feature),] ##Use the 'feature' column to determine whether the CpG matched any annotation
  output <- output[,-c(1,2)]
  check <- summaryBy(total_coverage + count_c ~ chr + feature + start + end, data=output, FUN=sum) 
  check$weightedMeth <- (check$count_c.sum)/(check$total_coverage.sum)
  myfile <- file.path(output_dir, paste0(names(samples[i]),"_","weighted_meth_per_genomic_region.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
  cat("Written: ", myfile, "\n")
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  cat("Sample", names(samples[i]), "done in", elapsed, "seconds\n")
}

# STEP 4: Merge results into one dataframe
message("STEP 4: Merging all sample results into one dataframe...")
# List all individual sample result files
file.list = list.files(path = output_dir, pattern="*weighted_meth_per_genomic_region.txt", full.names = TRUE)

# Extract sample names from file names
sample_names <- gsub("_weighted_meth_per_genomic_region.txt", "", basename(file.list))
# Sort files and sample names in the same order to ensure correct pairing
ordering <- order(sample_names)
file.list <- file.list[ordering]
sample_names <- sample_names[ordering]
# Read in each file and assign corresponding sample name
samples_wm <- lapply(file.list, read_file1)
names(samples_wm) <- sample_names

# Add sample name as a column in each dataframe
samples_with_name <- lapply(names(samples_wm),
  function(current_name) transform(samples_wm[[current_name]], sample = current_name)
)

# Merge dataframes
all <- as.data.frame(bind_rows(samples_with_name))
# Write merged dataframe to file
output_file <- file.path(output_dir, "weighted_meth_per_genomic_region_all_samples.txt")
write.table(all, file=output_file, col.names = T,
            row.names = F, quote = F, sep = "\t")
message("Combined results written to: ", output_file)
