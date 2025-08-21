# -------------------------------------------------------------------------
# STEP 4: Merge weighted methylation results across all samples
# -------------------------------------------------------------------------

# Load required libraries
library(readr)
library(dplyr)

# Define helper function to read each file
read_file1 <- function(x) {
  read_delim(x, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
}

# Set your output directory containing the per-sample WML results
output_dir <- "Hollie_split"

# Message
message("STEP 4: Merging all sample results into one dataframe...")

# List all individual sample result files
file.list <- list.files(path = output_dir, pattern = "*weighted_meth_per_genomic_region.txt", full.names = TRUE)

# Extract sample names from file names
sample_names <- gsub("_weighted_meth_per_genomic_region.txt", "", basename(file.list))

# Ensure correct pairing
ordering <- order(sample_names)
file.list <- file.list[ordering]
sample_names <- sample_names[ordering]

# Read in each file
samples_wm <- lapply(file.list, read_file1)
names(samples_wm) <- sample_names

# Add sample name column to each dataframe
samples_with_name <- lapply(names(samples_wm), function(current_name) {
  transform(samples_wm[[current_name]], sample = current_name)
})

# Merge all dataframes into one
all <- bind_rows(samples_with_name)

# Output file path
output_file <- file.path(output_dir, "weighted_meth_per_genomic_region_all_samples.txt")

# Write to file
write.table(all, file = output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Completion message
message("Combined results written to: ", output_file)
