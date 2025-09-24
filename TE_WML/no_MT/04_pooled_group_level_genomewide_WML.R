# This script reads per‐sample CpG coverage files (filtered to coverage ≥ 10),
# pools all retained CpGs across samples within each groups (Adult Female, Adult Male, 3rd Instar Female, 3rd Instar Male), and then
# directly computes each group’s genome‐wide weighted methylation as (sum of all count_C across group) / (sum of all total_coverage across group).


library(readr)
library(dplyr)

# STEP 1: Load all sample files
file.list <- list.files("./", pattern = "*final_coverage.txt", full.names = TRUE)
# Extract sample names from filenames
sample_names <- gsub("_final_coverage.txt", "", basename(file.list))

# Sort files and names in consistent order
ordering <- order(sample_names)
file.list <- file.list[ordering]
sample_names <- sample_names[ordering]

# Load each file into memory
read_file1 <- function(x){
  read_delim(x, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
}

samples <- lapply(file.list, read_file1)
names(samples) <- sample_names

# Check file-sample mapping
cat("File check:\n")
for (i in seq_along(sample_names)) {
  cat(sample_names[i], " --> ", basename(file.list[i]), "\n")
}

# Assign each sample to a group
group_labels <- case_when(
  grepl("^adultfemale", sample_names)  ~ "Adult Female",
  grepl("^adultmale",   sample_names)  ~ "Adult Male",
  grepl("^3rd_female",  sample_names)  ~ "3rd Instar Female",
  grepl("^3rd_male",    sample_names)  ~ "3rd Instar Male",
  TRUE                                 ~ "Unknown"
)

# Create a dataframe mapping samples to groups and file paths
sample_group_map <- data.frame(
  sample = sample_names,
  group = group_labels,
  file = file.list,
  stringsAsFactors = FALSE
)

# Show the result
cat("\nSample grouping:\n")
print(sample_group_map)

# STEP 2: Compute pooled genome-wide methylation per group
group_meth_list <- sample_group_map %>%
  group_by(group) %>%
  group_split()

group_results <- lapply(group_meth_list, function(group_df) {
  group_name <- unique(group_df$group)
  
  # Initialize counters
  total_count_c <- 0
  total_coverage <- 0
  
  for (i in seq_len(nrow(group_df))) {
    df <- read_delim(group_df$file[i], "\t", col_names = TRUE, trim_ws = TRUE)
    df_filt <- df %>% filter(total_coverage >= 10)
    
    total_count_c <- total_count_c + sum(df_filt$count_c, na.rm = TRUE)
    total_coverage <- total_coverage + sum(df_filt$total_coverage, na.rm = TRUE)
  }
  
  genomewide_meth <- total_count_c / total_coverage
  data.frame(group = group_name,
             total_count_c = total_count_c,
             total_coverage = total_coverage,
             genomewide_meth = genomewide_meth)
})

group_meth_df <- bind_rows(group_results)

# STEP 3: Save result
output_dir <- "24_cores"

write.table(group_meth_df,
            file = file.path(output_dir, "pooled_genomewide_weighted_meth_by_group.csv"),
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

message("Pooled group-level genome-wide methylation saved to:\n",
        file.path(output_dir, "pooled_genomewide_weighted_meth_by_group.csv"))
