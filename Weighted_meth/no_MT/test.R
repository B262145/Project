# Load required libraries
library(readr)
library(dplyr)
library(fuzzyjoin)
library(doBy)

# Set test file name
sample_file <- "adultfemale_04_final_coverage.txt"
sample_name <- gsub("_final_coverage\\.txt$", "", sample_file)

# Read methylation data
message("Reading sample file...")
df <- read_delim(sample_file, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
df <- df %>% filter(total_coverage > 10)

# Read TE annotation
message("Reading TE annotations...")
annotation <- read_table("../TEAnnot.gff", col_names = FALSE)
colnames(annotation) <- c("chr", "feature", "start", "end")

# Join methylation data with annotation by position
message("Joining methylation data with TE annotations...")
output <- fuzzy_left_join(
  df,
  annotation,
  by = c("chr" = "chr",
         "cpg" = "start",   # match cpg >= start
         "cpg" = "end"),    # match cpg <= end
  match_fun = list(`==`, `>=`, `<=`)
) %>%
  filter(!is.na(feature)) %>%
  select(chr = chr.x, feature, start, end, total_coverage, count_c)

# Summarize: calculate weighted methylation
message("Calculating weighted methylation...")
check <- summaryBy(total_coverage + count_c ~ chr + feature + start + end, data = output, FUN = sum)
check$weightedMeth <- check$count_c.sum / check$total_coverage.sum

# Save result
output_dir <- "test"
myfile <- file.path(output_dir, paste0(sample_name, "_weighted_meth.txt"))
write.table(check, file = myfile, quote = FALSE, sep = "\t", row.names = FALSE)

message("Written: ", myfile)
