# This is a modified version of a script originally written by Hollie Marshall,
# available at: https://github.com/MooHoll/physalia-DNAm-EcoEvo/tree/main/5_Friday

library(readr)

#List all *_coverage.txt files in current directory
file.list <- list.files("./", pattern="^DNA.*_coverage\\.txt$")

#Read function
read_file1 <- function(x) {
  read_delim(x, "\t", col_names = FALSE)
}

#Read all files into a list
samples <- lapply(file.list, read_file1)

#Define sample IDs
sample.ids <- c(
  paste0("adultfemale_", sprintf("%02d", 1:6)),
  paste0("adultmale_", sprintf("%02d", 7:12)),
  paste0("3rd_female_", sprintf("%02d", 13:18)),
  paste0("3rd_male_", sprintf("%02d", 19:24))
)

#Sort the file list if needed to match DNA01 to DNA24
#This step ensures consistent ordering if list.files doesn't return them sorted
file.list <- sort(file.list)

#Assign names
names(samples) <- sample.ids

#Process and export each sample
for(i in seq_along(samples)) {
  colnames(samples[[i]]) <- c("chr", "cpg", "count_c", "count_t")
  samples[[i]]$total_coverage <- samples[[i]]$count_c + samples[[i]]$count_t
  samples[[i]] <- samples[[i]][, -4]  # remove 'count_t'
  final_file <- samples[[i]]
  myfile <- file.path("./", paste0(names(samples[i]), "_final_coverage.txt"))
  write.table(final_file, file = myfile, quote = FALSE, sep = "\t", row.names = FALSE)
}
