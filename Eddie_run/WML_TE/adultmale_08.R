library(sqldf)
library(readr)
library(doBy)
library(dplyr)

meth_dir <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Weighted_meth/no_MT/"
input_file <- file.path(meth_dir, "adultmale_08_final_coverage.txt")
sample_name <- "adultmale_08"

df <- read_delim(input_file, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)

annotation <- read.delim("TEAnnot.gff", header=FALSE, sep="\t")
colnames(annotation) <- c("chr","feature", "start", "end")

df <- subset(df, total_coverage > 10)

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

output <- output[!is.na(output$feature), ]
output <- output[, -c(1,2)]
check <- summaryBy(total_coverage + count_c ~ chr + feature + start + end, data=output, FUN=sum)
check$weightedMeth <- (check$count_c.sum)/(check$total_coverage.sum)

output_dir <- "./"
output_path <- file.path(output_dir, paste0(sample_name, "_weighted_meth_per_genomic_region.txt"))
write.table(check, file=output_path, quote=FALSE, sep="\t", row.names=FALSE)
