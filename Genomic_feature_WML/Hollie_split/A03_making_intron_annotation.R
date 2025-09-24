
# Note: Corrected coordinate system by converting intron.bed from 0-based to 1-based to match gene annotation.

library(readr)
library(dplyr)
library(sqldf)

genes <-read.csv.sql("genes_with_start_and_end.txt",
                                      sql = "select * from file", sep = "\t", header = TRUE) 

intron_annot<-read.csv.sql("intron.bed",
                                     sql = "select * from file", sep = "\t", header = FALSE)
colnames(intron_annot) <- c("chr", "start", "end")

# 0-based to 1-based coordinates
intron_annot$start <- intron_annot$start + 1

# Add intron information
output <- sqldf("SELECT intron.chr,
                intron.start,
                intron.end,
                genes.chr,
                genes.start,
                genes.end,
                genes.strand,
                genes.gene_id
                FROM intron_annot AS intron
                LEFT JOIN genes AS genes
                ON intron.chr = genes.chr
                AND (intron.start >= genes.start AND intron.start <= genes.end)")

output<-subset(output, !output$gene_id=="NA")
output<-output[!duplicated(output),]
output <- output[,-c(4,5,6)]

output$feature<-"intron"
write.table(output, file="introns.csv", 
            row.names = F, col.names = T, sep = '\t', quote = F)
