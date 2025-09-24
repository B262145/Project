# ----------------------------------------------------------------
### Adding putative promotor information to the annotation file
# ----------------------------------------------------------------
library(readr)

genes <- read_delim("genes_with_start_and_end.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)
head(genes)

genes$promotor_start <- ifelse(genes$strand =="+", genes$start - 2000, genes$end + 2000)
nrow(genes) 
genes <- genes[genes$promotor_start > 0,]
nrow(genes) #promotor overlaps scaffold start: removed

plus_strand_genes <- genes[genes$strand == "+",]
plus_strand_genes <- plus_strand_genes[,-3] # rm redundent end column
colnames(plus_strand_genes)[2] <- "end"
colnames(plus_strand_genes)[5] <- "start"
plus_strand_genes$feature <- "promotor"

minus_strand_genes <- genes[genes$strand == "-",]
minus_strand_genes <- minus_strand_genes[,-2] # rm redundent start column
colnames(minus_strand_genes)[2] <- "start"
colnames(minus_strand_genes)[5] <- "end"
minus_strand_genes$feature <- "promotor"

gene_promotors <- rbind(plus_strand_genes, minus_strand_genes)

write.table(gene_promotors, file="gene_promotors.csv", 
            row.names = F, col.names = T, sep = '\t', quote = F)
