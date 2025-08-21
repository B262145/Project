#--------------------------------------------------------------------
# Making intergenic regions or unannotated regions to quantify background methylation
#--------------------------------------------------------------------

library(readr)
library(dplyr)

# 1. Promoter (2000bp upstream)
promoters <- read.delim("gene_promotors.csv", header=T, sep="\t")
promoters <- promoters[,c("chr", "start", "end", "gene_id")]
promoters$feature <- "promoter_2000bp"

# 2. Whole gene
genes <- read.delim("genes_with_start_and_end.txt", header=T)
genes <- genes[,c("chr", "start", "end", "gene_id")]
genes$feature <- "whole_gene"

# 3. TE
TEs <- read.delim("/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Genomic_feature_WML/filtered_TEannot.gff", header=F)
colnames(TEs) <- c("chr","annotation_source","something","start","end","something1","strand","something2","gene_id") # Assign TE ID to "gene_id" to standardize column names across features for merging.

TEs <- TEs[,c("chr", "start", "end", "gene_id")]
TEs$feature <- "TE"

merged_annotations <- rbind(promoters, genes, TEs)

# Order the file and make in a format the python script wants
ordered <- merged_annotations %>% arrange(chr, start)
head(ordered)
colnames(ordered) <- c("chr","start","end","gene")
ordered <- ordered[,c(4,1,2,3)]

write.table(ordered, file="annotations_for_getting_intergenic_regions.gtf",
            col.names = T, row.names = F, quote = F, sep = '\t')


# Run getting_intergenic_regions.py as:
# python B02_getting_intergenic_regions.py annotations_for_getting_intergenic_regions.gtf > intergenic.txt

intergenic <- read_delim("intergenic.txt", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

# Remove rows which are incorrect becasue they were created when a TE is found
# within gene conordinates
head(intergenic)

intergenic$remove <- ifelse(intergenic$X2 > intergenic$X3, "yes","no")
intergenic_for_sure <- intergenic[intergenic$remove == "no",]

# Neaten the file up for later use
intergenic_for_sure$X4 <- "intergenic"
intergenic_for_sure <- intergenic_for_sure[,-5]
colnames(intergenic_for_sure) <- c("chr","start","end","feature")

write.table(intergenic_for_sure, file="Curated_intergenic.txt",
            col.names = T, row.names = F, quote = F, sep = '\t')