library(readr)
library(dplyr)

# 1. Promoter (2kb upstream)
promoters <- read.delim("gene_promotors.csv", header=T, sep="\t")
promoters <- promoters[,c("chr", "start", "end", "gene_id")]
promoters$feature <- "promoter_2000bp"

# 2. Introns
introns <- read.delim("introns.csv", header=T)
introns <- introns[,c("chr", "start", "end", "gene_id")]
introns$feature <- "intron"

# 3. Exons
exons <- read.delim("exons.txt", header=T)
exons <- exons[,c("chr", "start", "end", "gene_id")]
exons$feature <- "exon"

# 4. Intergenic regions

intergenic <- read.delim("Curated_intergenic.txt", header=T)
intergenic <- intergenic[,c("chr", "start", "end")]
intergenic$gene_id <- "nope"
intergenic$feature <- "intergenic"

# 5. TE
TEs <- read.delim("/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Genomic_feature_WML/filtered_TEannot.gff", header=F)
colnames(TEs) <- c("chr","annotation_source","something","start","end","something1","strand","something2","gene_id") # Assign TE ID to "gene_id" to standardize column names across features for merging.

TEs <- TEs[,c("chr", "start", "end", "gene_id")]
TEs$feature <- "TE"


# 6. Merge all the regions
all_annotations <- rbind(promoters,
                         exons,
                         introns,
                         TEs,
                         intergenic)

# Sort
all_annotations <- all_annotations %>% arrange(chr, start)

# Output
write.table(all_annotations, file="merged_annotations.txt",
            sep="\t", quote=F, col.names=T, row.names=F)

