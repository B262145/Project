library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

# Load scaffold lengths from FNA
fai <- scanFaIndex("/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Genomic_feature_WML/genome_filtered.fna")

# Get scaffold lengths
fna_scaffold_lengths <- seqlengths(fai)

# View result
fna_scaffold_lengths
print(sum(fna_scaffold_lengths))

# GFF
gene_annot <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.augustus.hints.gff3"
annot <- import(gene_annot)

gff_scaffold_ranges <- split(ranges(annot), seqnames(annot))
gff_scaffold_lengths <- sapply(gff_scaffold_ranges, function(x) max(end(x)))
gff_scaffold_lengths

print(sum(gff_scaffold_lengths))