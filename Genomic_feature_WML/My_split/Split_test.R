# Splits Planococcus GFF3 + TE annotation into mutually exclusive features
# Following the order: TE > exon > intron > promoter > intergenic

# 1. Load required libraries
library(GenomicRanges)
library(rtracklayer)

# 2. Input file paths
gene_annot <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.augustus.hints.gff3"
te_annot    <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Genomic_feature_WML/filtered_TEannot.gff"

# 3. Import GFF3 as a GRanges object
annot <- import(gene_annot)
te_raw  <- import(te_annot)

# 4. Subset by feature type
exons      <- annot[annot$type == "exon"]
mrnas      <- annot[annot$type == "mRNA"]
genes      <- annot[annot$type == "gene"]

# 5. Compute introns: mRNAs minus exons
mrna_union <- reduce(mrnas)
exon_union <- reduce(exons)
introns    <- setdiff(mrna_union, exon_union)
intron_union <- reduce(introns)

# 6. Define promoters (2 kb upstream) manually via flank()
#    flank(x, width, start=TRUE) gives upstream regions for both strands
prom_raw <- flank(genes, width=2000, start=TRUE)
# Remove any that start before base 1
promoters <- prom_raw[start(prom_raw) >= 1]
promoters <- promoters[width(promoters) > 0]
promoter_union <- reduce(promoters)

# 7. Define intergenic regions as gaps between gene ranges
gene_union <- reduce(genes)
gaps_all    <- gaps(gene_union)
intergenic  <- gaps_all[width(gaps_all) > 0]
intergenic_union <- reduce(intergenic)

# 8. Partition by priority to avoid overlap: exon > intron > promoter > intergenic
exon_clean <- exon_union
intron_clean <- setdiff(intron_union, exon_clean)
promoter_clean <- setdiff(promoter_union, c(exon_clean, intron_clean))
inter_clean <- setdiff(intergenic_union, c(exon_clean, intron_clean, promoter_clean))

# 9. Add context label
exon_clean$context    <- "exon"
intron_clean$context  <- "intron"
promoter_clean$context    <- "promoter"
inter_clean$context   <- "intergenic"

# 10. Classify TE copy by genomic context
# Exonic TE
te_exon <- intersect(te_raw, exon_clean)
te_exon$context <- "exonic_TE"

# This operation triggers a warning because te_raw and exon_clean have different seqlevels (scaffold names).
# - te_raw contains scaffold CATLOI010000001.1, but exon_clean does not.
# - exon_clean contains scaffold CATLOI010000002.1, but after filtering (simple_repeats removed),
#   te_raw has no corresponding TE fragments on CATLOI010000002.1 anymore.
# GenomicRanges::intersect() will proceed by ignoring the scaffolds that are not shared between the two objects, 
# and only perform intersections on the common scaffolds.
# To avoid this warning, it's recommended to restrict both te_raw and exon_clean to the same set of scaffolds 
# using keepSeqlevels() before performing the intersection.

# Output the scaffolds actually used in te_exon
message("Scaffolds included in te_exon: ", paste(seqlevelsInUse(te_exon), collapse = ", "))

# Intronic TE
te_rest1 <- setdiff(te_raw, te_exon)
te_intron <- intersect(te_rest1, intron_clean) # By default, intersect() considers strand info
te_intron$context <- "intronic_TE"

# Output the scaffolds actually used in te_intron
message("Scaffolds included in te_intron: ", paste(seqlevelsInUse(te_intron), collapse = ", "))

# Promoter TE
te_rest2 <- setdiff(te_rest1, te_intron)
te_promoter <- intersect(te_rest2, promoter_clean)
te_promoter$context <- "promoter_TE"

# Output the scaffolds actually used in te_promoter
message("Scaffolds included in te_promoter: ", paste(seqlevelsInUse(te_promoter), collapse = ", "))

# Intergenic TE
te_rest3 <- setdiff(te_rest2, te_promoter)
te_intergenic <- intersect(te_rest3, inter_clean)
te_intergenic$context <- "intergenic_TE"

# Output the scaffolds actually used in te_intergenic
message("Scaffolds included in te_intergenic: ", paste(seqlevelsInUse(te_intergenic), collapse = ", "))

# Combine all
all_te_classified <- c(te_exon,
                       te_intron, 
                       te_promoter,
                       te_intergenic)

# 11. Combine, sort, and export
all_features <- c(all_te_classified,
                  exon_clean,
                  intron_clean,
                  promoter_clean,
                  inter_clean)
all_features <- sortSeqlevels(all_features)
all_features <- sort(all_features)

# Convert to data.frame and write

write.table(as.data.frame(all_features),
            file = "v1_all_features.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

message("Partition complete: see all_features.tsv")
