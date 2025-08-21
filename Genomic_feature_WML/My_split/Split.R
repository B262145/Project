# 1. Load required libraries
library(GenomicRanges)
library(rtracklayer)
library(plyranges) # For join_overlap_intersect()

# 2. Input file paths
gene_annot <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/TE_annotation/Pcitri_EarlGrey_5.1.0/GCA_950023065.1_ihPlaCitr1.1.augustus.hints.gff3"
te_annot   <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Genomic_feature_WML/filtered_TEannot.gff"
#Removed mitochondrial chromosome (OX465514.1), comment lines (#), 
# and entries annotated as Low_complexity, Satellite, Simple_repeat, or Unknown from the TE annotation.

# 3. Import GFF3 as GRanges
annot <- import(gene_annot)
te_raw <- import(te_annot)

# 4. Subset by feature type
exons <- annot[annot$type == "exon"]
mrnas <- annot[annot$type == "mRNA"]
genes <- annot[annot$type == "gene"]

# 5. Compute introns (ignore strand)
mrna_union <- reduce(mrnas, ignore.strand=TRUE)
exon_union <- reduce(exons, ignore.strand=TRUE)
introns <- setdiff(mrna_union, exon_union, ignore.strand=TRUE)
intron_union <- reduce(introns, ignore.strand=TRUE)


# 6. Get scaffold max end from GFF
gff_scaffold_ranges <- split(ranges(annot), seqnames(annot))
gff_scaffold_max_end <- sapply(gff_scaffold_ranges, function(x) max(end(x)))

print(gff_scaffold_max_end) # Check the length of each scaffold
# Calculate the total genome size based on GFF annotation.
print(sum(gff_scaffold_max_end)) # The genome size in GFF is 403554726.

# 7. Define promoters (2kb upstream, strand-specific!)
prom_raw <- flank(genes, width=2000, start=TRUE) # strand is considered by default
promoters <- prom_raw[start(prom_raw) >= 1]
promoters <- promoters[width(promoters) > 0]

# Correct - strand promoters to avoid scaffold overflow
end_pos <- end(promoters)
scaffold_ends <- gff_scaffold_max_end[as.character(seqnames(promoters))]
neg_strand <- as.vector(strand(promoters)) == "-"


end_pos[neg_strand] <- pmin(end_pos[neg_strand], scaffold_ends[neg_strand])

end(promoters) <- end_pos

promoter_union <- reduce(promoters, ignore.strand=TRUE) # Reduce with ignore.strand=TRUE for later partitioning

# 7. Define intergenic regions (gaps between genes, ignore strand)
gene_union <- reduce(genes, ignore.strand=TRUE)
gaps_all <- gaps(gene_union, ignore.strand=TRUE)
intergenic <- gaps_all[width(gaps_all) > 0]
intergenic_union <- reduce(intergenic, ignore.strand=TRUE)

# 8. Partition by priority: exon > intron > promoter > intergenic
exon_clean <- exon_union
intron_clean <- setdiff(intron_union, exon_clean, ignore.strand=TRUE)
promoter_clean <- setdiff(promoter_union, c(exon_clean, intron_clean), ignore.strand=TRUE)
inter_clean <- setdiff(intergenic_union, c(exon_clean, intron_clean, promoter_clean), ignore.strand=TRUE)

# 9. Add feature label
exon_clean$feature <- "exon"
intron_clean$feature <- "intron"
promoter_clean$feature <- "promoter"
inter_clean$feature <- "intergenic"

# Calculate the total length of each feature
exon_length <- sum(width(exon_clean))
intron_length <- sum(width(intron_clean))
promoter_length <- sum(width(promoter_clean))
inter_length <- sum(width(inter_clean))

cat("Exon  total length:", exon_length, "\n")
cat("Intron  total length:", intron_length, "\n")
cat("Promoter  total length:", promoter_length, "\n")
cat("Intergenic total length:", inter_length, "\n")

# Function to check if two genomic feature sets have overlaps
check_overlap <- function(gr1, gr2, label1, label2) {
  ov <- findOverlaps(gr1, gr2)
  if (length(ov) == 0) {
    cat(label1, "and", label2, "have no overlaps.\n")
  } else {
    cat(label1, "and", label2, "have", length(ov), "overlaps.\n")
  }
}

# Pairwise comparison of overlaps between the four feature sets
check_overlap(exon_clean, intron_clean, "Exon", "Intron")
check_overlap(exon_clean, promoter_clean, "Exon", "Promoter")
check_overlap(exon_clean, inter_clean, "Exon", "Intergenic")

check_overlap(intron_clean, promoter_clean, "Intron", "Promoter")
check_overlap(intron_clean, inter_clean, "Intron", "Intergenic")

check_overlap(promoter_clean, inter_clean, "Promoter", "Intergenic")
# Those 4 genomic features are mutually exclusive

# 10. Classify TE copies by genomic feature (ignore strand, preserve metadata)

# Exonic TE
te_exon <- join_overlap_intersect(te_raw, exon_clean) 
# join_overlap_intersect() ignores strand information and retains TE metadata (e.g., superfamily).

# This operation triggers a warning because te_raw and exon_clean have different seqlevels (scaffold names).
# - te_raw contains scaffold CATLOI010000001.1, but exon_clean does not.
# - exon_clean contains scaffold CATLOI010000002.1, but after filtering (simple_repeats removed),
#   te_raw has no corresponding TE fragments on CATLOI010000002.1 anymore.
# join_overlap_intersect() will proceed by ignoring the scaffolds that are not shared between the two objects, 
# and only perform intersections on the common scaffolds.

#-------------------------------------------------------------------------------
# # Try this example test:
# # Simulate two GRanges objects
# 
# # TE data: only located on scaffold1 and scaffold3
# test_te_raw <- GRanges(seqnames = c("scaffold1", "scaffold3"),
#                   ranges = IRanges(start = c(100, 300), end = c(200, 400)))
# 
# # Exon data: only located on scaffold2 and scaffold3
# test_exon_clean <- GRanges(seqnames = c("scaffold2", "scaffold3"),
#                       ranges = IRanges(start = c(150, 300), end = c(250, 350)))
# 
# # Check seqlevels (scaffolds/chromosomes) of each object
# seqlevels(test_te_raw)      # "scaffold1", "scaffold3"
# seqlevels(test_exon_clean)  # "scaffold2", "scaffold3"
# 
# # Perform intersection using join_overlap_intersect()
# # Since scaffold1 only exists in te_raw and scaffold2 only exists in exon_clean,
# # join_overlap_intersect() will automatically ignore these scaffolds.
# # It will only perform the intersection on the common scaffold: scaffold3.
# 
# result <- join_overlap_intersect(test_te_raw, test_exon_clean)
# 
# # View the result
# result
# # GRanges object with 1 range and 0 metadata columns:
# # seqnames    ranges strand
# # <Rle> <IRanges>  <Rle>
# #   [1] scaffold3   300-350      *
# #   -------
# #   seqinfo: 2 sequences from an unspecified genome; no seqlengths

#-------------------------------------------------------------------------------

te_exon$feature <- "exonic_TE"

message("Scaffolds included in te_exon: ", paste(seqlevelsInUse(te_exon), collapse = ", "))

# Intronic TE
te_rest1 <- setdiff(te_raw, te_exon, ignore.strand=TRUE)
te_intron <- join_overlap_intersect(te_rest1, intron_clean)
te_intron$feature <- "intronic_TE"

message("Scaffolds included in te_intron: ", paste(seqlevelsInUse(te_intron), collapse = ", "))

# Promoter TE
te_rest2 <- setdiff(te_rest1, te_intron, ignore.strand=TRUE)
te_promoter <- join_overlap_intersect(te_rest2, promoter_clean)
te_promoter$feature <- "promoter_TE"

message("Scaffolds included in te_promoter: ", paste(seqlevelsInUse(te_promoter), collapse = ", "))

# Intergenic TE
te_rest3 <- setdiff(te_rest2, te_promoter, ignore.strand=TRUE)
te_intergenic <- join_overlap_intersect(te_rest3, inter_clean)
te_intergenic$feature <- "intergenic_TE"

message("Scaffolds included in te_intergenic: ", paste(seqlevelsInUse(te_intergenic), collapse = ", "))

# Pairwise comparison of overlaps between the four TE feature sets
check_overlap(te_exon, te_intron, "Exonic_TE", "Intronic_TE") 

check_overlap(te_exon, te_promoter, "Exonic_TE", "Promoter_TE") 

check_overlap(te_exon, te_intergenic, "Exonic_TE", "Intergenic_TE") 

check_overlap(te_intron, te_promoter, "Intronic_TE", "Promoter_TE") 

check_overlap(te_intron, te_intergenic, "Intronic_TE", "Intergenic_TE") 

check_overlap(te_promoter, te_intergenic, "Promoter_TE", "Intergenic_TE") 
# Those 4 TE categories are mutually exclusive


# Calculate the total length of each TE subcategory
te_exon_length <- sum(width(te_exon))
te_intron_length <- sum(width(te_intron))
te_promoter_length <- sum(width(te_promoter))
te_intergenic_length <- sum(width(te_intergenic))

cat("Exonic TE  total length:", te_exon_length, "\n")
cat("Intronic TE  total length:", te_intron_length, "\n")
cat("Promoter TE  total length:", te_promoter_length, "\n")
cat("Intergenic TE total length:", te_intergenic_length, "\n")

# Combine all TE annotations
all_te_classified <- c(te_exon,
                       te_intron, 
                       te_promoter,
                       te_intergenic)

# 11. Combine all features
all_features <- c(all_te_classified,
                  exon_clean,
                  intron_clean,
                  promoter_clean,
                  inter_clean)

# Warning message:
#   In .merge_two_Seqinfo_objects(x, y) :
#   Each of the 2 combined objects has sequence levels not in the other:
#   - in 'x': CATLOI010000001.1
# - in 'y': CATLOI010000002.1
# Make sure to always combine/compare objects based on the same reference
# genome (use suppressWarnings() to suppress this warning).

#-------------------------------------------------------------------------------
# #Test
# # Simulate 5 GRanges objects with different scaffolds
# 
# # TE (all_te_classified) has scaffold1 and scaffold3
# test_all_te_classified <- GRanges(seqnames = c("scaffold1", "scaffold3"),
#                              ranges = IRanges(start = c(100, 300), end = c(200, 400)),
#                              superfamily = c("LINE", "Gypsy"))
# 
# # Exon only has scaffold2 and scaffold3
# test_exon_clean <- GRanges(seqnames = c("scaffold2", "scaffold3"),
#                       ranges = IRanges(start = c(150, 300), end = c(250, 350)),
#                       feature = "exon")
# 
# # Intron only has scaffold2
# test_intron_clean <- GRanges(seqnames = "scaffold2",
#                         ranges = IRanges(start = 400, end = 500),
#                         feature = "intron")
# 
# # Promoter only has scaffold4
# test_promoter_clean <- GRanges(seqnames = "scaffold4",
#                           ranges = IRanges(start = 50, end = 150),
#                           feature = "promoter")
# 
# # Intergenic only has scaffold5
# test_inter_clean <- GRanges(seqnames = "scaffold5",
#                        ranges = IRanges(start = 10, end = 100),
#                        feature = "intergenic")
# 
# # Combine all features
# test_all_features <- c(test_all_te_classified,
#                        test_exon_clean,
#                        test_intron_clean,
#                        test_promoter_clean,
#                        test_inter_clean)
# 
# # Check seqlevels in the combined object
# seqlevels(test_all_features)
# 
# # Check the combined object content
# test_all_features
# # c() retains all scaffolds, even if a scaffold is present in only one GRanges object.
# # The result is:
# # GRanges object with 7 ranges and 2 metadata columns:
# #   seqnames    ranges strand | superfamily     feature
# # <Rle> <IRanges>  <Rle> | <character> <character>
# #   [1] scaffold1   100-200      * |        LINE        <NA>
# #   [2] scaffold3   300-400      * |       Gypsy        <NA>
# #   [3] scaffold2   150-250      * |        <NA>        exon
# # [4] scaffold3   300-350      * |        <NA>        exon
# # [5] scaffold2   400-500      * |        <NA>      intron
# # [6] scaffold4    50-150      * |        <NA>    promoter
# # [7] scaffold5    10-100      * |        <NA>  intergenic
# # -------
# #   seqinfo: 5 sequences from an unspecified genome; no seqlengths
#-------------------------------------------------------------------------------

all_features <- sortSeqlevels(all_features)
all_features <- sort(all_features)

# Convert to data.frame and write

df <- as.data.frame(all_features)
colnames(df)[colnames(df) == "seqnames"] <- "chr"

write.table(df,
            file = "all_features.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

message("Partition complete: see all_features.tsv")

# Double check:
# 1. Calculate the total length of the original genome partitioning
# The genome is divided into four mutually exclusive regions:
# exon, intron, promoter, and intergenic
total_clean <- sum(width(exon_clean)) +
  sum(width(intron_clean)) +
  sum(width(promoter_clean)) +
  sum(width(inter_clean))

# 2: Calculate the total length of TE and non-TE partitions
# For each region (exon, intron, promoter, intergenic), partition it into:
# - The part that overlaps with TEs (te_*)
# - The remaining non-TE part (setdiff(region, te_region))
total_partitioned <- sum(width(te_exon)) + 
  sum(width(setdiff(exon_clean, te_exon, ignore.strand=TRUE))) +
  sum(width(te_intron)) +
  sum(width(setdiff(intron_clean, te_intron, ignore.strand=TRUE))) +
  sum(width(te_promoter)) +
  sum(width(setdiff(promoter_clean, te_promoter, ignore.strand=TRUE))) +
  sum(width(te_intergenic)) +
  sum(width(setdiff(inter_clean, te_intergenic, ignore.strand=TRUE)))

# Compare the partitioned total to:
# - total_clean: the original mutually exclusive feature partition (exon, intron, promoter, intergenic)
# - sum(gff_scaffold_max_end): the total genome size covered by GFF annotation (scaffold max ends)
if (total_clean == total_partitioned && total_clean == sum(gff_scaffold_max_end)) {
  cat("Partitioning is consistent: all totals match the GFF genome size:", sum(gff_scaffold_max_end), "\n")
} else {
  cat("Partitioning inconsistency detected: the total lengths do not match.\n")
  cat("Original total:", total_clean, "\n")
  cat("Partitioned total:", total_partitioned, "\n")
  cat("GFF genome size:", sum(gff_scaffold_max_end), "\n")
}

# All totals match: total_clean, total_partitioned, and GFF genome size are 403554726