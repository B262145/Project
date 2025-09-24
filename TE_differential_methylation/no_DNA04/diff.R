# TE Copy Level Differential Methylation Analysis using edgeR
# === Key modifications ===
# 1) Remove DNA04
# 2) Eliminate hard-coded indices — select samples dynamically by group
# ==========================================

suppressPackageStartupMessages({
  library(edgeR)
  library(GenomicRanges)
  library(rtracklayer)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
})

# ---------------------------
# 1. Sample information (DNA04 removed)
# ---------------------------
cov_dir      <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Cov_files/"
all_samples  <- sprintf("DNA%02d", 1:24)

# Metadata in the original order
sex_vec   <- c(rep("F",6), rep("M",6), rep("F",6), rep("M",6))
stage_vec <- c(rep("Adult",12), rep("Third_instar",12))

samples_df <- data.frame(
  Sample = all_samples,
  Sex    = factor(sex_vec,   levels=c("F","M")),
  Stage  = factor(stage_vec, levels=c("Adult","Third_instar")),
  File   = file.path(cov_dir, paste0(all_samples, ".cov.txt")),
  stringsAsFactors = FALSE
)

# Remove DNA04
samples_df <- subset(samples_df, Sample != "DNA04")
row.names(samples_df) <- NULL
print(samples_df)
cat("\nGroup sizes (Stage x Sex):\n"); print(table(samples_df$Stage, samples_df$Sex))

# ---------------------------
# 2. Read Bismark coverage
# ---------------------------
yall <- readBismark2DGE(samples_df$File, sample.names=samples_df$Sample)
head(yall)

# ---------------------------
# 3. TE annotations
# ---------------------------
te_gr <- import("/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/TE_differential_methylation/filtered_no_unplaced_TE.gff")

# ---------------------------
# 4. Summarize CpG -> TE regions
# ---------------------------
cpg_gr <- GRanges(seqnames=yall$genes$Chr,
                  ranges=IRanges(start=yall$genes$Locus, width=1))

hits <- findOverlaps(cpg_gr, te_gr)

counts_te <- sapply(1:ncol(yall), function(i) {
  counts_i <- yall$counts[, i]
  matched_counts <- counts_i[queryHits(hits)]
  te_ids <- subjectHits(hits)
  tapply(matched_counts, te_ids, sum, na.rm=TRUE)
})

rownames(counts_te) <- unique(subjectHits(hits))
te_used <- te_gr[as.integer(rownames(counts_te))]
colnames(counts_te) <- colnames(yall)

te_yall <- DGEList(
  counts = counts_te,
  genes  = as.data.frame(te_used)[, c("seqnames","start","end","type")]
)

cat("Number of TE rows (pre-filter): ", nrow(te_yall), "\n")

# ---------------------------
# 5. Filtering (coverage + contains both Me and Un)
# ---------------------------
Methylation <- gl(2, 1, ncol(yall), labels=c("Me","Un"))
Me <- te_yall$counts[, Methylation=="Me"]
Un <- te_yall$counts[, Methylation=="Un"]
Cov <- Me + Un

HasCoverage <- rowSums(Cov >= 5) == ncol(Me)  # coverage >=5 in every sample
HasBoth     <- rowSums(Me) > 0 & rowSums(Un) > 0

print(table(HasCoverage, HasBoth))
te_yall <- te_yall[HasCoverage & HasBoth,, keep.lib.sizes=FALSE]
cat("Number of TE rows (post-filter): ", nrow(te_yall), "\n")

# ---------------------------
# 6. Update library sizes
# ---------------------------
TotalLibSize <- 0.5 * te_yall$samples$lib.size[Methylation=="Me"] +
                0.5 * te_yall$samples$lib.size[Methylation=="Un"]
te_yall$samples$lib.size <- rep(TotalLibSize, each=2)
te_yall$samples

# ---------------------------
# 7. Data exploration (MDS)
# ---------------------------
Me_DE <- te_yall$counts[, Methylation=="Me"]
Un_DE <- te_yall$counts[, Methylation=="Un"]
M_DE  <- log2(Me_DE + 2) - log2(Un_DE + 2)

# Align column names with samples_df
colnames(M_DE) <- samples_df$Sample

# Dynamic groups
samples_df$Group <- factor(
  paste0(samples_df$Stage, ".", samples_df$Sex),
  levels=c("Adult.F","Adult.M","Third_instar.F","Third_instar.M")
)

group_map <- c("Adult.F"="Adult Female", "Adult.M"="Adult Male",
               "Third_instar.F"="3rd Instar Female", "Third_instar.M"="3rd Instar Male")
group.labels <- group_map[as.character(samples_df$Group)]

group.colors <- c(
  "Adult Female"      = "#00bec4",
  "Adult Male"        = "#c67bff",
  "3rd Instar Female" = "#f8766c",
  "3rd Instar Male"   = "#7bad00"
)
sample.colors <- group.colors[group.labels]

pdf("TE_copy_level_MDS_dim1_2.pdf", width=7, height=7)
plotMDS(M_DE, col=sample.colors, main="TE copy level M-values MDS: dim 1&2")
legend("bottomright", fill=group.colors, legend=names(group.colors), cex=0.8)
dev.off()

pdf("TE_copy_level_MDS_dim3_4.pdf", width=7, height=7)
plotMDS(M_DE, dim=c(3,4), col=sample.colors, main="TE copy level M-values MDS: dim 3&4")
legend("topleft", fill=group.colors, legend=names(group.colors), cex=0.8)
dev.off()

# Helpers: select sample IDs by pretty group label; append "-Me"
samples_by_label <- function(df, pretty_label) {
  inv_map <- setNames(names(group_map), group_map)
  machine_group <- inv_map[pretty_label]
  df$Sample[df$Group == machine_group]
}
samples_me <- function(x) paste0(x, "-Me")

# ---------------------------
# 8. Design matrix (unbalanced design supported)
# ---------------------------
designSL <- model.matrix(~0 + Group, data=samples_df)
colnames(designSL) <- sub("^Group", "", colnames(designSL))  # Adult.F, Adult.M, ...
design <- modelMatrixMeth(designSL)

fit <- glmQLFit(te_yall, design)
cat("Design columns:\n"); print(colnames(design))

# Contrasts
my.contrasts <- makeContrasts(
  AdultM_vs_adultF = Adult.M - Adult.F,
  ThirdF_vs_AdultF = Third_instar.F - Adult.F,
  ThirdM_vs_adultM = Third_instar.M - Adult.M,
  ThirdM_vs_ThirdF = Third_instar.M - Third_instar.F,
  M_vs_F           = 0.5*(Adult.M + Third_instar.M) - 0.5*(Adult.F + Third_instar.F),
  Third_vs_Adult   = 0.5*(Third_instar.M + Third_instar.F) - 0.5*(Adult.M + Adult.F),
  levels = design
)

# ---------------------------
# 9. Heatmap function
# ---------------------------
plot_TE_heatmap <- function(contrast_table, te_yall, output_pdf, sample_order, sample_labels,
                            class_order = c("DNA", "LINE", "LTR", "RC", "SINE")) {
  contrast_table$TE_ID <- paste0(contrast_table$seqnames, ":", contrast_table$start, "-", contrast_table$end)
  contrast_table <- contrast_table[!is.na(contrast_table$FDR), ]
  if (nrow(contrast_table) == 0) {
    warning("No rows in contrast_table for heatmap: ", output_pdf)
    return(invisible(NULL))
  }
  topN <- min(100, nrow(contrast_table))
  top100 <- contrast_table[order(contrast_table$FDR), ][1:topN, ]

  Methylation <- gl(2, 1, ncol(te_yall$counts), labels=c("Me","Un"))
  Me_TE <- te_yall$counts[, Methylation=="Me"]
  Un_TE <- te_yall$counts[, Methylation=="Un"]
  M_TE <- log2(Me_TE + 2) - log2(Un_TE + 2)

  rownames(M_TE) <- paste0(te_yall$genes$seqnames, ":", te_yall$genes$start, "-", te_yall$genes$end)

  keep_ids <- intersect(top100$TE_ID, rownames(M_TE))
  if (length(keep_ids) == 0) {
    warning("No matching TE IDs for heatmap: ", output_pdf)
    return(invisible(NULL))
  }

  M_top <- M_TE[keep_ids, , drop=FALSE]
  M_top <- M_top[, sample_order, drop=FALSE]
  colnames(M_top) <- sample_labels

  superfamily_info <- te_yall$genes$type
  names(superfamily_info) <- rownames(M_TE)
  superfamily_top <- superfamily_info[rownames(M_top)]

  rownames(M_top) <- paste0(rownames(M_top), " [", superfamily_top, "]")
  annotation_row <- data.frame(Superfamily = superfamily_top)
  rownames(annotation_row) <- rownames(M_top)

  # Build superfamily color map with ordered TE classes
  superfamily_char <- as.character(superfamily_top)
  split_superfamily <- strsplit(superfamily_char, "/")
  superfamily_df <- do.call(rbind, lapply(split_superfamily, function(x) {
    if (length(x) == 2) data.frame(class=x[1], superfamily=x[2], full_name=paste(x, collapse="/"))
    else data.frame(class=NA, superfamily=NA, full_name=NA)
  }))
  superfamily_df$class_rank <- match(superfamily_df$class, class_order)
  sort_df <- superfamily_df; sort_df$full_name <- superfamily_char
  sorted_superfamilies <- unique(sort_df[order(sort_df$class_rank, sort_df$superfamily), "full_name"])

  superfamily_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(sorted_superfamilies))
  names(superfamily_colors) <- sorted_superfamilies
  ann_colors <- list(Superfamily = superfamily_colors)

  comparison <- gsub("TE_top100_heatmap_(.*).pdf", "\\1", basename(output_pdf))
  comparison <- gsub("_", " ", comparison)
  main_title <- paste("Top", topN, "Differentially Methylated TE Copies:", comparison)

  pdf(output_pdf, width=10, height=10)
  pheatmap(M_top,
           annotation_row = annotation_row,
           annotation_colors = ann_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           fontsize_row = 5,
           angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(n=9, "RdBu")))(100),
           main = main_title)
  dev.off()
}

# ---------------------------
# 10. Superfamily enrichment functions
# ---------------------------
run_TE_superfamily_enrichment <- function(input_csv, output_prefix) {
  all_TE <- read.csv(input_csv, stringsAsFactors = FALSE)
  all_TE <- all_TE[!is.na(all_TE$FDR), ]
  de_TE  <- subset(all_TE, FDR < 0.05)

  all_super <- all_TE$type
  de_super  <- de_TE$type
  all_classes <- unique(all_super)

  enrich_results <- data.frame(Superfamily=character(),
                               DE_Count=integer(),
                               Background_Count=integer(),
                               p_value=numeric(),
                               stringsAsFactors=FALSE)

  for (sf in all_classes) {
    a <- sum(de_super == sf)
    b <- sum(de_super != sf)
    c <- sum(all_super == sf) - a
    d <- sum(all_super != sf) - b
    mat <- matrix(c(a, c, b, d), nrow=2, byrow=TRUE)
    fisher_res <- fisher.test(mat, alternative="greater")
    enrich_results <- rbind(enrich_results,
                            data.frame(Superfamily=sf,
                                       DE_Count=a,
                                       Background_Count=a+c,
                                       p_value=fisher_res$p.value))
  }

  enrich_results$FDR <- p.adjust(enrich_results$p_value, method="BH")
  write.csv(enrich_results, paste0(output_prefix, "_TE_superfamily_enrichment.csv"), row.names=FALSE)

  sig_enrich <- subset(enrich_results, FDR < 0.05)
  sig_enrich <- sig_enrich[order(sig_enrich$FDR), ]
  write.csv(sig_enrich, paste0(output_prefix, "_TE_superfamily_enrichment_sig.csv"), row.names=FALSE)

  message("Enrichment analysis completed for ", output_prefix)
}

run_TE_superfamily_enrichment_hyper <- function(input_csv, output_prefix) {
  all_TE <- read.csv(input_csv, stringsAsFactors = FALSE)
  all_TE <- all_TE[!is.na(all_TE$FDR), ]
  hyper_TE <- subset(all_TE, logFC > 0 & FDR < 0.05)

  all_super  <- all_TE$type
  hyper_super <- hyper_TE$type
  all_classes <- unique(all_super)

  enrich_results <- data.frame(Superfamily=character(),
                               Hyper_Count=integer(),
                               Background_Count=integer(),
                               p_value=numeric(),
                               stringsAsFactors=FALSE)

  for (sf in all_classes) {
    a <- sum(hyper_super == sf)
    b <- sum(hyper_super != sf)
    c <- sum(all_super == sf) - a
    d <- sum(all_super != sf) - b
    mat <- matrix(c(a, c, b, d), nrow=2, byrow=TRUE)
    fisher_res <- fisher.test(mat, alternative="greater")
    enrich_results <- rbind(enrich_results,
                            data.frame(Superfamily=sf,
                                       Hyper_Count=a,
                                       Background_Count=a+c,
                                       p_value=fisher_res$p.value))
  }

  enrich_results$FDR <- p.adjust(enrich_results$p_value, method="BH")
  write.csv(enrich_results, paste0(output_prefix, "_TE_superfamily_enrichment_hyper.csv"), row.names=FALSE)

  sig <- subset(enrich_results, FDR < 0.05)
  sig <- sig[order(sig$FDR), ]
  write.csv(sig, paste0(output_prefix, "_TE_superfamily_enrichment_hyper_sig.csv"), row.names=FALSE)

  message("Hyper enrichment completed for ", output_prefix)
}

run_TE_superfamily_enrichment_hypo <- function(input_csv, output_prefix) {
  all_TE <- read.csv(input_csv, stringsAsFactors = FALSE)
  all_TE <- all_TE[!is.na(all_TE$FDR), ]
  hypo_TE <- subset(all_TE, logFC < 0 & FDR < 0.05)

  all_super  <- all_TE$type
  hypo_super <- hypo_TE$type
  all_classes <- unique(all_super)

  enrich_results <- data.frame(Superfamily=character(),
                               Hypo_Count=integer(),
                               Background_Count=integer(),
                               p_value=numeric(),
                               stringsAsFactors=FALSE)

  for (sf in all_classes) {
    a <- sum(hypo_super == sf)
    b <- sum(hypo_super != sf)
    c <- sum(all_super == sf) - a
    d <- sum(all_super != sf) - b
    mat <- matrix(c(a, c, b, d), nrow=2, byrow=TRUE)
    fisher_res <- fisher.test(mat, alternative="greater")
    enrich_results <- rbind(enrich_results,
                            data.frame(Superfamily=sf,
                                       Hypo_Count=a,
                                       Background_Count=a+c,
                                       p_value=fisher_res$p.value))
  }

  enrich_results$FDR <- p.adjust(enrich_results$p_value, method="BH")
  write.csv(enrich_results, paste0(output_prefix, "_TE_superfamily_enrichment_hypo.csv"), row.names=FALSE)

  sig <- subset(enrich_results, FDR < 0.05)
  sig <- sig[order(sig$FDR), ]
  write.csv(sig, paste0(output_prefix, "_TE_superfamily_enrichment_hypo_sig.csv"), row.names=FALSE)

  message("Hypo enrichment completed for ", output_prefix)
}

# ---------------------------
# 11. A helper to run a full contrast block
# ---------------------------
run_full_contrast <- function(qlf, contrast_name_pretty, out_prefix,
                              sample_groups_for_heatmap) {
  # qlf -> CSV
  out_csv <- paste0("TE_DMA_", out_prefix, ".csv")
  write.csv(topTags(qlf, n=Inf)$table, out_csv, row.names=FALSE)

  # summary (at 5% FDR)
  sum_txt <- paste0("TE_DMA_", out_prefix, "_summary.txt")
  sum_decide <- summary(decideTests(qlf))
  capture.output(
    list(
      header  = paste0("=== ", contrast_name_pretty, " at FDR 5% ==="),
      summary = sum_decide
    ),
    file = sum_txt
  )

  # MD plot
  pdf(paste0("TE_DMA_", out_prefix, "_MDplot.pdf"), width=6, height=6)
  plotMD(qlf, main = contrast_name_pretty)
  dev.off()

  # Heatmap (top 100 by FDR) on Me channels of selected groups
  contrast_table <- read.csv(out_csv, stringsAsFactors = FALSE)

  # Build dynamic sample order and labels
  sample_order <- character(0)
  sample_labels <- character(0)
  for (pretty_label in sample_groups_for_heatmap) {
    ids <- samples_by_label(samples_df, pretty_label)
    sample_order  <- c(sample_order,  samples_me(ids))
    sample_labels <- c(sample_labels, paste0(pretty_label, " ", ids))
  }

  plot_TE_heatmap(contrast_table, te_yall,
                  paste0("TE_top100_heatmap_", out_prefix, ".pdf"),
                  sample_order, sample_labels)

  # Volcano
  df <- contrast_table
  pdf(paste0("Volcano_", out_prefix, ".pdf"), width=8, height=8)
  EnhancedVolcano(df,
                  lab = rep("", nrow(df)),
                  x = 'logFC',
                  y = 'FDR',
                  title = contrast_name_pretty,
                  xlab = bquote(Log[2]~Fold~Change),
                  ylab = bquote(-Log[10]~FDR),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 2.5,
                  labSize = 0,
                  drawConnectors = FALSE,
                  legendPosition = 'right',
                  legendLabSize = 12,
                  legendIconSize = 4.0)
  dev.off()

  # Enrichment analyses
  run_TE_superfamily_enrichment(out_csv, out_prefix)
  run_TE_superfamily_enrichment_hyper(out_csv, out_prefix)
  run_TE_superfamily_enrichment_hypo (out_csv, out_prefix)
}

# ---------------------------
# 12. Run ALL contrasts
# ---------------------------

# Adult Male vs Adult Female
qlf_AdultM_vs_adultF <- glmQLFTest(fit, contrast=my.contrasts[,"AdultM_vs_adultF"])
run_full_contrast(
  qlf_AdultM_vs_adultF,
  "Adult Male vs Adult Female",
  "AdultM_vs_AdultF",
  c("Adult Female", "Adult Male")
)

# 3rd Instar Female vs Adult Female
qlf_3rdF_vs_AdultF <- glmQLFTest(fit, contrast=my.contrasts[,"ThirdF_vs_AdultF"])
run_full_contrast(
  qlf_3rdF_vs_AdultF,
  "3rd Instar Female vs Adult Female",
  "3rdF_vs_AdultF",
  c("Adult Female", "3rd Instar Female")
)

# 3rd Instar Male vs Adult Male
qlf_3rdM_vs_adultM <- glmQLFTest(fit, contrast=my.contrasts[,"ThirdM_vs_adultM"])
run_full_contrast(
  qlf_3rdM_vs_adultM,
  "3rd Instar Male vs Adult Male",
  "3rdM_vs_AdultM",
  c("Adult Male", "3rd Instar Male")
)

# 3rd Instar Male vs 3rd Instar Female
qlf_3rdM_vs_3rdF <- glmQLFTest(fit, contrast=my.contrasts[,"ThirdM_vs_ThirdF"])
run_full_contrast(
  qlf_3rdM_vs_3rdF,
  "3rd Instar Male vs 3rd Instar Female",
  "3rdM_vs_3rdF",
  c("3rd Instar Female", "3rd Instar Male")
)

# Average Male vs Average Female
qlf_M_vs_F <- glmQLFTest(fit, contrast=my.contrasts[,"M_vs_F"])
run_full_contrast(
  qlf_M_vs_F,
  "Average Male vs Average Female",
  "M_vs_F",
  c("Adult Female", "Adult Male", "3rd Instar Female", "3rd Instar Male")
)

# Average 3rd instar vs Average Adult
qlf_Third_vs_Adult <- glmQLFTest(fit, contrast=my.contrasts[,"Third_vs_Adult"])
run_full_contrast(
  qlf_Third_vs_Adult,
  "Average 3rd Instar vs Average Adult",
  "Third_vs_Adult",
  c("Adult Female", "Adult Male", "3rd Instar Female", "3rd Instar Male")
)

# ---------------------------
# 13. Done
# ---------------------------
message("All differential methylation analysis results have been saved as CSV/PDF files!")
