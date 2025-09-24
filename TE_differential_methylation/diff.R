# TE Copy Level Differential Methylation Analysis using edgeR

# Load required packages
library(edgeR)
library(GenomicRanges)
library(rtracklayer)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

# 1. Prepare file paths, sample names, and sample information
cov_dir      <- "/mnt/loki/ross/scales/pseudococcidae/Planococcus_citri/Yuhan/WGBS_newassembly/Cov_files/"
sample_names <- sprintf("DNA%02d", 1:24)
files        <- file.path(cov_dir, paste0(sample_names, ".cov.txt"))
sex          <- factor(c(rep("F",6), rep("M",6), rep("F",6), rep("M",6)))
stage        <- factor(c(rep("Adult",12), rep("Third_instar",12)))
samples_df   <- data.frame(Sample=sample_names, Sex=sex, Stage=stage, File=files, stringsAsFactors=FALSE)

# Print the sample information table
print(samples_df)

# 2. Read Bismark coverage data and extract CpG counts
yall <- readBismark2DGE(samples_df$File, sample.names=samples_df$Sample)
head(yall)

# 3. Load TE annotations as GRanges object
te_gr <- import("filtered_no_unplaced_TE.gff")

# 4. Summarize CpG methylation counts over TE regions
# Create CpG positions as GRanges
cpg_gr <- GRanges(seqnames=yall$genes$Chr,
                  ranges=IRanges(start=yall$genes$Locus,width=1))
                  
# Find overlaps between CpG sites and TEs
hits <- findOverlaps(cpg_gr, te_gr)

# Sum methylation counts for CpGs overlapping each TE for each sample
counts_te <- sapply(1:ncol(yall), function(i) {
  counts_i <- yall$counts[, i]
  matched_counts <- counts_i[queryHits(hits)]         # CpG count
  te_ids         <- subjectHits(hits)                 # TE index
  tapply(matched_counts, te_ids, sum, na.rm = TRUE)   
})
# The counts_te is a matrix where each row corresponds to a TE copy and each column corresponds to a sample. 
# Each value in the matrix represents the total CpG count for that TE in that sample (either methylated or unmethylated reads).
head(counts_te)

# Set row names to TE indices
# This line sets the row names of the counts_te matrix to the unique indices of the TEs in te_gr 
# that have at least one overlapping CpG. This way, each row of counts_te can be directly matched 
# to the corresponding TE region in the annotation.
rownames(counts_te) <- unique(subjectHits(hits))

# Extract the TEs actually used (overlapped)
# Extract from te_gr the TE copies whose indices match the row names of counts_te
te_used <- te_gr[as.integer(rownames(counts_te))]

colnames(counts_te) <- colnames(yall)

# Create DGEList for edgeR using TE counts
te_yall <- DGEList(counts = counts_te,
                   genes = as.data.frame(te_used)[, c('seqnames','start','end', 'type')]
                   )

nrow(te_yall)                 

# 5. Filter TE copies for sufficient coverage
Methylation <- gl(2, 1, ncol(yall), labels = c("Me", "Un"))
Me   <- te_yall$counts[,Methylation=="Me"]
Un   <- te_yall$counts[,Methylation=="Un"]
Cov  <- Me + Un
HasCoverage <- rowSums(Cov >= 5) == ncol(Me) # Set the threshold as 5
HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0
table(HasCoverage, HasBoth)
te_yall   <- te_yall[HasCoverage & HasBoth,,keep.lib.sizes=FALSE]

nrow(te_yall) 

# 6. Update library sizes for filtered TEs
TotalLibSize <- 0.5 * te_yall$samples$lib.size[Methylation=="Me"] +
0.5 * te_yall$samples$lib.size[Methylation=="Un"]

te_yall$samples$lib.size <- rep(TotalLibSize, each=2)
te_yall$samples

# 7. Data exploration
Me_DE <- te_yall$counts[, Methylation=="Me"]
Un_DE <- te_yall$counts[, Methylation=="Un"]
M_DE <- log2(Me_DE + 2) - log2(Un_DE + 2)
colnames(M_DE) <- samples_df$Sample

# Define colors and groups
group.labels <- rep(c("Adult Female", "Adult Male", "3rd Instar Female", "3rd Instar Male"), each=6)

group.colors <- c(
  "Adult Female" = "#00bec4",
  "Adult Male" = "#c67bff",
  "3rd Instar Female" = "#f8766c",
  "3rd Instar Male" = "#7bad00"
)

sample.colors <- group.colors[group.labels]

# Plot dim 1&2
pdf("TE_copy_level_MDS_dim1_2.pdf", width=7, height=7)
plotMDS(M_DE, col=sample.colors, main="TE copy level M-values MDS: dim 1&2")
legend("bottomright", fill=group.colors, legend=names(group.colors), cex=0.8)
dev.off()

# Plot dim 3&4
pdf("TE_copy_level_MDS_dim3_4.pdf", width=7, height=7)
plotMDS(M_DE, dim=c(3,4), col=sample.colors, main="TE copy level M-values MDS: dim 3&4")
legend("topleft", fill=group.colors, legend=names(group.colors), cex=0.8)
dev.off()

# Define the heatmap function
plot_TE_heatmap <- function(contrast_table, te_yall, output_pdf, sample_order, sample_labels, class_order = c("DNA", "LINE", "LTR", "RC", "SINE")) {
  
  # Construct TE ID as seqnames:start-end
  contrast_table$TE_ID <- paste0(contrast_table$seqnames, ":", contrast_table$start, "-", contrast_table$end)
  
  # Select top 100 TE copies by FDR
  top100 <- contrast_table[order(contrast_table$FDR), ][1:100, ]
  
  # Compute M-values (log2 methylated/unmethylated counts)
  Methylation <- gl(2, 1, ncol(te_yall$counts), labels = c("Me", "Un"))
  Me_TE <- te_yall$counts[, Methylation=="Me"]
  Un_TE <- te_yall$counts[, Methylation=="Un"]
  M_TE <- log2(Me_TE + 2) - log2(Un_TE + 2)
  
  # Assign TE IDs to rows of M_TE
  rownames(M_TE) <- paste0(te_yall$genes$seqnames, ":", te_yall$genes$start, "-", te_yall$genes$end)
  
  # Subset top 100 TE copies
  M_top100 <- M_TE[top100$TE_ID, ]
  
  # Select samples and reformat sample names
  M_top100 <- M_top100[, sample_order]
  colnames(M_top100) <- sample_labels
  
  # Extract superfamily annotation
  superfamily_info <- te_yall$genes$type
  names(superfamily_info) <- rownames(M_TE)
  superfamily_top100 <- superfamily_info[top100$TE_ID]
  
  # Generate new rownames with TE ID and superfamily
  rownames(M_top100) <- paste0(top100$TE_ID, " [", superfamily_top100, "]")
  
  # Create annotation data frame for color bar (superfamily)
  annotation_row <- data.frame(Superfamily = superfamily_top100)
  rownames(annotation_row) <- rownames(M_top100)
  
  # Define TE class hierarchy for better color ordering
  superfamily_char <- as.character(superfamily_top100)
  split_superfamily <- strsplit(superfamily_char, "/")
  
  superfamily_df <- do.call(rbind, lapply(split_superfamily, function(x) {
    if (length(x) == 2) {
      return(data.frame(class = x[1], superfamily = x[2], full_name = paste(x[1], x[2], sep="/")))
    } else {
      return(data.frame(class = NA, superfamily = NA, full_name = NA))
    }
  }))
  
  # Assign rank to TE classes for sorting
  superfamily_df$class_rank <- match(superfamily_df$class, class_order)
  
  # Sort by class rank, then alphabetically within class
  sort_df <- superfamily_df
  sort_df$full_name <- superfamily_char
  
  sorted_superfamilies <- sort_df[order(sort_df$class_rank, sort_df$superfamily), "full_name"]
  sorted_superfamilies <- unique(sorted_superfamilies)
  
  # Assign colors based on sorted order
  superfamily_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(sorted_superfamilies))
  names(superfamily_colors) <- sorted_superfamilies
  
  # Prepare annotation color list
  ann_colors <- list(Superfamily = superfamily_colors)
  
  # Extract comparison name from file name
  comparison <- gsub("TE_top100_heatmap_(.*).pdf", "\\1", basename(output_pdf))
  comparison <- gsub("_", " ", comparison)

  # Define main title
  main_title <- paste("Top 100 Differentially Methylated TE Copies:", comparison)
  # Plot heatmap
  pdf(output_pdf, width=10, height=10)
  
  pheatmap(M_top100,
           annotation_row = annotation_row,
           annotation_colors = ann_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",                      
           fontsize_row = 5,
           angle_col = 45,                     
           color = colorRampPalette(rev(brewer.pal(n=9, "RdBu")))(100),
           main = main_title
  )
  
  dev.off()
}

# Define the Fisher's Exact Test Function - no direction
run_TE_superfamily_enrichment <- function(input_csv, output_prefix) {
  
  all_TE <- read.csv(input_csv, stringsAsFactors = FALSE)
  
  # Select the DM-TEs (FDR < 0.05)
  de_TE <- subset(all_TE, FDR < 0.05)
  
  # superfamily
  all_super <- all_TE$type
  de_super <- de_TE$type
  
  all_classes <- unique(all_super)
  
  enrich_results <- data.frame(Superfamily=character(),
                               DE_Count=integer(),
                               Background_Count=integer(),
                               p_value=numeric(),
                               stringsAsFactors=FALSE)
  
  # Fisher Exact Test (greater)
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
  
  # FDR correction
  enrich_results$FDR <- p.adjust(enrich_results$p_value, method="BH")
  write.csv(enrich_results,
            paste0(output_prefix, "_TE_superfamily_enrichment.csv"),
            row.names=FALSE)
  
  # FDR < 0.05
  sig_enrich <- subset(enrich_results, FDR < 0.05)
  sig_enrich <- sig_enrich[order(sig_enrich$FDR), ]
  
  write.csv(sig_enrich,
            paste0(output_prefix, "_TE_superfamily_enrichment_sig.csv"),
            row.names=FALSE)
  
  message("Enrichment analysis completed for ", output_prefix)
}



# Function: Superfamily enrichment for hypermethylated DM-TEs
run_TE_superfamily_enrichment_hyper <- function(input_csv, output_prefix) {
  all_TE <- read.csv(input_csv, stringsAsFactors = FALSE)
  
  # Filter hypermethylated DM-TEs (logFC > 0 & FDR < 0.05)
  hyper_TE <- subset(all_TE, logFC > 0 & FDR < 0.05)
  
  all_super <- all_TE$type
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
                                       Background_Count=a + c,
                                       p_value=fisher_res$p.value))
  }
  
  enrich_results$FDR <- p.adjust(enrich_results$p_value, method="BH")
  write.csv(enrich_results,
            paste0(output_prefix, "_TE_superfamily_enrichment_hyper.csv"),
            row.names=FALSE)
  
  sig <- subset(enrich_results, FDR < 0.05)
  sig <- sig[order(sig$FDR), ]
  write.csv(sig,
            paste0(output_prefix, "_TE_superfamily_enrichment_hyper_sig.csv"),
            row.names=FALSE)
  message("Hyper enrichment completed for ", output_prefix)
}

# Function: Superfamily enrichment for hypomethylated DM-TEs
run_TE_superfamily_enrichment_hypo <- function(input_csv, output_prefix) {
  all_TE <- read.csv(input_csv, stringsAsFactors = FALSE)
  
  # Filter hypomethylated DM-TEs (logFC < 0 & FDR < 0.05)
  hypo_TE <- subset(all_TE, logFC < 0 & FDR < 0.05)
  
  all_super <- all_TE$type
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
                                       Background_Count=a + c,
                                       p_value=fisher_res$p.value))
  }
  
  enrich_results$FDR <- p.adjust(enrich_results$p_value, method="BH")
  write.csv(enrich_results,
            paste0(output_prefix, "_TE_superfamily_enrichment_hypo.csv"),
            row.names=FALSE)
  
  sig <- subset(enrich_results, FDR < 0.05)
  sig <- sig[order(sig$FDR), ]
  write.csv(sig,
            paste0(output_prefix, "_TE_superfamily_enrichment_hypo_sig.csv"),
            row.names=FALSE)
  message("Hypo enrichment completed for ", output_prefix)
}


# 8. Construct design matrix for Sex and Stage (and their interaction)
# --- redefine groups as the four combinations ---
samples_df$Group <- factor(
  paste0(samples_df$Stage, ".", samples_df$Sex),
  levels = c("Adult.F","Adult.M","Third_instar.F","Third_instar.M")
)

# quick check
table(samples_df$Group)

designSL <- model.matrix(~0+Group,data=samples_df)
designSL
# design expands designSL to model both methylated and unmethylated counts for each sample
design <- modelMatrixMeth(designSL)
design

fit <- glmQLFit(te_yall, design)
# In the factorial model, the coefficient names are:
colnames(design)

# Make the contrasts
my.contrasts <- makeContrasts(
                              AdultM_vs_adultF = GroupAdult.M - GroupAdult.F,
                              ThirdF_vs_AdultF = GroupThird_instar.F - GroupAdult.F,
                              ThirdM_vs_adultM = GroupThird_instar.M - GroupAdult.M,
                              ThirdM_vs_ThirdF = GroupThird_instar.M - GroupThird_instar.F,
                              M_vs_F = 0.5*(GroupAdult.M + GroupThird_instar.M) - 0.5*(GroupAdult.F + GroupThird_instar.F),
                              Third_vs_Adult = 0.5*(GroupThird_instar.M + GroupThird_instar.F) - 0.5*(GroupAdult.M + GroupAdult.F),
                              levels=design
                              )




# 9. Test for Adult Male - Adult Female
qlf_AdultM_vs_adultF <- glmQLFTest(fit, contrast=my.contrasts[,"AdultM_vs_adultF"])
topTags(qlf_AdultM_vs_adultF)
sum_AdultM_vs_adultF <- summary(decideTests(qlf_AdultM_vs_adultF))
capture.output(
  list(
    header  = "=== Adult Male vs Female at FDR 5% ===",
    summary = sum_AdultM_vs_adultF
  ),
  file = "TE_DMA_AdultM_vs_AdultF_summary.txt"
)

pdf("TE_DMA_AdultM_vs_AdultF_MDplot.pdf", width=6, height=6)
plotMD(qlf_AdultM_vs_adultF, main = "Adult M vs Adult F")
dev.off()
write.csv(topTags(qlf_AdultM_vs_adultF, n=Inf)$table, "TE_DMA_AdultM_vs_AdultF.csv", row.names=FALSE)
#################################################################
# Heatmap
contrast_table <- read.csv("TE_DMA_AdultM_vs_AdultF.csv", stringsAsFactors = FALSE)

sample_order <- c(paste0("DNA", sprintf("%02d", 1:6), "-Me"),   # Adult Female
                  paste0("DNA", sprintf("%02d", 7:12), "-Me"))  # Adult Male

sample_labels <- c(paste0("Adult Female ", sprintf("%02d", 1:6)),
                   paste0("Adult Male ", sprintf("%02d", 7:12)))

plot_TE_heatmap(contrast_table, te_yall, "TE_top100_heatmap_AdultM_vs_AdultF.pdf", sample_order, sample_labels)

# Volcano plot
df1 <- read.csv("TE_DMA_AdultM_vs_AdultF.csv", stringsAsFactors = FALSE)

pdf("Volcano_AdultM_vs_AdultF.pdf", width=8, height=8)
EnhancedVolcano(df1,
                lab = rep("", nrow(df1)),      
                x = 'logFC',
                y = 'FDR',
                title = "Adult Male vs Adult Female",
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
# Enrichment analysis
run_TE_superfamily_enrichment("TE_DMA_AdultM_vs_AdultF.csv", "AdultM_vs_AdultF")
run_TE_superfamily_enrichment_hyper("TE_DMA_AdultM_vs_AdultF.csv", "AdultM_vs_AdultF")
run_TE_superfamily_enrichment_hypo("TE_DMA_AdultM_vs_AdultF.csv", "AdultM_vs_AdultF")
#-------------------------------------------------------------------------------------------------------------
# Test for 3rd Female - Adult Female
qlf_3rdF_vs_AdultF <- glmQLFTest(fit, contrast=my.contrasts[,"ThirdF_vs_AdultF"])
topTags(qlf_3rdF_vs_AdultF)
sum_3rdF_vs_AdultF <- summary(decideTests(qlf_3rdF_vs_AdultF))
capture.output(
  list(
    header  = "=== 3rd-instar F vs Adult F at FDR 5% ===",
    summary = sum_3rdF_vs_AdultF
  ),
  file = "TE_DMA_3rdF_vs_AdultF_summary.txt"
)

pdf("TE_DMA_3rdF_vs_AdultF_MDplot.pdf", width=6, height=6)
plotMD(qlf_3rdF_vs_AdultF, main = "3rd instar F vs Adult F")
dev.off()
write.csv(topTags(qlf_3rdF_vs_AdultF, n=Inf)$table, "TE_DMA_3rdF_vs_AdultF.csv", row.names=FALSE)
#################################################################
# Heatmap
contrast_table <- read.csv("TE_DMA_3rdF_vs_AdultF.csv", stringsAsFactors = FALSE)

sample_order <- c(paste0("DNA", sprintf("%02d", 1:6), "-Me"),   # Adult Female
                  paste0("DNA", sprintf("%02d", 13:18), "-Me")) # 3rd Instar Female

sample_labels <- c(paste0("Adult Female ", sprintf("%02d", 1:6)),
                   paste0("3rd Instar Female ", sprintf("%02d", 13:18)))

plot_TE_heatmap(contrast_table, te_yall, "TE_top100_heatmap_3rdF_vs_AdultF.pdf", sample_order, sample_labels)

# Volcano plot
df2 <- read.csv("TE_DMA_3rdF_vs_AdultF.csv", stringsAsFactors = FALSE)

pdf("Volcano_3rdF_vs_AdultF.pdf", width=8, height=8)
EnhancedVolcano(df2,
                lab = rep("", nrow(df2)),    
                x = 'logFC',
                y = 'FDR',
                title = "3rd Instar Female vs Adult Female",
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

# Enrichment analysis
run_TE_superfamily_enrichment("TE_DMA_3rdF_vs_AdultF.csv", "3rdF_vs_AdultF")
run_TE_superfamily_enrichment_hyper("TE_DMA_3rdF_vs_AdultF.csv", "3rdF_vs_AdultF")
run_TE_superfamily_enrichment_hypo("TE_DMA_3rdF_vs_AdultF.csv", "3rdF_vs_AdultF")
#------------------------------------------------------------------------------------------------------
##################################################################################################
# Contrasts: 3rd Male - Adult Male
qlf_3rdM_vs_adultM <- glmQLFTest(fit, contrast=my.contrasts[,"ThirdM_vs_adultM"])
topTags(qlf_3rdM_vs_adultM)
sum_3rdM_vs_adultM <- summary(decideTests(qlf_3rdM_vs_adultM))
capture.output(
  list(
    header  = "=== 3rd-instar M vs Adult M at FDR 5% ===",
    summary = sum_3rdM_vs_adultM
  ),
  file = "TE_DMA_3rdM_vs_AdultM_summary.txt"
)

pdf("TE_DMA_3rdM_vs_AdultM_MDplot.pdf", width=6, height=6)
plotMD(qlf_3rdM_vs_adultM, main = "3rd instar M vs Adult M")
dev.off()
write.csv(topTags(qlf_3rdM_vs_adultM, n=Inf)$table, "TE_DMA_3rdM_vs_AdultM.csv", row.names=FALSE)
#################################################################
# Heatmap
contrast_table <- read.csv("TE_DMA_3rdM_vs_AdultM.csv", stringsAsFactors = FALSE)

sample_order <- c(paste0("DNA", sprintf("%02d", 7:12), "-Me"),   # Adult Male
                  paste0("DNA", sprintf("%02d", 19:24), "-Me"))  # 3rd Instar Male

sample_labels <- c(paste0("Adult Male ", sprintf("%02d", 7:12)),
                   paste0("3rd Instar Male ", sprintf("%02d", 19:24)))

plot_TE_heatmap(contrast_table, te_yall, "TE_top100_heatmap_3rdM_vs_AdultM.pdf", sample_order, sample_labels)

# Volcano plot
df3 <- read.csv("TE_DMA_3rdM_vs_AdultM.csv", stringsAsFactors = FALSE)

pdf("Volcano_3rdM_vs_AdultM.pdf", width=8, height=8)
EnhancedVolcano(df3,
                lab = rep("", nrow(df3)),   
                x = 'logFC',
                y = 'FDR',
                title = "3rd Instar Male vs Adult Male",
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

# Enrichment analysis
run_TE_superfamily_enrichment("TE_DMA_3rdM_vs_AdultM.csv", "3rdM_vs_AdultM")
run_TE_superfamily_enrichment_hyper("TE_DMA_3rdM_vs_AdultM.csv", "3rdM_vs_AdultM")
run_TE_superfamily_enrichment_hypo("TE_DMA_3rdM_vs_AdultM.csv", "3rdM_vs_AdultM")
#-------------------------------------------------------------------------------------------------------------
########################################################################################################
# Contrasts: 3rd Male - 3rd Female
qlf_3rdM_vs_3rdF <- glmQLFTest(fit, contrast=my.contrasts[,"ThirdM_vs_ThirdF"])
topTags(qlf_3rdM_vs_3rdF)
sum_3rdM_vs_3rdF <- summary(decideTests(qlf_3rdM_vs_3rdF))
capture.output(
  list(
    header  = "=== 3rd-instar M vs 3rd-instar F at FDR 5% ===",
    summary = sum_3rdM_vs_3rdF
  ),
  file = "TE_DMA_3rdM_vs_3rdF_summary.txt"
)

pdf("TE_DMA_3rdM_vs_3rdF_MDplot.pdf", width=6, height=6)
plotMD(qlf_3rdM_vs_3rdF, main = "3rd instar M vs 3rd instar F")
dev.off()

write.csv(topTags(qlf_3rdM_vs_3rdF, n=Inf)$table, "TE_DMA_3rdM_vs_3rdF.csv", row.names=FALSE)
#################################################################
# Heatmap
contrast_table <- read.csv("TE_DMA_3rdM_vs_3rdF.csv", stringsAsFactors = FALSE)

sample_order <- c(paste0("DNA", sprintf("%02d", 13:18), "-Me"),   # 3rd Instar Female
                  paste0("DNA", sprintf("%02d", 19:24), "-Me"))   # 3rd Instar Male

sample_labels <- c(paste0("3rd Instar Female ", sprintf("%02d", 13:18)),
                   paste0("3rd Instar Male ", sprintf("%02d", 19:24)))

plot_TE_heatmap(contrast_table, te_yall, "TE_top100_heatmap_3rdM_vs_3rdF.pdf", sample_order, sample_labels)

# Volcano plot
df4 <- read.csv("TE_DMA_3rdM_vs_3rdF.csv", stringsAsFactors = FALSE)

pdf("Volcano_3rdM_vs_3rdF.pdf", width=8, height=8)
EnhancedVolcano(df4,
                lab = rep("", nrow(df4)),   
                x = 'logFC',
                y = 'FDR',
                title = "3rd Instar Male vs 3rd Instar Female",
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

# Enrichment analysis
run_TE_superfamily_enrichment("TE_DMA_3rdM_vs_3rdF.csv", "3rdM_vs_3rdF")
run_TE_superfamily_enrichment_hyper("TE_DMA_3rdM_vs_3rdF.csv", "3rdM_vs_3rdF")
run_TE_superfamily_enrichment_hypo("TE_DMA_3rdM_vs_3rdF.csv", "3rdM_vs_3rdF")

#--------------------------------------------------------------------------------------------------
################################################################################################
# Contrasts: Average Male - Average Female
qlf_M_vs_F <- glmQLFTest(fit, contrast=my.contrasts[,"M_vs_F"])
topTags(qlf_M_vs_F)
sum_M_vs_F <- summary(decideTests(qlf_M_vs_F))
capture.output(
  list(
    header  = "=== Average Male vs Average Female at FDR 5% ===",
    summary = sum_M_vs_F
  ),
  file = "TE_DMA_M_vs_F_summary.txt"
)

pdf("TE_DMA_M_vs_F_MDplot.pdf", width=6, height=6)
plotMD(qlf_M_vs_F, main = "Average M vs Average F")
dev.off()

write.csv(topTags(qlf_M_vs_F, n=Inf)$table, "TE_DMA_M_vs_F.csv", row.names=FALSE)
#################################################################
# Heatmap
contrast_table <- read.csv("TE_DMA_M_vs_F.csv", stringsAsFactors = FALSE)

sample_order <- c(paste0("DNA", sprintf("%02d", 1:6), "-Me"),    # Adult Female
                  paste0("DNA", sprintf("%02d", 7:12), "-Me"),   # Adult Male
                  paste0("DNA", sprintf("%02d", 13:18), "-Me"),  # 3rd Instar Female
                  paste0("DNA", sprintf("%02d", 19:24), "-Me"))  # 3rd Instar Male

sample_labels <- c(paste0("Adult Female ", sprintf("%02d", 1:6)),
                   paste0("Adult Male ", sprintf("%02d", 7:12)),
                   paste0("3rd Instar Female ", sprintf("%02d", 13:18)),
                   paste0("3rd Instar Male ", sprintf("%02d", 19:24)))

plot_TE_heatmap(contrast_table, te_yall, "TE_top100_heatmap_M_vs_F.pdf", sample_order, sample_labels)
# Volcano plot
df5 <- read.csv("TE_DMA_M_vs_F.csv", stringsAsFactors = FALSE)

pdf("Volcano_M_vs_F.pdf", width=8, height=8)
EnhancedVolcano(df5,
                lab = rep("", nrow(df5)),  
                x = 'logFC',
                y = 'FDR',
                title = "Average Male vs Average Female",
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
# Enrichment analysis
run_TE_superfamily_enrichment("TE_DMA_M_vs_F.csv", "M_vs_F")
run_TE_superfamily_enrichment_hyper("TE_DMA_M_vs_F.csv", "M_vs_F")
run_TE_superfamily_enrichment_hypo("TE_DMA_M_vs_F.csv", "M_vs_F")
################################################################################################
# Contrasts: Average Third instar - Average Adult
qlf_Third_vs_Adult <- glmQLFTest(fit, contrast=my.contrasts[,"Third_vs_Adult"])
topTags(qlf_Third_vs_Adult)
sum_Third_vs_Adult <- summary(decideTests(qlf_Third_vs_Adult))
capture.output(
  list(
    header  = "=== Average 3rd instar vs Average Adult at FDR 5% ===",
    summary = sum_Third_vs_Adult
  ),
  file = "TE_DMA_Third_vs_Adult_summary.txt"
)

pdf("TE_DMA_Third_vs_Adult_MDplot.pdf", width=6, height=6)
plotMD(qlf_Third_vs_Adult, main = "Average 3rd instar vs Average Adult")
dev.off()

write.csv(topTags(qlf_Third_vs_Adult, n=Inf)$table, "TE_DMA_Third_vs_Adult.csv", row.names=FALSE)

#################################################################
# Heatmap
contrast_table <- read.csv("TE_DMA_Third_vs_Adult.csv", stringsAsFactors = FALSE)

sample_order <- c(paste0("DNA", sprintf("%02d", 1:6), "-Me"),    # Adult Female
                  paste0("DNA", sprintf("%02d", 7:12), "-Me"),   # Adult Male
                  paste0("DNA", sprintf("%02d", 13:18), "-Me"),  # 3rd Instar Female
                  paste0("DNA", sprintf("%02d", 19:24), "-Me"))  # 3rd Instar Male

sample_labels <- c(paste0("Adult Female ", sprintf("%02d", 1:6)),
                   paste0("Adult Male ", sprintf("%02d", 7:12)),
                   paste0("3rd Instar Female ", sprintf("%02d", 13:18)),
                   paste0("3rd Instar Male ", sprintf("%02d", 19:24)))

plot_TE_heatmap(contrast_table, te_yall, "TE_top100_heatmap_3rd_vs_Adult.pdf", sample_order, sample_labels)

# Volcano plot
df6 <- read.csv("TE_DMA_Third_vs_Adult.csv", stringsAsFactors = FALSE)

pdf("Volcano_3rd_vs_Adult.pdf", width=8, height=8)
EnhancedVolcano(df6,
                lab = rep("", nrow(df6)),   
                x = 'logFC',
                y = 'FDR',
                title = "Average 3rd Instar vs Average Adult",
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

# Enrichment analysis
run_TE_superfamily_enrichment("TE_DMA_Third_vs_Adult.csv", "Third_vs_Adult")
run_TE_superfamily_enrichment_hyper("TE_DMA_Third_vs_Adult.csv", "Third_vs_Adult")
run_TE_superfamily_enrichment_hypo("TE_DMA_Third_vs_Adult.csv", "Third_vs_Adult")
# 11. Output message -------------------------------------------------------------
message("All differential methylation analysis results have been saved as CSV files!")
















