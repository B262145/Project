# Purpose: Calculate a,b,c,d counts, Odds Ratio (OR), p-values, and BH-FDR 
#          for LTR/Pao superfamily across multiple contrasts

compute_LTRPao_abcd_with_FDR <- function(input_csv, contrast_name,
                                         modes = c("all","hyper","hypo"),
                                         fdr_cutoff = 0.05,
                                         alternative = "greater",
                                         target_sf = "LTR/Pao") {
  all_TE <- read.csv(input_csv, stringsAsFactors = FALSE)
  all_sf <- all_TE$type
  BG_total <- length(all_sf)                  # Total number of tested TE copies
  SF_total_BG <- sum(all_sf == target_sf)     # Total LTR/Pao copies in the background set

  one_mode <- function(mode) {
    # Select DE set based on mode
    if (mode == "all") {
      de_TE <- subset(all_TE, FDR < fdr_cutoff)
    } else if (mode == "hyper") {
      de_TE <- subset(all_TE, FDR < fdr_cutoff & logFC > 0)
    } else if (mode == "hypo") {
      de_TE <- subset(all_TE, FDR < fdr_cutoff & logFC < 0)
    } else stop("mode must be all/hyper/hypo")

    DE_total <- nrow(de_TE)
    de_sf    <- de_TE$type
    sflist   <- sort(unique(all_sf))

    # Compute 2x2 contingency tables + Fisher's exact test for all superfamilies
    by_sf <- lapply(sflist, function(sf) {
      sf_total <- sum(all_sf == sf)
      a <- sum(de_sf == sf)                   # DE & superfamily
      c <- sf_total - a                       # NotDE & superfamily
      b <- DE_total - a                       # DE & not-superfamily
      d <- (BG_total - sf_total) - b          # NotDE & not-superfamily

      mat <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE,
                    dimnames = list(c("SF","NonSF"), c("DE","NotDE")))
      ft  <- fisher.test(mat, alternative = alternative)

      data.frame(
        Superfamily       = sf,
        Mode              = mode,
        a_DE_and_SF       = a,
        b_DE_and_nonSF    = b,
        c_NotDE_and_SF    = c,
        d_NotDE_and_nonSF = d,
        DE_total          = DE_total,
        BG_total          = BG_total,
        SF_total_in_BG    = sf_total,
        Prop_in_DE        = ifelse(DE_total > 0, a / DE_total, NA_real_),
        Prop_in_BG        = sf_total / BG_total,
        OR                = unname(ft$estimate),
        OR_CI_low         = ft$conf.int[1],
        OR_CI_high        = ft$conf.int[2],
        p_value           = ft$p.value,
        stringsAsFactors  = FALSE
      )
    })

    df_all <- do.call(rbind, by_sf)
    # Adjust p-values within this mode across all superfamilies
    df_all$FDR <- p.adjust(df_all$p_value, method = "BH")
    # Keep only the LTR/Pao row (or NA if absent)
    out <- subset(df_all, Superfamily == target_sf)
    if (nrow(out) == 0) {
      out <- data.frame(
        Superfamily       = target_sf, Mode = mode,
        a_DE_and_SF       = NA, b_DE_and_nonSF = NA,
        c_NotDE_and_SF    = NA, d_NotDE_and_nonSF = NA,
        DE_total          = DE_total, BG_total = BG_total,
        SF_total_in_BG    = SF_total_BG,
        Prop_in_DE        = NA, Prop_in_BG = SF_total_BG / BG_total,
        OR = NA, OR_CI_low = NA, OR_CI_high = NA,
        p_value = NA, FDR = NA, stringsAsFactors = FALSE
      )
    }
    out
  }

  res <- do.call(rbind, lapply(modes, one_mode))
  res$Contrast <- contrast_name
  # Reorder columns
  res <- res[, c("Contrast","Mode","Superfamily",
                 "a_DE_and_SF","b_DE_and_nonSF","c_NotDE_and_SF","d_NotDE_and_nonSF",
                 "DE_total","BG_total","SF_total_in_BG",
                 "Prop_in_DE","Prop_in_BG","OR","OR_CI_low","OR_CI_high","p_value","FDR")]
  rownames(res) <- NULL
  res
}

# === Run for six contrasts and save results ===
LTRPao_abcd_FDR_all <- do.call(rbind, list(
  compute_LTRPao_abcd_with_FDR("TE_DMA_3rdF_vs_AdultF.csv",  "3rdF_vs_AdultF"),
  compute_LTRPao_abcd_with_FDR("TE_DMA_Third_vs_Adult.csv",  "Third_vs_Adult"),
  compute_LTRPao_abcd_with_FDR("TE_DMA_AdultM_vs_AdultF.csv","AdultM_vs_AdultF"),
  compute_LTRPao_abcd_with_FDR("TE_DMA_3rdM_vs_AdultM.csv",  "3rdM_vs_AdultM"),
  compute_LTRPao_abcd_with_FDR("TE_DMA_3rdM_vs_3rdF.csv",    "3rdM_vs_3rdF"),
  compute_LTRPao_abcd_with_FDR("TE_DMA_M_vs_F.csv",          "M_vs_F")
))

write.csv(LTRPao_abcd_FDR_all, "LTRPao_abcd_OR_FDR.csv", row.names = FALSE)
