# R script: Extract MD plot data and re-plot with ggplot2 for six contrasts and by superfamily (PDF output)

# 1. Load needed libraries
library(edgeR)
library(ggplot2)
library(dplyr)

# 2. Assume `te_yall` and all qlf objects (e.g. qlf_AdultM_vs_adultF, qlf_3rdF_vs_AdultF, etc.) are in your workspace
contrast_objs <- list(
  AdultM_vs_adultF = qlf_AdultM_vs_adultF,
  ThirdF_vs_AdultF = qlf_3rdF_vs_AdultF,
  ThirdM_vs_AdultM = qlf_3rdM_vs_adultM,
  ThirdM_vs_ThirdF = qlf_3rdM_vs_3rdF,
  M_vs_F           = qlf_M_vs_F,
  Third_vs_Adult   = qlf_Third_vs_Adult
)

# 3. Extract MD data for each contrast
data_list <- lapply(names(contrast_objs), function(cname) {
  df <- topTags(contrast_objs[[cname]], n=Inf)$table
  df <- df %>%
    mutate(
      Comparison = cname,
      status = case_when(
        FDR < 0.05 & logFC <  0 ~ "Down",
        FDR < 0.05 & logFC >  0 ~ "Up",
        TRUE                    ~ "NotSig"
      ),
      logCPM  = logCPM,
      logFC   = logFC,
      Superfamily = te_yall$genes$type
    )
  df
})
md_all <- bind_rows(data_list)

# Set Comparison as factor with desired plotting order for facets
desired_order <- c(
  "AdultM_vs_adultF",
  "ThirdM_vs_ThirdF",
  "M_vs_F",
  "ThirdM_vs_AdultM",
  "ThirdF_vs_AdultF",
  "Third_vs_Adult"
)
md_all$Comparison <- factor(md_all$Comparison, levels = desired_order)

# Custom PDF save wrapper to enforce consistent background
save_pdf <- function(plot_obj, filename, width, height) {
  ggsave(
    filename = filename,
    plot     = plot_obj,
    width    = width,
    height   = height,
    bg       = "white"
  )
}

# Function to plot ordered layers: NotSig first, then Down, then Up
# Adds horizontal dashed line at y=0 and explicit white backgrounds
gg_md_plot <- function(df, title, filename, width=6, height=5) {
  p <- ggplot() +
    geom_point(data = filter(df, status == "NotSig"),
               aes(x = logCPM, y = logFC), color = "grey50", alpha = 0.6) +
    geom_point(data = filter(df, status == "Down"),
               aes(x = logCPM, y = logFC), color = "blue", alpha = 0.6) +
    geom_point(data = filter(df, status == "Up"),
               aes(x = logCPM, y = logFC), color = "red", alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = title,
         x = "Average logCPM", y = "log2 Fold Change") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.border     = element_rect(fill = NA, color = "black"),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  save_pdf(p, filename, width, height)
}

# 4. Generate MD plots for six contrasts in specified order saved as PDF
custom_order <- c(
  "AdultM_vs_adultF",
  "ThirdM_vs_ThirdF",
  "M_vs_F",
  "ThirdM_vs_AdultM",
  "ThirdF_vs_AdultF",
  "Third_vs_Adult"
)
for(cname in custom_order) {
  df_sub <- filter(md_all, Comparison == cname)
  if(nrow(df_sub) == 0) next
  fname  <- paste0("MDplot_", cname, ".pdf")
  ttl    <- paste0("MD Plot: ", cname)
  gg_md_plot(df_sub, ttl, fname)
}

# 5. List of 29 superfamilies
superfams <- c(
  "DNA", "DNA/Dada", "DNA/hAT-Ac", "DNA/hAT-Tag1", "DNA/hAT-Tip100",
  "DNA/Maverick", "DNA/Merlin", "DNA/MULE-MuDR", "DNA/MULE-NOF",
  "DNA/PIF-Harbinger", "DNA/PiggyBac", "DNA/Sola2", "DNA/TcMar-Mariner",
  "DNA/TcMar-Pogo", "DNA/TcMar-Tc2", "LINE/I", "LINE/I-Jockey", "LINE/L2",
  "LINE/Penelope", "LINE/R1", "LINE/RTE-BovB", "LINE/RTE-RTE", "LTR",
  "LTR/Copia", "LTR/Gypsy", "SINE/5S", "SINE/MIR", "SINE/tRNA-RTE", "LTR/Pao"
)

# 6. MD plots per superfamily with facets, saved as PDF
for(sf in superfams) {
  df_sf <- filter(md_all, Superfamily == sf)
  if(nrow(df_sf) == 0) next
  p_sf <- ggplot() +
    geom_point(data = filter(df_sf, status == "NotSig"),
               aes(x = logCPM, y = logFC), color = "grey50", alpha = 0.6) +
    geom_point(data = filter(df_sf, status == "Down"),
               aes(x = logCPM, y = logFC), color = "blue", alpha = 0.6) +
    geom_point(data = filter(df_sf, status == "Up"),
               aes(x = logCPM, y = logFC), color = "red", alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ Comparison, ncol = 3) +
    labs(title = paste0("MD Plot: Superfamily = ", sf),
         x = "Average logCPM", y = "log2 Fold Change") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.border     = element_rect(fill = NA, color = "black"),
      strip.text       = element_text(face = "bold"),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  fname <- gsub("[\\/ ]", "_", sf)
  save_pdf(p_sf, paste0("MDplot_superfam_", fname, ".pdf"), 12, 8)
}

# 7. Optional: Save superfamily plots as PNG too (transparent background removed)
for(sf in superfams) {
  df_sf <- filter(md_all, Superfamily == sf)
  if(nrow(df_sf) == 0) next
  p_sf <- ggplot() +
    geom_point(data = filter(df_sf, status == "NotSig"),
               aes(x = logCPM, y = logFC), color = "grey50", alpha = 0.6) +
    geom_point(data = filter(df_sf, status == "Down"),
               aes(x = logCPM, y = logFC), color = "blue", alpha = 0.6) +
    geom_point(data = filter(df_sf, status == "Up"),
               aes(x = logCPM, y = logFC), color = "red", alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ Comparison, ncol = 3) +
    labs(title = paste0("MD Plot: Superfamily = ", sf),
         x = "Average logCPM", y = "log2 Fold Change") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.border     = element_rect(fill = NA, color = "black"),
      strip.text       = element_text(face = "bold"),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  fname <- gsub("[\\/ ]", "_", sf)
  ggsave(paste0("MDplot_superfam_", fname, ".png"), p_sf, width = 12, height = 8, bg = "white")
}
