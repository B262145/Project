# This R script reads weighted methylation data for all samples,
# classifies features by TE superfamily, and computes group-level summary statistics. 
# It then fits simple linear models comparing male vs female methylation levels for both adult and 3rd instar stages. 
# Finally, the script outputs diagnostic plots, an XY scatter plot with regression lines, 
# and density plots for each specified TE superfamily.



library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tune)

# Set directory path
output_dir <- "WML_TE"
input_file  <- file.path(output_dir, "weighted_meth_per_genomic_region_all_samples.txt")

# Read combined methylation data
message("Reading combined methylation file...")
all <- read_delim(input_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Extract TE class from feature
# - If the feature contains '/', take the substring before '/'
# - If the feature exactly matches a known class, use it directly
# - Otherwise, assign "Other"
known_classes <- c("DNA", "LINE", "LTR", "RC", "SINE")
all <- all %>%
  mutate(class = case_when(
    str_detect(feature, "/") ~ str_extract(feature, "^[^/]+"),
    feature %in% known_classes  ~ feature,
    TRUE                        ~ "Other"
  ))

# Reorder class so that "Other" appears last
all$class <- factor(
  all$class,
  levels = c(setdiff(sort(unique(all$class)), "Other"), "Other")
)

# Save the mapping from feature to inferred class
write_csv(
  distinct(all, feature, class),
  file.path(output_dir, "TE_feature_class_map.csv")
)

# Define a function to compute summary statistics (mean, SD, SE, CI) by group
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                      conf.interval = .95, .drop = TRUE) {
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x)) else length(x)
  }
  
  datac <- plyr::ddply(
    data, groupvars, .drop = .drop,
    .fun = function(xx, col) {
      c(
        N    = length2(xx[[col]], na.rm = na.rm),
        mean = mean(xx[[col]], na.rm = na.rm),
        sd   = sd(xx[[col]], na.rm = na.rm)
      )
    },
    measurevar
  )
  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)
  ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

# Assign each sample to a group based on its name prefix
all <- all %>%
  mutate(group = case_when(
    str_detect(sample, "^adultfemale") ~ "Adult Female",
    str_detect(sample, "^adultmale")   ~ "Adult Male",
    str_detect(sample, "^3rd_female")  ~ "3rd Instar Female",
    str_detect(sample, "^3rd_male")    ~ "3rd Instar Male",
    TRUE                               ~ "Unknown"
  ))

# Compute summary statistics (mean, SD, SE, CI) for weightedMeth grouped by class, feature, and group
summary_grouped <- summarySE(
  all,
  measurevar = "weightedMeth",
  groupvars  = c("class", "feature", "group")
)

# Determine ordering of features: first all non-"Other", then all "Other"
all_feats <- summary_grouped %>% select(class, feature) %>% distinct()
non_other <- all_feats %>% filter(class != "Other") %>% pull(feature)
other     <- all_feats %>% filter(class == "Other") %>% pull(feature)
new_levels <- c(non_other, other)

# Reassign feature as a factor with the new ordering
summary_grouped <- summary_grouped %>%
  mutate(feature = factor(feature, levels = new_levels))

# -----------------------------------------------------------------------------
# STEP: Build an XY dataset comparing Male vs Female per stage
# -----------------------------------------------------------------------------

# 1) Extract adult female vs adult male
adult_xy <- summary_grouped %>%
  filter(group %in% c("Adult Female", "Adult Male")) %>%
  select(class, feature, group, weightedMeth) %>%
  pivot_wider(
    names_from  = group,
    values_from = weightedMeth
  ) %>%
  rename(
    female_WML = `Adult Female`,
    male_WML   = `Adult Male`
  ) %>%
  mutate(stage = "Adult")

# 2) Extract 3rd instar female vs 3rd instar male
juvenile_xy <- summary_grouped %>%
  filter(group %in% c("3rd Instar Female", "3rd Instar Male")) %>%
  select(class, feature, group, weightedMeth) %>%
  pivot_wider(
    names_from  = group,
    values_from = weightedMeth
  ) %>%
  rename(
    female_WML = `3rd Instar Female`,
    male_WML   = `3rd Instar Male`
  ) %>%
  mutate(stage = "3rd Instar")

# 3) Combine adult and juvenile for plotting
xy_df <- bind_rows(adult_xy, juvenile_xy)

# STEP: Fit models and save diagnostic plots to PDF
# -----------------------------------------------------------------------------
adult_mod    <- lm(male_WML ~ female_WML, data = filter(xy_df, stage == "Adult"))
juvenile_mod <- lm(male_WML ~ female_WML, data = filter(xy_df, stage == "3rd Instar"))

# Adult model diagnostics with main title
pdf(file = file.path(output_dir, "model_diagnostics_adult.pdf"), width = 8, height = 8)
par(mfrow = c(2, 2),       # arrange 4 plots in 2×2
    oma   = c(0, 0, 3, 0))  # add space at top for an overall title
plot(adult_mod, which = c(1, 2, 3, 5))
# Draw an overall title in the outer margin
mtext("Adult Model Diagnostics", side = 3, outer = TRUE, line = 1, cex = 1.5)
dev.off()

# 3rd Instar model diagnostics with main title
pdf(file = file.path(output_dir, "model_diagnostics_3rdInstar.pdf"), width = 8, height = 8)
par(mfrow = c(2, 2),
    oma   = c(0, 0, 3, 0))
plot(juvenile_mod, which = c(1, 2, 3, 5))
mtext("3rd Instar Model Diagnostics", side = 3, outer = TRUE, line = 1, cex = 1.5)
dev.off()

# Restore single-plot layout
par(mfrow = c(1, 1))

# Extract coefficients for drawing with geom_abline()
coef_adult <- coef(adult_mod)
coef_juv   <- coef(juvenile_mod)

# build the equation strings
eq_adult <- paste0("Adult:   y = ",
                   round(coef_adult[2], 3), " x + ",
                   round(coef_adult[1], 3))
eq_juv   <- paste0("3rd Instar:   y = ",
                   round(coef_juv[2], 3), " x + ",
                   round(coef_juv[1], 3))

# pick positions for the text (e.g. top-right corner of each facet)
# adjust x and y to suit the data range
xpos <- max(xy_df$female_WML, na.rm = TRUE) * 0.6
ypos_adult <- max(xy_df$male_WML[xy_df$stage=="Adult"], na.rm=TRUE) * 0.9
ypos_juv  <- max(xy_df$male_WML[xy_df$stage=="3rd Instar"], na.rm=TRUE) * 0.9



# -----------------------------------------------------------------------------
# STEP: Plot Male vs Female WML with regression lines via geom_abline()
# -----------------------------------------------------------------------------
message("Generating Male vs Female XY plot with regression lines...")

p_xy <- ggplot(xy_df, aes(x = female_WML, y = male_WML, color = class, shape = stage)) +
  geom_point(size = 5, alpha = 0.8) +
  
  # Add a diagonal reference line (y = x)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  # Adult regression line from pre-fitted model
  geom_abline(intercept = coef_adult[1], slope = coef_adult[2],
              color = "darkred", linetype = "solid", linewidth = 1) +

  # 3rd Instar regression line from pre-fitted model
  geom_abline(intercept = coef_juv[1], slope = coef_juv[2],
              color = "darkblue", linetype = "solid", linewidth = 1) +
  coord_obs_pred(ratio = 1) +    # Assign the same length to x and y axis  
  theme_bw() +
  xlab("Female Weighted Methylation Level (WML)") +
  ylab("Male Weighted Methylation Level (WML)") +
  ggtitle("Male vs Female TE WML, by Developmental Stage and TE Order") +
  theme(
    legend.title      = element_text(size = 20),
    legend.text       = element_text(size = 15),
    axis.text         = element_text(size = 15),
    axis.title        = element_text(size = 20),
    plot.title        = element_text(size = 20, hjust = 0),
    strip.text        = element_text(size = 16)
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17)) +
  guides(
    color = guide_legend(title = "TE Order"),
    shape = guide_legend(title = "Developmental Stage")
  ) +
  # annotate adult equation (only in the Adult facet)
  annotate("text", 
           x    = xpos, 
           y    = ypos_adult, 
           label= eq_adult,
           size = 5,
           color= "darkred",
           hjust= 0) +
  # annotate juvenile equation (only in the 3rd Instar facet)
  annotate("text", 
           x    = xpos, 
           y    = ypos_juv, 
           label= eq_juv,
           size = 5,
           color= "darkblue",
           hjust= 0)

# Save the XY plot
plot_xy_file <- file.path(output_dir, "TE_WML_Male_vs_Female_XY_withLM.pdf")
ggsave(plot_xy_file, plot = p_xy, width = 10, height = 8)

message("XY plot with regression lines saved to: ", plot_xy_file)


# ------------------------------------------------------------
# Function to plot per–copy Female vs Male WML density, faceted by stage
# ------------------------------------------------------------
library(ggpointdensity)
library(viridis)

plot_TE_copy_density <- function(feature_name) {
  # 1) subset to that TE superfamily and annotate stage + sex
  df_feat <- all %>%
    filter(feature == feature_name) %>%
    mutate(
      stage = case_when(
        str_detect(sample, "^adult")   ~ "Adult",
        str_detect(sample, "^3rd")     ~ "3rd Instar",
        TRUE                           ~ NA_character_
      ),
      sex = case_when(
        str_detect(sample, "female")   ~ "Female",
        str_detect(sample, "male")     ~ "Male",
        TRUE                           ~ NA_character_
      )
    ) %>%
    drop_na(stage, sex)
  
  # 2) compute mean WML per copy (chr, start, end) × stage × sex
  df_mean <- df_feat %>%
    group_by(chr, start, end, stage, sex) %>%
    summarise(mean_WML = mean(weightedMeth, na.rm = TRUE), .groups = "drop") %>%
    mutate(length = end - start + 1)  # Calculate the TE length
     
  # 3) pivot into Female vs Male columns, keep stage
  df_xy <- df_mean %>%
    pivot_wider(
      id_cols     = c(chr, start, end, stage, length),
      names_from  = sex,
      values_from = mean_WML
    ) %>%
    rename(Female = Female, Male = Male)
    
  # 3.5) save filtered-out cases (i.e. only present in one sex)
  df_filtered_out <- df_xy %>%
    filter(is.na(Female) | is.na(Male))
  
  if(nrow(df_filtered_out) > 0) {
    out_filtered_file <- file.path(
      output_dir,
      paste0("filtered_out_", gsub("[^A-Za-z0-9]", "_", feature_name), ".csv")
    )
    write_csv(df_filtered_out, out_filtered_file)
    message("Filtered out entries saved to:\n  ", out_filtered_file)
  }
  
  # remove rows with NA in either Female or Male before plotting
  df_xy <- df_xy %>% drop_na(Female, Male)
            
  # 3.6) calculate the mean
  mean_df <- df_xy %>%
    group_by(stage) %>%
    summarise(
      meanF = mean(Female, na.rm = TRUE),
      meanM = mean(Male,   na.rm = TRUE),
      label = paste0(
        "Mean:\n(",
        round(meanF, 3), ", ",
        round(meanM, 3), ")"
      ),
      .groups = "drop"
    )
    
  # Before plotting, sort by length from smallest to largest
  df_xy <- df_xy %>% arrange(length)

  # 4) build and save a faceted scatter plot colored by length
  p <- ggplot(df_xy, aes(x = Female, y = Male, color = length)) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_viridis_c(name = "TE Length (bp)", option = "viridis") +
    geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "black") +
    facet_wrap(~ stage, ncol = 2) +
    coord_obs_pred(ratio = 1) +    # Assign the same length to x and y axis 
    theme_bw() +
    labs(
      title = paste0("Per-copy WML density for ", feature_name),
      x     = "Female Weighted Methylation Level (WML)",
      y     = "Male Weighted Methylation Level (WML)"
    ) +
    theme(
      panel.spacing = unit(2, "lines"),
      plot.title   = element_text(size = 20, hjust = 0.5),
      strip.text   = element_text(size = 15),
      axis.text    = element_text(size = 15),
      axis.title   = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.text  = element_text(size = 12)
    ) +
    
    # 5) Add the mean and label for each facet
    geom_point(
      data        = mean_df,
      aes(x = meanF, y = meanM),
      inherit.aes = FALSE,
      shape       = 16,
      size        = 2,
      color       = "red"
    ) +
    
    geom_text(
      data        = mean_df,
      aes(x = meanF, y = meanM, label = "Mean"),
      inherit.aes = FALSE,
      color       = "red",
      vjust       = 0.5,    
      hjust       = 0,  
      nudge_x     = 0.05,
      size        = 5,
      fontface    = "plain"
    )

  # 6) Save to PDF
  out_pdf <- file.path(
    output_dir,
    paste0("density_", gsub("[^A-Za-z0-9]", "_", feature_name), "_faceted.pdf")
  )
  ggsave(out_pdf, p, width = 10, height = 5)
  message("Saved density plot for ", feature_name, " to:\n  ", out_pdf)
  
  return(invisible(p))
}

# ------------------------------------------------------------
# Now call it for the enriched TE superfamilies:
# ------------------------------------------------------------
for(feat in c(
  "DNA", "DNA/Dada", "DNA/hAT-Ac", "DNA/hAT-Tag1", "DNA/hAT-Tip100", 
  "DNA/Maverick", "DNA/Merlin", "DNA/MULE-MuDR", "DNA/MULE-NOF", "DNA/PIF-Harbinger", 
  "DNA/PiggyBac", "DNA/Sola2", "DNA/TcMar-Mariner", "DNA/TcMar-Pogo", "DNA/TcMar-Tc2", 
  "LINE/I", "LINE/I-Jockey", "LINE/L2", "LINE/Penelope", "LINE/R1", "LINE/RTE-BovB", 
  "LINE/RTE-RTE", "LTR", "LTR/Copia", "LTR/Gypsy", "SINE/5S", "SINE/MIR", "SINE/tRNA-RTE"
)) {
  plot_TE_copy_density(feat)
}
