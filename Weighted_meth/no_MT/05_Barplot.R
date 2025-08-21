# This R script reads weighted methylation data and computes summary statistics per TE class, TE superfamily, and sample group.
# It then generates a single-row faceted bar plot showing group-wise weighted methylation levels with confidence intervals for each TE superfamily,
# including group-specific average genome-wide weighted methylation as reference lines.


library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

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
  # Helper to compute length while optionally removing NAs
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x)) else length(x)
  }
  
  # Use ddply to calculate N, mean, and SD for each combination of grouping variables
  datac <- ddply(
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
  
  # Rename the 'mean' column to the actual measurement variable name
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  # Calculate standard error
  datac$se <- datac$sd / sqrt(datac$N)
  
  # Compute t-multiplier for the confidence interval
  ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
  
  # Calculate the margin of error
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

print(colnames(summary_grouped))


# Determine ordering of features: first all non-"Other", then all "Other"
all_feats <- summary_grouped %>%
  select(class, feature) %>%
  distinct()

non_other <- all_feats %>%
  filter(class != "Other") %>%
  pull(feature)
other <- all_feats %>%
  filter(class == "Other") %>%
  pull(feature)

new_levels <- c(non_other, other)

# Reassign feature as a factor with the new ordering
summary_grouped <- summary_grouped %>%
  mutate(feature = factor(feature, levels = new_levels))

# Define genomewide methylation levels (rounded to three decimal places)
genomewide_meth <- data.frame(
  group = c("3rd Instar Female", "3rd Instar Male", "Adult Female", "Adult Male"),
  meth  = c(0.061, 0.097, 0.068, 0.094)
)

# -----------------------------------------------------------------------------
# NEW: Identify TE superfamilies whose mean WML > genome‐wide WML in ALL groups
# 1) Extract feature–group mean
feat_means <- summary_grouped %>%
  select(feature, group, weightedMeth)

# 2) Join genomewide and filter
feat_comp <- feat_means %>%
  left_join(genomewide_meth, by = "group") %>%
  filter(weightedMeth > meth)

# 3) Keep only features that pass in all 4 groups
good_feats <- feat_comp %>%
  count(feature) %>%
  filter(n == 4) %>%
  pull(feature)

# 4) Build output table: one row per feature, columns = group mean WML
out_table <- feat_means %>%
  filter(feature %in% good_feats) %>%
  pivot_wider(
    names_from  = group,
    values_from = weightedMeth
  ) %>%
  arrange(feature)

# 5) Write to CSV and print
out_path <- file.path(output_dir, "TE_superfamilies_above_genomewide_WML.csv")
write.csv(out_table, file = out_path, row.names = FALSE, quote = FALSE)
message("TE superfamilies above genome‐wide WML saved to:\n", out_path)
print(out_table)
# -----------------------------------------------------------------------------

# Specify custom colors for each group
colors <- c(
  "3rd Instar Female" = "#f8766c",
  "3rd Instar Male"   = "#7bad00",
  "Adult Female"      = "#00bec4",
  "Adult Male"        = "#c67bff"
)

# Generate the plot using a single row of facets, one facet per class
message("Generating single-row faceted plot by TE class...")

p <- ggplot(summary_grouped, aes(x = feature, y = weightedMeth, fill = group)) +
  # 1) Draw bar chart with dodge positioning for different groups
  geom_bar(position = position_dodge(), stat = "identity") +
  
  # 2) Add error bars for the confidence intervals
  geom_errorbar(
    aes(ymin = weightedMeth - ci, ymax = weightedMeth + ci),
    width = 0.2, position = position_dodge(0.9)
  ) +
  
  # 3) Add genomewide methylation level as dotted horizontal lines
  geom_hline(
    data = genomewide_meth,
    aes(yintercept = meth, color = group),
    linetype = "dotted", linewidth = 1
  ) +
  
  # 4) Facet by TE class in a single row
  facet_grid(
    cols   = vars(class),
    scales = "free_x",    # allow each facet to have its own x-axis scale
    space  = "free_x",    # allow each facet’s width to vary based on number of features
    labeller = labeller(class = c("DNA" = "DNA transposon"))
  ) +
  
  # 5) Apply a clean theme and set axis labels
  theme_bw() +
  xlab("TE") +
  ylab("Weighted Methylation Level (WML)") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y  = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text   = element_text(size = 16),  # top facet strip text size
    legend.text  = element_text(size = 15),
    legend.title = element_blank()
  ) +
  
  # 6) Use the custom fill and line colors for groups
  scale_fill_manual(
    values = colors
  ) +
  scale_color_manual(
    name   = "Average Genome-wide WML",
    values = colors,
    labels = c(
      "3rd Instar Female Genome-wide WML",
      "3rd Instar Male Genome-wide WML",
      "Adult Female Genome-wide WML",
      "Adult Male Genome-wide WML")
  ) +
  
  # 7) Ensure that bars and dotted lines appear in the same legend
  guides(
    fill  = guide_legend(override.aes = list(linetype = 0)),
    color = guide_legend(override.aes = list(linetype = "dotted"))
  )

# Save the plot to a PDF file
plot_file <- file.path(output_dir, "TE_methylation.pdf")
ggsave(plot_file, plot = p, width = 30, height = 10, limitsize = FALSE)

message("Faceted plot saved to: ", plot_file)
