# This R script reads weighted methylation data and computes summary statistics per genomic feature and sample group.
# It then generates a single-row faceted bar plot showing group-wise weighted methylation levels with confidence intervals for each feature,
# including group-specific genome-wide weighted methylation as reference lines.


library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)

# Set directory path
dir <- "My_split"
input_file  <- file.path(dir, "weighted_meth_per_genomic_region_all_samples.txt")

# Read combined methylation data
message("Reading combined methylation file...")
#all <- read_delim(input_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
all <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


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

# Compute summary statistics (mean, SD, SE, CI) for weightedMeth grouped by genomic feature and sample group
summary_grouped <- summarySE(
  all,
  measurevar = "weightedMeth",
  groupvars  = c("feature", "group")
)

# Output the summary statistics
write.csv(summary_grouped, file = file.path(dir, "summary_grouped_stats.csv"), row.names = FALSE)

# Define the genomic features
new_levels <- c("promoter", "exon", "intron", "intergenic",
                "promoter_TE", "exonic_TE", "intronic_TE", "intergenic_TE")

# Reassign feature as a factor with the new ordering
summary_grouped <- summary_grouped %>%
  mutate(feature = factor(feature, levels = new_levels))


# Define genomewide methylation levels (rounded to three decimal places)
genomewide_meth <- data.frame(
  group = c("3rd Instar Female", "3rd Instar Male", "Adult Female", "Adult Male"),
  meth  = c(0.061, 0.097, 0.068, 0.094)
)

# Specify custom colors for each group
colors <- c(
  "3rd Instar Female" = "#f8766c",
  "3rd Instar Male"   = "#7bad00",
  "Adult Female"      = "#00bec4",
  "Adult Male"        = "#c67bff"
)

### ---- Compute region counts per genomic feature ---- ###
region_counts_grouped <- all %>%
  dplyr::group_by(feature, group) %>%
  dplyr::summarise(region_number = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(feature = factor(feature, levels = new_levels))

# Generate the plot
message("Generating grouped bar plot across genomic features...")

### ---- Upper panel: WML bar plot ---- ###
p1 <- ggplot(summary_grouped, aes(x = feature, y = weightedMeth, fill = group)) +
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
   
  # 4) Apply a clean theme and set axis labels
  theme_bw() +
  xlab("") +
  ylab("Weighted Methylation Level (WML)") +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1, size = 15),
    axis.text.y  = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text  = element_text(size = 15),
    legend.title = element_blank()
  ) +
  
  # 5) Use the custom fill and line colors for groups
  scale_fill_manual(
    values = colors
  ) +
  scale_color_manual(
    name   = "Genome-wide WML",
    values = colors,
    labels = c(
      "3rd Instar Female Genome-wide WML",
      "3rd Instar Male Genome-wide WML",
      "Adult Female Genome-wide WML",
      "Adult Male Genome-wide WML")
  ) +
  
  # 6) Ensure that bars and dotted lines appear in the same legend
  guides(
    fill  = guide_legend(override.aes = list(linetype = 0)),
    color = guide_legend(override.aes = list(linetype = "dotted"))
  )

### ---- Lower panel: Region count bar plot ---- ###
p2 <- ggplot(region_counts_grouped, aes(x = feature, y = region_number, fill = group)) + 
  geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
  geom_text(aes(label = region_number),
            position = position_dodge(width = 0.9),
            angle = 60, hjust = 0, vjust = -0.1,
            size = 4, fontface = "italic") +
  scale_fill_manual(values = colors) +
  #scale_y_log10() +  # log scale on y-axis
  coord_cartesian(clip = "off") +
  theme_bw() +
  xlab("Genomic Feature") +
  ylab("Region Count") +  #(not log10 scale) # label updated accordingly
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y  = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.position = "none"
  )

### ---- Combine the two plots ---- ###
final_plot <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(4, 2))

### ---- Save the final plot ---- ###
plot_file <- file.path(dir, "LATEST_WML_with_feature_counts_angle60.pdf")
ggsave(plot_file, plot = final_plot, width = 15, height = 12, limitsize = FALSE)

message("Final plot saved to: ", plot_file)
