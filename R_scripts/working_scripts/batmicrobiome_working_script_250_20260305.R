# Ternary Plot for Kruger Wild Mammal Data
#install.packages("ggtern")
#install.packages("rlang")
library(ggtern)
library(dplyr)
library(readr)
library(RColorBrewer)
library(tidyverse)
library(caret)
library(randomForest)
library(viridis)
library(ggpubr)
library(patchwork)
library(cowplot)

# data preparation
setwd("C:/Users/mmayy/OneDrive - University of Florida/Projects/Bat Microbiome/2_Analysis/variable_importance/250")
#setwd("C:/Users/thongthum.t/OneDrive - University of Florida/Projects/Bat Microbiome/2_Analysis/variable_importance/250")
taxonomy_df <- read.delim("data/taxonomy.tsv", row.names = 1, check.names = FALSE)
otu_table_original <- read.delim("data/otu_table_corrected_bioproj_primer.txt", row.names = 1, check.names = FALSE)
#otu_table_original <- read.delim("data/otu_table_corrected_primer_bioproj.txt", row.names = 1, check.names = FALSE)
otu_table <- otu_table_original

## ------------------------- summary of statistics ----------------------------

library(tidyverse)

# Function to calculate OTU table statistics
calculate_otu_stats <- function(otu_table, table_name = "OTU_Table") {
  
  # Calculate sample sums (columns)
  sample_sums <- colSums(otu_table)
  
  # Calculate feature sums (rows)
  feature_sums <- rowSums(otu_table)
  
  # Create statistics dataframe
  stats_df <- data.frame(
    Metric = c(
      "Number of samples",
      "Number of features",
      "Total frequency",
      "Median frequency per sample",
      "3rd quartile frequency per sample",
      "Mean frequency per sample",
      "Median frequency per feature",
      "3rd quartile frequency per feature",
      "Mean frequency per feature"
    ),
    Value = c(
      ncol(otu_table),                           # Number of samples
      nrow(otu_table),                           # Number of features
      sum(otu_table),                            # Total frequency
      median(sample_sums),                       # Median per sample
      quantile(sample_sums, 0.75),              # 3rd quartile per sample
      mean(sample_sums),                         # Mean per sample
      median(feature_sums),                      # Median per feature
      quantile(feature_sums, 0.75),             # 3rd quartile per feature
      mean(feature_sums)                         # Mean per feature
    )
  )
  
  # Add table name column
  stats_df$Table <- table_name
  
  return(stats_df)
}

stats <- calculate_otu_stats(otu_table, "OTU_Table")

write.csv(stats, "output/otu_table_statistics.csv", row.names = FALSE)

# Print summary
print(stats)

## ----------------------------------------------------------------------------

# Create aggregated genus-level table for most analyses
taxonomic <- taxonomy_df[, 1]
otu_table$bac_tax_group <- taxonomic[match(rownames(otu_table), rownames(taxonomy_df))]
otu_tax <- aggregate(. ~ bac_tax_group, data = otu_table, FUN = sum)
rownames(otu_tax) <- otu_tax$bac_tax_group
otu_tax$bac_tax_group <- NULL
otu_tax <- otu_tax[!grepl("^d__Archaea;", rownames(otu_tax)), ]
otu_tax <- t(otu_tax)
otu_tax <- as.data.frame(otu_tax)

################################################################################

# --------------- 1. LOAD AND MERGE
metadata <- read.delim("data/metadata_250_20260216.txt", row.names = 1, check.names = FALSE)

if ("reads" %in% colnames(metadata)) {
  metadata <- metadata[metadata$reads != "single", ]
  message(paste("Filtered metadata for 'single' reads. Remaining rows:", nrow(metadata)))
} else {
  warning("'reads' column not found in metadata.")
}

df <- merge(otu_tax, metadata, by = "row.names")

# Fix row names after merge
rownames(df) <- df$Row.names
df$Row.names <- NULL

# --------------- 2. APPLY QUALITY FILTERS (Duplicates & Read Depth)

# --- A. REMOVE DUPLICATE SAMPLES ---
# Extract only the bacteria abundance columns to check for duplicates
abundance_cols <- grep("^d__", colnames(df), value = TRUE)
abundance_matrix <- df[, abundance_cols]

# Check for duplicates based on abundance profile
is_dup <- duplicated(abundance_matrix)

if (sum(is_dup) > 0) {
  print(paste("WARNING: Removed", sum(is_dup), "samples with identical bacterial counts."))
  # Keep only unique rows
  df <- df[!is_dup, ]
} else {
  print("No duplicate samples found.")
}

# --- B. REMOVE LOW READ COUNT SAMPLES ---
# Recalculate abundance matrix on the filtered dataframe
abundance_matrix <- df[, abundance_cols]

# Calculate total reads per sample
sample_sums <- rowSums(abundance_matrix)

# Define Threshold (e.g., 5000 reads)
min_reads <- 5000
keep_samples <- sample_sums >= min_reads

# Print summary of what will be removed
print("Read Depth Summary:")
print(summary(sample_sums))
print(paste("Samples with <", min_reads, "reads removed:", sum(!keep_samples)))
print(paste("Samples remaining:", sum(keep_samples)))

# Apply Filter
df <- df[keep_samples, ]

# Check if we lost everything
if (nrow(df) == 0) stop("Error: All samples were filtered out! Check your data or lower 'min_reads'.")

# ------------------ 3. METADATA CLEANING & FORMATTING

# Clean Diet Labels
df <- df %>%
  mutate(diet_1 = trimws(diet_1)) %>%
  mutate(diet_1 = case_when(
    diet_1 %in% c("Fruit", "Fruits.", "Fruits") ~ "Fruit",
    diet_1 %in% c("Fish", "Fishes.", "Fishes") ~ "Fish",
    diet_1 %in% c("Insect", "Insects", "Insectivore") ~ "Insect",
    diet_1 %in% c("Flower", "Flowers", "Nectar", "Pollen", "Flower. Nectar") ~ "Fruit",
    TRUE ~ diet_1
  )) %>%
  mutate(cw_diet1 = paste(captive_wild, diet_1, sep = "_")) %>%
  mutate(cw_diet1_sampletype = paste(captive_wild, diet_1, sample_type, sep = "_"))

# Rename specific taxonomy columns if they exist
if("genus" %in% colnames(df)) colnames(df)[colnames(df) == "genus"] <- "taxonomic_genus"
if("species" %in% colnames(df)) colnames(df)[colnames(df) == "species"] <- "taxonomic_species"

df$Diet_type <- df$diet_1

# Quick sanity check
print(table(df$Diet_type))

# ---------------- 4. PREPARE OUTPUT (METADATA + TAXA)

df_metadata <- df %>%
  rownames_to_column("sample_id") 

df_metadata <- df_metadata %>%
  rename(
    index = sample_id,
    Common.name_Species = common_name,
    Diet_type = Diet_type)

print(paste("Final rows:", nrow(df_metadata)))

write.csv(df_metadata, "output/df_taxa_metadata_250.csv", row.names = FALSE)

################################################################################
## ----------------------- rarefaction curve -------------------------------- ##

library(vegan)

abundance_cols <- grep("^d__", colnames(df), value = TRUE)
otu_matrix <- df[, abundance_cols]
otu_matrix <- as.matrix(otu_matrix)
class(otu_matrix) <- "numeric"

# Define your custom colors (matching your previous script)
diet_colors <- c(
  "Insect" = "#ff5e2b",
  "Blood"  = "#ba0970",
  "Fruit"  = "#2E8B57",
  "Fish"   = "#303F9F"
)

rare_data <- rarecurve(otu_matrix, step = 100, label = FALSE, tidy = TRUE) 

rare_data$Site <- as.character(rare_data$Site)

# Create a temporary lookup table from your metadata
meta_lookup <- df %>% 
  rownames_to_column("Site") %>% 
  select(Site, Diet_type)

# Merge rarefaction data with Diet info
rare_plot_data <- left_join(rare_data, meta_lookup, by = "Site")

# Create the plot object (assign it to a variable, don't just print it)
alpha_rarefaction_plot <- ggplot(rare_plot_data, aes(x = Sample, y = Species, group = Site, color = Diet_type)) +
  geom_line(alpha = 0.7, size = 0.5) + 
  scale_color_manual(values = diet_colors) +
  scale_x_continuous(
    limits = c(0, 10000),       # Set max to 10000
    breaks = seq(0, 10000, by = 1000) # Set interval to 1000
  ) +
  theme_bw() +
  labs(
    title = "Alpha Rarefaction Curves by Diet (250bp)",
    x = "Sequencing Depth (Reads)",
    y = "Observed Features (ASVs)",
    color = "Diet Group"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(alpha_rarefaction_plot)

out_dir <- "output_bioproj_prim"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("alpha_rarefaction_by_diet.", ext)),
    plot = alpha_rarefaction_plot,
    width = 12,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  print(paste("Saved:", ext))
}

################################################################################
### --------------------- alpha and beta diversity ------------------------- ###

#install.packages(c("vegan", "ggplot2", "dplyr", "tidyr", "ggpubr"))
#install.packages("picante")
# Load necessary libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(readr)
library(picante) 
library(ape) 

# ------------------ 0. DATA LOADING AND PREPARATION ---------------------------

#setwd("C:/Users/thongthum.t/OneDrive - University of Florida/Projects/Bat Microbiome/2_Analysis/variable_importance/250")
setwd("C:/Users/mmayy/OneDrive - University of Florida/Projects/Bat Microbiome/2_Analysis/variable_importance/250")
df_input <- read_csv("output/df_taxa_metadata_250.csv")

# We need to separate the community matrix (abundance data) from the metadata
# Based on your previous script, taxa columns start with "d__Bacteria"
# We extract columns that contain taxa abundance
community_matrix <- df_input %>%
  select(starts_with("d__Bacteria"))

# Extract metadata (everything that isn't a taxa column)
metadata_sub <- df_input %>%
  select(-starts_with("d__Bacteria"))

# Ensure the row names match for vegan
rownames(community_matrix) <- df_input$index # Assuming 'index' is the sample_id column created earlier
rownames(metadata_sub) <- df_input$index

# Optional: Rarefaction (Normalizing sequencing depth)
set.seed(1000)
min_depth <- min(rowSums(community_matrix))
community_matrix <- rrarefy(community_matrix, min_depth)

table(metadata_sub$Diet_type)

diet_colors <- c(
  "Insect" = "#ff5e2b",
  "Blood"  = "#ba0970",
  "Fruit"  = "#2E8B57",
  "Fish"   = "#303F9F")


table(metadata_sub$cw_diet1)

cw_diet_colors <- c(
  "Wild_Insect" = "#ff5e2b", 
  "Captive_Insect" = "#e37f05", 
  "Wild_probably_Insect" = "#fa8246", 
  "Captive_probably_Insect" = "#fcb644",
  "Wild_Fruit" = "#2E8B57", 
  "Wild_probably_Fruit" = "#1b5e3a", 
  "Captive_probably_Fruit" = "#80e35f",
  "Wild_Blood" = "#ba0970", 
  "Wild_probably_Blood" = "#e329da", 
  "Wild_probably_Fish" = "#303F9F"
)

#------------------------- 1. DATA PREPARATION ---------------------------------

# (Assuming df_input is already loaded or created as per previous steps)
df_input <- read_csv("output/df_taxa_metadata_250.csv")

# Ensure cw_diet1 exists
if(!"cw_diet1" %in% colnames(df_input)) {
  df_input <- df_input %>%
    mutate(cw_diet1 = paste(captive_wild, Diet_type, sep = "_"))
}

community_matrix <- df_input %>% select(starts_with("d__Bacteria"))
metadata_sub <- df_input %>% select(-starts_with("d__Bacteria"))

################################################################################
## ----------------- ALPHA, BETA, CONFUSION MATRIX: DIET1 ------------------- ##
################################################################################

library(tidyverse)
library(vegan)
library(picante)
library(rstatix)
library(ggpubr)
library(ape) 
library(rstatix)

set.seed(1000)
tree_file <- rtree(ncol(community_matrix), tip.label = colnames(community_matrix))

shannon  <- diversity(community_matrix, index = "shannon")
observed <- specnumber(community_matrix)
simpson <- diversity(community_matrix, index = "simpson")
#pd_result <- pd(community_matrix, tree_file, include.root = FALSE)

# Calculate Indices
alpha_div <- data.frame(
  Shannon = shannon,
  Simpson = simpson, 
  Observed = observed
)

# Merge back with metadata
alpha_plot_data <- cbind(alpha_div, metadata_sub)

# ============================================================================ #
# ------------------ NOT INCLUDED: ALHA DIVERSITY PANEL ---------------------- #
# ------------------------------ Predictor: DIET1 ---------------------------- #
# ============================================================================ #

library(tidyverse)
library(rstatix)
library(ggpubr)

# Note: Ensure you have defined 'diet_colors' before running this function.
# Example: diet_colors <- c("red", "blue", "green", "orange") 

plot_alpha_diet <- function(data, y_metric, title) {
  
  # 1. PREPARE DATA & CALCULATE DUNN'S TEST
  # Create a clean dataframe for calculation
  stats_data <- data %>%
    select(Diet_type, all_of(y_metric)) %>%
    filter(!is.na(.data[[y_metric]]), !is.na(Diet_type))
  
  # Run Dunn's Test (No positions yet)
  dunn_stats <- stats_data %>%
    dunn_test(as.formula(paste(y_metric, "~ Diet_type")), p.adjust.method = "fdr") %>%
    add_significance()
  
  # 2. MANUALLY CALCULATE Y-POSITIONS FOR NEAT STACKING
  # This guarantees consistent spacing for ALL brackets (Sig + NS)
  if(nrow(dunn_stats) > 0) {
    # Find the maximum value in the plot to know where to start drawing
    max_y <- max(stats_data[[y_metric]], na.rm = TRUE)
    
    # Define a step size (e.g., 10% of the data range)
    y_range <- max_y - min(stats_data[[y_metric]], na.rm = TRUE)
    step_size <- y_range * 0.05
    
    # Assign y.position incrementally for every comparison in the list
    # The first bracket goes at max_y + step, the second at max_y + 2*step, etc.
    dunn_stats$y.position <- max_y + (1:nrow(dunn_stats)) * step_size
  }
  
  # 3. Calculate global p-value (Kruskal-Wallis)
  kruskal_p <- kruskal_test(stats_data, as.formula(paste(y_metric, "~ Diet_type")))$p
  subtitle_text <- paste0("Kruskal-Wallis, p = ", sprintf("%.3f", kruskal_p))
  
  # 4. Plot Base
  p <- ggplot(data, aes(x = Diet_type, y = .data[[y_metric]], fill = Diet_type)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.4, color = "black") +
    theme_bw() +
    labs(title = title, subtitle = subtitle_text, x = "Diet Group", y = y_metric) +
    scale_fill_manual(values = diet_colors) +
    
    # --- LAYER A: NON-SIGNIFICANT (RED) ---
    stat_pvalue_manual(
      dunn_stats %>% filter(p.adj > 0.05),
      label = "p.adj.signif",
      color = "red",
      tip.length = 0.005,
      vjust = 0.5,
      remove.bracket = FALSE
    ) +
    
    # --- LAYER B: SIGNIFICANT (BLACK) ---
    stat_pvalue_manual(
      dunn_stats %>% filter(p.adj <= 0.05),
      label = "p.adj.signif",
      color = "black",
      tip.length = 0.005,
      vjust = 0.5,
      remove.bracket = FALSE
    ) +
    
    # Theme adjustments
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10)
    )
  
  return(p)
}

p1 <- plot_alpha_diet(alpha_plot_data, "Shannon", "Shannon Entropy Index")
p2 <- plot_alpha_diet(alpha_plot_data, "Simpson", "Simpson's Diversity Index")
p3 <- plot_alpha_diet(alpha_plot_data, "Observed", "Observed Richness")

# Combine and Save
alpha_plots_diet1 <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
print(alpha_plots_diet1)

out_dir <- "output_bioproj_prim"
dir.create("output_bioproj_prim", showWarnings = FALSE, recursive = TRUE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("alpha_diet1_sig_nonsig.", ext)),
    plot = alpha_plots_diet1,
    width = 13,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================ #
# ------------------ 1.1) ALHA DIVERSITY PANEL (SIG ONLY) -------------------- #
# ------------------------------ Predictor: DIET1 ---------------------------- #
# ============================================================================ #

plot_alpha_sig_only <- function(data, y_metric, title) {
  
  # 1. CALCULATE DUNN'S TEST
  # We calculate the statistics first (without positions)
  stat.test <- data %>%
    dunn_test(as.formula(paste(y_metric, "~ Diet_type")), p.adjust.method = "fdr") %>%
    add_significance()
  
  # 2. FILTER FOR SIGNIFICANT PAIRS
  sig_stats <- stat.test %>% filter(p.adj <= 0.05)
  
  # 3. MANUALLY CALCULATE Y-POSITIONS FOR NEAT STACKING
  # This replaces add_xy_position to ensure tight, gap-free stacking
  if(nrow(sig_stats) > 0) {
    
    # Find the maximum value in the plot to know where to start drawing
    max_y <- max(data[[y_metric]], na.rm = TRUE)
    
    # Define a step size (e.g., 10% of the data range)
    y_range <- max_y - min(data[[y_metric]], na.rm = TRUE)
    step_size <- y_range * 0.05
    
    # Assign y.position incrementally
    # The first significant bracket goes at max_y + step
    # The second goes at max_y + 2*step, etc.
    sig_stats$y.position <- max_y + (1:nrow(sig_stats)) * step_size
  }
  
  # 4. GENERATE PLOT
  p <- ggplot(data, aes(x = Diet_type, y = .data[[y_metric]], fill = Diet_type)) +
    
    # Boxplot
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    
    # Add points
    geom_jitter(width = 0.2, size = 1.0, alpha = 0.6, color = "black") +
    
    theme_bw() +
    labs(title = title, x = "Diet Group", y = y_metric) +
    scale_fill_manual(values = diet_colors) +
    
    # 5. ADD BRACKETS
    # We pass the manually calculated sig_stats
    stat_pvalue_manual(
      sig_stats,
      label = "p.adj.signif",    
      tip.length = 0.01,
      hide.ns = TRUE,
      vjust = 0.5
    ) +
    
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10)
    )
  
  return(p)
}

# Generate Plots
p1_sig <- plot_alpha_sig_only(alpha_plot_data, "Shannon", "Shannon Diversity")
p2_sig <- plot_alpha_sig_only(alpha_plot_data, "Simpson", "Simpson's Diversity Index")
p3_sig <- plot_alpha_sig_only(alpha_plot_data, "Observed", "Observed Richness")

# Combine
alpha_plots_diet1_sig <- ggarrange(p1_sig, p2_sig, p3_sig, ncol = 3, nrow = 1)
print(alpha_plots_diet1_sig)

out_dir <- "output_bioproj_prim"
dir.create("output_bioproj_prim", showWarnings = FALSE, recursive = TRUE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("alpha_plots_diet1_sig.", ext)),
    plot = alpha_plots_diet1_sig,
    width = 14, 
    height = 10, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================ #
# ----------------------------- STATISTICS ----------------------------------- #
# ============================================================================ #

#install.packages("BiocManager")
#remotes::install_github("joey711/phyloseq", force = TRUE)

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(phyloseq)
library(ape)

# A. OTU Table (Taxa must be columns, Samples must be rows based on your previous code)
otu_mat <- as.matrix(community_matrix)

# B. Taxonomy Table (Optional for distance calc, but good practice)
# For 'taxonomy_df', use: tax_mat <- as.matrix(taxonomy_df)
tax_mat <- matrix(colnames(otu_mat), ncol = 1)
rownames(tax_mat) <- colnames(otu_mat)
colnames(tax_mat) <- "Species"

# C. Metadata
meta_df <- as.data.frame(metadata_sub)
rownames(meta_df) <- rownames(metadata_sub)

# D. Phylogenetic Tree
# IF USE A TREE FILE, UNCOMMENT AND USE THIS LINE:
# tree_file <- read_tree("path/to/your/tree.nwk")

set.seed(1000)
tree_file <- rtree(ncol(otu_mat), tip.label = colnames(otu_mat))

# Create the Phyloseq Object
physeq <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = FALSE),
  tax_table(tax_mat),
  sample_data(meta_df),
  phy_tree(tree_file)
)

# --------------------------- 1) Kruskal-Wallis ------------------------------ #

library(rstatix)
library(dplyr)
library(tibble)

kruskal_shannon <- alpha_plot_data %>% kruskal_test(Shannon ~ Diet_type)
kruskal_simpson <- alpha_plot_data %>% kruskal_test(Simpson ~ Diet_type)
kruskal_observed <- alpha_plot_data %>% kruskal_test(Observed ~ Diet_type)

# Combine the global test results into a single data frame
all_kruskal_results <- bind_rows(kruskal_shannon, kruskal_simpson, kruskal_observed)

# Define the output directory and filename for Kruskal-Wallis results
out_dir <- "output_bioproj_prim/statistics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
kruskal_output_file <- file.path(out_dir, "alpha_diet1_kruskal_results.csv")

# Save the global test results to CSV
write.csv(all_kruskal_results, file = kruskal_output_file, row.names = FALSE)
message("Global Kruskal-Wallis test results saved to: ", kruskal_output_file)


# --- 2. Perform Pairwise Wilcoxon Tests ---
# This replicates the post-hoc part of your original script.

# Use wilcox_test() which is the rstatix equivalent of pairwise.wilcox.test
# It returns a tidy data frame. We specify the same 'BH' adjustment method.
pairwise_shannon <- alpha_plot_data %>%
  wilcox_test(Shannon ~ Diet_type, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(metric = "Shannon")

pairwise_simpson <- alpha_plot_data %>%
  wilcox_test(Simpson ~ Diet_type, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(metric = "Simpson")

pairwise_observed <- alpha_plot_data %>%
  wilcox_test(Observed ~ Diet_type, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(metric = "Observed")

# Combine all pairwise results into a single data frame
all_pairwise_results <- bind_rows(pairwise_shannon, pairwise_simpson, pairwise_observed)

# Define the output filename for the pairwise results
pairwise_output_file <- file.path(out_dir, "alpha_diet1_pairwise_wilcox_results.csv")

# Save the pairwise test results to CSV
write.csv(all_pairwise_results, file = pairwise_output_file, row.names = FALSE)
message("Pairwise Wilcoxon test results saved to: ", pairwise_output_file)

# -------------------------- 2) Dunn's Test ---------------------------------- #

library(phyloseq)
library(vegan)
library(picante)
library(rstatix)
library(dplyr)
library(tibble)

# --- Assumptions ---
# This script assumes you have a phyloseq object named 'physeq'
# with 'Diet_type' in its sample data and a phylogenetic tree.

# --- 1. Calculate Alpha Diversity Metrics ---
alpha_div_richness <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson"))
otu_table_for_pd <- t(as(otu_table(physeq), "matrix"))
phylo_tree_for_pd <- phy_tree(physeq)
alpha_div_faith <- pd(otu_table_for_pd, phylo_tree_for_pd, include.root = FALSE)

# --- 2. Combine All Data into a Single Data Frame ---
metadata <- data.frame(sample_data(physeq))
all_alpha_data <- metadata %>%
  rownames_to_column("SampleID") %>%
  left_join(rownames_to_column(alpha_div_richness, "SampleID"))

# --- 3. Perform Pairwise Tests Independently for Each Metric ---
metrics_to_test <- c("Observed", "Shannon", "Simpson")
all_results_list <- list() 

for (metric_name in metrics_to_test) {
  
  # Create a formula for the test, e.g., Observed ~ Diet_type
  formula <- as.formula(paste(metric_name, "~ Diet_type"))
  
  # Prepare data specifically for this metric, removing any NAs
  current_data <- all_alpha_data %>%
    select(Diet_type, all_of(metric_name)) %>%
    filter(!is.na(.data[[metric_name]]))
  
  # Check for sufficient group sizes for this specific metric
  group_counts <- current_data %>%
    count(Diet_type) %>%
    filter(n > 1)
  
  if (nrow(group_counts) >= 2) {
    # If we have at least two groups with more than one sample, run the test
    message(paste("Performing Dunn's test for:", metric_name))
    
    dunn_result <- dunn_test(data = current_data, formula = formula, p.adjust.method = "fdr") %>%
      mutate(metric = metric_name)
    
    all_results_list[[metric_name]] <- dunn_result # Add result to our list
  } else {
    # Otherwise, skip this metric and inform the user
    message(paste("Skipping Dunn's test for:", metric_name, "- insufficient data after removing NAs or only one group with n>1."))
  }
}

# --- 4. Combine Successful Results and Export ---
if (length(all_results_list) > 0) {
  # Combine all the data frames from the list into one
  final_dunn_results <- bind_rows(all_results_list)
  
  # Define output directory and filename
  out_dir <- "output_bioproj_prim/statistics"
  output_filename <- "alpha_diet1_dunn_test.csv"
  
  # Ensure the output directory exists before trying to save the file
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Construct the full path to the output file
  full_path <- file.path(out_dir, output_filename)
  
  # Write the results to the specified path
  write.csv(final_dunn_results, file = full_path, row.names = FALSE)
  
  # Notify the user with the full, correct path
  message("Successfully saved pairwise test results to: ", full_path)
  
} else {
  message("Analysis complete, but no metrics had sufficient data to perform pairwise tests.")
}

# ============================================================================ #
# ---------------- 1.2) BETA DIVERSITY (Visuals + PERMANOVA) ----------------- #
# ------------------------------ Predictor: DIET1 ---------------------------- #
# ============================================================================ #

# A. Calculate Distance Matrix (Bray-Curtis)
# Bray-Curtis is standard for abundance data
dist_bray <- vegdist(community_matrix, method = "bray")

# B. Ordination (PCoA/MDS)
pcoa_res <- cmdscale(dist_bray, eig = TRUE, k = 2) # k=2 dimensions

# Extract coordinates for plotting
pcoa_df <- data.frame(
  PC1 = pcoa_res$points[,1],
  PC2 = pcoa_res$points[,2]
)

diet_shapes <- c(
  "Blood"  = 16,
  "Insect" = 3,
  "Fruit"  = 15,
  "Fish"   = 17) 

# Calculate percent variance explained by axes
eig_percent <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)

# Merge with metadata for plotting
pcoa_plot_data <- cbind(pcoa_df, metadata_sub)

# C. Visualization (PCoA Scatterplot)
beta_plot <- ggplot(pcoa_plot_data, aes(x = PC1, y = PC2, color = Diet_type, shape = Diet_type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Diet_type), geom = "polygon", alpha = 0.1, level = 0.95) +
  theme_bw() +
  labs(
    title = "PCoA of Bray-Curtis Dissimilarity",
    x = paste0("PCoA 1 (", eig_percent[1], "%)"),
    y = paste0("PCoA 2 (", eig_percent[2], "%)"),
    color = "Diet Group",
    shape = "Diet Group",
    fill  = "Diet Group"
  ) +
  scale_color_manual(values = diet_colors) +
  scale_fill_manual(values = diet_colors) +
  scale_shape_manual(values = diet_shapes) +  # This applies your specific shapes
  theme(
    plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10)
  )

print(beta_plot)

out_dir <- "output_bioproj_prim"
dir.create("output_bioproj_prim", showWarnings = FALSE, recursive = TRUE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("beta_bray_diet1.", ext)),
    plot = beta_plot,
    width = 10, 
    height = 8, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}

## ------------- 1.3 BETA DIVERSITY PANEL (bray, jaccard, wunifrac, unifrac)

#install.packages("BiocManager")
#remotes::install_github("joey711/phyloseq", force = TRUE)

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(phyloseq)
library(ape)

# A. OTU Table (Taxa must be columns, Samples must be rows based on your previous code)
otu_mat <- as.matrix(community_matrix)

# B. Taxonomy Table (Optional for distance calc, but good practice)
# For 'taxonomy_df', use: tax_mat <- as.matrix(taxonomy_df)
tax_mat <- matrix(colnames(otu_mat), ncol = 1)
rownames(tax_mat) <- colnames(otu_mat)
colnames(tax_mat) <- "Species"

# C. Metadata
meta_df <- as.data.frame(metadata_sub)
rownames(meta_df) <- rownames(metadata_sub)

# D. Phylogenetic Tree
# IF USE A TREE FILE, UNCOMMENT AND USE THIS LINE:
# tree_file <- read_tree("path/to/your/tree.nwk")

set.seed(1000)
tree_file <- rtree(ncol(otu_mat), tip.label = colnames(otu_mat))

# Create the Phyloseq Object
physeq <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = FALSE),
  tax_table(tax_mat),
  sample_data(meta_df),
  phy_tree(tree_file)
)

# -------------- A. FUNCTION TO GENERATE PCoA PLOTS
plot_beta_metric <- function(physeq_obj, distance_method, title, weighted = FALSE) {
  
  # Calculate Distance
  # UniFrac requires special handling in the ordinate function
  if (distance_method == "wunifrac") {
    ord <- ordinate(physeq_obj, method = "PCoA", distance = "unifrac", weighted = TRUE)
  } else if (distance_method == "unifrac") {
    ord <- ordinate(physeq_obj, method = "PCoA", distance = "unifrac", weighted = FALSE)
  } else {
    ord <- ordinate(physeq_obj, method = "PCoA", distance = distance_method)
  }
  
  # Plot using phyloseq's built-in plotter (easier wrapper for ggplot)
  p <- plot_ordination(physeq_obj, ord, color = "Diet_type", shape = "Diet_type") +
    geom_point(size = 1.5, alpha = 0.8) +
    stat_ellipse(aes(fill = Diet_type), geom = "polygon", alpha = 0.1, level = 0.95) +
    theme_bw() +
    labs(title = title) +
    
    # Custom Colors and Shapes
    scale_color_manual(values = diet_colors) +
    scale_fill_manual(values = diet_colors) +
    scale_shape_manual(values = diet_shapes) + 
    
    # Theme settings
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      legend.position = "none"
    )
  
  return(p)
}

# ------------ B. GENERATE THE FOUR PLOTS

# 1. Bray-Curtis (Abundance based)
p_bray <- plot_beta_metric(physeq, "bray", "Bray-Curtis (Abundance)")

# 2. Jaccard (Presence/Absence based)
# binary = TRUE is handled automatically by phyloseq for "jaccard" usually, 
# but "jaccard" in vegan is inherently quantitative. 
# For true binary Jaccard in phyloseq, use distance="jaccard", binary=TRUE
# Simulated here by converting data to presence/absence first for this plot

physeq_binary <- transform_sample_counts(physeq, function(x) ifelse(x > 0, 1, 0))
#p_jaccard <- plot_beta_metric(physeq_binary, "jaccard", "Jaccard (Presence/Absence)")

# 3. Weighted UniFrac (Phylogeny + Abundance)
p_wunifrac <- plot_beta_metric(physeq, "wunifrac", "Weighted UniFrac (Phylogeny + Abundance)")

# 4. Unweighted UniFrac (Phylogeny + Presence/Absence)
p_unifrac <- plot_beta_metric(physeq, "unifrac", "Unweighted UniFrac (Phylogeny + Presence/Absence)")

# ------------ C. COMBINE AND SAVE
# Extract the legend from one of the plots to share it
p_bray <- p_bray + 
  scale_color_manual(name = "Diet Group", values = diet_colors) +
  scale_fill_manual(name = "Diet Group", values = diet_colors) +
  scale_shape_manual(name = "Diet Group", values = diet_shapes) +   
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10)
  )

# 2. Extract the legend again with the new title
legend_bray <- get_legend(p_bray + theme(legend.position = "right"))

beta_plots_diet1 <- ggarrange(
  p_bray, p_wunifrac, p_unifrac,
  ncol = 3, nrow = 1,
  common.legend = TRUE,
  legend = "right",
  legend.grob = legend_bray
)

print(beta_plots_diet1)

out_dir <- "output_bioproj_prim"
dir.create("output_bioproj_prim", showWarnings = FALSE, recursive = TRUE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("beta_plots_diet1.", ext)),
    plot = beta_plots_diet1,
    width = 16, 
    height = 6, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}

## ----------------------------- STATISTICS --------------------------------- ##

library(phyloseq)
library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(tibble)
library(broom) 

# --- 1. Pre-computation of Distances ---

metadata <- data.frame(sample_data(physeq))
dist_bray <- distance(physeq, method = "bray")
dist_wunifrac <- distance(physeq, method = "wunifrac")
dist_uunifrac <- distance(physeq, method = "unifrac")

# Define output directory
out_dir <- "output_bioproj_prim/statistics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. Homogeneity of Dispersion (PERMDISP) to CSV ---

permadisp_bray <- tidy(anova(betadisper(dist_bray, metadata$Diet_type))) %>% mutate(metric = "Bray-Curtis")
permadisp_wunifrac <- tidy(anova(betadisper(dist_wunifrac, metadata$Diet_type))) %>% mutate(metric = "Weighted UniFrac")
permadisp_uunifrac <- tidy(anova(betadisper(dist_uunifrac, metadata$Diet_type))) %>% mutate(metric = "Unweighted UniFrac")

all_permadisp_results <- bind_rows(permadisp_bray, permadisp_wunifrac, permadisp_uunifrac)
permadisp_output_file <- file.path(out_dir, "beta_diet1_permadisp.csv")
write.csv(all_permadisp_results, file = permadisp_output_file, row.names = FALSE)
message("PERMDISP results saved to: ", permadisp_output_file)

# --- 3. Global PERMANOVA (adonis2) to CSV ---

permanova_bray <- tidy(adonis2(dist_bray ~ Diet_type, data = metadata, permutations = 999)) %>% mutate(metric = "Bray-Curtis")
permanova_wunifrac <- tidy(adonis2(dist_wunifrac ~ Diet_type, data = metadata, permutations = 999)) %>% mutate(metric = "Weighted UniFrac")
permanova_uunifrac <- tidy(adonis2(dist_uunifrac ~ Diet_type, data = metadata, permutations = 999)) %>% mutate(metric = "Unweighted UniFrac")

all_permanova_results <- bind_rows(permanova_bray, permanova_wunifrac, permanova_uunifrac)
permanova_output_file <- file.path(out_dir, "beta_diet1_permanova_global.csv")
write.csv(all_permanova_results, file = permanova_output_file, row.names = FALSE)
message("Global PERMANOVA results saved to: ", permanova_output_file)

# --- 4. Pairwise PERMANOVA to TXT ---
metadata <- data.frame(sample_data(physeq))

# --- 1. Pre-computation of Distances ---
metadata <- data.frame(sample_data(physeq))
dist_bray <- distance(physeq, method = "bray")
physeq_binary <- transform_sample_counts(physeq, function(x) ifelse(x > 0, 1, 0))
dist_wunifrac <- distance(physeq, method = "wunifrac")
dist_uunifrac <- distance(physeq, method = "unifrac")

# --- 2. Run Statistical Analyses and Generate Report ---
output_file <- "output_bioproj_prim/statistics/beta_diet1_pairwise.txt"
sink(output_file)
# --- Header ---
cat("==========================================================\n")
cat("          BETA DIVERSITY STATISTICAL ANALYSIS\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Grouping Variable: Diet_type\\\n")
cat("==========================================================\n\n")
# --- A. BRAY-CURTIS ANALYSIS ---
cat("A. BRAY-CURTIS (Abundance)\n")
cat("==========================================================\n\n")
cat("1. Homogeneity of Dispersion (PERMDISP):\n")
print(anova(betadisper(dist_bray, metadata$Diet_type)))
cat("\n2. Overall PERMANOVA (adonis2):\n")
print(adonis2(dist_bray ~ Diet_type, data = metadata, permutations = 999))
cat("\n3. Pairwise PERMANOVA (p-values adjusted with FDR):\n")
pairwise_bray <- pairwise.adonis2(dist_bray ~ Diet_type, data = metadata, p.adjust.m = "fdr")
print(pairwise_bray)
cat("\n\n")

# --- C. WEIGHTED UNIFRAC ANALYSIS ---
cat("C. WEIGHTED UNIFRAC (Phylogeny + Abundance)\n")
cat("==========================================================\n\n")
cat("1. Homogeneity of Dispersion (PERMDISP):\n")
print(anova(betadisper(dist_wunifrac, metadata$Diet_type)))
cat("\n2. Overall PERMANOVA (adonis2):\n")
print(adonis2(dist_wunifrac ~ Diet_type, data = metadata, permutations = 999))
cat("\n3. Pairwise PERMANOVA (p-values adjusted with FDR):\n")
pairwise_wunifrac <- pairwise.adonis2(dist_wunifrac ~ Diet_type, data = metadata, p.adjust.m = "fdr")
print(pairwise_wunifrac)
cat("\n\n")

# --- D. UNWEIGHTED UNIFRAC ANALYSIS ---
cat("D. UNWEIGHTED UNIFRAC (Phylogeny + Presence/Absence)\n")
cat("==========================================================\n\n")
cat("1. Homogeneity of Dispersion (PERMDISP):\n")
print(anova(betadisper(dist_uunifrac, metadata$Diet_type)))
cat("\n2. Overall PERMANOVA (adonis2):\n")
print(adonis2(dist_uunifrac ~ Diet_type, data = metadata, permutations = 999))
cat("\n3. Pairwise PERMANOVA (p-values adjusted with FDR):\n")
pairwise_uunifrac <- pairwise.adonis2(dist_uunifrac ~ Diet_type, data = metadata, p.adjust.m = "fdr")
print(pairwise_uunifrac)
cat("\n\n")

# --- End of Report ---
cat("==========================================================\n")
cat("                      End of Report\n")
cat("==========================================================\n")
sink()
message("Complete statistical results with pairwise tests have been saved to: ", output_file)

# ------------------------ 1.3) CONFUSION MATRIX  ---------------------------- #
# ------------------------------ Predictor: DIET1 ---------------------------- #
# ============================================================================ #

# Random Forest Classification Analysis
library(tidyverse)
library(caret)
library(patchwork)

# Function to generate confusion matrix plot
make_confusion_plot <- function(data, remove_outliers = TRUE, plot_title = "") {
  if (remove_outliers) {
    outliers <- c("K058683-R", "K058706-R", "K058789-R", "K058689-R", "K058814-R", "K058785-R")
    data <- data[!(data$index %in% outliers), ]
  }
  
  # Prepare metadata
  metadata <- data.frame(
    SampleID = data$index,
    Diet_type = factor(make.names(data$Diet_type)),
    stringsAsFactors = FALSE
  )
  
  # Extract abundance data
  abundance_table <- data[, grep("d__Bacteria", colnames(data))]
  rownames(abundance_table) <- metadata$SampleID
  
  # Extract family names
  extract_family <- function(tax_string) {
    if (grepl("g__", tax_string)) {
      gsub(".*;g__([^;]+).*", "\\1", tax_string)
    } else if (grepl("f__", tax_string)) {
      gsub(".*;f__([^;]+).*", "\\1", tax_string)
    } else {
      NA
    }
  }
  
  family_names <- sapply(colnames(abundance_table), extract_family)
  abundance_table_family <- rowsum(t(abundance_table), group = family_names, na.rm = TRUE)
  abundance_table_family <- t(abundance_table_family)
  abundance_table_family <- abundance_table_family[, !is.na(colnames(abundance_table_family))]
  
  # CLR transformation
  clr_transform <- function(x) {
    x <- x + 1
    log(x / exp(mean(log(x))))
  }
  abundance_table_clr <- t(apply(abundance_table_family, 1, clr_transform))
  
  # Prepare data for random forest
  rf_data <- cbind(Diet_type = metadata$Diet_type, as.data.frame(abundance_table_clr))
  rf_data$Diet_type <- as.factor(rf_data$Diet_type)
  
  # Train random forest model
  set.seed(1000)
  rf_model <- train(
    Diet_type ~ ., data = rf_data,
    method = "rf",
    trControl = trainControl(method = "cv", number = 5, classProbs = TRUE, savePredictions = "final"),
    ntree = 500, importance = TRUE, metric = "Accuracy"
  )
  
  # Get predictions
  cv_predictions <- rf_model$pred
  cv_predictions$obs <- factor(cv_predictions$obs, levels = levels(rf_data$Diet_type))
  cv_predictions$pred <- factor(cv_predictions$pred, levels = levels(rf_data$Diet_type))
  
  # Create confusion matrix
  conf_matrix_cv <- confusionMatrix(cv_predictions$pred, cv_predictions$obs)
  cm_df <- as.data.frame(conf_matrix_cv$table)
  
  # --- Calculate Ratios ---
  cm_df <- cm_df %>%
    group_by(Reference) %>%
    mutate(
      Total = sum(Freq),
      Ratio = Freq / Total,
      Label = sprintf("%.2f", Ratio)
    ) %>%
    ungroup()
  # --- MODIFICATION END ---
  
  # Create plot
  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Ratio)) + # Fill mapped to Ratio
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = Label),
              color = ifelse(cm_df$Ratio > 0.5, "black", "white"), 
              size = 5, fontface = "bold") +
    # MODIFIED: Changed labels to format as 2 decimal places
    scale_fill_viridis_c(name = "Ratio", labels = function(x) sprintf("%.2f", x)) +
    labs(x = "Expected Diet", y = "Predicted Diet", title = plot_title) +
    theme_minimal(base_size = 14) +
    theme(
      plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

library(readr)
library(ggplot2)
library(patchwork)

# Load data
# Ensure this path matches your local environment
data <- read_csv("output/df_taxa_metadata_250.csv")

# Create plot
conf_mat_plot_diet1 <- make_confusion_plot(data, remove_outliers = FALSE, plot_title = "")

# Display plot
conf_mat_plot_diet1

# Save the plot in multiple formats
out_dir <- "output_bioproj_prim"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("confusion_matrix_diet1.", ext)),
    plot = conf_mat_plot_diet1,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
}

## ------------------- confusion matrix counts ---------------------------------

make_confusion_plot_counts <- function(data, remove_outliers = TRUE, plot_title = "") {
  
  if (remove_outliers) {
    outliers <- c("K058683-R", "K058706-R", "K058789-R", "K058689-R", "K058814-R", "K058785-R")
    data <- data[!(data$index %in% outliers), ]
  }
  
  # Prepare metadata
  metadata <- data.frame(
    SampleID = data$index,
    Diet_type = factor(make.names(data$Diet_type)),
    stringsAsFactors = FALSE
  )
  
  # Extract abundance data
  abundance_table <- data[, grep("d__Bacteria", colnames(data))]
  rownames(abundance_table) <- metadata$SampleID
  
  # Extract family names
  extract_family <- function(tax_string) {
    if (grepl("g__", tax_string)) {
      gsub(".*;g__([^;]+).*", "\\1", tax_string)
    } else if (grepl("f__", tax_string)) {
      gsub(".*;f__([^;]+).*", "\\1", tax_string)
    } else {
      NA
    }
  }
  
  family_names <- sapply(colnames(abundance_table), extract_family)
  abundance_table_family <- rowsum(t(abundance_table), group = family_names, na.rm = TRUE)
  abundance_table_family <- t(abundance_table_family)
  abundance_table_family <- abundance_table_family[, !is.na(colnames(abundance_table_family))]
  
  # CLR transformation
  clr_transform <- function(x) {
    x <- x + 1
    log(x / exp(mean(log(x))))
  }
  
  abundance_table_clr <- t(apply(abundance_table_family, 1, clr_transform))
  
  # Prepare data for random forest
  rf_data <- cbind(Diet_type = metadata$Diet_type, as.data.frame(abundance_table_clr))
  rf_data$Diet_type <- as.factor(rf_data$Diet_type)
  
  # Train random forest model
  set.seed(123)
  rf_model <- train(
    Diet_type ~ ., data = rf_data,
    method = "rf",
    trControl = trainControl(method = "cv", number = 5, classProbs = TRUE, savePredictions = "final"),
    ntree = 500, importance = TRUE, metric = "Accuracy"
  )
  
  # Get predictions
  cv_predictions <- rf_model$pred
  cv_predictions$obs <- factor(cv_predictions$obs, levels = levels(rf_data$Diet_type))
  cv_predictions$pred <- factor(cv_predictions$pred, levels = levels(rf_data$Diet_type))
  
  # Create confusion matrix
  conf_matrix_cv <- confusionMatrix(cv_predictions$pred, cv_predictions$obs)
  cm_df <- as.data.frame(conf_matrix_cv$table)
  
  # Create plot
  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = Freq),
              color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "black", "white"),
              size = 5, fontface = "bold") +
    scale_fill_viridis_c(name = "Count") +
    labs(x = "Expected Diet", y = "Predicted Diet", title = plot_title) +
    theme_minimal(base_size = 14) +
    theme(
      plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

data <- read_csv("output/df_taxa_metadata_250.csv")

# Create plot
conf_mat_plot_diet1_counts <- make_confusion_plot_counts(data, remove_outliers = FALSE, plot_title = "")

# Display plot
conf_mat_plot_diet1_counts

# Save the plot in multiple formats
out_dir <- "output_bioproj_prim"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("confusion_matrix_diet1_counts.", ext)),
    plot = conf_mat_plot_diet1_counts,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
}

## ------------ 1.4) identifying misclassified samples

library(tidyverse)
library(caret)

export_misclassified_diet1 <- function(data, 
                                       remove_outliers = TRUE, 
                                       out_dir = "output_bioproj_prim", 
                                       filename = "misclassified_samples_diet1.csv") {
  
  # --- 1. SETUP OUTPUT DIRECTORY ---
  # Create the directory if it does not exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message(paste("Created directory:", out_dir))
  }
  
  # --- 2. DATA PREPROCESSING ---
  if (remove_outliers) {
    outliers <- c("K058683-R", "K058706-R", "K058789-R", "K058689-R", "K058814-R", "K058785-R")
    data <- data[!(data$index %in% outliers), ]
  }
  
  # Prepare metadata (Ensuring we keep the tracking columns)
  metadata <- data.frame(
    SampleID = data$index,
    Diet_type = factor(make.names(data$Diet_type)),
    cw_diet1 = data$cw_diet1,
    bioproject = data$bioproject,
    sample_type = data$sample_type,
    stringsAsFactors = FALSE
  )
  
  # Extract abundance data
  abundance_table <- data[, grep("d__Bacteria", colnames(data))]
  rownames(abundance_table) <- metadata$SampleID
  
  # Extract family names
  extract_family <- function(tax_string) {
    if (grepl("g__", tax_string)) {
      gsub(".*;g__([^;]+).*", "\\1", tax_string)
    } else if (grepl("f__", tax_string)) {
      gsub(".*;f__([^;]+).*", "\\1", tax_string)
    } else {
      NA
    }
  }
  
  family_names <- sapply(colnames(abundance_table), extract_family)
  abundance_table_family <- rowsum(t(abundance_table), group = family_names, na.rm = TRUE)
  abundance_table_family <- t(abundance_table_family)
  abundance_table_family <- abundance_table_family[, !is.na(colnames(abundance_table_family))]
  
  # CLR transformation
  clr_transform <- function(x) {
    x <- x + 1
    log(x / exp(mean(log(x))))
  }
  abundance_table_clr <- t(apply(abundance_table_family, 1, clr_transform))
  
  # Prepare data for random forest
  rf_data <- cbind(Diet_type = metadata$Diet_type, as.data.frame(abundance_table_clr))
  rf_data$Diet_type <- as.factor(rf_data$Diet_type)
  
  # --- 3. TRAIN MODEL ---
  set.seed(1000)
  rf_model <- train(
    Diet_type ~ ., data = rf_data,
    method = "rf",
    # savePredictions = "final" allows us to map back to specific rows later
    trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
    ntree = 500, importance = TRUE
  )
  
  # --- 4. IDENTIFY MISCLASSIFIED SAMPLES ---
  
  # Extract predictions from the final model
  cv_results <- rf_model$pred
  
  # Filter for rows where Prediction (pred) does not match Observation (obs)
  # rowIndex tells us which row in 'metadata' this prediction corresponds to
  misclassified_rows <- cv_results[cv_results$pred != cv_results$obs, ]
  
  if (nrow(misclassified_rows) == 0) {
    message("Great news! The model achieved 100% accuracy. No misclassified samples to export.")
    return(NULL)
  }
  
  # Extract the metadata for these specific rows using rowIndex
  misclassified_info <- metadata[misclassified_rows$rowIndex, ]
  
  # Add the Predicted value to this dataframe
  misclassified_info$Predicted_Diet <- misclassified_rows$pred
  
  # Organize columns
  final_output <- misclassified_info %>%
    rename(Observed_Diet = Diet_type) %>%
    select(SampleID, Observed_Diet, Predicted_Diet, cw_diet1, bioproject, sample_type)
  
  # --- 5. EXPORT ---
  # Combine directory and filename
  full_path <- file.path(out_dir, filename)
  
  write_csv(final_output, full_path)
  
  message(paste("Export complete. Found", nrow(final_output), "misclassified samples."))
  message(paste("Saved to:", full_path))
  
  return(final_output)
}

data <- read_csv("output/df_taxa_metadata_250.csv")

wrong_preds <- export_misclassified_diet1(data, 
                                          remove_outliers = TRUE,
                                          out_dir = "output_bioproj_prim", 
                                          filename = "misclassified_samples_diet1.csv")

## ========================================================================== ##
## ------------------------ FIG 1 PANEL ARRANGEMENT ------------------------- ##
# ------------------------------ Predictor: DIET1 ---------------------------- #

library(ggpubr)

combined_plots_diet1 <- ggarrange(
  alpha_plots_diet1_sig,
  beta_plots_diet1,
  conf_mat_plot_diet1,
  ncol = 1, 
  nrow = 3,
  labels = c("A", "B", "C")
)

# Display the combined plot
print(combined_plots_diet1)

# Save the plot in multiple formats
out_dir <- "output_bioproj_prim/_FIGURES"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("fig1_diet1_panel.", ext)),
    plot = combined_plots_diet1,
    width = 14,
    height = 18,
    units = "in",
    dpi = 300
  )
}

################################################################################
## ----------------- ALPHA, BETA, CONFUSION MATRIX: CW_DIET1 ---------------- ##
################################################################################

# ============================================================================ #
# -----------  2.1 ALPHA DIVERSITY (Visuals + Kruskal-Wallis) ---------------- #
# --------------------------- Predictor: CW_DIET1 ---------------------------- #
# ============================================================================ #

set.seed(1000)
tree_file <- rtree(ncol(community_matrix), tip.label = colnames(community_matrix))

shannon  <- diversity(community_matrix, index = "shannon")
observed <- specnumber(community_matrix)
simpson <- diversity(community_matrix, index = "simpson")
#pd_result <- pd(community_matrix, tree_file, include.root = FALSE)

# Calculate Indices
alpha_div <- data.frame(
  Shannon = shannon,
  Simpson = simpson,
  Observed = observed
)

# Merge back with metadata
alpha_plot_data <- cbind(alpha_div, metadata_sub)

# ============================================================================ #
# ---------- NOT INCLUDED: BOXPLOT FUNCTION WITH SIGNIFICANCE BRACKETS ------- #
# ---------------------------- Predictor: CW_DIET1 --------------------------- #
# ============================================================================ #

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix) # Helper for statistical pipes

plot_alpha <- function(data, y_metric, title) {
  
  # 1. CALCULATE DUNN'S TEST (No positions yet)
  # Ensure the grouping variable is a factor
  data$cw_diet1 <- as.factor(data$cw_diet1)
  
  # Run the test
  dunn_stats <- data %>%
    dunn_test(as.formula(paste(y_metric, "~ cw_diet1")), p.adjust.method = "fdr") %>%
    add_significance()
  
  # 2. MANUALLY CALCULATE Y-POSITIONS FOR NEAT STACKING
  # This guarantees consistent spacing for ALL brackets (Sig + NS)
  if(nrow(dunn_stats) > 0) {
    # Find the maximum value in the plot to know where to start drawing
    max_y <- max(data[[y_metric]], na.rm = TRUE)
    
    # Define a step size (e.g., 5% of the data range to avoid overlap)
    y_range <- max_y - min(data[[y_metric]], na.rm = TRUE)
    step_size <- y_range * 0.05
    
    # Assign y.position incrementally for every comparison in the list
    # This mimics the "slot" behavior of your original script
    dunn_stats$y.position <- max_y + (1:nrow(dunn_stats)) * step_size
  }
  
  # 3. Calculate global p-value (Kruskal-Wallis) for subtitle
  kruskal_p <- kruskal_test(data, as.formula(paste(y_metric, "~ cw_diet1")))$p
  subtitle_text <- paste0("Kruskal-Wallis, p = ", sprintf("%.3f", kruskal_p))
  
  # 4. Plot Base
  p <- ggplot(data, aes(x = cw_diet1, y = .data[[y_metric]], fill = cw_diet1)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.4, color = "black") +
    theme_bw() +
    labs(title = title, subtitle = subtitle_text, x = "Captive/Wild Status + Diet Group", y = y_metric) +
    scale_fill_manual(values = cw_diet_colors) +
    
    # --- LAYER A: NON-SIGNIFICANT (RED) ---
    # Filter for NS, but use the manually calculated y.positions
    stat_pvalue_manual(
      dunn_stats %>% filter(p.adj > 0.05),
      label = "p.adj.signif",
      color = "red",
      tip.length = 0.005,
      vjust = 0.5,
      remove.bracket = FALSE
    ) +
    
    # --- LAYER B: SIGNIFICANT (BLACK) ---
    # Filter for Sig, using the same manual position logic
    stat_pvalue_manual(
      dunn_stats %>% filter(p.adj <= 0.05),
      label = "p.adj.signif",
      color = "black",
      tip.length = 0.005,
      vjust = 0.5,
      remove.bracket = FALSE
    ) +
    
    # Theme adjustments
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10)
    )
  
  return(p)
}

# Generate Plots
p1 <- plot_alpha(alpha_plot_data, "Shannon", "Shannon Entropy Index")
p2 <- plot_alpha(alpha_plot_data, "Simpson", "Simpson's Diversity Index")
p3 <- plot_alpha(alpha_plot_data, "Observed", "Observed Richness")

# Combine and Save
alpha_plots_cw_diet1 <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
print(alpha_plots_cw_diet1)

out_dir <- "output_bioproj_prim"
dir.create("output_bioproj_prim", showWarnings = FALSE, recursive = TRUE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("alpha_cw_diet1_sig_nonsig.", ext)),
    plot = alpha_plots_cw_diet1,
    width = 16, 
    height = 10, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================ #
# -------------- 2.1) BOXPLOT FUNCTION WITH SIG BRACKETS ONLY ---------------- #
# --------------------------- Predictor: CW_DIET1 ---------------------------- #
# ============================================================================ #

plot_alpha_sig_only <- function(data, y_metric, title) {
  
  # 1. CALCULATE DUNN'S TEST
  stat.test <- data %>%
    dunn_test(as.formula(paste(y_metric, "~ cw_diet1")), p.adjust.method = "fdr") %>%
    add_significance()
  
  # 2. FILTER FOR SIGNIFICANT PAIRS
  sig_stats <- stat.test %>% filter(p.adj <= 0.05)
  
  # 3. MANUALLY CALCULATE Y-POSITIONS FOR NEAT STACKING
  # This mimics what stat_compare_means does internally.
  if(nrow(sig_stats) > 0) {
    
    # Find the maximum value in the plot to know where to start drawing
    max_y <- max(data[[y_metric]], na.rm = TRUE)
    
    # Define a step size (e.g., 5% of the range)
    y_range <- max_y - min(data[[y_metric]], na.rm = TRUE)
    step_size <- y_range * 0.05
    
    # Assign y.position incrementally
    # The first bracket goes at max_y + step, the second at max_y + 2*step, etc.
    sig_stats$y.position <- max_y + (1:nrow(sig_stats)) * step_size
  }
  
  # 4. GENERATE PLOT
  p <- ggplot(data, aes(x = cw_diet1, y = .data[[y_metric]], fill = cw_diet1)) +
    
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.0, alpha = 0.6, color = "black") +
    
    theme_bw() +
    labs(title = title, x = "Captive/Wild Status + Diet Group", y = y_metric) +
    scale_fill_manual(values = cw_diet_colors) +
    
    # 5. ADD BRACKETS
    stat_pvalue_manual(
      sig_stats, 
      label = "p.adj.signif",    
      tip.length = 0.01,
      hide.ns = TRUE,
      vjust = 0.7
    ) +
    
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10)
    )
  
  return(p)
}
# Generate Plots
p1_sig <- plot_alpha_sig_only(alpha_plot_data, "Shannon", "Shannon Diversity")
p2_sig <- plot_alpha_sig_only(alpha_plot_data, "Simpson", "Simpson's Diversity Index")
p3_sig <- plot_alpha_sig_only(alpha_plot_data, "Observed", "Observed Richness")

# Combine
alpha_plots_cw_diet1_sig <- ggarrange(p1_sig, p2_sig, p3_sig, ncol = 3, nrow = 1)
print(alpha_plots_cw_diet1_sig)

out_dir <- "output_bioproj_prim"
dir.create("output_bioproj_prim", showWarnings = FALSE, recursive = TRUE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("alpha_cw_diet1_sig.", ext)),
    plot = alpha_plots_cw_diet1_sig,
    width = 16, 
    height = 10, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================ #
# ----------------------------- STATISTICS: CW_DIET1 ------------------------- #
# ============================================================================ #

library(rstatix)
library(dplyr)
library(tibble)

# --- Assumption ---
# This script assumes you have a data frame named 'alpha_plot_data'
# which contains the columns 'Shannon', 'Simpson', 'Observed', and 'cw_diet1'.

# --- 1. Perform Global Kruskal-Wallis Tests ---
# The formula is updated to use 'cw_diet1' as the grouping variable.

kruskal_shannon <- alpha_plot_data %>% kruskal_test(Shannon ~ cw_diet1) 
kruskal_simpson <- alpha_plot_data %>% kruskal_test(Simpson ~ cw_diet1) 
kruskal_observed <- alpha_plot_data %>% kruskal_test(Observed ~ cw_diet1) 

# Combine the global test results into a single data frame
all_kruskal_results <- bind_rows(kruskal_shannon, kruskal_simpson, kruskal_observed)

# Define the output directory and filename for Kruskal-Wallis results
out_dir <- "output_bioproj_prim/statistics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
kruskal_output_file <- file.path(out_dir, "alpha_cw_diet1_kruskal_results.csv")

# Save the global test results to CSV
write.csv(all_kruskal_results, file = kruskal_output_file, row.names = FALSE)
message("Global Kruskal-Wallis test results saved to: ", kruskal_output_file)


# --- 2. Perform Pairwise Wilcoxon Tests ---
# The formula is updated here as well to use 'cw_diet1'.

pairwise_shannon <- alpha_plot_data %>%
  wilcox_test(Shannon ~ cw_diet1, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(metric = "Shannon")

pairwise_simpson <- alpha_plot_data %>%
  wilcox_test(Simpson ~ cw_diet1, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>%
  mutate(metric = "Simpson")

pairwise_observed <- alpha_plot_data %>%
  wilcox_test(Observed ~ cw_diet1, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(metric = "Observed")

# Combine all pairwise results into a single data frame
all_pairwise_results <- bind_rows(pairwise_shannon, pairwise_simpson, pairwise_observed)

# Define the output filename for the pairwise results
pairwise_output_file <- file.path(out_dir, "alpha_cw_diet1_pairwise_wilcox_results.csv")

# Save the pairwise test results to CSV
write.csv(all_pairwise_results, file = pairwise_output_file, row.names = FALSE)
message("Pairwise Wilcoxon test results saved to: ", pairwise_output_file)

# -------------------------- 2) Dunn's Test to CSV ---------------------------------- #

# Load necessary libraries
library(phyloseq)
library(vegan)
library(rstatix)
library(dplyr)
library(tibble)

# --- Assumptions ---
# This script assumes you have a phyloseq object named 'physeq'
# with a 'cw_diet1' column in its sample data.

# --- 1. Calculate Alpha Diversity Metrics ---
alpha_div_richness <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson"))

# --- 2. Combine All Data into a Single Data Frame ---
metadata <- data.frame(sample_data(physeq))
all_alpha_data <- metadata %>%
  rownames_to_column("SampleID") %>%
  left_join(rownames_to_column(alpha_div_richness, "SampleID"))

# --- 3. Perform Pairwise Dunn's Tests using 'cw_diet1' ---
metrics_to_test <- c("Observed", "Shannon", "Simpson")
all_results_list <- list()

for (metric_name in metrics_to_test) {
  
  formula <- as.formula(paste(metric_name, "~ cw_diet1"))
  
  current_data <- all_alpha_data %>%
    select(cw_diet1, all_of(metric_name)) %>%
    filter(!is.na(.data[[metric_name]]))
  
  group_counts <- current_data %>%
    count(cw_diet1) %>% ### CHANGED ###
    filter(n > 1)
  
  if (nrow(group_counts) >= 2) {
    message(paste("Performing Dunn's test for:", metric_name))
    dunn_result <- dunn_test(data = current_data, formula = formula, p.adjust.method = "fdr") %>%
      mutate(metric = metric_name)
    all_results_list[[metric_name]] <- dunn_result
  } else {
    message(paste("Skipping Dunn's test for:", metric_name, "- insufficient data."))
  }
}

if (length(all_results_list) > 0) {
  final_dunn_results <- bind_rows(all_results_list)
  
  out_dir <- "output_bioproj_prim/statistics"
  output_filename <- "alpha_cw_diet1_dunn_test.csv"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  full_path <- file.path(out_dir, output_filename)
  
  write.csv(final_dunn_results, file = full_path, row.names = FALSE)
  message("Successfully saved pairwise Dunn's test results to: ", full_path)
} else {
  message("Analysis complete, but no metrics had sufficient data to perform pairwise tests.")
}

# ============================================================================ #
## ---------------- NOT INCLUDED: Beta diversity PCOA ------------------------ #
# ============================================================================ #

# Merge with metadata for plotting
pcoa_plot_data <- cbind(pcoa_df, metadata_sub)

# 1. DEFINE CUSTOM SHAPES
diet_shapes <- c(
  "Blood"  = 16,
  "Insect" = 3,
  "Fruit"  = 15,
  "Fish"   = 17) 

# 16 = Circle, 15 = Square, 17 = Triangle, 18 = Diamond
cw_diet_shapes <- c(
  # --- INSECT (diamond) ---
  "Wild_Insect"              = 18,
  "Captive_Insect"           = 3,
  "Wild_probably_Insect"     = 18,
  "Captive_probably_Insect"  = 3,
  
  # --- FRUIT (Squares) ---
  "Wild_Fruit"               = 15,
  "Wild_probably_Fruit"      = 15,
  "Captive_probably_Fruit"   = 15,
  
  # --- BLOOD (cirlce) ---
  "Wild_Blood"               = 16,
  "Wild_probably_Blood"      = 16,
  
  # --- FISH (triangle) ---
  "Wild_probably_Fish"       = 17
)

# 2. PLOT WITH MANUAL SHAPES
beta_plot_cwdiet <- ggplot(pcoa_plot_data, aes(x = PC1, y = PC2, color = cw_diet1, shape = cw_diet1)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = cw_diet1), geom = "polygon", alpha = 0.1, level = 0.95) +
  theme_bw() +
  labs(
    title = "PCoA of Bray-Curtis Dissimilarity",
    x = paste0("PCoA 1 (", eig_percent[1], "%)"),
    y = paste0("PCoA 2 (", eig_percent[2], "%)"),
    
    color = "Captive/Wild Status + Diet Group",
    shape = "Captive/Wild Status + Diet Group",
    fill  = "Captive/Wild Status + Diet Group"
  ) +
  
  
  # --- COLORS ---
  scale_color_manual(values = cw_diet_colors) +  
  scale_fill_manual(values = cw_diet_colors) +
  
  # --- SHAPES (NEW) ---
  scale_shape_manual(values = cw_diet_shapes) +  
  
  # --- THEME ---
  theme(
    plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10), 
    axis.title = element_text(size = 8, face = "bold"),              
    axis.text = element_text(size = 8),                              
    legend.title = element_text(size = 10, face = "bold"),            
    legend.text = element_text(size = 8)                             
  )

print(beta_plot_cwdiet)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("beta_bray_cwdiet1.", ext)),
    plot = beta_plot_cwdiet,
    width = 10, 
    height = 8, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================ #
## ----------------- 2.2) BETA DIVERSITY MULTIPLE PCOA PLOTS ---------------- ##
## --------------------------- Predictor: CW_DIET1 -------------------------- ##
# ============================================================================ #

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(phyloseq)
library(ape) 

# ---------- 0. SETUP COLORS & SHAPES (From previous steps)

cw_diet_colors <- c(
  "Wild_Insect" = "#ff5e2b", 
  "Captive_Insect" = "#e37f05", 
  "Wild_probably_Insect" = "#fa8246", 
  "Captive_probably_Insect" = "#fcb644",
  "Wild_Fruit" = "#2E8B57", 
  "Wild_probably_Fruit" = "#1b5e3a", 
  "Captive_probably_Fruit" = "#80e35f",
  "Wild_Blood" = "#ba0970", 
  "Wild_probably_Blood" = "#e329da", 
  "Wild_probably_Fish" = "#303F9F"
)

cw_diet_shapes <- c(
  # --- INSECT (diamond) ---
  "Wild_Insect"              = 18,
  "Captive_Insect"           = 3,
  "Wild_probably_Insect"     = 18,
  "Captive_probably_Insect"  = 3,
  
  # --- FRUIT (Squares) ---
  "Wild_Fruit"               = 15,
  "Wild_probably_Fruit"      = 15,
  "Captive_probably_Fruit"   = 15,
  
  # --- BLOOD (cirlce) ---
  "Wild_Blood"               = 16,
  "Wild_probably_Blood"      = 16,
  
  # --- FISH (triangle) ---
  "Wild_probably_Fish"       = 17
)

# ------------ A. PREPARE PHYLOSEQ OBJECT

# A. OTU Table (Taxa must be columns, Samples must be rows based on your previous code)
otu_mat <- as.matrix(community_matrix)

# B. Taxonomy Table (Optional for distance calc, but good practice)
# create a dummy tax table if you don't have the original 'taxonomy_df' handy
# For 'taxonomy_df', use: tax_mat <- as.matrix(taxonomy_df)
tax_mat <- matrix(colnames(otu_mat), ncol = 1)
rownames(tax_mat) <- colnames(otu_mat)
colnames(tax_mat) <- "Species"

# C. Metadata
meta_df <- as.data.frame(metadata_sub)
rownames(meta_df) <- rownames(metadata_sub)

# D. Phylogenetic Tree (CRITICAL FOR UNIFRAC)
# ------------------------------------------------------------------
# IF USE A TREE FILE, UNCOMMENT AND USE THIS LINE:
# tree_file <- read_tree("path/to/your/tree.nwk")

set.seed(1000)
tree_file <- rtree(ncol(otu_mat), tip.label = colnames(otu_mat))
# ------------------------------------------------------------------

# Create the Phyloseq Object
physeq <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = FALSE),
  tax_table(tax_mat),
  sample_data(meta_df),
  phy_tree(tree_file)
)

# -------------- B. FUNCTION TO GENERATE PCoA PLOTS
plot_beta_metric <- function(physeq_obj, distance_method, title, weighted = FALSE) {
  
  # Calculate Distance
  # UniFrac requires special handling in the ordinate function
  if (distance_method == "wunifrac") {
    ord <- ordinate(physeq_obj, method = "PCoA", distance = "unifrac", weighted = TRUE)
  } else if (distance_method == "unifrac") {
    ord <- ordinate(physeq_obj, method = "PCoA", distance = "unifrac", weighted = FALSE)
  } else {
    ord <- ordinate(physeq_obj, method = "PCoA", distance = distance_method)
  }
  
  # Plot using phyloseq's built-in plotter (easier wrapper for ggplot)
  p <- plot_ordination(physeq_obj, ord, color = "cw_diet1", shape = "cw_diet1") +
    geom_point(size = 1.5, alpha = 0.8) +
    stat_ellipse(aes(fill = cw_diet1), geom = "polygon", alpha = 0.1, level = 0.95) +
    theme_bw() +
    labs(title = title) +
    
    # Custom Colors and Shapes
    scale_color_manual(values = cw_diet_colors) +
    scale_fill_manual(values = cw_diet_colors) +
    scale_shape_manual(values = cw_diet_shapes) +
    
    # Theme settings
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      legend.position = "none" # Hide legend for individual plots
    )
  
  return(p)
}

# ------------ C. GENERATE THE THREE PCOA PLOTS

# 1. Bray-Curtis (Abundance based)
p_bray <- plot_beta_metric(physeq, "bray", "Bray-Curtis (Abundance)")

# 2. Jaccard (Presence/Absence based)
# binary = TRUE is handled automatically by phyloseq for "jaccard" usually, 
# but "jaccard" in vegan is inherently quantitative. 
# For true binary Jaccard in phyloseq, use distance="jaccard", binary=TRUE
# Simulated here by converting data to presence/absence first for this plot

physeq_binary <- transform_sample_counts(physeq, function(x) ifelse(x > 0, 1, 0))
#p_jaccard <- plot_beta_metric(physeq_binary, "jaccard", "Jaccard (Presence/Absence)")

# 3. Weighted UniFrac (Phylogeny + Abundance)
p_wunifrac <- plot_beta_metric(physeq, "wunifrac", "Weighted UniFrac (Phylogeny + Abundance)")

# 4. Unweighted UniFrac (Phylogeny + Presence/Absence)
p_unifrac <- plot_beta_metric(physeq, "unifrac", "Unweighted UniFrac (Phylogeny + Presence/Absence)")

# ------------ D. COMBINE AND SAVE
# Extract the legend from one of the plots to share it
p_bray <- p_bray + 
  scale_color_manual(name = "Captive/Wild Status + Diet Group", values = cw_diet_colors) +
  scale_fill_manual(name = "Captive/Wild Status + Diet Group", values = cw_diet_colors) +
  scale_shape_manual(name = "Captive/Wild Status + Diet Group", values = cw_diet_shapes) +   
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10)
  )

# 2. Extract the legend again with the new title
legend_bray <- get_legend(p_bray + theme(legend.position = "right"))

beta_plots_cw_diet1 <- ggarrange(
  p_bray, p_wunifrac, p_unifrac,
  ncol = 3, nrow = 1,
  common.legend = TRUE,
  legend = "right",
  legend.grob = legend_bray
)

print(beta_plots_cw_diet1)

# Define output directory
out_dir <- "output_bioproj_prim"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("beta_cwdiet1_panel.", ext)),
    plot = beta_plots_cw_diet1,
    width = 16, 
    height = 6, 
    units = "in", 
    dpi = 300,
    bg = "white"
  )
}

# ----------------------------- STATISTICS ----------------------------------- #

library(phyloseq)
library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(tibble)
library(broom) 

# --- 1. Pre-computation of Distances ---

metadata <- data.frame(sample_data(physeq))
dist_bray <- distance(physeq, method = "bray")
dist_wunifrac <- distance(physeq, method = "wunifrac")
dist_uunifrac <- distance(physeq, method = "unifrac")

# Define output directory
out_dir <- "output_bioproj_prim/statistics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. Homogeneity of Dispersion (PERMDISP) to CSV ---

permadisp_bray <- tidy(anova(betadisper(dist_bray, metadata$cw_diet1))) %>% mutate(metric = "Bray-Curtis")
permadisp_wunifrac <- tidy(anova(betadisper(dist_wunifrac, metadata$cw_diet1))) %>% mutate(metric = "Weighted UniFrac")
permadisp_uunifrac <- tidy(anova(betadisper(dist_uunifrac, metadata$cw_diet1))) %>% mutate(metric = "Unweighted UniFrac")

all_permadisp_results <- bind_rows(permadisp_bray, permadisp_wunifrac, permadisp_uunifrac)
permadisp_output_file <- file.path(out_dir, "beta_cw_diet1_permadisp.csv")
write.csv(all_permadisp_results, file = permadisp_output_file, row.names = FALSE)
message("PERMDISP results saved to: ", permadisp_output_file)

# --- 3. Global PERMANOVA (adonis2) to CSV ---

permanova_bray <- tidy(adonis2(dist_bray ~ cw_diet1, data = metadata, permutations = 999)) %>% mutate(metric = "Bray-Curtis")
permanova_wunifrac <- tidy(adonis2(dist_wunifrac ~ cw_diet1, data = metadata, permutations = 999)) %>% mutate(metric = "Weighted UniFrac")
permanova_uunifrac <- tidy(adonis2(dist_uunifrac ~ cw_diet1, data = metadata, permutations = 999)) %>% mutate(metric = "Unweighted UniFrac")

all_permanova_results <- bind_rows(permanova_bray, permanova_wunifrac, permanova_uunifrac)
permanova_output_file <- file.path(out_dir, "beta_cw_diet1_permanova_global.csv")
write.csv(all_permanova_results, file = permanova_output_file, row.names = FALSE)
message("Global PERMANOVA results saved to: ", permanova_output_file)

# --- 4. Pairwise PERMANOVA to TXT ---
metadata <- data.frame(sample_data(physeq))

# --- 1. Pre-computation of Distances ---
metadata <- data.frame(sample_data(physeq))
dist_bray <- distance(physeq, method = "bray")
physeq_binary <- transform_sample_counts(physeq, function(x) ifelse(x > 0, 1, 0))
dist_wunifrac <- distance(physeq, method = "wunifrac")
dist_uunifrac <- distance(physeq, method = "unifrac")

# --- 2. Run Statistical Analyses and Generate Report ---
output_file <- "output_bioproj_prim/statistics/beta_cw_diet1_pairwise.txt"
sink(output_file)
# --- Header ---
cat("==========================================================\n")
cat("          BETA DIVERSITY STATISTICAL ANALYSIS\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Grouping Variable: cw_diet1\\\n")
cat("==========================================================\n\n")
# --- A. BRAY-CURTIS ANALYSIS ---
cat("A. BRAY-CURTIS (Abundance)\n")
cat("==========================================================\n\n")
cat("1. Homogeneity of Dispersion (PERMDISP):\n")
print(anova(betadisper(dist_bray, metadata$cw_diet1)))
cat("\n2. Overall PERMANOVA (adonis2):\n")
print(adonis2(dist_bray ~ cw_diet1, data = metadata, permutations = 999))
cat("\n3. Pairwise PERMANOVA (p-values adjusted with FDR):\n")
pairwise_bray <- pairwise.adonis2(dist_bray ~ cw_diet1, data = metadata, p.adjust.m = "fdr")
print(pairwise_bray)
cat("\n\n")

# --- C. WEIGHTED UNIFRAC ANALYSIS ---
cat("C. WEIGHTED UNIFRAC (Phylogeny + Abundance)\n")
cat("==========================================================\n\n")
cat("1. Homogeneity of Dispersion (PERMDISP):\n")
print(anova(betadisper(dist_wunifrac, metadata$cw_diet1)))
cat("\n2. Overall PERMANOVA (adonis2):\n")
print(adonis2(dist_wunifrac ~ cw_diet1, data = metadata, permutations = 999))
cat("\n3. Pairwise PERMANOVA (p-values adjusted with FDR):\n")
pairwise_wunifrac <- pairwise.adonis2(dist_wunifrac ~ cw_diet1, data = metadata, p.adjust.m = "fdr")
print(pairwise_wunifrac)
cat("\n\n")

# --- D. UNWEIGHTED UNIFRAC ANALYSIS ---
cat("D. UNWEIGHTED UNIFRAC (Phylogeny + Presence/Absence)\n")
cat("==========================================================\n\n")
cat("1. Homogeneity of Dispersion (PERMDISP):\n")
print(anova(betadisper(dist_uunifrac, metadata$cw_diet1)))
cat("\n2. Overall PERMANOVA (adonis2):\n")
print(adonis2(dist_uunifrac ~ cw_diet1, data = metadata, permutations = 999))
cat("\n3. Pairwise PERMANOVA (p-values adjusted with FDR):\n")
pairwise_uunifrac <- pairwise.adonis2(dist_uunifrac ~ cw_diet1, data = metadata, p.adjust.m = "fdr")
print(pairwise_uunifrac)
cat("\n\n")

# --- End of Report ---
cat("==========================================================\n")
cat("                      End of Report\n")
cat("==========================================================\n")
sink()
message("Complete statistical results with pairwise tests have been saved to: ", output_file)

# ============================================================================ #
## ----------------------- 2.3) CONFUSION MATRIX ---------------------------- ##
## --------------------------- Predictor: CW_DIET1 -------------------------- ##
# ============================================================================ #

# Random Forest Classification Analysis
library(tidyverse)
library(caret)
library(patchwork)
library(caret)
library(dplyr)
library(ggplot2)
library(tibble)

# Function to generate confusion matrix plot
make_confusion_plot <- function(data, remove_outliers = TRUE, plot_title = "") {
  
  # 1. Clean Data & Outliers
  if (remove_outliers) {
    outliers <- c("K058683-R", "K058706-R", "K058789-R", "K058689-R", "K058814-R", "K058785-R")
    data <- data %>% filter(!index %in% outliers)
  }
  
  # 2. Prepare Metadata & Response Variable
  # make.names is crucial here because Random Forest formulas hate spaces/special chars
  clean_diet <- make.names(factor(data$cw_diet1))
  
  metadata <- data.frame(
    SampleID = data$index,
    cw_diet1 = factor(clean_diet)
  )
  
  # 3. Extract Abundance Data
  # We assume abundance columns contain "d__Bacteria"
  # Check if we have columns to work with
  abund_cols <- grep("d__Bacteria", colnames(data))
  if(length(abund_cols) == 0) stop("No columns containing 'd__Bacteria' found.")
  
  abundance_table <- data[, abund_cols]
  rownames(abundance_table) <- metadata$SampleID
  
  # 4. Extract Family Names & Aggregate
  extract_family <- function(tax_string) {
    # Looks for g__ first, if not then f__
    if (grepl("f__", tax_string)) {
      # Extract text between f__ and the next ; or end of string
      sub(".*;f__([^;]+).*", "\\1", tax_string)
    } else {
      "Unclassified"
    }
  }
  
  # Get family names for every column
  family_names <- sapply(colnames(abundance_table), extract_family)
  
  # Transpose so taxa are rows, sum by family, then transpose back
  # rowsum requires the aggregation group (family_names) to match the number of rows
  t_abundance <- t(abundance_table)
  abundance_table_family <- rowsum(t_abundance, group = family_names, na.rm = TRUE)
  
  # Transpose back so Samples are Rows
  abundance_table_family <- t(abundance_table_family)
  
  # 5. CLR Transformation
  # CLR requires handling zeros. Adding 1 is a simple pseudo-count method.
  clr_transform <- function(x) {
    x <- x + 1
    log(x) - mean(log(x)) # Standard CLR formula
  }
  
  abundance_table_clr <- as.data.frame(t(apply(abundance_table_family, 1, clr_transform)))
  
  # 6. Prepare Data for Random Forest
  # Combine the response variable (cw_diet1) with the features
  rf_data <- cbind(cw_diet1 = metadata$cw_diet1, abundance_table_clr)
  
  # 7. Train Random Forest
  set.seed(1000)
  
  # Setup Cross Validation
  ctrl <- trainControl(
    method = "cv", 
    number = 5, 
    classProbs = TRUE, 
    savePredictions = "final"
  )
  
  rf_model <- train(
    cw_diet1 ~ ., 
    data = rf_data,
    method = "rf",
    trControl = ctrl,
    ntree = 500, 
    metric = "Accuracy"
  )
  
  # 8. Generate Confusion Matrix
  # Extract predictions from the 'final' hold-out fold
  cv_preds <- rf_model$pred
  
  # Ensure the factor levels match the model's levels exactly
  # (This fixes the "overlap" error)
  obs_levels <- levels(rf_data$cw_diet1)
  
  # Create the matrix
  conf_matrix_cv <- confusionMatrix(
    data = factor(cv_preds$pred, levels = obs_levels), 
    reference = factor(cv_preds$obs, levels = obs_levels)
  )
  
  # Convert to dataframe for plotting
  cm_df <- as.data.frame(conf_matrix_cv$table)
  
  # 9. Calculate Ratios for Heatmap
  cm_df <- cm_df %>%
    group_by(Reference) %>%
    mutate(
      Total = sum(Freq),
      Ratio = Freq / Total,
      Ratio = ifelse(is.nan(Ratio), 0, Ratio), 
      Label = sprintf("%.2f", Ratio)
    ) %>%
    ungroup()
  
  # 10. Create Plot
  p <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Ratio)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = Label),
              color = ifelse(cm_df$Ratio > 0.5, "black", "white"),
              size = 5, fontface = "bold") +
    scale_fill_viridis_c(name = "Ratio", labels = function(x) sprintf("%.2f", x)) +
    labs(x = "Expected Diet", y = "Predicted Diet", title = plot_title) +
    theme_minimal(base_size = 14) +
    theme(
      plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  return(p)
}

# ---------

library(readr)
library(ggplot2)
library(patchwork)

data <- read_csv("output/df_taxa_metadata_250.csv")

conf_mat_plot_cw_diet1 <- make_confusion_plot(data, remove_outliers = TRUE, plot_title = "")

print(conf_mat_plot_cw_diet1)

# Save the plot in multiple formats
out_dir <- "output_bioproj_prim"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("confusion_matrix_cw_diet1.", ext)),
    plot = conf_mat_plot_cw_diet1,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
}

## ------------------- confusion matrix counts ---------------------------------

make_confusion_plot_counts <- function(data, remove_outliers = TRUE, plot_title = "") {
  
  if (remove_outliers) {
    outliers <- c("K058683-R", "K058706-R", "K058789-R", "K058689-R", "K058814-R", "K058785-R")
    # Ensure 'index' exists in your data, otherwise change to the correct ID column
    data <- data[!(data$index %in% outliers), ]
  }
  
  metadata <- data.frame(
    SampleID = data$index,
    cw_diet1 = factor(make.names(data$cw_diet1)), 
    stringsAsFactors = FALSE
  )
  
  # Extract abundance data
  abundance_table <- data[, grep("d__Bacteria", colnames(data))]
  rownames(abundance_table) <- metadata$SampleID
  
  # Extract family names (Logic unchanged)
  extract_family <- function(tax_string) {
    if (grepl("g__", tax_string)) {
      gsub(".*;g__([^;]+).*", "\\1", tax_string)
    } else if (grepl("f__", tax_string)) {
      gsub(".*;f__([^;]+).*", "\\1", tax_string)
    } else {
      NA
    }
  }
  
  family_names <- sapply(colnames(abundance_table), extract_family)
  
  # Aggregate by extracted taxonomy
  abundance_table_family <- rowsum(t(abundance_table), group = family_names, na.rm = TRUE)
  abundance_table_family <- t(abundance_table_family)
  abundance_table_family <- abundance_table_family[, !is.na(colnames(abundance_table_family))]
  
  # CLR transformation
  clr_transform <- function(x) {
    x <- x + 1
    log(x / exp(mean(log(x))))
  }
  
  abundance_table_clr <- t(apply(abundance_table_family, 1, clr_transform))
  
  # Prepare data for random forest (Changed Diet_type to cw_diet1)
  rf_data <- cbind(cw_diet1 = metadata$cw_diet1, as.data.frame(abundance_table_clr))
  rf_data$cw_diet1 <- as.factor(rf_data$cw_diet1)
  
  # Train random forest model (Formula changed to cw_diet1 ~ .)
  set.seed(123)
  rf_model <- train(
    cw_diet1 ~ ., data = rf_data,
    method = "rf",
    trControl = trainControl(method = "cv", number = 5, classProbs = TRUE, savePredictions = "final"),
    ntree = 500, importance = TRUE, metric = "Accuracy"
  )
  
  # Get predictions
  cv_predictions <- rf_model$pred
  
  # Align levels based on cw_diet1
  cv_predictions$obs <- factor(cv_predictions$obs, levels = levels(rf_data$cw_diet1))
  cv_predictions$pred <- factor(cv_predictions$pred, levels = levels(rf_data$cw_diet1))
  
  # Create confusion matrix
  conf_matrix_cv <- confusionMatrix(cv_predictions$pred, cv_predictions$obs)
  cm_df <- as.data.frame(conf_matrix_cv$table)
  
  # Create plot
  ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = Freq),
              color = ifelse(cm_df$Freq > max(cm_df$Freq) * 0.5, "black", "white"),
              size = 5, fontface = "bold") +
    scale_fill_viridis_c(name = "Count") +
    labs(x = "Expected cw_diet1", y = "Predicted cw_diet1", title = plot_title) +
    theme_minimal(base_size = 14) +
    theme(
      plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

data <- read_csv("output/df_taxa_metadata_250.csv")

# Create plot
conf_mat_plot_cw_diet1_counts <- make_confusion_plot_counts(data, remove_outliers = FALSE, plot_title = "")

# Display plot
conf_mat_plot_cw_diet1_counts

# Save the plot in multiple formats
out_dir <- "output_bioproj_prim"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("confusion_matrix_cw_diet1_counts.", ext)),
    plot = conf_mat_plot_cw_diet1_counts,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
}

## ------------ 1.4) identifying misclassified samples

library(tidyverse)
library(caret)

export_misclassified_cw_diet1 <- function(data,
                                          remove_outliers = TRUE,
                                          out_dir = "output_bioproj_prim",
                                          filename = "misclassified_samples_cw_diet1.csv") {
  
  # --- 1. SETUP OUTPUT DIRECTORY ---
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message(paste("Created directory:", out_dir))
  }
  
  # --- 2. DATA PREPROCESSING ---
  if (remove_outliers) {
    outliers <- c("K058683-R", "K058706-R", "K058789-R", "K058689-R", "K058814-R", "K058785-R")
    data <- data[!(data$index %in% outliers), ]
  }
  
  # Prepare metadata
  # CHANGED: cw_diet1 is now the primary factor target
  metadata <- data.frame(
    SampleID = data$index,
    cw_diet1 = factor(make.names(data$cw_diet1)), 
    bioproject = data$bioproject,
    sample_type = data$sample_type,
    stringsAsFactors = FALSE
  )
  
  # Extract abundance data
  abundance_table <- data[, grep("d__Bacteria", colnames(data))]
  rownames(abundance_table) <- metadata$SampleID
  
  # Extract family names
  extract_family <- function(tax_string) {
    if (grepl("g__", tax_string)) {
      gsub(".*;g__([^;]+).*", "\\1", tax_string)
    } else if (grepl("f__", tax_string)) {
      gsub(".*;f__([^;]+).*", "\\1", tax_string)
    } else {
      NA
    }
  }
  
  family_names <- sapply(colnames(abundance_table), extract_family)
  abundance_table_family <- rowsum(t(abundance_table), group = family_names, na.rm = TRUE)
  abundance_table_family <- t(abundance_table_family)
  abundance_table_family <- abundance_table_family[, !is.na(colnames(abundance_table_family))]
  
  # CLR transformation
  clr_transform <- function(x) {
    x <- x + 1
    log(x / exp(mean(log(x))))
  }
  abundance_table_clr <- t(apply(abundance_table_family, 1, clr_transform))
  
  # Prepare data for random forest
  # CHANGED: Renamed the target column to cw_diet1
  rf_data <- cbind(cw_diet1 = metadata$cw_diet1, as.data.frame(abundance_table_clr))
  
  # --- 3. TRAIN MODEL ---
  set.seed(1000)
  
  # CHANGED: Updated formula to predict cw_diet1
  rf_model <- train(
    cw_diet1 ~ ., data = rf_data,
    method = "rf",
    trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
    ntree = 500, importance = TRUE
  )
  
  # --- 4. IDENTIFY MISCLASSIFIED SAMPLES ---
  cv_results <- rf_model$pred
  
  # Filter for rows where Prediction (pred) does not match Observation (obs)
  misclassified_rows <- cv_results[cv_results$pred != cv_results$obs, ]
  
  if (nrow(misclassified_rows) == 0) {
    message("Great news! The model achieved 100% accuracy. No misclassified samples to export.")
    return(NULL)
  }
  
  # Extract the metadata for these specific rows
  misclassified_info <- metadata[misclassified_rows$rowIndex, ]
  
  # Add the Predicted value to this dataframe
  misclassified_info$Predicted_Diet <- misclassified_rows$pred
  
  # Organize columns
  # CHANGED: Updated select statement to match the new variable names
  final_output <- misclassified_info %>%
    rename(Observed_Diet = cw_diet1) %>%
    select(SampleID, Observed_Diet, Predicted_Diet, bioproject, sample_type)
  
  # --- 5. EXPORT ---
  full_path <- file.path(out_dir, filename)
  write_csv(final_output, full_path)
  
  message(paste("Export complete. Found", nrow(final_output), "misclassified samples."))
  message(paste("Saved to:", full_path))
  
  return(final_output)
}

data <- read_csv("output/df_taxa_metadata_250.csv")

wrong_preds <- export_misclassified_cw_diet1(data, 
                                             remove_outliers = TRUE,
                                             out_dir = "output_bioproj_prim", 
                                             filename = "misclassified_samples_cw_diet1.csv")

## ------------------------FIG 2: PANEL ARRANGEMENT ------------------------- ##
## --------------------------- Predictor: CW_DIET1 -------------------------- ##

library(ggpubr)

combined_plots_cw_diet1 <- ggarrange(
  alpha_plots_cw_diet1_sig,
  beta_plots_cw_diet1,
  conf_mat_plot_cw_diet1,
  ncol = 1, 
  nrow = 3,
  labels = c("A", "B", "C"),
  heights = c(2, 1, 1) 
)

# Display the combined plot
print(combined_plots_cw_diet1)

# Save the plot in multiple formats
out_dir <- "output_bioproj_prim/_FIGURES"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("fig2_cw_diet1_panel.", ext)),
    plot = combined_plots_cw_diet1,
    width = 14,
    height = 18,
    units = "in",
    dpi = 300
  )
}

################################################################################
### --------------- COMPOSITE FIGURE 3: RELATIVE ABUNDANCE ----------------- ###
################################################################################

library(tidyverse)
library(RColorBrewer)

# Load data
#setwd("C:/Users/thongthum.t/OneDrive - University of Florida/Projects/Bat Microbiome/2_Analysis/variable_importance/250")
setwd("C:/Users/mmayy/OneDrive - University of Florida/Projects/Bat Microbiome/2_Analysis/variable_importance/250")
data <- read_csv("output/df_taxa_metadata_250.csv")

# ------------------- A. CLEANING AND FILTERING
# Extract abundance columns for checking
tax_cols_check <- grep("^d__", names(data), value = TRUE)
abundance_check <- data[, tax_cols_check]

# Filter 1: Remove exact duplicate rows (Technical error check)
is_not_dup <- !duplicated(abundance_check)
data <- data[is_not_dup, ]
print(paste("Removed", sum(!is_not_dup), "exact duplicate samples."))

# Re-extract abundance and metadata after duplicate removal
tax_cols <- grep("^d__", names(data), value = TRUE)
abundance <- data[, tax_cols]
sample_ids <- data$index
diet_type <- data$Diet_type
species_name <- data$taxonomic_species
animals <- data$Common.name_Species
bioproject <- data$bioproject

# Extract family names (Corrected for spaces)
family_names <- sapply(tax_cols, function(col) {
  parts <- unlist(strsplit(col, ";"))
  parts <- trimws(parts)
  f_part <- parts[grep("^f__", parts)]
  if (length(f_part) == 0) return("Unclassified")
  family <- sub("^f__", "", f_part)
  if (family == "") return("Unclassified")
  return(family)
})

# Aggregate to family level
abundance_t <- t(abundance)
abundance_family <- rowsum(abundance_t, group = family_names, reorder = FALSE)

# Filter 2: Remove Contaminants (Chloroplast/Mitochondria)
abundance_family <- abundance_family[rownames(abundance_family) != "Chloroplast", ]
abundance_family <- abundance_family[rownames(abundance_family) != "Mitochondria", ]
abundance_family <- t(abundance_family)
rownames(abundance_family) <- sample_ids

# Convert to dataframe
abundance_family_df <- as.data.frame(abundance_family)

# Remove "Unclassified" column if present
if ("Unclassified" %in% colnames(abundance_family_df)) {
  abundance_family_df$Unclassified <- NULL 
}

# Filter 3: CRITICAL READ DEPTH FILTER (Fixes the "Flat Block" issue)

# Calculate total reads per sample
row_sums <- rowSums(abundance_family_df)

# Print summary to see distribution before filtering
print("Read count summary before filtering:")
print(summary(row_sums))

min_reads <- 5000
valid_rows <- row_sums >= min_reads

# Print report
print(paste("Total samples:", length(row_sums)))
print(paste("Valid samples (>= ", min_reads, " reads):", sum(valid_rows)))
print(paste("Discarded samples (< ", min_reads, " reads):", sum(!valid_rows)))

if(sum(valid_rows) == 0) {
  stop("Error: All samples were filtered out! Try lowering 'min_reads' or check your data.")
}

# Apply the filter to Data and Metadata
abundance_family_df <- abundance_family_df[valid_rows, ]
diet_type <- diet_type[valid_rows]
species_name <- species_name[valid_rows]
animals <- animals[valid_rows]
sample_ids <- sample_ids[valid_rows] 
bioproject <- bioproject[valid_rows]

# ----------- B. NORMALIZATION AND PREPARATION

# Normalize to relative abundance (%)
abundance_family_rel <- as.data.frame(abundance_family_df / rowSums(abundance_family_df) * 100)

# Add metadata back to the relative abundance dataframe
abundance_family_rel$Sample <- sample_ids
abundance_family_rel$Diet_type <- diet_type
abundance_family_rel$Species <- species_name
abundance_family_rel$Animal <- animals
abundance_family_rel$BioProject <- bioproject

# Identify numeric columns (the bacteria)
bacteria_cols <- sapply(abundance_family_rel, is.numeric)

# Get top 20 families based on overall abundance
overall_top <- colSums(abundance_family_rel[, bacteria_cols])
top_families <- names(sort(overall_top, decreasing = TRUE)[1:20])
print(top_families)

# Calculate 'Other' category
top_data <- abundance_family_rel[, top_families, drop = FALSE]
other <- 100 - rowSums(top_data)
final_data <- cbind(top_data, Other = other)

# Add metadata to final object
final_data$Sample <- abundance_family_rel$Sample
final_data$Diet_type <- factor(abundance_family_rel$Diet_type,
                               levels = c("Fish", "Insect", "Fruit", "Blood"))
final_data$Species <- abundance_family_rel$Species
final_data$Animal <- abundance_family_rel$Animal
final_data$BioProject <- bioproject 

# ------------------- C. EXPORT FULL SUMMARY (ALL FAMILIES) -------------------
library(dplyr)
library(tidyr)
library(readr)

# 1. Identify all bacterial family columns in the full relative abundance dataset
# We exclude the metadata columns we added earlier
metadata_cols <- c("Sample", "Diet_type", "Species", "Animal", "BioProject")
all_bacteria_cols <- setdiff(names(abundance_family_rel), metadata_cols)

# 2. Group by Diet_type and calculate the mean for ALL families (not just top 20)
full_summary_table <- abundance_family_rel %>%
  select(Diet_type, all_of(all_bacteria_cols)) %>%
  group_by(Diet_type) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

# 3. Pivot to long format (Diet_type, Bacterial_Family, Average_Relative_Abundance)
full_long_summary <- full_summary_table %>%
  pivot_longer(
    cols = -Diet_type, 
    names_to = "Bacterial_Family", 
    values_to = "Average_Relative_Abundance"
  )

# 4. Filter for > 1% abundance
# This will capture Top 20 families AND any "Other" families that might exceed 1% 
# in specific diet groups even if they weren't in the global Top 20.
full_long_summary_filtered <- full_long_summary %>%
  filter(Average_Relative_Abundance > 1) %>%
  arrange(Diet_type, desc(Average_Relative_Abundance))

# 5. Export the table
write_csv(full_long_summary_filtered, "output_bioproj_prim/Average_Rel_Abundance_Families.csv")

# Print check
print(head(full_long_summary_filtered))

# ------------------------- D. PLOTTING SETUP ----------------------------------

# Pivot to long format

long_data <- final_data %>%
  pivot_longer(cols = -c(Sample, Diet_type, Species, Animal, BioProject),
               names_to = "Family", values_to = "Abundance") %>%
  mutate(Family = factor(Family, levels = c(top_families, "Other")))

# Define Palettes
bright_families <- c("Aeromonadaceae", "Bacillaceae", "Bacteroidaceae", 
                     "Burkholderiaceae", "Chlamydiaceae", "Clostridiaceae", 
                     "Enterobacteriaceae", "Enterococcaceae", 
                     "Fusobacteriaceae", "Hafniaceae", "Helicobacteraceae", 
                     "Moraxellaceae", "Morganellaceae", "Mycoplasmataceae", 
                     "Pasteurellaceae", "Peptostreptococcaceae", "Staphylococcaceae", 
                     "Streptococcaceae", "Yersiniaceae")

n_bright <- sum(top_families %in% bright_families)
n_subtle <- length(top_families) - n_bright
base_colors <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set1"))
bright_palette <- colorRampPalette(base_colors)(n_bright)
subtle_palette <- colorRampPalette(c("#CDFADB", "#F0C1E1", "#A5B68D", "#F4DEB3",
                                     "#C5D3E8", "#CBE2B5", "#E8C5E5", "#D6DAC8"))(n_subtle)

# Assign colors
family_colors <- c()
bright_idx <- 1
subtle_idx <- 1
for (fam in top_families) {
  if (fam %in% bright_families) {
    family_colors[fam] <- bright_palette[bright_idx]
    bright_idx <- bright_idx + 1
  } else {
    family_colors[fam] <- subtle_palette[subtle_idx]
    subtle_idx <- subtle_idx + 1
  }
}

# ----------------- color assignment
family_colors["Other"] <- "grey70"

# Genus: Escherichia-Shigella (#D62728) -> Family: Enterobacteriaceae
if("Enterobacteriaceae" %in% names(family_colors)) family_colors["Enterobacteriaceae"] <- "#D62728"
if("Clostridiaceae" %in% names(family_colors)) family_colors["Clostridiaceae"] <- "#FF7F0E"
if("Streptococcaceae" %in% names(family_colors)) family_colors["Streptococcaceae"] <- "#8760d6"
if("Enterococcaceae" %in% names(family_colors)) family_colors["Enterococcaceae"] <- "#17BECF"
if("Helicobacteraceae" %in% names(family_colors)) family_colors["Helicobacteraceae"] <- "#2b7ced"
if("Mycoplasmataceae" %in% names(family_colors)) family_colors["Mycoplasmataceae"] <- "#20a10a"
if("Hafniaceae" %in% names(family_colors)) family_colors["Hafniaceae"] <- "#8C564B"
if("Aeromonadaceae" %in% names(family_colors)) family_colors["Aeromonadaceae"] <- "#AEC7E8"
if("Morganellaceae" %in% names(family_colors)) family_colors["Morganellaceae"] <- "#FF9896"
if("Bacillaceae" %in% names(family_colors)) family_colors["Bacillaceae"] <- "#ebeb36"
if("Staphylococcaceae" %in% names(family_colors)) family_colors["Staphylococcaceae"] <- "#ad4b75"
if("Fusobacteriaceae" %in% names(family_colors)) family_colors["Fusobacteriaceae"] <- "#E377C2"
if("Pasteurellaceae" %in% names(family_colors)) family_colors["Pasteurellaceae"] <- "#f5bb1d"
if("Corynebacteriaceae" %in% names(family_colors)) family_colors["Corynebacteriaceae"] <- "#FF4500"
if("Lachnospiraceae" %in% names(family_colors)) family_colors["Lachnospiraceae"] <- "#32CD32"
if("Leuconostocaceae" %in% names(family_colors)) family_colors["Leuconostocaceae"] <- "#C71585"
if("Orbaceae" %in% names(family_colors)) family_colors["Orbaceae"] <- "#b6afed"
if("Yersiniaceae" %in% names(family_colors)) family_colors["Yersiniaceae"] <- "#acad0c"
if("Peptostreptococcaceae" %in% names(family_colors)) family_colors["Peptostreptococcaceae"] <- "#edd9af"
if("Bacteroidaceae" %in% names(family_colors)) family_colors["Bacteroidaceae"] <- "#24bfa3"
if("Burkholderiaceae" %in% names(family_colors)) family_colors["Burkholderiaceae"] <- "#9467bd"
if("Chlamydiaceae" %in% names(family_colors)) family_colors["Chlamydiaceae"] <- "#8c564b"
if("Moraxellaceae" %in% names(family_colors)) family_colors["Moraxellaceae"] <- "#e377c2" 

# ----------------- D. GENERATE PLOT

abundance_fam_all <- ggplot(long_data, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_col(width = 1) +
  facet_grid(. ~ Diet_type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = family_colors) +
  theme_bw(base_size = 12) +
  theme(
    plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(y = "Relative Abundance (%)", fill = "Family") +
  guides(fill = guide_legend(ncol = 1))

print(abundance_fam_all)

# Save
# (Make sure 'out_dir' is defined in your environment, or replace with a path string)
if (!exists("out_dir")) { out_dir <- "output" } 
if (!dir.exists(out_dir)) { dir.create(out_dir) }

out_dir <- "output_bioproj_prim"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("Relative_abundance_family_all.", ext)),
    plot = abundance_fam_all,
    width = 20,
    height = 10,
    units = "in",
    dpi = 300
  )
}

## ---------------- SUPPLEMENTAL FIG: AVERAGE ABUNDANCE (FAMILY) ---------------

long_data_summary <- long_data %>%
  group_by(Diet_type, BioProject, Family) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Create plot
abundance_fam_avg <- ggplot(long_data_summary, aes(x = BioProject, y = Abundance, fill = Family)) +
  geom_col(width = 0.9) + # Slight gap between bars looks cleaner
  # Facet by Diet_type (the "Big Column")
  facet_grid(. ~ Diet_type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = family_colors) +
  theme_bw(base_size = 12) +
  theme(
    plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Turned axis text back on so you can see the BioProject IDs
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    strip.background = element_rect(fill = "grey90", color = NA), # Optional: makes headers distinct
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(y = "Mean Relative Abundance (%)", fill = "Family") +
  guides(fill = guide_legend(ncol = 1))

print(abundance_fam_avg)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("Relative_abundance_family_avg.", ext)),
    plot = abundance_fam_avg,
    width = 10, 
    height = 8, 
    units = "in", 
    dpi = 300
  )
}

## ========================================================================== ##
## ---------------- SUPPLEMENTAL: REL ABUNDANCE GENUS ALL ------------------- ##
## ========================================================================== ##

#setwd("C:/Users/thongthum.t/OneDrive - University of Florida/Projects/Bat Microbiome/2_Analysis/variable_importance/250")
setwd("C:/Users/mmayy/OneDrive - University of Florida/Projects/Bat Microbiome/2_Analysis/variable_importance/250")
data <- read_csv("output/df_taxa_metadata_250.csv")

library(tidyverse)
library(RColorBrewer)

#------------- 1. SETUP & DATA LOADING
if (!exists("out_dir")) { out_dir <- "output" } 
if (!dir.exists(out_dir)) { dir.create(out_dir) }

# Extract abundance columns for checking
tax_cols_check <- grep("^d__", names(data), value = TRUE)
abundance_check <- data[, tax_cols_check]

# Filter 1: Remove exact duplicate rows (Technical error check)
is_not_dup <- !duplicated(abundance_check)
data <- data[is_not_dup, ]
print(paste("Removed", sum(!is_not_dup), "exact duplicate samples."))

# Re-extract abundance and metadata
tax_cols <- grep("^d__", names(data), value = TRUE)
abundance <- data[, tax_cols]
sample_ids <- data$index
diet_type <- data$Diet_type
species_name <- data$taxonomic_species
animals <- data$Common.name_Species
bioproject <- data$bioproject

## -------------- 2. EXTRACT GENUS NAMES (User Provided Logic)

genus_names <- sapply(tax_cols, function(col) {
  parts <- unlist(strsplit(col, ";"))
  parts <- trimws(parts)
  g_part <- parts[grep("^g__", parts)]
  if (length(g_part) == 0) return("Unclassified")
  genus <- sub("^g__", "", g_part)
  if (genus == "" || genus == "uncultured") return("Unclassified")
  return(genus)
})

# Check the extraction
print(head(table(genus_names)))

## -------------------- 3. AGGREGATION & FILTERING

# Aggregate to GENUS level
abundance_t <- t(abundance)
abundance_genus <- rowsum(abundance_t, group = genus_names, reorder = FALSE)

# Filter 2: Remove Contaminants (Chloroplast/Mitochondria)

abundance_genus <- abundance_genus[!rownames(abundance_genus) %in% c("Chloroplast", "Mitochondria"), ]

abundance_genus <- t(abundance_genus)
rownames(abundance_genus) <- sample_ids

abundance_genus_df <- as.data.frame(abundance_genus)

if ("Unclassified" %in% colnames(abundance_genus_df)) {
  abundance_genus_df$Unclassified <- NULL 
}

##--- Filter 3: CRITICAL READ DEPTH FILTER
row_sums <- rowSums(abundance_genus_df)

print("Read count summary before filtering:")
print(summary(row_sums))

# Set threshold (Keep consistent with Family level analysis)
min_reads <- 5000
valid_rows <- row_sums >= min_reads

print(paste("Total samples:", length(row_sums)))
print(paste("Valid samples (>= ", min_reads, " reads):", sum(valid_rows)))

if(sum(valid_rows) == 0) {
  stop("Error: All samples were filtered out! Try lowering 'min_reads'.")
}

# Apply the filter
abundance_genus_df <- abundance_genus_df[valid_rows, ]
diet_type <- diet_type[valid_rows]
species_name <- species_name[valid_rows]
animals <- animals[valid_rows]
sample_ids <- sample_ids[valid_rows]
bioproject <- bioproject[valid_rows]

## ----------------- 4. NORMALIZATION & PREPARATION

# Normalize to relative abundance (%)
abundance_genus_rel <- as.data.frame(abundance_genus_df / rowSums(abundance_genus_df) * 100)

# Identify numeric columns (the bacteria)
bacteria_cols <- sapply(abundance_genus_rel, is.numeric)

# Get top 20 Genera
overall_top <- colSums(abundance_genus_rel[, bacteria_cols])
top_genera <- names(sort(overall_top, decreasing = TRUE)[1:20])

print("Top 20 Genera found:")
print(top_genera)

# Calculate 'Other' category
top_data <- abundance_genus_rel[, top_genera, drop = FALSE]
other <- 100 - rowSums(top_data)
final_data <- cbind(top_data, Other = other)

# Add metadata
final_data$Sample <- sample_ids
final_data$Diet_type <- factor(diet_type, levels = c("Fish", "Insect", "Fruit", "Blood"))
final_data$Species <- species_name
final_data$Animal <- animals
final_data$BioProject <- bioproject

# Pivot to long format
long_data <- final_data %>%
  pivot_longer(cols = -c(Sample, Diet_type, Species, Animal, BioProject),
               names_to = "Genus", values_to = "Abundance") %>%
  mutate(Genus = factor(Genus, levels = c(top_genera, "Other")))

# ------------------- 5. EXPORT GENUS SUMMARY TABLE ----------------------------
library(dplyr)
library(tidyr)
library(readr)

# 1. Prepare the full relative abundance dataframe with metadata
# (We need to attach Diet_type to abundance_genus_rel for grouping)
abundance_genus_rel_meta <- abundance_genus_rel
abundance_genus_rel_meta$Diet_type <- factor(diet_type, levels = c("Fish", "Insect", "Fruit", "Blood"))

# 2. Group by Diet_type and calculate the MEAN relative abundance for EVERY genus
summary_genus <- abundance_genus_rel_meta %>%
  group_by(Diet_type) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()

# 3. Pivot to long format (Diet_type, Genus, Average_Abundance)
long_summary_genus <- summary_genus %>%
  pivot_longer(
    cols = -Diet_type, 
    names_to = "Bacterial_Genus", 
    values_to = "Average_Relative_Abundance"
  )

# 4. Filter for genera with > 1% abundance
# This keeps the "Top 20" plus any "Other" genera that are locally abundant >1%
long_summary_genus_filtered <- long_summary_genus %>%
  filter(Average_Relative_Abundance > 1) %>%
  arrange(Diet_type, desc(Average_Relative_Abundance))

# 5. Export the table
write_csv(long_summary_genus_filtered, "output_bioproj_prim/Average_Rel_Abundance_Genera.csv")

print(head(long_summary_genus_filtered))

## --------------------------- 6. PLOTTING ---------------------------------- ##

# Dynamic Color Palette (Genera names change, so we generate colors automatically)
# We use a large palette to handle 21 distinct colors
n_taxa <- length(top_genera) + 1 # +1 for Other
# Create a palette that interpolates distinct colors
genus_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_taxa)
names(genus_colors) <- c(top_genera, "Other")
genus_colors["Other"] <- "grey70"

# --- Genus Color Assignments based on Family ---
# Family: Aeromonadaceae (#AEC7E8)
if("Aeromonas" %in% names(genus_colors)) genus_colors["Aeromonas"] <- "#AEC7E8"

# Family: Bacillaceae (#ebeb36)
if("Bacillus" %in% names(genus_colors)) genus_colors["Bacillus"] <- "#ebeb36"

# Family: Bacteroidaceae (#24bfa3)
if("Bacteroides" %in% names(genus_colors)) genus_colors["Bacteroides"] <- "#24bfa3"

# Family: Burkholderiaceae (#9467bd)
if("Burkholderia-Caballeronia-Paraburkholderia" %in% names(genus_colors)) genus_colors["Burkholderia-Caballeronia-Paraburkholderia"] <- "#d94ce6"
if("Cupriavidus" %in% names(genus_colors)) genus_colors["Cupriavidus"] <- "#8a1184"

# Family: Clostridiaceae (#FF7F0E)
if("Clostridium_sensu_stricto_1"            %in% names(genus_colors)) genus_colors["Clostridium_sensu_stricto_1"] <- "#FF7F0E"

# Family: Corynebacteriaceae (#FF4500)
if("Corynebacterium"                        %in% names(genus_colors)) genus_colors["Corynebacterium"] <- "#FF4500"

# Family: Enterobacteriaceae (#D62728)
if("Escherichia-Shigella"                   %in% names(genus_colors)) genus_colors["Escherichia-Shigella"] <- "#D62728"

# Family: Enterococcaceae (#17BECF)
if("Enterococcus"                           %in% names(genus_colors)) genus_colors["Enterococcus"] <- "#54ecf7"

# Family: Fusobacteriaceae (#E377C2)
if("Cetobacterium"                          %in% names(genus_colors)) genus_colors["Cetobacterium"] <- "#ff63e2"

# Family: Hafniaceae (#8C564B)
if("Hafnia-Obesumbacterium"                 %in% names(genus_colors)) genus_colors["Hafnia-Obesumbacterium"] <- "#8C564B"
if("Edwardsiella"                           %in% names(genus_colors)) genus_colors["Edwardsiella"] <- "#e39b8d"

# Family: Helicobacteraceae (#2b7ced)
if("Helicobacter"                           %in% names(genus_colors)) genus_colors["Helicobacter"] <- "#2b7ced"

# Family: Morganellaceae (#FF9896)
if("Morganella"                             %in% names(genus_colors)) genus_colors["Morganella"] <- "#FF9896"
if("Proteus"                                %in% names(genus_colors)) genus_colors["Proteus"] <- "#FF9896"

# Family: Mycoplasmataceae (#20a10a)
if("Mycoplasma"                            %in% names(genus_colors)) genus_colors["Mycoplasma"] <- "#20a10a"

# Family: Pasteurellaceae (#f5bb1d)
if("Actinobacillus"                        %in% names(genus_colors)) genus_colors["Actinobacillus"] <- "#f5bb1d"

# Family: Staphylococcaceae (#ad4b75)
if("Staphylococcus"                        %in% names(genus_colors)) genus_colors["Staphylococcus"] <- "#ad4b75"

# Family: Streptococcaceae (#8760d6)
if("Lactococcus"                           %in% names(genus_colors)) genus_colors["Lactococcus"] <- "#8760d6"
if("Streptococcus"                         %in% names(genus_colors)) genus_colors["Streptococcus"] <- "#4d1991"


abundance_genus_plot <- ggplot(long_data, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(width = 1) +
  facet_grid(. ~ Diet_type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 12) +
  theme(
    plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(), # Hide X labels to avoid clutter
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "italic"), # Italic for Genus
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(y = "Relative Abundance (%)", fill = "Genus") +
  guides(fill = guide_legend(ncol = 1))

print(abundance_genus_plot)

# Save
for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("Relative_abundance_genus_all.", ext)),
    plot = abundance_genus_plot,
    width = 20,
    height = 10,
    units = "in",
    dpi = 300
  )
}

# ============================================================================ #
# ----------------------- 3.1 ) AVERAGE ABUNDANCE (GENUS) -------------------- #
# ============================================================================ #

# Summarize data by Diet_type AND BioProject
long_data_summary <- long_data %>%
  group_by(Diet_type, BioProject, Genus) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Create plot
abundance_genus_avg <- ggplot(long_data_summary, aes(x = BioProject, y = Abundance, fill = Genus)) +
  geom_col(width = 0.9) + # Slight gap between bars looks cleaner
  facet_grid(. ~ Diet_type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 14) +
  theme(
    plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Turn x-axis text back on so you can see the BioProject IDs
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(y = "Mean Relative Abundance (%)", fill = "Genus") +
  guides(fill = guide_legend(ncol = 1))

print(abundance_genus_avg)

out_dir = "output_bioproj_prim"
#out_dir = "output_prim_bioproj"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("Relative_abundance_genus_avg.", ext)),
    plot = abundance_genus_avg,
    width = 12, 
    height = 8, 
    units = "in", 
    dpi = 300
  )
}

# ============================================================================ #
# ----------- 3.2) Variable Importance CLR Abundance (GENUS level) ----------- #
## ------------------------ Predictor: Diet_type ---------------------------- ##
# ============================================================================ #

data <- df

# Prepare metadata with CLEAN diet labels (no make.names)
metadata_df <- data.frame(
  SampleID = rownames(df),
  Diet_type = factor(df$diet_1),
  stringsAsFactors = FALSE
)

# Extract abundance data (taxonomy columns)
abundance_table <- data[, grep("d__Bacteria", colnames(data))]
rownames(abundance_table) <- metadata_df$SampleID

# genus extraction
extract_genus <- function(t) {
  if (!grepl("g__", t)) return(NA_character_)
  g <- sub(".*g__([^;]+).*", "\\1", t)
  if (g == "" || is.na(g)) return(NA_character_) else return(g)
}
genus_names <- sapply(colnames(abundance_table), extract_genus)
keep <- !is.na(genus_names)

# Aggregate to GENUS level (drop NA genera)
abundance_table_genus <- rowsum(t(abundance_table[, keep, drop = FALSE]),
                                group = genus_names[keep], na.rm = TRUE)
abundance_table_genus <- t(abundance_table_genus)

# EXCLUDE organelle genera from calculations
abundance_table_genus <- abundance_table_genus[, !(colnames(abundance_table_genus) %in% c("Chloroplast", "Mitochondria")), drop = FALSE]

# Sanity checks
print(table(metadata_df$Diet_type))

# CLR transformation
clr_transform <- function(x) {
  x <- x + 1
  log(x / exp(mean(log(x))))
}
abundance_table_clr <- t(apply(abundance_table_genus, 1, clr_transform))

# Combine with metadata
rf_data <- cbind(Diet_type = metadata_df$Diet_type, as.data.frame(abundance_table_clr))
rf_data$Diet_type <- as.factor(rf_data$Diet_type)

# Ensure unique predictor names
colnames(rf_data) <- make.unique(colnames(rf_data))

# RF Classification
set.seed(1000)
train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  savePredictions = "final"
)
rf_model <- train(
  Diet_type ~ .,
  data = rf_data,
  method = "rf",
  trControl = train_control,
  ntree = 500,
  importance = TRUE,
  metric = "Accuracy"
)
print(rf_model)

# Variable Importance (one-vs-rest)
diet_classes <- levels(rf_data$Diet_type)
importance_list <- list()

for (diet in diet_classes) {
  rf_data_binary <- rf_data
  rf_data_binary$Diet_binary <- ifelse(rf_data_binary$Diet_type == diet, diet, paste0("Not_", diet))
  rf_data_binary$Diet_binary <- factor(rf_data_binary$Diet_binary)
  x <- rf_data_binary %>% select(-Diet_type, -Diet_binary)
  y <- rf_data_binary$Diet_binary
  rf_bin <- randomForest(x = x, y = y, ntree = 500, importance = TRUE)
  imp <- importance(rf_bin, type = 1, scale = TRUE)
  imp_df <- data.frame(
    Taxa = rownames(imp),
    Importance = imp[, 1],
    Diet = diet
  )
  importance_list[[diet]] <- imp_df
}
importance_long <- bind_rows(importance_list)

# Identify Diet with Highest CLR Abundance per genus
clr_df <- as.data.frame(abundance_table_clr)
clr_df$Diet_type <- rf_data$Diet_type
clr_long <- clr_df %>%
  pivot_longer(cols = -Diet_type, names_to = "Taxa", values_to = "CLR")
mean_abund <- clr_long %>%
  group_by(Taxa, Diet_type) %>%
  summarise(MeanCLR = mean(CLR), .groups = "drop")
winner_diet <- mean_abund %>%
  group_by(Taxa) %>%
  slice_max(order_by = MeanCLR, n = 1) %>%
  select(Taxa, Winner = Diet_type)

# Merge with importance data
importance_long <- importance_long %>%
  left_join(winner_diet, by = "Taxa") %>%
  mutate(
    Diet = case_when(
      Diet == "Fishes" ~ "Fish",
      Diet == "Fruits" ~ "Fruit",
      Diet == "Insects" ~ "Insect",
      Diet == "Blood" ~ "Blood",
      TRUE ~ Diet
    ),
    Winner = case_when(
      Winner == "Fishes" ~ "Fish",
      Winner == "Fruits" ~ "Fruit",
      Winner == "Insects" ~ "Insect",
      Winner == "Blood" ~ "Blood",
      TRUE ~ Winner
    )
  )

# Get Top 20 genera and Set Order
top_taxa <- importance_long %>%
  group_by(Taxa) %>%
  summarise(MaxImportance = max(Importance), .groups = "drop") %>%
  slice_max(order_by = MaxImportance, n = 20, with_ties = FALSE) %>%
  pull(Taxa)

taxa_order <- importance_long %>%
  filter(Taxa %in% top_taxa) %>%
  group_by(Taxa) %>%
  summarise(MaxImportance = max(Importance), .groups = "drop") %>%
  arrange(desc(MaxImportance)) %>%
  pull(Taxa)

# Convert to factors with the correct order
importance_long$Taxa <- factor(importance_long$Taxa, levels = taxa_order)

# Prepare CLR data for the distribution plot
clr_long_top <- clr_long %>%
  filter(Taxa %in% top_taxa) %>%
  mutate(
    Taxa = factor(Taxa, levels = taxa_order),
    Diet_type = case_when(
      Diet_type == "Fishes" ~ "Fish",
      Diet_type == "Fruits" ~ "Fruit",
      Diet_type == "Insects" ~ "Insect",
      Diet_type == "Blood" ~ "Blood",
      TRUE ~ Diet_type
    ),
    CLR = as.numeric(CLR)
  )

# Create Background Stripes Data
importance_long_char <- importance_long %>%
  filter(Taxa %in% top_taxa) %>%
  mutate(Taxa_char = as.character(Taxa))
clr_long_top_char <- clr_long_top %>%
  mutate(Taxa_char = as.character(Taxa))

taxa_order_bg <- clr_long_top_char %>%
  group_by(Taxa_char) %>%
  summarise(mean_clr = mean(CLR, na.rm = TRUE)) %>%
  arrange(mean_clr) %>%
  pull(Taxa_char)

bg_stripes <- data.frame(Taxa_char = taxa_order_bg) %>%
  mutate(stripe = (row_number() %% 2 == 1))

importance_long_char$Taxa_char <- factor(importance_long_char$Taxa_char, levels = taxa_order_bg)
clr_long_top_char$Taxa_char <- factor(clr_long_top_char$Taxa_char, levels = taxa_order_bg)

# ---- Plot 1: Variable importance (Genus) ----
custom_plot1 <- ggplot() +
  geom_tile(
    data = bg_stripes %>% filter(stripe == TRUE),
    aes(y = Taxa_char, x = 0, height = 0.9, width = Inf),
    fill = "grey90", alpha = 0.6
  ) +
  geom_col(
    data = importance_long_char,
    aes(x = Importance, y = Taxa_char, fill = Diet),
    position = position_dodge(0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "Insect" = "#ff5e2b",
      "Blood"  = "#ba0970",
      "Fruit"  = "#2E8B57",
      "Fish"   = "#303F9F")
    
  ) +
  labs(x = "Mean Decrease Accuracy", y = "Genus") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "italic", size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

# ---- Plot 2: CLR abundance (Genus) ----
custom_plot2 <- ggplot() +
  geom_tile(
    data = bg_stripes %>% filter(stripe == TRUE),
    aes(y = Taxa_char, x = 0, height = 0.9, width = Inf),
    fill = "grey90", alpha = 0.6
  ) +
  geom_boxplot(
    data = clr_long_top_char,
    aes(x = CLR, y = Taxa_char, fill = Diet_type),
    outlier.size = 0.8,
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "Insect" = "#ff5e2b",
      "Blood"  = "#ba0970",
      "Fruit"  = "#2E8B57",
      "Fish"   = "#303F9F")
    
  ) +
  labs(x = "CLR Abundance", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

# Legend
legend_plot <- ggplot(importance_long_char,
                      aes(x = Importance, y = Taxa_char, fill = Diet)) +
  geom_col(position = position_dodge(0.8)) +
  scale_fill_manual(
    values = c(
      "Insect" = "#ff5e2b",
      "Blood"  = "#ba0970",
      "Fruit"  = "#2E8B57",
      "Fish"   = "#303F9F")
    ,
    name = "Diet Group",
    breaks = c("Insect", "Fruit", "Fish", "Blood"),
    labels = c("Insect", "Fruit", "Fish", "Blood")
  ) +
  theme(legend.position = "bottom",     
        legend.text = element_text(size = 12),              # Size of "Insect", "Fruit", etc.
        legend.title = element_text(size = 12, face = "bold"), # Size of "Diet Group"
        legend.key.size = unit(0.5, "cm")                     # Optional: Makes the color boxes bigger
  )

plot_legend <- cowplot::get_legend(legend_plot)

# Combine plots with legend
CLR_combined_genus <- (custom_plot1 + custom_plot2) /
  plot_legend +
  plot_layout(heights = c(10, 1))

print(CLR_combined_genus)

out_dir = "output_bioproj_prim"
#out_dir = "output_prim_bioproj"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("CLR_genus.", ext)),
    plot = CLR_combined_genus,
    width = 16, 
    height = 10, 
    units = "in", 
    dpi = 300
  )
}

# ============================================================================ #
## ---------------------------- ternary plot -------------------------------- ##
## ------------------------ Predictor: Diet_type ---------------------------- ##
# ============================================================================ #

library(ggplot2)
library(ggtern)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(patchwork) # For combining plots

# 0. Define Diet Colors
diet_colors <- c(
  "Insect" = "#ff5e2b",
  "Blood"  = "#ba0970",
  "Fruit"  = "#2E8B57",
  "Fish"   = "#303F9F"
)

# 1. Prepare Data (df_long)
df_long <- df %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(
    cols = starts_with("d__Bacteria"),
    names_to = "taxon",
    values_to = "abundance"
  ) %>%
  filter(abundance > 0) %>%
  mutate(
    genus   = stringr::str_extract(taxon, "g__[^;|]+") %>% sub("^g__", "", .),
    species = stringr::str_extract(taxon, "s__[^;|]+") %>% sub("^s__", "", .)
  ) %>%
  filter(!is.na(genus), genus != "") %>%
  filter(!genus %in% c("Chloroplast", "Mitochondria")) %>%
  filter(is.na(species) | !grepl("^uncultured", species, ignore.case = TRUE)) %>%
  rename(taxonomic_group = Diet_type) %>%
  select(sample_id, genus, species, abundance, taxonomic_group)

# 2. Aggregate Data into Two Sets
genus_summary_fish <- df_long %>%
  filter(taxonomic_group %in% c("Insect", "Fruit", "Fish")) %>%
  group_by(genus, taxonomic_group) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = taxonomic_group, values_from = total_abundance, values_fill = 0) %>%
  mutate(total = Insect + Fruit + Fish, Fish_prop = Fish/total, Insect_prop = Insect/total, Fruit_prop = Fruit/total)

genus_summary_blood <- df_long %>%
  filter(taxonomic_group %in% c("Insect", "Fruit", "Blood")) %>%
  group_by(genus, taxonomic_group) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = taxonomic_group, values_from = total_abundance, values_fill = 0) %>%
  mutate(total = Insect + Fruit + Blood, Blood_prop = Blood/total, Insect_prop = Insect/total, Fruit_prop = Fruit/total)

# 3. Identify Top Genera across BOTH sets for a consistent legend
all_genera <- bind_rows(
  select(genus_summary_fish, genus, total),
  select(genus_summary_blood, genus, total)
) %>%
  group_by(genus) %>%
  summarise(grand_total = sum(total)) %>%
  arrange(desc(grand_total))

top_genera_names <- all_genera %>%
  slice_head(n = 20) %>%
  pull(genus)

# 4. Filter data to only top genera
plot_data_fish <- genus_summary_fish %>% filter(genus %in% top_genera_names)
plot_data_blood <- genus_summary_blood %>% filter(genus %in% top_genera_names)

# 5. Generate Genus Color Palette
n_taxa <- length(top_genera_names)
genus_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_taxa)
names(genus_colors) <- top_genera_names

# --- Genus Color Assignments based on Family ---
# Family: Aeromonadaceae (#AEC7E8)
if("Aeromonas" %in% names(genus_colors)) genus_colors["Aeromonas"] <- "#AEC7E8"

# Family: Bacillaceae (#ebeb36)
if("Bacillus" %in% names(genus_colors)) genus_colors["Bacillus"] <- "#ebeb36"

# Family: Bacteroidaceae (#24bfa3)
if("Bacteroides" %in% names(genus_colors)) genus_colors["Bacteroides"] <- "#24bfa3"

# Family: Burkholderiaceae (#9467bd)
if("Burkholderia-Caballeronia-Paraburkholderia" %in% names(genus_colors)) genus_colors["Burkholderia-Caballeronia-Paraburkholderia"] <- "#d94ce6"
if("Cupriavidus" %in% names(genus_colors)) genus_colors["Cupriavidus"] <- "#8a1184"

# Family: Clostridiaceae (#FF7F0E)
if("Clostridium_sensu_stricto_1"            %in% names(genus_colors)) genus_colors["Clostridium_sensu_stricto_1"] <- "#FF7F0E"

# Family: Corynebacteriaceae (#FF4500)
if("Corynebacterium"                        %in% names(genus_colors)) genus_colors["Corynebacterium"] <- "#FF4500"

# Family: Enterobacteriaceae (#D62728)
if("Escherichia-Shigella"                   %in% names(genus_colors)) genus_colors["Escherichia-Shigella"] <- "#D62728"

# Family: Enterococcaceae (#17BECF)
if("Enterococcus"                           %in% names(genus_colors)) genus_colors["Enterococcus"] <- "#54ecf7"

# Family: Fusobacteriaceae (#E377C2)
if("Cetobacterium"                          %in% names(genus_colors)) genus_colors["Cetobacterium"] <- "#ff63e2"

# Family: Hafniaceae (#8C564B)
if("Hafnia-Obesumbacterium"                 %in% names(genus_colors)) genus_colors["Hafnia-Obesumbacterium"] <- "#8C564B"
if("Edwardsiella"                           %in% names(genus_colors)) genus_colors["Edwardsiella"] <- "#e39b8d"

# Family: Helicobacteraceae (#2b7ced)
if("Helicobacter"                           %in% names(genus_colors)) genus_colors["Helicobacter"] <- "#2b7ced"

# Family: Morganellaceae (#FF9896)
if("Morganella"                             %in% names(genus_colors)) genus_colors["Morganella"] <- "#FF9896"
if("Proteus"                                %in% names(genus_colors)) genus_colors["Proteus"] <- "#FF9896"

# Family: Mycoplasmataceae (#20a10a)
if("Mycoplasma"                            %in% names(genus_colors)) genus_colors["Mycoplasma"] <- "#20a10a"

# Family: Pasteurellaceae (#f5bb1d)
if("Actinobacillus"                        %in% names(genus_colors)) genus_colors["Actinobacillus"] <- "#f5bb1d"

# Family: Staphylococcaceae (#ad4b75)
if("Staphylococcus"                        %in% names(genus_colors)) genus_colors["Staphylococcus"] <- "#ad4b75"

# Family: Streptococcaceae (#8760d6)
if("Lactococcus"                           %in% names(genus_colors)) genus_colors["Lactococcus"] <- "#8760d6"
if("Streptococcus"                         %in% names(genus_colors)) genus_colors["Streptococcus"] <- "#4d1991"

# 6. Create Plot 1: Fish
ternary_fish <- ggtern(plot_data_fish,
                       aes(x = Insect_prop, y = Fruit_prop, z = Fish_prop, fill = genus)) +
  geom_point(size = 6, alpha = 0.8, shape = 21, color = "black") +
  labs(x = "Insect", y = "Fruit", z = "Fish", fill = "Genus") +
  theme_rgbw() +
  theme(
    axis.title = element_blank(), 
    tern.axis.arrow.show = FALSE,
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    tern.axis.title.L = element_text(color = diet_colors["Insect"], face="bold"),
    tern.axis.title.T = element_text(color = diet_colors["Fruit"], face="bold"),
    tern.axis.title.R = element_text(color = diet_colors["Fish"], face="bold")
  ) +
  scale_fill_manual(values = genus_colors, name = "Genus")

print(ternary_fish)

# 7. Create Plot 2: Blood
ternary_blood <- ggtern(plot_data_blood,
                        aes(x = Insect_prop, y = Fruit_prop, z = Blood_prop, fill = genus)) +
  geom_point(size = 6, alpha = 0.8, shape = 21, color = "black") +
  labs(x = "Insect", y = "Fruit", z = "Blood", fill = "Genus") +
  theme_rgbw() +
  theme(
    axis.title = element_blank(), 
    tern.axis.arrow.show = FALSE,
    legend.position = "none", # Remove legend from this plot
    tern.axis.title.L = element_text(color = diet_colors["Insect"], face="bold"),
    tern.axis.title.T = element_text(color = diet_colors["Fruit"], face="bold"),
    tern.axis.title.R = element_text(color = diet_colors["Blood"], face="bold")
  ) +
  scale_fill_manual(values = genus_colors, name = "Genus")

print(ternary_blood)

# 8. Combine plots with a shared legend using patchwork
combined_ternary <- ternary_fish + ternary_blood + plot_layout(guides = 'collect')

print(combined_ternary)

# 9. Save
out_dir = "output_bioproj_prim"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("ternary_combined_plots.", ext)),
    plot = combined_ternary,
    width = 14,
    height = 6,
    units = "in",
    dpi = 300
  )
}

# -------------------------FIGURE 3 PANEL ARRANGEMENT ------------------------ #

library(ggpubr)

# 2. Stack everything vertically (A, then B, then the bottom_row)
abun_clr_ternary_panel <- ggarrange(
  abundance_genus_avg, 
  CLR_combined_genus, 
  combined_ternary,
  labels = c("A", "B", "C"),
  ncol = 1,
  nrow = 3
)

print(abun_clr_ternary_panel)

out_dir = "output_bioproj_prim/_FIGURES"
#out_dir = "output_prim_bioproj"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("fig3_abun_clr_ternary_genus.", ext)),
    plot = abun_clr_ternary_panel,
    width = 14, 
    height = 18, 
    units = "in", 
    dpi = 300
  )
}

################################################################################
## ------------- SUPPLEMENTAL: Variable Importance (FAMILY level) ------------ ##

# Packages
library(tidyverse)
library(caret)
library(randomForest)
library(patchwork)
library(cowplot)

# Use df already in memory
data <- df

# Prepare metadata with CLEAN diet labels (no make.names)
metadata_df <- data.frame(
  SampleID  = rownames(data),
  Diet_type = factor(data$diet_1, levels = c("Insect","Fruit","Fish", "Blood")),
  stringsAsFactors = FALSE
)

# Sanity check (Fruits should be present)
print(table(metadata_df$Diet_type))

# Extract abundance data (taxonomy columns)
abundance_table <- data[, grep("d__Bacteria", colnames(data)), drop = FALSE]
rownames(abundance_table) <- metadata_df$SampleID

# Family extraction
extract_family <- function(t) {
  if (!grepl("f__", t)) return(NA_character_)
  f <- sub(".*f__([^;]+).*", "\\1", t)
  if (f == "" || is.na(f)) NA_character_ else f
}
family_names <- sapply(colnames(abundance_table), extract_family)
keep <- !is.na(family_names)

# Aggregate to FAMILY level (drop NA families)
abundance_table_family <- rowsum(t(abundance_table[, keep, drop = FALSE]),
                                 group = family_names[keep], na.rm = TRUE)
abundance_table_family <- t(abundance_table_family)

# EXCLUDE organelle families from calculations
abundance_table_family <- abundance_table_family[, !(colnames(abundance_table_family) %in% c("Chloroplast", "Mitochondria")), drop = FALSE]

# Sanity check
print(table(metadata_df$Diet_type))

# CLR transformation
clr_transform <- function(x) { x <- x + 1; log(x / exp(mean(log(x)))) }
abundance_table_clr <- t(apply(abundance_table_family, 1, clr_transform))

# Combine with metadata
rf_data <- cbind(Diet_type = metadata_df$Diet_type, as.data.frame(abundance_table_clr))
rf_data$Diet_type <- as.factor(rf_data$Diet_type)

# Ensure unique predictor names
colnames(rf_data) <- make.unique(colnames(rf_data))

# Sanity: class counts (Fruits should appear)
print(table(rf_data$Diet_type))

# RF Classification (optional overall multiclass summary)
set.seed(1000)
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, savePredictions = "final")
rf_model <- train(Diet_type ~ ., data = rf_data, method = "rf",
                  trControl = train_control, ntree = 500, importance = TRUE, metric = "Accuracy")
print(rf_model)

# Variable Importance (one-vs-rest)
diet_classes <- levels(rf_data$Diet_type)
importance_list <- list()
for (diet in diet_classes) {
  rf_data_binary <- rf_data
  rf_data_binary$Diet_binary <- factor(ifelse(rf_data_binary$Diet_type == diet, diet, paste0("Not_", diet)))
  x <- rf_data_binary %>% select(-Diet_type, -Diet_binary)
  y <- rf_data_binary$Diet_binary
  
  rf_bin <- randomForest(x = x, y = y, ntree = 500, importance = TRUE)
  imp <- importance(rf_bin, type = 1, scale = TRUE)
  
  imp_df <- data.frame(Taxa = rownames(imp), Importance = imp[, 1], Diet = diet)
  importance_list[[diet]] <- imp_df
}
importance_long <- bind_rows(importance_list)

# Identify Diet with Highest CLR Abundance per family
clr_df <- as.data.frame(abundance_table_clr)
clr_df$Diet_type <- rf_data$Diet_type
clr_long <- clr_df %>% pivot_longer(cols = -Diet_type, names_to = "Taxa", values_to = "CLR")

mean_abund <- clr_long %>%
  group_by(Taxa, Diet_type) %>%
  summarise(MeanCLR = mean(CLR), .groups = "drop")

winner_diet <- mean_abund %>%
  group_by(Taxa) %>%
  slice_max(order_by = MeanCLR, n = 1) %>%
  select(Taxa, Winner = Diet_type)

# Merge with importance data
importance_long <- importance_long %>%
  left_join(winner_diet, by = "Taxa") %>%
  mutate(
    Diet = case_when(
      Diet == "Fishes" ~ "Fish",
      Diet == "Fruits" ~ "Fruit",
      Diet == "Insects" ~ "Insect",
      Diet == "Blood" ~ "Blood",
      TRUE ~ Diet
    ),
    Winner = case_when(
      Winner == "Fishes" ~ "Fish",
      Winner == "Fruits" ~ "Frui",
      Winner == "Insects" ~ "Insect",
      Winner == "Blood" ~ "Blood",
      TRUE ~ Winner
    )
  )

# Get Top 20 families and Set Order
top_taxa <- importance_long %>%
  group_by(Taxa) %>%
  summarise(MaxImportance = max(Importance), .groups = "drop") %>%
  slice_max(order_by = MaxImportance, n = 20, with_ties = FALSE) %>%
  pull(Taxa)

taxa_order <- importance_long %>%
  filter(Taxa %in% top_taxa) %>%
  group_by(Taxa) %>%
  summarise(MaxImportance = max(Importance), .groups = "drop") %>%
  arrange(desc(MaxImportance)) %>%
  pull(Taxa)

# Convert to factors with the correct order
importance_long$Taxa <- factor(importance_long$Taxa, levels = taxa_order)

# Prepare CLR data for the distribution plot
clr_long_top <- clr_long %>%
  filter(Taxa %in% top_taxa) %>%
  mutate(
    Taxa = factor(Taxa, levels = taxa_order),
    Diet_type = case_when(
      Diet_type == "Fishes" ~ "Fish",
      Diet_type == "Fruits" ~ "Fruit",
      Diet_type == "Insects" ~ "Insect",
      Diet_type == "Blood" ~ "Blood",
      TRUE ~ Diet_type
    ),
    CLR = as.numeric(CLR)
  )

# Create Background Stripes Data
importance_long_char <- importance_long %>%
  filter(Taxa %in% top_taxa) %>%
  mutate(Taxa_char = as.character(Taxa))

clr_long_top_char <- clr_long_top %>%
  mutate(Taxa_char = as.character(Taxa))

taxa_order_bg <- clr_long_top_char %>%
  group_by(Taxa_char) %>%
  summarise(mean_clr = mean(CLR, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_clr) %>%
  pull(Taxa_char)

bg_stripes <- data.frame(Taxa_char = taxa_order_bg) %>%
  mutate(stripe = (row_number() %% 2 == 1))

importance_long_char$Taxa_char <- factor(importance_long_char$Taxa_char, levels = taxa_order_bg)
clr_long_top_char$Taxa_char <- factor(clr_long_top_char$Taxa_char, levels = taxa_order_bg)

# ---- Plot 1: Variable importance (Family) ----
custom_plot1 <- ggplot() +
  geom_tile(
    data = bg_stripes %>% filter(stripe == TRUE),
    aes(y = Taxa_char, x = 0, height = 0.9, width = Inf),
    fill = "grey90", alpha = 0.6
  ) +
  geom_col(
    data = importance_long_char,
    aes(x = Importance, y = Taxa_char, fill = Diet),
    position = position_dodge(0.8),
    width = 0.7
  ) +
  scale_fill_manual(values = c(    "Insect" = "#ff5e2b",
                                   "Blood"  = "#ba0970",
                                   "Fruit"  = "#2E8B57",
                                   "Fish"   = "#303F9F")
  ) +
  labs(x = "Mean Decrease Accuracy", y = "Family") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10),
        axis.title = element_text(face = "bold"),
        legend.position = "none")

# ---- Plot 2: CLR abundance (Family) ----
custom_plot2 <- ggplot() +
  geom_tile(
    data = bg_stripes %>% filter(stripe == TRUE),
    aes(y = Taxa_char, x = 0, height = 0.9, width = Inf),
    fill = "grey90", alpha = 0.6
  ) +
  geom_boxplot(
    data = clr_long_top_char,
    aes(x = CLR, y = Taxa_char, fill = Diet_type),
    outlier.size = 0.8,
    width = 0.7
  ) +
  scale_fill_manual(values = c(    "Insect" = "#ff5e2b",
                                   "Blood"  = "#ba0970",
                                   "Fruit"  = "#2E8B57",
                                   "Fish"   = "#303F9F")
  ) +
  labs(x = "CLR Abundance", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.position = "none")

# Legend
legend_plot <- ggplot(importance_long_char,
                      aes(x = Importance, y = Taxa_char, fill = Diet)) +
  geom_col(position = position_dodge(0.8)) +
  scale_fill_manual(
    values = c(
      "Insect" = "#ff5e2b",
      "Blood"  = "#ba0970",
      "Fruit"  = "#2E8B57",
      "Fish"   = "#303F9F")
    ,
    name = "Diet Group",
    breaks = c("Insect", "Fruit", "Fish", "Blood"),
    labels = c("Insect", "Fruit", "Fish", "Blood")
  ) +
  theme(legend.position = "bottom",     
        legend.text = element_text(size = 12),              # Size of "Insect", "Fruit", etc.
        legend.title = element_text(size = 12, face = "bold"), # Size of "Diet Group"
        legend.key.size = unit(0.5, "cm")                     # Optional: Makes the color boxes bigger
  )

plot_legend <- cowplot::get_legend(legend_plot)

# Combine plots with legend
final_combined_family <- (custom_plot1 + custom_plot2) / plot_legend + plot_layout(heights = c(10, 1))

print(final_combined_family)

# Print and save final plot
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ext in c("jpg", "pdf", "tiff")) {
  ggsave(
    filename = file.path(out_dir, paste0("CLR_family.", ext)),
    plot = final_combined_family,
    width = 16, 
    height = 10, 
    units = "in", 
    dpi = 300
  )
}
