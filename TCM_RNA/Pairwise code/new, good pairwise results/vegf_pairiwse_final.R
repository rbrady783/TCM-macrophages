# =========================================
# Flexible Quartile-Based DEG Analysis
# =========================================

# 0) Load libraries ----
library(dplyr)
library(DESeq2)
library(ggplot2)
library(tibble)
library(ashr)
library(ggrepel)

# =========================================
# CONFIGURATION - Change these parameters
# =========================================
CYTOKINE <- "VEGF"  # Change to: VEGF, IL.8, KC.like, IL.10, TNF.a, TGF.b, CCL2
count_threshold <- 3
min_samples <- 6
ycap <- 6

# 1) Read data ----
raw <- read.csv("rawcountdata.csv", stringsAsFactors = FALSE)
cytokine_data <- read.csv("means_modz.csv", stringsAsFactors = FALSE)

# 2) Create sample mapping between files ----
# Handle naming differences (X1771 <-> 1771, D.17 <-> D-17, etc.)
create_sample_mapping <- function(raw_samples, cytokine_treatments) {
  cytokine_clean <- trimws(as.character(cytokine_treatments))
  
  message("Attempting to match sample names...")
  message("Raw count samples (first 10):", paste(head(raw_samples, 10), collapse = ", "))
  message("Cytokine data samples (first 10):", paste(head(cytokine_clean, 10), collapse = ", "))
  
  mapping <- sapply(raw_samples, function(sample) {
    # Skip the first column if it's gene IDs
    if (sample %in% c("Hugo_Symbol", "Ensembl_Gene_Id", "gene", "Gene", "GENE")) {
      return(NA)
    }
    
    # Try exact match
    if (sample %in% cytokine_clean) {
      message(sprintf("  %s -> %s (exact match)", sample, sample))
      return(sample)
    }
    
    # Try with/without X prefix
    sample_no_x <- sub("^X", "", sample)
    if (sample_no_x %in% cytokine_clean) {
      message(sprintf("  %s -> %s (removed X prefix)", sample, sample_no_x))
      return(sample_no_x)
    }
    
    # Try adding X prefix
    sample_with_x <- paste0("X", sample)
    if (sample_with_x %in% cytokine_clean) {
      message(sprintf("  %s -> %s (added X prefix)", sample, sample_with_x))
      return(sample_with_x)
    }
    
    # Try dots vs dashes
    sample_dash <- gsub("\\.", "-", sample)
    sample_dot <- gsub("-", "\\.", sample)
    if (sample_dash %in% cytokine_clean) {
      message(sprintf("  %s -> %s (dot to dash)", sample, sample_dash))
      return(sample_dash)
    }
    if (sample_dot %in% cytokine_clean) {
      message(sprintf("  %s -> %s (dash to dot)", sample, sample_dot))
      return(sample_dot)
    }
    
    # Try case-insensitive matching
    case_match <- cytokine_clean[tolower(cytokine_clean) == tolower(sample)]
    if (length(case_match) > 0) {
      message(sprintf("  %s -> %s (case insensitive)", sample, case_match[1]))
      return(case_match[1])
    }
    
    # Try partial matching (for cases like spaces or special characters)
    clean_sample <- gsub("[^a-zA-Z0-9]", "", tolower(sample))
    clean_cytokines <- gsub("[^a-zA-Z0-9]", "", tolower(cytokine_clean))
    partial_match <- cytokine_clean[clean_cytokines == clean_sample]
    if (length(partial_match) > 0) {
      message(sprintf("  %s -> %s (partial/clean match)", sample, partial_match[1]))
      return(partial_match[1])
    }
    
    message(sprintf("  %s -> NO MATCH FOUND", sample))
    return(NA)
  })
  
  matched_count <- sum(!is.na(mapping))
  message(sprintf("\nSample matching summary: %d/%d samples matched", 
                  matched_count, length(mapping)))
  
  if (matched_count < 4) {
    warning("Too few samples matched for quartile analysis. Need at least 4 samples.")
    message("Unmatched raw samples:", paste(names(mapping)[is.na(mapping)], collapse = ", "))
    message("Available cytokine samples:", paste(cytokine_clean, collapse = ", "))
  }
  
  return(mapping)
}

# 3) Generate quartile groups from matched samples only ----
# Get samples that exist in raw data (excluding gene ID column)
raw_sample_names <- colnames(raw)[-1]
sample_mapping <- create_sample_mapping(raw_sample_names, cytokine_data$treatment)

# Get cytokine values for ONLY the matched samples
matched_data <- data.frame()
for (sample in names(sample_mapping)) {
  mapped_name <- sample_mapping[sample]
  if (!is.na(mapped_name)) {
    cytokine_value <- cytokine_data[[CYTOKINE]][cytokine_data$treatment == mapped_name]
    if (length(cytokine_value) > 0 && !is.na(cytokine_value)) {
      matched_data <- rbind(matched_data, 
                            data.frame(raw_sample = sample, 
                                       cytokine_sample = mapped_name,
                                       cytokine_value = cytokine_value))
    }
  }
}

# Sort by cytokine value to determine quartiles
matched_data <- matched_data[order(matched_data$cytokine_value), ]

message(sprintf("Matched samples for %s analysis (n=%d):", CYTOKINE, nrow(matched_data)))
for (i in 1:nrow(matched_data)) {
  message(sprintf("  %d. %s: %.3f", i, matched_data$raw_sample[i], matched_data$cytokine_value[i]))
}

# Calculate quartile positions for the matched samples
n_matched <- nrow(matched_data)
n_per_quartile <- round(n_matched * 0.25)  # Use round instead of floor for proper quartiles

message(sprintf("Quartile calculation: %d samples total, %d samples per quartile (25%%)", 
                n_matched, n_per_quartile))

# Bottom quartile (lowest values)
bottom_quartile_samples <- matched_data$raw_sample[1:n_per_quartile]
# Top quartile (highest values) 
top_quartile_samples <- matched_data$raw_sample[(n_matched - n_per_quartile + 1):n_matched]

# Create groups dataframe
groups_df <- data.frame(
  id = c(bottom_quartile_samples, top_quartile_samples),
  group = factor(c(rep("low", length(bottom_quartile_samples)), 
                   rep("high", length(top_quartile_samples))),
                 levels = c("low", "high"))
)

message(sprintf("Generated groups for %s:", CYTOKINE))
message(sprintf("  Low group (n=%d): %s", 
                sum(groups_df$group == "low"),
                paste(groups_df$id[groups_df$group == "low"], collapse = ", ")))
message(sprintf("  High group (n=%d): %s", 
                sum(groups_df$group == "high"),
                paste(groups_df$id[groups_df$group == "high"], collapse = ", ")))

# 4) Subset samples based on generated groups ----
selected_cols <- c(1, which(colnames(raw) %in% groups_df$id))
raw_subset <- raw[, selected_cols]

# 5) Move gene IDs into rownames & drop that column ----
rownames(raw_subset) <- raw_subset[[1]]
raw_subset <- raw_subset[, -1]

# 6) Check sample names against metadata ----
if (!all(colnames(raw_subset) %in% groups_df$id)) {
  stop("Sample names in raw_subset do not match groups_df$id!")
} else {
  message("All sample names match.")
}

# 7) Prepare metadata ----
rownames(groups_df) <- groups_df$id

# 8) Prefilter low-count genes ----
keep <- rowSums(raw_subset >= count_threshold) >= min_samples
filtered_raw <- raw_subset[keep, ]

# 9) Build & run DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = filtered_raw,
  colData   = groups_df[colnames(filtered_raw), ],
  design    = ~ group
)

dds <- DESeq(dds)

res_shrunk <- lfcShrink(
  dds,
  coef = "group_high_vs_low",
  type = "ashr"
)

# 10) Export results ----
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# Write results
output_file <- paste0(tolower(CYTOKINE), "_quartile_results.csv")
write.csv(res_df, output_file, row.names = FALSE)

message(sprintf("Results exported to: %s", output_file))
message(sprintf("Total genes tested: %d", nrow(res_df)))
message(sprintf("Significant genes (padj < 0.05): %d", sum(res_df$padj < 0.05, na.rm = TRUE)))
message(sprintf("Upregulated (padj < 0.05, log2FC >= 1): %d", 
                sum(res_df$padj < 0.05 & res_df$log2FoldChange >= 1, na.rm = TRUE)))
message(sprintf("Downregulated (padj < 0.05, log2FC <= -1): %d", 
                sum(res_df$padj < 0.05 & res_df$log2FoldChange <= -1, na.rm = TRUE)))

