# Enhanced Volcano Plot Generator with Top 10 Up/Down Gene Annotations
# Author: Enhanced from original code
# Purpose: Generate publication-quality volcano plots for differential expression data

# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# Diagnostic function to check file structure
check_file_structure <- function(file_path) {
  cat(paste("Checking file:", file_path, "\n"))
  
  tryCatch({
    df <- read_csv(file_path, show_col_types = FALSE)
    cat(paste("  Columns:", paste(colnames(df), collapse = ", "), "\n"))
    cat(paste("  Rows:", nrow(df), "\n"))
    
    # Check for Hugo_Symbol column
    if("Hugo_Symbol" %in% colnames(df)) {
      cat("  Hugo_Symbol column: Present\n")
      hugo_na <- sum(is.na(df$Hugo_Symbol))
      hugo_empty <- sum(df$Hugo_Symbol == "", na.rm = TRUE)
      cat(paste("  Hugo_Symbol missing/NA:", hugo_na, "\n"))
      cat(paste("  Hugo_Symbol empty strings:", hugo_empty, "\n"))
    } else if("Hugo" %in% colnames(df)) {
      cat("  Hugo column: Present\n")
      hugo_na <- sum(is.na(df$Hugo))
      hugo_empty <- sum(df$Hugo == "", na.rm = TRUE)
      cat(paste("  Hugo missing/NA:", hugo_na, "\n"))
      cat(paste("  Hugo empty strings:", hugo_empty, "\n"))
    } else {
      cat("  Hugo/Hugo_Symbol column: MISSING\n")
    }
    
    # Check critical columns
    critical_cols <- c("gene", "log2FoldChange", "padj")
    for(col in critical_cols) {
      if(col %in% colnames(df)) {
        missing_count <- sum(is.na(df[[col]]))
        cat(paste("  ", col, "missing:", missing_count, "\n"))
      } else {
        cat(paste("  ", col, ": MISSING COLUMN\n"))
      }
    }
    
  }, error = function(e) {
    cat(paste("  Error reading file:", e$message, "\n"))
  })
  
  cat("\n")
}

# Function to create high-quality volcano plot with top gene annotations
create_volcano_plot <- function(res_df, cytokine_name, 
                                padj_threshold = 0.05, 
                                lfc_threshold = 1,
                                top_n = 10,
                                ycap = 8,
                                output_width = 8, 
                                output_height = 6) {
  
  cat(paste("Processing", cytokine_name, "...\n"))
  
  # Data preprocessing and classification
  volcano_df <- res_df %>%
    # Remove rows with missing critical values
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    # Handle Hugo_Symbol column
    mutate(
      Hugo_Symbol = ifelse(is.na(Hugo_Symbol) | Hugo_Symbol == "" | Hugo_Symbol == "NA", gene, Hugo_Symbol)
    ) %>%
    filter(!is.na(Hugo_Symbol), Hugo_Symbol != "", Hugo_Symbol != "NA") %>%
    mutate(
      negLog10padj = -log10(padj),
      direction = case_when(
        padj < padj_threshold & log2FoldChange >= lfc_threshold ~ "Up",
        padj < padj_threshold & log2FoldChange <= -lfc_threshold ~ "Down",
        TRUE ~ "NS"
      ),
      # Make NS points semi-transparent
      alpha_pt = ifelse(direction == "NS", 0.3, 1)
    )
  
  # Cap extreme -log10(padj) values
  volcano_df <- volcano_df %>%
    mutate(
      above_cap = negLog10padj > ycap,
      y_plot = pmin(negLog10padj, ycap)  # Cap at ycap
    )
  
  # Find most extreme outlier for asterisk annotation
  outlier <- volcano_df %>%
    filter(above_cap) %>%
    slice_max(negLog10padj, n = 1)
  
  # Get top N up and down regulated genes for annotation
  top_up <- volcano_df %>%
    filter(direction == "Up") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top_n)
  
  top_down <- volcano_df %>%
    filter(direction == "Down") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top_n)
  
  # Combine top genes for annotation
  genes_to_label <- bind_rows(top_up, top_down)
  
  # Print summary statistics
  cat(paste("  Total genes:", nrow(volcano_df), "\n"))
  cat(paste("  Significant up-regulated:", nrow(top_up), "\n"))
  cat(paste("  Significant down-regulated:", nrow(top_down), "\n"))
  cat(paste("  Genes to annotate:", nrow(genes_to_label), "\n"))
  
  # Create the volcano plot
  p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = y_plot)) +
    
    # Main scatter plot with appropriate alpha
    geom_point(aes(color = direction, alpha = alpha_pt),
               size = 1.8, show.legend = FALSE) +
    
    # Significance threshold lines
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.5, alpha = 0.8) +
    geom_hline(yintercept = -log10(padj_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.5, alpha = 0.8) +
    
    # Cytokine label in upper-left corner
    annotate("text", x = -Inf, y = Inf, label = cytokine_name,
             hjust = -0.1, vjust = 1.3, size = 6, fontface = "bold",
             color = "black") +
    
    # Asterisk for capped values (if any)
    {if(nrow(outlier) > 0) {
      geom_text(data = outlier, aes(label = "*"), 
                vjust = -0.3, size = 6, color = "black", fontface = "bold")
    } else {
      NULL
    }} +
    
    # Highlight top genes with enhanced visibility
    {if(nrow(genes_to_label) > 0) {
      geom_point(data = genes_to_label,
                 aes(color = direction),
                 size = 3.5, stroke = 1.2, shape = 21, fill = "white",
                 show.legend = FALSE)
    } else {
      NULL
    }} +
    
    # Gene labels with improved repulsion
    {if(nrow(genes_to_label) > 0) {
      geom_text_repel(data = genes_to_label, 
                      aes(label = Hugo_Symbol, color = direction),
                      size = 3.5, fontface = "bold",
                      box.padding = 0.4, point.padding = 0.3,
                      segment.color = "grey30", segment.size = 0.4,
                      max.overlaps = 20, min.segment.length = 0.1,
                      show.legend = FALSE)
    } else {
      NULL
    }} +
    
    # Enhanced color palette (colorblind-friendly)
    scale_color_manual(values = c(
      NS = "#CCCCCC",     # Light grey for NS
      Down = "#3366CC",   # Blue for down-regulated
      Up = "#DC3912"      # Red for up-regulated
    )) +
    scale_alpha_identity() +
    
    # Axis labels
    labs(x = expression(log[2]*"(Fold Change)"),
         y = expression(-log[10]*"(adjusted p-value)")) +
    
    # Set plot limits with some padding
    coord_cartesian(ylim = c(0, ycap + 0.5)) +
    
    # Publication-ready theme
    theme_classic(base_size = 14) +
    theme(
      # Panel and axes
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),  # Remove axis lines since we have panel border
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 14, face = "bold"),
      
      # Grid lines for better readability
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      # Plot margins
      plot.margin = margin(20, 20, 20, 20),
      
      # Remove legend since we're using direct labeling
      legend.position = "none"
    )
  
  # Create output filename
  safe_name <- gsub("[^A-Za-z0-9_-]", "", cytokine_name)
  
  # Save high-quality outputs
  # Vector PDF for manuscripts
  ggsave(paste0("volcano_", safe_name, "_top", top_n, ".pdf"), 
         plot = p_volcano, 
         width = output_width, height = output_height, 
         device = cairo_pdf, dpi = 300)
  
  # High-resolution PNG
  ggsave(paste0("volcano_", safe_name, "_top", top_n, ".png"), 
         plot = p_volcano, 
         width = output_width, height = output_height, 
         dpi = 600)
  
  # Return the plot object and summary info
  return(list(
    plot = p_volcano,
    summary = list(
      cytokine = cytokine_name,
      total_genes = nrow(volcano_df),
      up_regulated = nrow(filter(volcano_df, direction == "Up")),
      down_regulated = nrow(filter(volcano_df, direction == "Down")),
      top_up_genes = top_up$Hugo_Symbol,
      top_down_genes = top_down$Hugo_Symbol,
      genes_labeled = nrow(genes_to_label)
    )
  ))
}

# Function to process multiple cytokines
process_multiple_cytokines <- function(file_paths, cytokine_names) {
  
  if(length(file_paths) != length(cytokine_names)) {
    stop("Number of file paths must match number of cytokine names")
  }
  
  results <- list()
  
  for(i in seq_along(file_paths)) {
    tryCatch({
      # Read the data
      cat(paste("\n--- Processing file", i, "of", length(file_paths), "---\n"))
      res_df <- read_csv(file_paths[i], show_col_types = FALSE)
      
      # Generate the volcano plot
      result <- create_volcano_plot(res_df, cytokine_names[i])
      results[[cytokine_names[i]]] <- result
      
      cat(paste("Successfully processed", cytokine_names[i], "\n"))
      
    }, error = function(e) {
      cat(paste("Error processing", cytokine_names[i], ":", e$message, "\n"))
    })
  }
  
  return(results)
}

# Example usage for single cytokine (using your provided data)
# Uncomment and modify paths as needed:

# For CCL2 example:
# res_df <- read_csv("ccl2_quartile_results.csv")
# ccl2_result <- create_volcano_plot(res_df, "CCL2")
# print(ccl2_result$plot)

# Example for processing your specific cytokines:
# Define your file paths and cytokine display names
file_paths <- c(
  "vegf_quartile_results.csv",
  "il.8_quartile_results.csv",
  "il.10_quartile_results.csv",
  "kc.like_quartile_results.csv",
  "tnfa_quartile_results.csv",
  "tgfb_quartile_results.csv",
  "ccl2_quartile_results.csv"
)

cytokine_names <- c("VEGF", "IL-8", "IL-10", "KC-like", "TNF-α", "TGF-β", "CCL2")

# Process all cytokines
all_results <- process_multiple_cytokines(file_paths, cytokine_names)

# Print summary for all cytokines
for(cytokine in names(all_results)) {
  cat(paste("\n=== Summary for", cytokine, "===\n"))
  summary_info <- all_results[[cytokine]]$summary
  cat(paste("Total genes:", summary_info$total_genes, "\n"))
  cat(paste("Up-regulated:", summary_info$up_regulated, "\n"))
  cat(paste("Down-regulated:", summary_info$down_regulated, "\n"))
  cat(paste("Top up genes:", paste(summary_info$top_up_genes, collapse = ", "), "\n"))
  cat(paste("Top down genes:", paste(summary_info$top_down_genes, collapse = ", "), "\n"))
}

# Print session info for reproducibility
cat("\n=== Session Information ===\n")
print(sessionInfo())