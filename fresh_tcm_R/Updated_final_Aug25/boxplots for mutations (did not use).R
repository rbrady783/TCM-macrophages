# =========================================
# Exploratory Findings: Visualization Script
# =========================================

# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# Read the results and data
exploratory_dir <- read_csv("exploratory_cytokine_mutation_associations.csv", 
                            show_col_types = FALSE)
df <- read_csv("means_modz_for_mutations.csv", show_col_types = FALSE) %>%
  rename(treatment = 1) %>%
  inner_join(
    read_csv("mutations_updated.csv", show_col_types = FALSE) %>%
      rename(treatment = 1),
    by = "treatment"
  )

# Convert mutation columns to factors
mutations <- setdiff(names(read_csv("mutations_updated.csv", 
                                    show_col_types = FALSE)), 
                     names(read_csv("exploratory_cytokine_mutation_associations.csv", 
                                    show_col_types = FALSE)))
mutations <- mutations[mutations != "treatment"]
df <- df %>% mutate(across(all_of(mutations), ~ as.factor(.x)))

print(paste("Plotting", nrow(exploratory_dir), "exploratory associations"))
print("These represent potential associations that warrant further investigation.")

# Enhanced boxplot function for exploratory findings
plot_exploratory <- function(df, mut_col, cyto_col, pval, out_file) {
  dfp <- df %>%
    select(all_of(c(mut_col, cyto_col))) %>%
    rename(mut = all_of(mut_col), cyto = all_of(cyto_col)) %>%
    mutate(mut = factor(mut))
  
  # Remove empty factor levels
  dfp$mut <- droplevels(dfp$mut)
  
  # Compute vertical position for p-value
  y_max <- max(dfp$cyto, na.rm = TRUE)
  y_min <- min(dfp$cyto, na.rm = TRUE)
  y_range <- y_max - y_min
  y_pos <- y_max + 0.08 * y_range
  
  # Calculate sample sizes for each group
  n_per_group <- dfp %>% 
    group_by(mut) %>% 
    summarise(n = n(), .groups = "drop")
  
  # Create subtitle with sample sizes
  n_text <- paste0("n = ", paste(n_per_group$n, collapse = ", "))
  
  plt <- ggplot(dfp, aes(x = mut, y = cyto, fill = mut)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.8, size = 2.5, color = "black") +
    annotate("text",
             x = Inf, y = y_pos,
             label = paste0("p = ", sprintf("%.3f", pval), "\n(exploratory)"),
             hjust = 1.1, vjust = 0,
             size = 3.5, lineheight = 0.9,
             color = "darkred") +
    labs(
      x        = mut_col,
      y        = paste0(cyto_col, " (z-score)"),
      title    = paste0(cyto_col, " by ", mut_col),
      subtitle = n_text
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "gray50"),
      plot.margin     = margin(5, 5, 20, 5),
      panel.grid.major.y = element_line(color = "gray90", size = 0.3)
    ) +
    coord_cartesian(clip = "off")
  
  ggsave(out_file, plot = plt, width = 5, height = 4.5, dpi = 300)
  invisible(plt)
}

# Create directory and generate plots
dir.create("exploratory_boxplots", showWarnings = FALSE)

print("\nGenerating exploratory boxplots...")
for(i in 1:nrow(exploratory_dir)) {
  mut_col  <- exploratory_dir$mut[i]
  cyto_col <- exploratory_dir$cytokine[i]
  pval     <- exploratory_dir$p.value[i]
  
  fn <- file.path("exploratory_boxplots",
                  paste0(sprintf("%02d", i), "_", mut_col, "_vs_", cyto_col, ".png"))
  
  plot_exploratory(df, mut_col, cyto_col, pval, fn)
  
  if(i %% 5 == 0) print(paste("Generated", i, "of", nrow(exploratory_dir), "plots"))
}

# Create summary visualization of all exploratory findings
summary_plot <- exploratory_dir %>%
  arrange(p.value) %>%
  mutate(
    association = paste0(cytokine, " ~ ", mut),
    rank = row_number()
  ) %>%
  ggplot(aes(x = rank, y = -log10(p.value))) +
  geom_point(size = 3, alpha = 0.7, color = "darkred") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  annotate("text", x = Inf, y = -log10(0.05) + 0.1, 
           label = "p = 0.05", hjust = 1.1, color = "gray50") +
  labs(
    title = "Exploratory Cytokine-Mutation Associations",
    subtitle = paste0(nrow(exploratory_dir), " associations with p < 0.05 (before multiple testing correction)"),
    x = "Rank (by p-value)",
    y = "-log10(p-value)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray60")
  )

ggsave("exploratory_boxplots/00_summary_volcano.png", summary_plot, 
       width = 8, height = 5, dpi = 300)

# Create a summary table sorted by p-value
summary_table <- exploratory_dir %>%
  arrange(p.value) %>%
  mutate(
    rank = row_number(),
    p_formatted = sprintf("%.4f", p.value),
    fdr_formatted = sprintf("%.4f", p.adj)
  ) %>%
  select(rank, cytokine, mut, p_formatted, fdr_formatted, direction) %>%
  rename(
    "Rank" = rank,
    "Cytokine" = cytokine,
    "Mutation/Feature" = mut,
    "P-value" = p_formatted,
    "FDR" = fdr_formatted,
    "Direction" = direction
  )

write_csv(summary_table, "exploratory_boxplots/exploratory_summary_table.csv")
