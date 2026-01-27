# ===============================================================
# Streamlined Cytokine-Histology Analysis
# Statistical analysis of cytokine secretion by tumor histology
# PLOS ONE COMPLIANT VERSION
# ===============================================================
install.packages("lme4")

# Load required libraries ----
library(tidyverse)
library(lmerTest)      # For LMMs with p-values
library(patchwork)     # For combining plots

# Load and prepare data ----
df <- read_csv("raw_no_ctrls_with_histo.csv", show_col_types = FALSE)

# Create long format with both histology systems
df_long <- df %>%
  pivot_longer(
    cols = starts_with("donor_"),
    names_to = "donor",
    values_to = "concentration"
  ) %>%
  mutate(
    donor = factor(str_remove(donor, "donor_")),
    log_conc = log10(concentration + 1),  # log10 for easier interpretation
    # Clean up cytokine names for plotting
    cytokine_clean = case_when(
      cytokine == "vegf" ~ "VEGF",
      cytokine == "il.8" ~ "IL-8", 
      cytokine == "kc.like" ~ "KC-like",
      cytokine == "il.10" ~ "IL-10",
      cytokine == "ccl2" ~ "CCL2",
      cytokine == "tnf.a" ~ "TNF-α",
      cytokine == "tgf.b" ~ "TGF-β",
      TRUE ~ cytokine
    ),
    # Reorder histology factors by frequency (most common first)
    histo_4 = fct_infreq(factor(histo_4)),
    histo_5 = fct_infreq(factor(histo_5))
  )

# Function to perform statistical analysis ----
analyze_histology <- function(data, histology_col) {
  # Fit linear mixed models and perform ANOVA
  model_results <- data %>%
    group_by(cytokine, cytokine_clean) %>%
    nest() %>%
    mutate(
      # Fit LMM with donor as random effect
      model = map(data, ~ {
        formula_str <- paste0("log_conc ~ ", histology_col, " + (1|donor)")
        lmer(as.formula(formula_str), data = .x)
      }),
      # Perform ANOVA on the model
      anova_result = map(model, ~ {
        anova_res <- anova(.x)
        tibble(
          f_statistic = anova_res$`F value`[1],
          p_value_anova = anova_res$`Pr(>F)`[1]
        )
      })
    )
  
  # Extract and format ANOVA results
  anova_results <- model_results %>%
    select(cytokine, cytokine_clean, anova_result) %>%
    unnest(anova_result) %>%
    mutate(
      p_adj_anova = p.adjust(p_value_anova, method = "BH"),
      significance = case_when(
        p_adj_anova < 0.001 ~ "***",
        p_adj_anova < 0.01 ~ "**", 
        p_adj_anova < 0.05 ~ "*",
        p_adj_anova < 0.1 ~ "†",
        TRUE ~ "ns"
      )
    )
  
  return(list(
    models = model_results,
    anova = anova_results
  ))
}

# Analyze both histology systems ----
cat("Analyzing 4-category histology system...\n")
results_4cat <- analyze_histology(df_long, "histo_4")

cat("Analyzing 5-category histology system...\n") 
results_5cat <- analyze_histology(df_long, "histo_5")

# Display results ----
cat("\n=== ANOVA RESULTS (4-category system) ===\n")
print(results_4cat$anova %>% arrange(p_adj_anova))

cat("\n=== ANOVA RESULTS (5-category system) ===\n")
print(results_5cat$anova %>% arrange(p_adj_anova))

# PLOS ONE compliant visualization function ----
create_cytokine_plot <- function(data, histology_col, include_title = TRUE) {
  p <- data %>%
    ggplot(aes(x = !!sym(histology_col), y = log_conc)) +
    geom_boxplot(fill = "lightgrey", color = "black", alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 1, color = "black") +
    facet_wrap(~ cytokine_clean, scales = "free_y", ncol = 3) +
    labs(
      title = if(include_title) "Cytokine Secretion by Tumor Type" else NULL,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 9, base_family = "sans") +  # PLOS: sans font, 9pt base
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", color = "black", size = 9),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5, color = "black"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      panel.border = element_rect(color = "black", fill = NA),
      plot.margin = unit(c(2, 2, 2, 2), "pt"),  # 2pt margins
      text = element_text(family = "sans")
    )
  
  return(p)
}

# Generate plots ----
plot_5cat <- create_cytokine_plot(df_long, "histo_5", include_title = TRUE)
plot_4cat <- create_cytokine_plot(df_long, "histo_4", include_title = FALSE)

# Create stacked plot with y-axis label ----
final_plot <- (plot_5cat / plot_4cat) + 
  plot_layout(heights = c(1, 1))

# Add y-axis label
y_label <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Log10(Concentration + 1)", 
           angle = 90, size = 3.5, fontface = "bold", family = "sans") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# Combine with y-axis label
final_plot_with_label <- wrap_plots(y_label, final_plot, widths = c(0.05, 1))

# Display the final plot
print(final_plot_with_label)

# Save as TIFF with LZW compression (PLOS ONE requirement)
ggsave("cytokine_histology_plot.tiff", 
       plot = final_plot_with_label,
       device = "tiff",
       width = 7,        # Within PLOS 2.63-7.5" width limit
       height = 8.5,     # Within PLOS 8.75" height limit
       units = "in",
       dpi = 600,        # PLOS: 300-600 dpi
       compression = "lzw",  # PLOS requirement
       bg = "white")

# Optional: Also save as EPS (vector format, also accepted by PLOS)
ggsave("cytokine_histology_plot.eps", 
       plot = final_plot_with_label,
       device = cairo_ps,
       width = 7, 
       height = 8.5, 
       units = "in",
       bg = "white")

# Summary statistics ----
cat("\n=== SUMMARY STATISTICS ===\n")
summary_stats <- df_long %>%
  group_by(cytokine_clean, histo_4) %>%
  summarise(
    n_treatments = n_distinct(treatment),
    mean_log_conc = mean(log_conc, na.rm = TRUE),
    sd_log_conc = sd(log_conc, na.rm = TRUE),
    median_conc = median(concentration, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

cat("\n=== PLOS ONE Compliance Checklist ===\n")
cat("✓ Format: TIFF with LZW compression\n")
cat("✓ Resolution: 600 dpi (300-600 required)\n")
cat("✓ Dimensions: 7 x 8.5 inches (within limits)\n")
cat("✓ Font: Sans-serif 8-11 pt (8-12 pt required)\n")
cat("✓ Color mode: RGB (default)\n")
cat("\nAnalysis complete!\n")
