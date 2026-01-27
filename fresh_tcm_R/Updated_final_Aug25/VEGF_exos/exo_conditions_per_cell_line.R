# Load required libraries
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(emmeans)

# Use your existing full_data
# Analyze each cell line separately
analyze_cell_line <- function(cell_line_name) {
  # Filter data for specific cell line
  cell_data <- full_data %>% filter(Cell_Line == cell_line_name)
  
  cat(sprintf("\n=== %s ANALYSIS ===\n", cell_line_name))
  
  # Mixed-effects model with donor as random effect
  model <- lmer(VEGF ~ Condition + (1|Donor), data = cell_data)
  
  # ANOVA
  anova_result <- anova(model)
  print(anova_result)
  
  # Post-hoc comparisons
  emm <- emmeans(model, ~ Condition)
  pairwise <- pairs(emm, adjust = "tukey")
  
  cat("\nPairwise comparisons:\n")
  print(pairwise)
  
  # Extract p-values
  pairwise_summary <- summary(pairwise)
  cat("\nP-VALUES:\n")
  for(i in 1:nrow(pairwise_summary)) {
    comparison <- pairwise_summary$contrast[i]
    p_value <- pairwise_summary$p.value[i]
    significance <- ifelse(p_value < 0.001, "***", 
                           ifelse(p_value < 0.01, "**",
                                  ifelse(p_value < 0.05, "*", "ns")))
    cat(sprintf("%-30s p = %.4f %s\n", comparison, p_value, significance))
  }
  
  # Calculate means
  means <- cell_data %>% 
    group_by(Condition) %>% 
    summarise(Mean = mean(VEGF), SD = sd(VEGF), .groups = 'drop')
  
  cat("\nMeans:\n")
  print(means)
  
  cat(paste0("\n", paste(rep("=", 50), collapse = ""), "\n"))
  
  return(list(model = model, pairwise = pairwise_summary, means = means))
}

# Run analysis for all cell lines
cell_lines <- c("Jones", "MacKinley", "Nike", "OS2.4", "OSA8", "Parks", "STS-1")
results <- list()

for(cell_line in cell_lines) {
  results[[cell_line]] <- analyze_cell_line(cell_line)
}

# Create a comprehensive plot showing all cell lines
p1 <- ggplot(full_data, aes(x = Condition, y = VEGF, fill = Condition)) +
  # Bars
  stat_summary(fun = "mean", geom = "col", alpha = 0.8, width = 0.7, 
               color = "black", linewidth = 0.3) +
  # Individual points
  geom_point(position = position_jitter(width = 0.2, seed = 123), 
             size = 1, alpha = 0.7, color = "black") +
  # Error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.2, linewidth = 0.5, color = "black") +
  # Colors
  scale_fill_manual(values = c("Whole TCM" = "#E91E63", 
                               "Exosome Depleted" = "#696969", 
                               "Exosome Only" = "#00BCD4")) +
  # Facet by cell line
  facet_wrap(~ Cell_Line, scales = "free_y", ncol = 4) +
  # Clean labels
  scale_x_discrete(labels = c("Whole", "Depleted", "Exo Only")) +
  labs(y = "VEGF (pg/ml)", x = NULL) +
  # Clean theme
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(color = "black", face = "bold"),
    strip.text = element_text(face = "bold", size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # Y-axis
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

print(p1)

# Save plot
ggsave("cellline_exosome_analysis.pdf", plot = p1, 
       width = 10, height = 6, units = "in", dpi = 300)

# Summary across all cell lines
cat("\n=== SUMMARY ACROSS ALL CELL LINES ===\n")
cat("Examining which cell lines show significant exosome effects...\n")

# Check which comparisons are significant for each cell line
for(cell_line in cell_lines) {
  pairwise_data <- results[[cell_line]]$pairwise
  
  whole_vs_exo <- pairwise_data$p.value[grep("Whole.*Exosome Only", pairwise_data$contrast)]
  depleted_vs_exo <- pairwise_data$p.value[grep("Depleted.*Exosome Only", pairwise_data$contrast)]
  
  cat(sprintf("%-10s: Whole vs Exo-Only p=%.3f %s, Depleted vs Exo-Only p=%.3f %s\n", 
              cell_line,
              whole_vs_exo, ifelse(whole_vs_exo < 0.05, "*", "ns"),
              depleted_vs_exo, ifelse(depleted_vs_exo < 0.05, "*", "ns")))
}

# ==========================================
# CREATE S8 TABLE WITH 95% CIs
# ==========================================

# Extract complete statistics for whole vs exosome-only comparison
s8_table_data <- list()

for(cell_line in cell_lines) {
  # Get model and emmeans
  cell_data <- full_data %>% filter(Cell_Line == cell_line)
  model <- lmer(VEGF ~ Condition + (1|Donor), data = cell_data)
  emm <- emmeans(model, ~ Condition)
  
  # Get pairwise with CIs
  pairwise <- pairs(emm, adjust = "tukey")
  pairwise_ci <- confint(pairwise)  # Get confidence intervals
  
  # Extract whole vs exosome-only comparison
  comparison_idx <- grep("Whole.*Exosome Only", pairwise_ci$contrast)
  
  s8_table_data[[cell_line]] <- data.frame(
    Cell_Line = cell_line,
    Estimate = pairwise_ci$estimate[comparison_idx],
    SE = pairwise_ci$SE[comparison_idx],
    CI_lower = pairwise_ci$lower.CL[comparison_idx],
    CI_upper = pairwise_ci$upper.CL[comparison_idx],
    df = pairwise_ci$df[comparison_idx],
    t = pairwise_ci$estimate[comparison_idx] / pairwise_ci$SE[comparison_idx],
    p = summary(pairwise)$p.value[comparison_idx]
  )
}

# Combine into single table
s8_table <- bind_rows(s8_table_data) %>%
  mutate(
    `Cell Line` = Cell_Line,
    `Estimate` = round(Estimate, 1),
    `SE` = round(SE, 1),
    `95% CI` = sprintf("[%.0f, %.0f]", CI_lower, CI_upper),
    `df` = round(df, 0),
    `t` = round(t, 3),
    `p` = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  ) %>%
  select(`Cell Line`, Estimate, SE, `95% CI`, df, t, p)

# Print table
cat("\n=== S8 TABLE: Complete Statistics ===\n")
print(s8_table, row.names = FALSE)

# Save as CSV
write.csv(s8_table, "S8_Table_Individual_CellLines.csv", row.names = FALSE)

cat("\nTable saved as S8_Table_Individual_CellLines.csv\n")

####################################

# ==========================================
# CREATE COMPLETE S8 TABLE - ALL PAIRWISE COMPARISONS
# ==========================================

s8_complete_data <- list()

for(cell_line in cell_lines) {
  # Get model and emmeans
  cell_data <- full_data %>% filter(Cell_Line == cell_line)
  model <- lmer(VEGF ~ Condition + (1|Donor), data = cell_data)
  emm <- emmeans(model, ~ Condition)
  
  # Get all pairwise comparisons with CIs
  pairwise <- pairs(emm, adjust = "tukey")
  pairwise_ci <- confint(pairwise)
  pairwise_summary <- summary(pairwise)
  
  # Extract all three comparisons
  for(i in 1:nrow(pairwise_ci)) {
    s8_complete_data[[paste(cell_line, i)]] <- data.frame(
      Cell_Line = cell_line,
      Comparison = as.character(pairwise_ci$contrast[i]),
      Estimate = pairwise_ci$estimate[i],
      SE = pairwise_ci$SE[i],
      CI_lower = pairwise_ci$lower.CL[i],
      CI_upper = pairwise_ci$upper.CL[i],
      df = pairwise_ci$df[i],
      t = pairwise_summary$t.ratio[i],
      p = pairwise_summary$p.value[i]
    )
  }
}

# Combine and format
s8_complete <- bind_rows(s8_complete_data) %>%
  mutate(
    `Cell Line` = Cell_Line,
    `Comparison` = Comparison,
    `Estimate` = round(Estimate, 1),
    `SE` = round(SE, 1),
    `95% CI` = sprintf("[%.0f, %.0f]", CI_lower, CI_upper),
    `df` = round(df, 0),
    `t` = round(t, 3),
    `p-value` = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
    `Sig` = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(`Cell Line`, Comparison, Estimate, SE, `95% CI`, df, t, `p-value`, Sig)

# Print
cat("\n=== COMPLETE S8 TABLE ===\n")
print(s8_complete, row.names = FALSE)

# Save
write.csv(s8_complete, "S8_Table_AllComparisons.csv", row.names = FALSE)

cat("\nTable saved!\n")
