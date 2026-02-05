# Load required libraries ----
library(tidyverse)
library(lmerTest)
library(MuMIn)  # For R²
library(knitr)

# Load and prepare data ----
df <- read_csv("raw_no_ctrls_with_histo.csv", show_col_types = FALSE)

df_long <- df %>%
  pivot_longer(
    cols = starts_with("donor_"),
    names_to = "donor",
    values_to = "concentration"
  ) %>%
  mutate(
    donor = factor(str_remove(donor, "donor_")),
    log_conc = log10(concentration + 1),
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
    histo_4 = fct_infreq(factor(histo_4)),
    histo_5 = fct_infreq(factor(histo_5))
  )

# Helper functions
safe_confint <- function(fit) {
  if (is.null(fit)) {
    tibble(term = character(0), conf.low = double(0), conf.high = double(0))
  } else {
    tryCatch({
      ci <- confint(fit, method = "Wald", parm = "beta_")
      tibble(
        term = rownames(ci),
        conf.low = ci[, 1],
        conf.high = ci[, 2]
      )
    }, error = function(e) {
      tibble(term = character(0), conf.low = double(0), conf.high = double(0))
    })
  }
}

safe_r2 <- function(fit) {
  if (is.null(fit)) {
    tibble(marginal = NA_real_, conditional = NA_real_)
  } else {
    tryCatch({
      r2 <- r.squaredGLMM(fit)
      tibble(marginal = r2[1, "R2m"], conditional = r2[1, "R2c"])
    }, error = function(e) {
      tibble(marginal = NA_real_, conditional = NA_real_)
    })
  }
}

# Complete analysis function
analyze_histology_complete <- function(data, histology_col, system_name) {
  model_results <- data %>%
    group_by(cytokine, cytokine_clean) %>%
    nest() %>%
    mutate(
      # Fit LMM
      model = map(data, ~ {
        formula_str <- paste0("log_conc ~ ", histology_col, " + (1|donor)")
        lmer(as.formula(formula_str), data = .x)
      }),
      # Get ALL fixed effects (intercept + all contrasts)
      fixed_effects = map(model, ~ {
        broom.mixed::tidy(.x, effects = "fixed") %>%
          select(term, estimate, std.error, statistic, p.value)
      }),
      # Get confidence intervals
      ci = map(model, safe_confint),
      # Get R²
      r2 = map(model, safe_r2),
      # Get overall F-test
      anova_result = map(model, ~ {
        anova_res <- anova(.x)
        tibble(
          f_statistic = anova_res$`F value`[1],
          df_num = anova_res$NumDF[1],
          df_denom = anova_res$DenDF[1],
          p_value_overall = anova_res$`Pr(>F)`[1]
        )
      })
    )
  
  # Unnest fixed effects
  fixed_df <- model_results %>%
    select(cytokine, cytokine_clean, fixed_effects) %>%
    unnest(fixed_effects)
  
  # Unnest CIs
  ci_df <- model_results %>%
    select(cytokine, ci) %>%
    unnest(ci)
  
  # Unnest R² and ANOVA
  r2_anova_df <- model_results %>%
    select(cytokine, r2, anova_result) %>%
    unnest(c(r2, anova_result))
  
  # Combine everything
  complete_results <- fixed_df %>%
    left_join(ci_df, by = c("cytokine", "term")) %>%
    left_join(r2_anova_df, by = "cytokine") %>%
    mutate(
      system = system_name,
      # Clean up term names
      term_clean = str_remove(term, paste0(histology_col))
    )
  
  return(complete_results)
}

# Run analysis for both systems
cat("Analyzing histology effects - this may take a moment...\n")
complete_4cat <- analyze_histology_complete(df_long, "histo_4", "4-category")
complete_5cat <- analyze_histology_complete(df_long, "histo_5", "5-category")

# Combine
all_results <- bind_rows(complete_4cat, complete_5cat)

# Apply FDR correction to individual contrasts WITHIN each system
all_results <- all_results %>%
  group_by(system) %>%
  mutate(p_adj_contrast = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  # Also FDR correction for overall F-tests
  group_by(system) %>%
  mutate(p_adj_overall = p.adjust(p_value_overall, method = "BH")) %>%
  ungroup()

# Format for publication
final_table <- all_results %>%
  mutate(
    Cytokine = cytokine_clean,
    System = system,
    Term = term_clean,
    Estimate = sprintf("%.3f", estimate),
    SE = sprintf("%.3f", std.error),
    `95% CI` = sprintf("[%.3f, %.3f]", conf.low, conf.high),
    t = sprintf("%.2f", statistic),
    `p-value` = if_else(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    `Adj. p` = if_else(p_adj_contrast < 0.001, "<0.001", sprintf("%.3f", p_adj_contrast)),
    `Overall F` = sprintf("%.2f", f_statistic),
    `Overall F df` = sprintf("(%d, %.1f)", df_num, df_denom),
    `Overall p` = if_else(p_value_overall < 0.001, "<0.001", sprintf("%.3f", p_value_overall)),
    `Overall Adj. p` = if_else(p_adj_overall < 0.001, "<0.001", sprintf("%.3f", p_adj_overall)),
    `Marginal R²` = sprintf("%.3f", marginal),
    `Conditional R²` = sprintf("%.3f", conditional)
  ) %>%
  select(Cytokine, System, Term, Estimate, SE, `95% CI`, t, `p-value`, `Adj. p`,
         `Overall F`, `Overall F df`, `Overall p`, `Overall Adj. p`, 
         `Marginal R²`, `Conditional R²`) %>%
  arrange(System, Cytokine, Term)

# Print summary
cat("\n=== COMPLETE S4 TABLE ===\n")
print(kable(final_table, caption = "Linear mixed-effects model results for tumor histology effects on cytokine secretion"))

# Save
write_csv(final_table, "S4_Table_histology_COMPLETE.csv")

cat("\nComplete table saved as S4_Table_histology_COMPLETE.csv\n")
cat("Total rows:", nrow(final_table), "\n")
cat("This includes all coefficient estimates for each tumor type contrast\n")
