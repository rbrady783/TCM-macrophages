############################################################
# Analysis: Do cell counts predict cytokine secretion?
# Model: Linear mixed-effects with donor as random intercept
# Transformation selection based on residual normality
# Workflow:
#   1) Read & reshape raw cytokine data
#   2) Join cell counts by treatment
#   3) Standardize predictor (cell_count per million)
#   4) Fit RAW + LOG models; choose scale by residual normality
#   5) Tidy fixed effects (broom.mixed), FDR-correct, export table
############################################################

# ---- Libraries ----
library(tidyverse)
library(lmerTest)        # lmer() with p-values
library(broom.mixed)     # tidy() for mixed models
library(knitr)           # kable()

# ---- 1) Read and reshape RAW cytokine data ----
# Expected columns in raw_no_ctrls.csv:
#   cytokine | treatment | donor_1 | donor_2 | donor_3
raw_df <- read_csv("raw_no_ctrls.csv", show_col_types = FALSE)

raw_long <- raw_df %>%
  pivot_longer(
    cols         = starts_with("donor_"),
    names_to     = "donor",
    names_prefix = "donor_",
    values_to    = "raw_val"
  ) %>%
  mutate(
    donor     = str_trim(donor),
    treatment = str_trim(treatment),
    cytokine  = str_trim(cytokine)
  )

# ---- 2) Read cell counts and align naming ----
# Expected columns in cell_counts.csv:
#   cell_line | cell_count
cell_counts <- read_csv("cell_counts.csv", show_col_types = FALSE) %>%
  rename(treatment = cell_line) %>%
  mutate(treatment = str_trim(treatment))

# ---- 3) Join by treatment ----
df_cc_raw <- raw_long %>%
  left_join(cell_counts, by = "treatment") %>%
  mutate(
    donor    = factor(donor),
    cytokine = factor(cytokine)
  )

# ---- 4) Quick sanity checks ----
df_cc_raw %>%
  summarise(
    total_rows    = n(),
    missing_cells = sum(is.na(cell_count)),
    missing_raw   = sum(is.na(raw_val))
  ) %>% print()

df_cc_raw %>%
  select(cytokine, treatment, donor, raw_val, cell_count) %>%
  slice(1:10) %>%
  print(n = 10)

# ---- 5) Prepare predictor: scale cell count PER MILLION ----
# Standardize globally (shared scale across cytokines)
df_cc_scaled <- df_cc_raw %>%
  mutate(
    cell_millions   = cell_count / 1e6,
    cell_millions_z = as.numeric(scale(cell_millions))  # mean 0, sd 1
  )

# ---- Helpers: robust model fit + robust tidy ----
safe_lmer <- function(formula, data) {
  tryCatch(lmer(formula, data = data), error = function(e) NULL)
}

# Always return a tibble with expected columns; empty on failure
safe_tidy_fixed <- function(fit) {
  if (is.null(fit)) {
    tibble(
      term      = character(0),
      estimate  = double(0),
      std.error = double(0),
      statistic = double(0),
      p.value   = double(0)
    )
  } else {
    broom.mixed::tidy(fit, effects = "fixed") %>%
      select(term, estimate, std.error, statistic, p.value)
  }
}

# Extract confidence intervals
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

# Extract R-squared values
safe_r2 <- function(fit) {
  if (is.null(fit)) {
    tibble(marginal = NA_real_, conditional = NA_real_)
  } else {
    tryCatch({
      r2 <- MuMIn::r.squaredGLMM(fit)
      tibble(marginal = r2[1, "R2m"], conditional = r2[1, "R2c"])
    }, error = function(e) {
      tibble(marginal = NA_real_, conditional = NA_real_)
    })
  }
}

# Check residual normality
check_residuals <- function(fit) {
  if(is.null(fit)) return(NA)
  resid <- residuals(fit)
  if(length(resid) >= 3 && length(resid) <= 5000) shapiro.test(resid)$p.value else NA
}

# ---- 6) Per-cytokine analysis: COMPLETE models ----
library(MuMIn)  # For r.squaredGLMM

cell_effects_analysis <- df_cc_scaled %>%
  filter(!is.na(raw_val), !is.na(cell_millions_z)) %>%
  group_by(cytokine) %>%
  nest() %>%
  mutate(
    # Fit both models
    fit_raw = map(data, ~ safe_lmer(raw_val ~ cell_millions_z + (1 | donor), data = .x)),
    fit_log = map(data, ~ safe_lmer(log1p(raw_val) ~ cell_millions_z + (1 | donor), data = .x)),
    
    # Check residual normality
    shapiro_resid_raw = map_dbl(fit_raw, check_residuals),
    shapiro_resid_log = map_dbl(fit_log, check_residuals),
    
    # Tidy fixed effects (ALL terms, not just slope)
    tab_raw = map(fit_raw, safe_tidy_fixed),
    tab_log = map(fit_log, safe_tidy_fixed),
    
    # Get confidence intervals
    ci_raw = map(fit_raw, safe_confint),
    ci_log = map(fit_log, safe_confint),
    
    # Get R-squared
    r2_raw = map(fit_raw, safe_r2),
    r2_log = map(fit_log, safe_r2)
  ) %>%
  ungroup()

# ---- 7) Decide transformation ----
normality_check <- cell_effects_analysis %>%
  transmute(
    cytokine,
    shapiro_resid_raw,
    shapiro_resid_log,
    use_log = case_when(
      is.na(shapiro_resid_raw) & is.na(shapiro_resid_log) ~ FALSE,
      is.na(shapiro_resid_raw) ~ TRUE,
      is.na(shapiro_resid_log) ~ FALSE,
      shapiro_resid_log > shapiro_resid_raw ~ TRUE,
      TRUE ~ FALSE
    )
  )

# ---- 8) Select appropriate results and combine ALL info (FIXED) ----
cell_effects_complete <- cell_effects_analysis %>%
  select(cytokine, tab_raw, tab_log, ci_raw, ci_log, r2_raw, r2_log) %>%
  left_join(normality_check, by = "cytokine") %>%
  mutate(
    selected_tab = if_else(use_log, tab_log, tab_raw),
    selected_ci = if_else(use_log, ci_log, ci_raw),
    selected_r2 = if_else(use_log, r2_log, r2_raw)
  ) %>%
  select(cytokine, selected_tab, selected_ci, selected_r2, use_log)

# Unnest tab and r2 first
cell_effects_tab <- cell_effects_complete %>%
  unnest(c(selected_tab, selected_r2))

# Unnest CI separately
cell_effects_ci <- cell_effects_complete %>%
  select(cytokine, selected_ci) %>%
  unnest(selected_ci)

# Join them together
cell_effects_complete <- cell_effects_tab %>%
  left_join(cell_effects_ci, by = c("cytokine", "term"))

# ---- 9) FDR correction on SLOPE only ----
slope_results <- cell_effects_complete %>%
  filter(term == "cell_millions_z") %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr"))

# ---- 10) COMPLETE publication table ----
complete_table <- slope_results %>%
  left_join(normality_check %>% select(cytokine, shapiro_resid_raw, shapiro_resid_log), 
            by = "cytokine") %>%
  mutate(
    Transformation = if_else(use_log, "log(x+1)", "raw"),
    `Shapiro p` = if_else(use_log, 
                          if_else(is.na(shapiro_resid_log), "NA", 
                                  if_else(shapiro_resid_log < 0.001, "<0.001", sprintf("%.3f", shapiro_resid_log))),
                          if_else(is.na(shapiro_resid_raw), "NA", 
                                  if_else(shapiro_resid_raw < 0.001, "<0.001", sprintf("%.3f", shapiro_resid_raw)))),
    `95% CI` = sprintf("[%.3f, %.3f]", conf.low, conf.high),
    `p-value` = if_else(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    `Adj. p` = if_else(p.adj < 0.001, "<0.001", sprintf("%.3f", p.adj))
  ) %>%
  select(
    Cytokine = cytokine,
    Estimate = estimate,
    SE = std.error,
    `95% CI`,
    t = statistic,
    `p-value`,
    `Adj. p`,
    `Marginal R²` = marginal,
    `Conditional R²` = conditional,
    `Shapiro p`,
    Transformation
  ) %>%
  arrange(Cytokine)

# Print and save
kable(complete_table, digits = 3)
write_csv(complete_table, "S3_Table_cell_count_COMPLETE.csv")
