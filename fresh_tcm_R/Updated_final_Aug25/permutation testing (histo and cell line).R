# =========================================
# Histology-Level Top-Quartile Proportion Test (4-Category System)
# =========================================

# 1) INSTALL / LOAD PACKAGES ----
# install.packages(c("readr","dplyr","tidyr"))
library(readr)
library(dplyr)
library(tidyr)

# 2) READ DATA & PREPARE MATRIX ----
df <- read_csv("for_perm_testing_updated.csv", show_col_types = FALSE)

# rows = cell lines (treatments), cols = cytokines
mat       <- df %>% select(vegf:tgfb) %>% as.matrix()
histology_4 <- df$histo_4  # Using 4-category system
histology_5 <- df$histo_5  # Using 5-category system
n_cytok   <- ncol(mat)

# 3) COMPUTE 1st & 3rd QUARTILES FOR EACH CYTOKINE ----
Q1 <- apply(mat, 2, quantile, probs = 0.25, na.rm = TRUE)
Q3 <- apply(mat, 2, quantile, probs = 0.75, na.rm = TRUE)

# =========================================
# ANALYSIS WITH 4-CATEGORY HISTOLOGY SYSTEM
# =========================================

# 4) FLAG TOP-QUARTILE MEMBERSHIP ----
is_top <- sweep(mat, 2, Q3, FUN = ">=")

# 5) OBSERVED PROPORTIONS BY HISTOLOGY (4-category) ----
df_h4 <- df %>%
  mutate(top_count = rowSums(is_top, na.rm = TRUE)) %>%
  group_by(histo_4) %>%
  summarise(
    n_samples   = n(),
    total_calls = n_samples * n_cytok,
    obs_prop    = sum(top_count, na.rm = TRUE) / total_calls,
    .groups = 'drop'
  )

# 6) PERMUTATION TEST FOR PROPORTIONS (4-category) ----
set.seed(2025)
n_perm   <- 5000L
levels_h4 <- df_h4$histo_4

# build null distribution of proportions
perm_props_4 <- replicate(n_perm, {
  shuf <- sample(histology_4)
  sapply(levels_h4, function(h) {
    # count hits for this permuted histology group
    cnt <- sum(is_top[shuf == h, ], na.rm = TRUE)
    cnt / (sum(shuf == h) * n_cytok)
  })
})

# 7) EMPIRICAL P-VALUES & MULTIPLE TEST CORRECTION (4-category) ----
obs_props_4 <- df_h4$obs_prop
p_emp_4     <- (rowSums(perm_props_4 >= obs_props_4) + 1) / (n_perm + 1)

df_h4 <- df_h4 %>%
  mutate(
    p_emp      = p_emp_4,
    p_emp_adj  = p.adjust(p_emp_4, method = "BH")
  )

# 8) VIEW RESULTS (4-category TOP quartile) ----
cat("=== 4-CATEGORY HISTOLOGY: TOP QUARTILE RESULTS ===\n")
print(df_h4)

# =========================================
# BOTTOM QUARTILE ANALYSIS (4-category)
# =========================================

# Flag bottom quartile membership
is_bottom <- sweep(mat, 2, Q1, FUN = "<=")

# Observed bottom-quartile proportions per histology (4-category)
df_h4_bottom <- df %>%
  mutate(bottom_count = rowSums(is_bottom, na.rm = TRUE)) %>%
  group_by(histo_4) %>%
  summarise(
    n_samples   = n(),
    total_calls = n_samples * n_cytok,
    obs_prop    = sum(bottom_count, na.rm = TRUE) / total_calls,
    .groups = 'drop'
  )

# Permutation test for bottom-quartile proportions (4-category)
set.seed(2025)
perm_props_bot_4 <- replicate(n_perm, {
  shuf <- sample(histology_4)
  sapply(levels_h4, function(h) {
    cnt <- sum(is_bottom[shuf == h, ], na.rm = TRUE)
    cnt / (sum(shuf == h) * n_cytok)
  })
})

obs_bot_props_4 <- df_h4_bottom$obs_prop
p_emp_bot_4     <- (rowSums(perm_props_bot_4 >= obs_bot_props_4) + 1) / (n_perm + 1)

df_h4_bottom <- df_h4_bottom %>%
  mutate(
    p_emp        = p_emp_bot_4,
    p_emp_adj    = p.adjust(p_emp_bot_4, method = "BH")
  )

cat("\n=== 4-CATEGORY HISTOLOGY: BOTTOM QUARTILE RESULTS ===\n")
print(df_h4_bottom)

# =========================================
# ANALYSIS WITH 5-CATEGORY HISTOLOGY SYSTEM
# =========================================

# TOP QUARTILE (5-category)
df_h5 <- df %>%
  mutate(top_count = rowSums(is_top, na.rm = TRUE)) %>%
  group_by(histo_5) %>%
  summarise(
    n_samples   = n(),
    total_calls = n_samples * n_cytok,
    obs_prop    = sum(top_count, na.rm = TRUE) / total_calls,
    .groups = 'drop'
  )

# Permutation test (5-category)
set.seed(2025)
levels_h5 <- df_h5$histo_5

perm_props_5 <- replicate(n_perm, {
  shuf <- sample(histology_5)
  sapply(levels_h5, function(h) {
    cnt <- sum(is_top[shuf == h, ], na.rm = TRUE)
    cnt / (sum(shuf == h) * n_cytok)
  })
})

obs_props_5 <- df_h5$obs_prop
p_emp_5     <- (rowSums(perm_props_5 >= obs_props_5) + 1) / (n_perm + 1)

df_h5 <- df_h5 %>%
  mutate(
    p_emp      = p_emp_5,
    p_emp_adj  = p.adjust(p_emp_5, method = "BH")
  )

cat("\n=== 5-CATEGORY HISTOLOGY: TOP QUARTILE RESULTS ===\n")
print(df_h5)

# BOTTOM QUARTILE (5-category)
df_h5_bottom <- df %>%
  mutate(bottom_count = rowSums(is_bottom, na.rm = TRUE)) %>%
  group_by(histo_5) %>%
  summarise(
    n_samples   = n(),
    total_calls = n_samples * n_cytok,
    obs_prop    = sum(bottom_count, na.rm = TRUE) / total_calls,
    .groups = 'drop'
  )

set.seed(2025)
perm_props_bot_5 <- replicate(n_perm, {
  shuf <- sample(histology_5)
  sapply(levels_h5, function(h) {
    cnt <- sum(is_bottom[shuf == h, ], na.rm = TRUE)
    cnt / (sum(shuf == h) * n_cytok)
  })
})

obs_bot_props_5 <- df_h5_bottom$obs_prop
p_emp_bot_5     <- (rowSums(perm_props_bot_5 >= obs_bot_props_5) + 1) / (n_perm + 1)

df_h5_bottom <- df_h5_bottom %>%
  mutate(
    p_emp        = p_emp_bot_5,
    p_emp_adj    = p.adjust(p_emp_bot_5, method = "BH")
  )

cat("\n=== 5-CATEGORY HISTOLOGY: BOTTOM QUARTILE RESULTS ===\n")
print(df_h5_bottom)

# =========================================
# Cell Line-Level Permutation: Top & Bottom Quartiles
# =========================================

# 1) READ DATA & SETUP ----
mat        <- df %>% select(vegf:tgfb) %>% as.matrix()
rownames(mat) <- df$treatment  # Cell line names
n_treat    <- nrow(mat)
n_cytok    <- ncol(mat)

# 2) CALCULATE QUARTILE CUTOFFS (already done above) ----
# Q1 and Q3 already calculated

# 3) OBSERVED COUNTS ----
obs_top    <- rowSums(mat >= matrix(Q3, n_treat, n_cytok, byrow = TRUE))
obs_bottom <- rowSums(mat <= matrix(Q1, n_treat, n_cytok, byrow = TRUE))

# Quick sanity check:
cat("\n=== OBSERVED COUNTS BY CELL LINE ===\n")
print(data.frame(cell_line = rownames(mat), obs_top, obs_bottom))

# 4) PERMUTATION PARAMETERS ----
set.seed(202559)
n_perm <- 20000L

# pre-allocate null matrices
null_top    <- matrix(0L, n_treat, n_perm)
null_bottom <- matrix(0L, n_treat, n_perm)

# 5) RUN PERMUTATIONS ----
for (i in seq_len(n_perm)) {
  # permute each column independently
  permuted <- apply(mat, 2, sample)
  
  # recount top/bottom
  null_top[, i]    <- rowSums(permuted >= matrix(Q3, n_treat, n_cytok, byrow = TRUE))
  null_bottom[, i] <- rowSums(permuted <= matrix(Q1, n_treat, n_cytok, byrow = TRUE))
}

# 6) EMPIRICAL P-VALUES ----
p_emp_top <- sapply(seq_len(n_treat), function(r) {
  (sum(null_top[r, ] >= obs_top[r]) + 1) / (n_perm + 1)
})

p_emp_bot <- sapply(seq_len(n_treat), function(r) {
  (sum(null_bottom[r, ] >= obs_bottom[r]) + 1) / (n_perm + 1)
})

# 7) ASSEMBLE & ADJUST ----
results <- data.frame(
  cell_line        = rownames(mat),
  top_count        = obs_top,
  p_emp_top        = p_emp_top,
  p_emp_top_adj    = p.adjust(p_emp_top, "BH"),
  bottom_count     = obs_bottom,
  p_emp_bottom     = p_emp_bot,
  p_emp_bottom_adj = p.adjust(p_emp_bot, "BH")
)

# 8) DISPLAY CELL LINE RESULTS ----
cat("\n=== CELL LINE PERMUTATION TEST RESULTS ===\n")
print(results)

# Summary of significant results
cat("\n=== SUMMARY OF SIGNIFICANT RESULTS ===\n")
cat("Cell lines with significant over-representation in TOP quartile (p_adj < 0.05):\n")
sig_top <- results[results$p_emp_top_adj < 0.05, ]
if(nrow(sig_top) > 0) {
  print(sig_top[, c("cell_line", "top_count", "p_emp_top", "p_emp_top_adj")])
} else {
  cat("None found.\n")
}

cat("\nCell lines with significant over-representation in BOTTOM quartile (p_adj < 0.05):\n")
sig_bottom <- results[results$p_emp_bottom_adj < 0.05, ]
if(nrow(sig_bottom) > 0) {
  print(sig_bottom[, c("cell_line", "bottom_count", "p_emp_bottom", "p_emp_bottom_adj")])
} else {
  cat("None found.\n")
}

# Create a comprehensive summary file
summary_stats <- data.frame(
  Analysis = c("4-cat Histology Top", "4-cat Histology Bottom", 
               "5-cat Histology Top", "5-cat Histology Bottom",
               "Cell Lines Top", "Cell Lines Bottom"),
  N_Groups = c(nrow(df_h4), nrow(df_h4_bottom), 
               nrow(df_h5), nrow(df_h5_bottom),
               nrow(results), nrow(results)),
  N_Significant_p05 = c(
    sum(df_h4$p_emp_adj < 0.05, na.rm = TRUE),
    sum(df_h4_bottom$p_emp_adj < 0.05, na.rm = TRUE),
    sum(df_h5$p_emp_adj < 0.05, na.rm = TRUE),
    sum(df_h5_bottom$p_emp_adj < 0.05, na.rm = TRUE),
    sum(results$p_emp_top_adj < 0.05, na.rm = TRUE),
    sum(results$p_emp_bottom_adj < 0.05, na.rm = TRUE)
  ),
  N_Marginal_p10 = c(
    sum(df_h4$p_emp_adj < 0.10 & df_h4$p_emp_adj >= 0.05, na.rm = TRUE),
    sum(df_h4_bottom$p_emp_adj < 0.10 & df_h4_bottom$p_emp_adj >= 0.05, na.rm = TRUE),
    sum(df_h5$p_emp_adj < 0.10 & df_h5$p_emp_adj >= 0.05, na.rm = TRUE),
    sum(df_h5_bottom$p_emp_adj < 0.10 & df_h5_bottom$p_emp_adj >= 0.05, na.rm = TRUE),
    sum(results$p_emp_top_adj < 0.10 & results$p_emp_top_adj >= 0.05, na.rm = TRUE),
    sum(results$p_emp_bottom_adj < 0.10 & results$p_emp_bottom_adj >= 0.05, na.rm = TRUE)
  ),
  Min_P_Value = c(
    min(df_h4$p_emp_adj, na.rm = TRUE),
    min(df_h4_bottom$p_emp_adj, na.rm = TRUE),
    min(df_h5$p_emp_adj, na.rm = TRUE),
    min(df_h5_bottom$p_emp_adj, na.rm = TRUE),
    min(results$p_emp_top_adj, na.rm = TRUE),
    min(results$p_emp_bottom_adj, na.rm = TRUE)
  )
)

write_csv(summary_stats, "permutation_test_summary.csv")


