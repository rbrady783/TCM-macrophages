# =========================================
# Histology-Level Top-Quartile Proportion Test
# =========================================

# 1) INSTALL / LOAD PACKAGES ----
# install.packages(c("readr","dplyr","tidyr"))
library(readr)
library(dplyr)
library(tidyr)

# 2) READ DATA & PREPARE MATRIX ----
df <- read_csv("for_perm_testing.csv", show_col_types = FALSE)

# rows = treatments, cols = cytokines
mat       <- df %>% select(vegf:tgfb) %>% as.matrix()
histology <- df$histology
n_cytok   <- ncol(mat)

# 3) COMPUTE 1st & 3rd QUARTILES FOR EACH CYTOKINE ----
Q1 <- apply(mat, 2, quantile, probs = 0.25, na.rm = TRUE)
Q3 <- apply(mat, 2, quantile, probs = 0.75, na.rm = TRUE)

# 4) FLAG TOP-QUARTILE MEMBERSHIP ----
is_top <- sweep(mat, 2, Q3, FUN = ">=")

# 5) OBSERVED PROPORTIONS BY HISTOLOGY ----
df_h <- df %>%
  mutate(top_count = rowSums(is_top, na.rm = TRUE)) %>%
  group_by(histology) %>%
  summarise(
    n_samples   = n(),
    total_calls = n_samples * n_cytok,
    obs_prop    = sum(top_count, na.rm = TRUE) / total_calls
  ) %>%
  ungroup()

# 6) PERMUTATION TEST FOR PROPORTIONS ----
set.seed(2025)
n_perm   <- 5000L
levels_h <- df_h$histology

# build null distribution of proportions
perm_props <- replicate(n_perm, {
  shuf <- sample(histology)
  sapply(levels_h, function(h) {
    # count hits for this permuted histology group
    cnt <- sum(is_top[shuf == h, ], na.rm = TRUE)
    cnt / (sum(shuf == h) * n_cytok)
  })
})

# 7) EMPIRICAL P-VALUES & MULTIPLE TEST CORRECTION ----
obs_props <- df_h$obs_prop
p_emp     <- (rowSums(perm_props >= obs_props) + 1) / (n_perm + 1)

df_h <- df_h %>%
  mutate(
    p_emp      = p_emp,
    p_emp_adj  = p.adjust(p_emp, method = "BH")
  )

# 8) VIEW & SAVE RESULTS ----
print(df_h)

###################################Bottom quartile:

# --- Flag bottom quartile membership ---
is_bottom <- sweep(mat, 2, Q1, FUN = "<=")

# --- Observed bottom-quartile proportions per histology ---
df_h_bottom <- df %>%
  mutate(bottom_count = rowSums(is_bottom, na.rm = TRUE)) %>%
  group_by(histology) %>%
  summarise(
    n_samples   = n(),
    total_calls = n_samples * n_cytok,
    obs_prop    = sum(bottom_count, na.rm = TRUE) / total_calls
  ) %>%
  ungroup()

# --- Permutation test for bottom-quartile proportions ---
set.seed(2025)
perm_props_bot <- replicate(n_perm, {
  shuf <- sample(histology)
  sapply(levels_h, function(h) {
    cnt <- sum(is_bottom[shuf == h, ], na.rm = TRUE)
    cnt / (sum(shuf == h) * n_cytok)
  })
})

obs_bot_props <- df_h_bottom$obs_prop
p_emp_bot     <- (rowSums(perm_props_bot >= obs_bot_props) + 1) / (n_perm + 1)

df_h_bottom <- df_h_bottom %>%
  mutate(
    p_emp        = p_emp_bot,
    p_emp_adj    = p.adjust(p_emp, method = "BH")
  )

print(df_h_bottom)

