# =========================================
# Robust Treatment-Level Permutation: Top & Bottom Quartiles
# =========================================

# 0) INSTALL / LOAD PACKAGES ----
# install.packages(c("readr","dplyr"))
library(readr)
library(dplyr)

# 1) READ DATA & SETUP ----
df         <- read_csv("for_perm_testing.csv", show_col_types = FALSE)
mat        <- df %>% select(vegf:tgfb) %>% as.matrix()
rownames(mat) <- df$treatment
n_treat    <- nrow(mat)
n_cytok    <- ncol(mat)

# 2) CALCULATE QUARTILE CUTOFFS ----
Q1 <- apply(mat, 2, quantile, probs = 0.25, na.rm = TRUE)
Q3 <- apply(mat, 2, quantile, probs = 0.75, na.rm = TRUE)

# 3) OBSERVED COUNTS ----
obs_top    <- rowSums(mat >= matrix(Q3, n_treat, n_cytok, byrow = TRUE))
obs_bottom <- rowSums(mat <= matrix(Q1, n_treat, n_cytok, byrow = TRUE))

# quick sanity:
print(data.frame(treatment=rownames(mat), obs_top, obs_bottom))

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
results <- tibble(
  treatment        = rownames(mat),
  top_count        = obs_top,
  p_emp_top        = p_emp_top,
  p_emp_top_adj    = p.adjust(p_emp_top, "BH"),
  bottom_count     = obs_bottom,
  p_emp_bottom     = p_emp_bot,
  p_emp_bottom_adj = p.adjust(p_emp_bot, "BH")
)

# 8) DISPLAY ----
print(results)

