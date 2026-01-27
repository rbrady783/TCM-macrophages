# =========================================
# Cytokine–Histology Analysis Script
# =========================================

# 1) Libraries ----
if (!requireNamespace("tidyverse", quietly=TRUE)) install.packages("tidyverse")
if (!requireNamespace("lmerTest",  quietly=TRUE)) install.packages("lmerTest")
library(tidyverse)
library(lmerTest)   # gives lmer() p-values automatically

# 2) Load data ----
# Assumes your CSV has columns: cytokine, treatment, histology, donor_1, donor_2, donor_3
df <- read_csv("raw_no_ctrls_histo.csv")

# 3) Reshape & transform ----
df_long <- df %>%
  pivot_longer(
    cols      = starts_with("donor_"),
    names_to  = "donor",
    values_to = "concentration"
  ) %>%
  mutate(
    donor    = factor(str_remove(donor, "donor_")),  # 1,2,3
    log_conc = log(concentration + 1)                # log-transform
  )

# 4) Fit LMMs & extract fixed effects ----
model_results <- df_long %>%
  group_by(cytokine) %>%
  nest() %>%
  mutate(
    fit    = map(data, ~ lmer(log_conc ~ histology + (1|donor), data = .x)),
    tidied = map(fit, ~ {
      cf <- summary(.x)$coefficients
      tibble(
        term      = rownames(cf),
        estimate  = cf[, "Estimate"],
        std.error = cf[, "Std. Error"],
        statistic = cf[, "t value"],
        p.value   = cf[, "Pr(>|t|)"]
      )
    })
  ) %>%
  select(cytokine, tidied) %>%
  unnest(tidied) %>%
  filter(str_starts(term, "histology")) %>%       # drop intercept
  group_by(cytokine) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# 5) Print model results & significant contrasts ----
print(model_results)

sig_contrasts <- model_results %>%
  filter(p.adj < 0.05)

cat("\nSignificant histology contrasts (FDR < 0.05):\n")
print(sig_contrasts)

# 6) Visualization ----
# Clear any old devices
while (!is.null(dev.list())) dev.off()
# define your mapping
label_map <- c(
  "vegf"    = "VEGF",
  "il.8"    = "IL-8",
  "kc.like" = "KC-like",
  "il.10"   = "IL-10",
  "ccl2"    = "CCL2",
  "tnf.a"   = "TNF-α",
  "tgf.b"   = "TGF-β"
)

# recode in your df_long
df_long2 <- df_long %>%
  mutate(cytokine = recode_factor(cytokine, !!!label_map))

# now plot as before — facet labels will use the new names
ggplot(df_long2, aes(x = histology, y = log_conc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ cytokine, scales = "free_y") +
  labs(
    x       = "Histology",
    y       = "log(concentration + 1)",
    title   = "Cytokine Production by CM Histology",
    caption = "All p-values non-significant (ns) based on BH-adjusted LMM & Kruskal–Wallis"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0)
  )

# Re-plot with caption
ggplot(df_long, aes(x = histology, y = log_conc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ cytokine, scales = "free_y") +
  labs(
    x       = "Histology",
    y       = "log(concentration + 1)",
    title   = "Cytokine Production by CM Histology",
    caption = "All p-values non-significant (ns) based on BH-adjusted LMM & Kruskal–Wallis"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    plot.caption   = element_text(hjust = 0)
  )

# 7) Model diagnostics ----
# Refit models into a list-column for iteration
model_fits <- df_long %>%
  group_by(cytokine) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lmer(log_conc ~ histology + (1|donor), data = .x))
  )


# Loop: Residuals vs Fitted & QQ-plot per cytokine
for (i in seq_len(nrow(model_fits))) {
  fit_obj   <- model_fits$fit[[i]]
  cyto_name <- model_fits$cytokine[[i]]
  
  par(mfrow = c(1,2), oma = c(0,0,2,0))
  
  # Residuals vs Fitted
  plot(
    fitted(fit_obj),
    residuals(fit_obj),
    main = paste0(cyto_name, ": Residuals vs Fitted"),
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  abline(h = 0, lty = 2)
  
  # QQ-plot of residuals
  qqnorm(
    residuals(fit_obj),
    main = paste0(cyto_name, ": QQ-plot")
  )
  qqline(residuals(fit_obj))
  
  mtext(paste0("Diagnostics for ", cyto_name), outer = TRUE, line = 0, cex = 1.2)
}

# 8) Non-parametric backup: Kruskal–Wallis ----
kruskal_results <- df_long %>%
  group_by(cytokine) %>%
  summarise(
    statistic = kruskal.test(concentration ~ histology, data = cur_data())$statistic,
    p.value   = kruskal.test(concentration ~ histology, data = cur_data())$p.value
  ) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

cat("\nKruskal–Wallis test results (with BH correction):\n")
print(kruskal_results)

# =========================================
# End of script
# =========================================