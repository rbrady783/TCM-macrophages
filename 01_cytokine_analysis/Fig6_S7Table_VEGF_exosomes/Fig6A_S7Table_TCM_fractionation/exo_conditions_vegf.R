# Load required libraries
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)

# Create the data manually based on the CSV structure
# First, let me create vectors for each condition

# Whole samples (no suffix)
whole_samples <- c(
  J1 = 54.19390253, M1 = 49.83719781, N1 = 109.0414019, 
  O2_1 = 76.91210763, O8_1 = 78.67316244, P1 = 47.68476734, S1 = 123.5945402,
  J2 = 128.0901019, M2 = 71.09661839, N2 = 283.6162878,
  O2_2 = 210.753584, O8_2 = 326.10064, P2 = 99.12516953, S2 = 312.3227859,
  J3 = 361.8116509, M3 = 278.3863215, N3 = 336.120923,
  O2_3 = 722.871425, O8_3 = 883.7620791, P3 = 421.489871, S3 = 641.3944948
)

# Exosome-depleted samples (x suffix)
depleted_samples <- c(
  J1 = 268.7229365, M1 = 151.6670368, N1 = 168.4621522,
  O2_1 = 232.1787045, O8_1 = 149.6741848, P1 = 193.9028531, S1 = 309.2764953,
  J2 = 256.1851066, M2 = 212.1696414, N2 = 252.5192645,
  O2_2 = 152.3324179, O8_2 = 168.4621522, P2 = 134.5638016, S2 = 340.7650239,
  J3 = 499.4592159, M3 = 594.477161, N3 = 535.5603246,
  O2_3 = 471.2640577, O8_3 = 539.7893778, P3 = 522.9118873, S3 = 445.8358181
)

# Exosome-only samples (o suffix)
exosome_samples <- c(
  J1 = 715.7118863, M1 = 826.3661872, N1 = 786.1006462,
  O2_1 = 705.8903097, O8_1 = 1446.536776, P1 = 781.5501243, S1 = 1020.855785,
  J2 = 701.4347295, M2 = 562.736211, N2 = 493.6286596,
  O2_2 = 882.8303601, O8_2 = 1282.954119, P2 = 1216.071781, S2 = 1310.108418,
  J3 = 1882.182783, M3 = 1657.195575, N3 = 1320.19455,
  O2_3 = 2268.599089, O8_3 = 2187.656742, P3 = 2412.862539, S3 = 2354.704069
)

# Create a tidy data frame
create_condition_df <- function(values, condition) {
  data.frame(
    Sample_ID = names(values),
    VEGF = as.numeric(values),
    Condition = condition,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Cell_Line = case_when(
        grepl("^J", Sample_ID) ~ "Jones",
        grepl("^M", Sample_ID) ~ "MacKinley", 
        grepl("^N", Sample_ID) ~ "Nike",
        grepl("^O2", Sample_ID) ~ "OS2.4",
        grepl("^O8", Sample_ID) ~ "OSA8",
        grepl("^P", Sample_ID) ~ "Parks",
        grepl("^S", Sample_ID) ~ "STS-1",
        TRUE ~ "Unknown"
      ),
      Donor = case_when(
        grepl("1$|1_", Sample_ID) ~ "Donor_1",
        grepl("2$|2_", Sample_ID) ~ "Donor_2", 
        grepl("3$|3_", Sample_ID) ~ "Donor_3",
        TRUE ~ "Unknown"
      )
    )
}

# Combine all conditions
full_data <- bind_rows(
  create_condition_df(whole_samples, "Whole TCM"),
  create_condition_df(depleted_samples, "Exosome Depleted"),
  create_condition_df(exosome_samples, "Exosome Only")
) %>%
  mutate(
    Condition = factor(Condition, levels = c("Whole TCM",
                                             "Exosome Depleted",
                                             "Exosome Only")),
    Cell_Line = factor(Cell_Line),
    Donor = factor(Donor)
  )

# Check the data structure
print("Data structure:")
head(full_data, 15)

print("\nSample sizes:")
table(full_data$Condition, full_data$Cell_Line)


# Mixed-effects model comparing conditions (donor as random effect)
model_donor <- lmer(VEGF ~ Condition + (1|Donor), data = full_data)

cat("\n=== MIXED-EFFECTS MODEL RESULTS ===\n")
summary(model_donor)
cat("\n=== ANOVA ===\n")
anova(model_donor)

# Post-hoc pairwise comparisons
library(emmeans)
cat("\n=== PAIRWISE COMPARISONS ===\n")
emm <- emmeans(model_donor, ~ Condition)
pairwise_results <- pairs(emm, adjust = "tukey")
print(pairwise_results)

# Extract p-values for easy reference
pairwise_summary <- summary(pairwise_results)
cat("\n=== P-VALUES SUMMARY ===\n")
for(i in 1:nrow(pairwise_summary)) {
  comparison <- pairwise_summary$contrast[i]
  p_value <- pairwise_summary$p.value[i]
  significance <- ifelse(p_value < 0.001, "***", 
                         ifelse(p_value < 0.01, "**",
                                ifelse(p_value < 0.05, "*", "ns")))
  cat(sprintf("%-30s p = %.4f %s\n", comparison, p_value, significance))
}

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggsignif)

# Use your existing full_data from the previous analysis
# Clean, minimal plot following good data viz principles
p1 <- ggplot(full_data, aes(x = Condition, y = VEGF, fill = Condition)) +
  # Clean bars
  stat_summary(fun = "mean", geom = "col", alpha = 0.8, width = 0.6, 
               color = "black", linewidth = 0.5) +
  # Individual points
  geom_point(position = position_jitter(width = 0.15, seed = 123), 
             size = 1.5, alpha = 0.6, color = "black") +
  # Error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.25, linewidth = 0.7, color = "black") +
  # Significance bars - both comparisons
  geom_signif(comparisons = list(c("Whole TCM", "Exosome Only")), 
              annotations = "***", y_position = 2400, tip_length = 0.02, 
              textsize = 4, vjust = 0.4) +
  geom_signif(comparisons = list(c("Exosome Depleted", "Exosome Only")), 
              annotations = "***", y_position = 2600, tip_length = 0.02, 
              textsize = 4, vjust = 0.4) +
  # Colors
  scale_fill_manual(values = c("Whole TCM" = "#E91E63", 
                               "Exosome Depleted" = "#696969", 
                               "Exosome Only" = "#00BCD4")) +
  # Clean labels
  scale_x_discrete(labels = c("Whole\nTCM", "Exosome\nDepleted", "Exosome\nOnly")) +
  labs(y = "VEGF (pg/ml)", x = NULL) +
  # Clean theme
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # Reasonable axis limits with room for significance bars
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), 
                     breaks = scales::pretty_breaks(n = 6))

print(p1)

# Save plot
ggsave("exosome_vegf_clean.pdf", plot = p1, 
       width = 3, height = 6, units = "in", dpi = 1200)

# Simple statistical annotation in caption or text instead of cluttered brackets
cat("\nStatistical results:")
cat("\nExosome Only vs Whole TCM: p < 0.001")
cat("\nExosome Only vs Exosome Depleted: p < 0.001") 
cat("\nWhole TCM vs Exosome Depleted: p = 0.87 (ns)")

print(pairwise_results)
