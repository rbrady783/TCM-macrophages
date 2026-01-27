# --- Libraries ---
library(ggplot2)
library(dplyr)
library(nlme)
has_MuMIn <- requireNamespace("MuMIn", quietly = TRUE)

# --- Data & model ---
data <- read.csv("recom_ccl3_tnfa.csv", stringsAsFactors = FALSE)
data$donor   <- as.factor(data$donor)
data$inv_tnfa <- 1 / data$tnfa
final_model <- lme(inv_tnfa ~ dose, random = ~ 1 | donor, data = data, method = "REML")

anova_result <- anova(final_model)
dose_p <- anova_result$`p-value`[2]

# Create p-value text with asterisks for significance
if (is.na(dose_p)) {
  p_text <- "p = NA"
} else if (dose_p < 0.001) {
  p_text <- "p < 0.001***"
} else if (dose_p < 0.01) {
  p_text <- paste0("p = ", signif(dose_p, 3), "**")
} else if (dose_p < 0.05) {
  p_text <- paste0("p = ", signif(dose_p, 3), "*")
} else {
  p_text <- paste0("p = ", signif(dose_p, 3))
}

r2_marginal <- r2_conditional <- NA_real_
if (has_MuMIn) {
  r2_vals <- suppressWarnings(MuMIn::r.squaredGLMM(final_model))
  r2_marginal    <- round(as.numeric(r2_vals[1]), 3)
  r2_conditional <- round(as.numeric(r2_vals[2]), 3)
}

# --- Build plot ---
base_plot <- ggplot(data, aes(x = dose, y = tnfa)) +
  geom_point(aes(color = donor), size = 3.5, alpha = 0.85, stroke = 0.3) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#2166AC", fill = "#4393C3", alpha = 0.3, linewidth = 1.2) +
  scale_y_log10(breaks = c(5, 10, 20, 50), labels = c("5", "10", "20", "50"), limits = c(4, 80)) +
  scale_x_continuous(breaks = c(0, 1.25, 2.5, 5), labels = c("0", "1.25", "2.5", "5")) +
  scale_color_manual(values = c("1" = "#D73027", "2" = "blue", "3" = "purple", "4" = "magenta"),
                     name = "Donor") +
  labs(x = "CCL3 (ng/mL)", y = "TNF-α (pg/mL, log)") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 6), colour = "black"),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 6), colour = "black"),
    axis.text    = element_text(size = 10, colour = "black"),
    
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 9),
    legend.position = c(0.72, 0.16),
    legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
    legend.key.size = unit(0.6, "lines"),
    legend.margin = margin(4, 4, 4, 4),
    
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    
    plot.margin = margin(8, 12, 8, 8)
  )

# --- Simple text annotations ---
p_export <- base_plot +
  annotate("text", x = 0.1, y = 72, 
           label = paste0("Dose effect: ", p_text),
           size = 2.8, hjust = 0, colour = "black") +
  annotate("text", x = 0.1, y = 58, 
           label = paste0("Marginal R² = ", r2_marginal, " | Conditional R² = ", r2_conditional),
           size = 2.6, hjust = 0, colour = "black")

# --- Export ---
fig_width  <- 5
fig_height <- 2.8
fig_dpi    <- 600

ggsave("ccl3_tnfa_professional_halfheight.tiff", p_export,
       width = fig_width, height = fig_height, units = "in",
       dpi = fig_dpi, compression = "lzw", bg = "white")

#ggsave("ccl3_tnfa_professional_halfheight.eps", p_export,
 #      width = fig_width, height = fig_height, units = "in",
  #     device = cairo_ps, bg = "white")
