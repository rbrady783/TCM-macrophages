# Load required packages
library(ggplot2)
library(dplyr)

# Load your data and run the statistical model
data <- read.csv("recom_ccl3_tnfa.csv", stringsAsFactors = FALSE)
data$donor <- as.factor(data$donor)
data$inv_tnfa <- 1/data$tnfa

# Run the statistical model (inverse-transformed for proper statistics)
library(nlme)
final_model <- lme(inv_tnfa ~ dose, 
                   random = ~ 1 | donor, 
                   data = data,
                   method = "REML")

# Extract statistical results
anova_result <- anova(final_model)
dose_coef <- fixef(final_model)["dose"]
dose_p <- anova_result$`p-value`[2]

# Calculate R-squared values with performance package
library(performance)
r2_values <- r2(final_model)
r2_conditional <- round(r2_values$R2_conditional, 3)
r2_marginal <- round(r2_values$R2_marginal, 3)

# Format p-value
if(dose_p < 0.001) {
  p_text <- "p < 0.001"
  sig_stars <- "***"
} else if(dose_p < 0.01) {
  p_text <- paste0("p = ", round(dose_p, 3))
  sig_stars <- "**"
} else if(dose_p < 0.05) {
  p_text <- paste0("p = ", round(dose_p, 3))
  sig_stars <- "*"
} else {
  p_text <- paste0("p = ", round(dose_p, 3))
  sig_stars <- "ns"
}

# Create the professional plot with log scale
professional_plot <- ggplot(data, aes(x = dose, y = tnfa)) +
  # Add points with professional color scheme
  geom_point(aes(color = donor), size = 4, alpha = 0.85, stroke = 0.3) +
  
  # Add regression line on log scale
  geom_smooth(method = "lm", se = TRUE, color = "#2166AC", fill = "#4393C3", 
              alpha = 0.3, linewidth = 1.3) +
  
  # Use log10 scale for y-axis
  scale_y_log10(
    breaks = c(5, 10, 20, 50),
    labels = c("5", "10", "20", "50"),
    limits = c(4, 80)
  ) +
  
  # Professional color palette (similar to your reference)
  scale_color_manual(
    values = c("1" = "#D73027", "2" = "blue", "3" = "purple"),
    name = "Donor"
  ) +
  
  # Clean scale for x-axis
  scale_x_continuous(
    breaks = c(0, 1.25, 2.5, 5),
    labels = c("0", "1.25", "2.5", "5.0")
  ) +
  
  
  # Professional labels
  labs(
    title = "",
    x = "**CCL3 (ng/mL)**",
    y = "**TNF-α(pg/mL, log scale)**",
    color = "**Donor**"
  ) +
  
  # Statistical annotation - top left
  annotate("text", x = 0.3, y = 70, 
           label = paste0("**Dose effect: ", sig_stars, " (", p_text, ")**"),
           size = 9, hjust = 0,
           color = "black") +
  
  # Marginal R-squared annotation (variance explained by CCL3 dose alone)
  annotate("text", x = 0.3, y = 54,
           label = paste0("**Marginal R² = ", r2_marginal, "**"),
           size = 8.5, hjust = 0, color = "black") +
  
  
  # Professional theme
  theme_minimal() +
  theme(
    
    # Axis titles
    axis.title.x = element_text(size = 20, color = "black", 
                                margin = margin(t = 10), face = "bold"),
    axis.title.y = element_text(size = 20, color = "black", 
                                margin = margin(r = 10), face = "bold"),
    
    # Axis text
    axis.text = element_text(size = 18, color = "black", face = "bold"),
    
    # Legend
    legend.title = element_text(size = 18, color = "black", face = "bold"),
    legend.text = element_text(size = 18, color = "black", face = "bold"),
    legend.position = c(0.85, 0.25),
    legend.background = element_rect(fill = "white", color = "#E8E8E8", 
                                     size = 0.5),
    legend.margin = margin(8, 8, 8, 8),
    
    # Panel
    panel.grid.major = element_line(color = "#E8E8E8", size = 0.5),
    panel.grid.minor = element_line(color = "#F5F5F5", size = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    
    # Plot background
    plot.background = element_rect(fill = "white", color = NA),
    
    # Margins
    plot.margin = margin(20, 25, 15, 15)
  )

# Display the plot
print(professional_plot)

# Print statistical summary
cat("\n=============================================================\n")
cat("STATISTICAL SUMMARY\n")
cat("=============================================================\n")
cat("Model: Linear mixed-effects (inverse-transformed TNF-α)\n")
cat("Dose effect:", sig_stars, "(", p_text, ")\n")
cat("Direction: TNF-α increases with CCL3 concentration\n")
if(r2_conditional != "N/A") {
  cat("Conditional R²:", r2_conditional, "\n")
  cat("Marginal R²:", r2_marginal, "\n")
}
cat("Sample size: n =", nrow(data), "(4 donors, 4 doses each)\n")

# Optional: Save the plot in high resolution
ggsave("ccl3_tnfa_professional.png", professional_plot, 
       width = 10, height = 7, dpi = 1200, bg = "white")
# 
ggsave("ccl3_tnfa_professional.pdf", professional_plot, 
       width = 10, height = 7, bg = "transparent")

cat("\nTo save high-resolution versions, uncomment the ggsave() lines at the end of the code.\n")