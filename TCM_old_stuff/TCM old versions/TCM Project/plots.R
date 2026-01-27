library(ggpubr)



# BEGF x NF1, p=0.08
p <- ggplot(joined_data, aes(x = NF1, y = VEGF, fill = NF1)) +
  geom_boxplot(outlier.shape = NA) +      # Hide outliers since jitter shows all points
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  labs(title = "VEGF by NF1 Status",
       x = "NF1 Status",
       y = "VEGF") +
  theme_minimal()
 # stat_compare_means(method = "t.test", label = "p.signif")

# Display the plot
print(p)

# Create the boxplot with jittered points
p <- ggplot(joined_data, aes(x = MED12, y = IL.10, fill = MED12)) +
  geom_boxplot(outlier.shape = NA) +      # Hide outliers since jitter shows all points
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  labs(title = "IL-10 by MED12 Status",
       x = "MED12",
       y = "IL_10") +
  theme_minimal() +
  stat_compare_means(method = "t.test", label = "p.signif")

# Display the plot
print(p)


# Install ggforce if needed:
if (!require(ggforce)) {
  install.packages("ggforce")
}
library(ggforce)

p_zoom <- ggplot(joined_data, aes(x = MED12, y = IL.10, fill = MED12)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  labs(title = "IL-10 by MED12 Status", x = "MED12", y = "IL-10") +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_minimal() +
  facet_zoom(ylim = c(0, 100))  # Zoom in on the main cluster

print(p_zoom)
