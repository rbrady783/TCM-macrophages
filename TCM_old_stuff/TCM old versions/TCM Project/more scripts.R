# Define the columns to be divided by 100
cols_to_divide <- c("VEGF", "IL.8", "KC.like", "IL.10")

# Divide each of these columns by 100
joined_data[cols_to_divide] <- lapply(joined_data[cols_to_divide], function(x) x / 100)

# Identify numeric columns in joined_data
numeric_vars <- names(joined_data)[sapply(joined_data, is.numeric)]
print(numeric_vars)

# Compute the correlation matrix using pairwise complete observations
cor_matrix <- cor(joined_data[numeric_vars], use = "pairwise.complete.obs",
                  method = "spearman")
print(cor_matrix)

# Install and load GGally if you haven't already:
# install.packages("GGally")
library(GGally)

# Create and display a scatterplot matrix of all numeric variables
ggpairs(joined_data[, numeric_vars])

my_custom_cor <- function(data, mapping, ...) {
  # Extract x and y data from the mapping
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # Compute correlation test
  ct <- cor.test(x, y, method = "spearman")
  corr_coef <- ct$estimate
  p_val <- ct$p.value
  
  # Determine significance stars
  stars <- if (p_val < 0.001) {
    "***"
  } else if (p_val < 0.01) {
    "**"
  } else if (p_val < 0.05) {
    "*"
  } else {
    ""
  }
  
  # Create a label with correlation coefficient and significance
  label <- paste0("r = ", round(corr_coef, 2), "\n", stars)
  
  # Use ggally_text to create the panel; center the text with xP, yP
  ggally_text(label = label, mapping = aes(), xP = 0.5, yP = 0.5, ...) +
    theme_void()
}

# Identify your numeric columns (assuming they are in a data frame called joined_data)
numeric_vars <- names(joined_data)[sapply(joined_data, is.numeric)]

# Create a pairwise plot using ggpairs
pair_plot <- ggpairs(
  joined_data[, numeric_vars],
  upper = list(continuous = my_custom_cor),    # custom correlation with significance
  lower = list(continuous = wrap("smooth", method = "lm", se = FALSE)),
  diag  = list(continuous = "densityDiag")
)

# Display the plot
print(pair_plot)

ggsave("pairwise_plot_large.pdf", pair_plot, width = 12, height = 12)

ggsave("pairwise_plot_powerpoint.png", pair_plot, width = 10, height = 7.5, 
       units = "in", dpi = 600)

ggsave("pairwise_plot_powerpoint.png", pair_plot, width = 10, height = 7.5, units = "in", dpi = 300)
