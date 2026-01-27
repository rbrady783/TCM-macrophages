# Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)

# Read the CSV file
df_wide <- read.csv("breakpoints.csv", header = TRUE, stringsAsFactors = FALSE)

# (Optional) Rename the first column if it has a space; here we change "Cell Line" to "CellLine"
colnames(df_wide)[1] <- "CellLine"

# Remove the "%" sign from all cytokine columns and convert them to numeric.
# This applies the transformation to all columns except the "CellLine" column.
df_wide <- df_wide %>%
  mutate(across(-CellLine, ~ as.numeric(gsub("%", "", .))))

# View the cleaned wide-format data
head(df_wide)

# Pivot the data from wide to long format:
df_long <- df_wide %>%
  pivot_longer(
    cols = -CellLine,
    names_to = "Cytokine",
    values_to = "Value"
  )

# View the first few rows of the long-format data
head(df_long)

# Perform k-means clustering (k = 2) for each cytokine.
# For each Cytokine group, kmeans is run on the "Value" column.
set.seed(123)  # For reproducibility
df_clustered <- df_long %>%
  group_by(Cytokine) %>%
  mutate(Cluster = kmeans(Value, centers = 2)$cluster) %>%
  ungroup()

# View the first few rows of the clustered data
head(df_clustered)

# Plot the clustering results:
# Here we use facet_wrap() so each cytokine is plotted separately.
# We reorder CellLine on the x-axis based on Value for better visualization.
ggplot(df_clustered, aes(x = reorder(CellLine, Value), y = Value, color = as.factor(Cluster))) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering (k = 2) of Cytokine Values by Cell Line",
       x = "Cell Line",
       y = "Cytokine Value",
       color = "Cluster") +
  facet_wrap(~ Cytokine, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Perform k-means clustering with k = 3 for each cytokine.
# We reassign cluster labels based on the ordering of the cluster centers.
set.seed(123)  # For reproducibility

df_clustered <- df_long %>%
  group_by(Cytokine) %>%
  mutate(Cluster = {
    km <- kmeans(Value, centers = 3)
    # 'km$centers' is a 3x1 matrix. Order the clusters by their center (mean)
    center_order <- order(km$centers[,1])
    # Map each original cluster label to its order (1 = lowest, 2 = middle, 3 = highest)
    new_labels <- match(km$cluster, center_order)
    new_labels
  }) %>%
  ungroup() %>%
  # Optionally, convert numeric cluster labels to descriptive factors
  mutate(Cluster = factor(Cluster, levels = c(1, 2, 3), labels = c("Low", "Middle", "High")))

# View the first few rows of the clustered data
head(df_clustered)

# If you only want to use the extreme groups (Low and High) for pairwise analysis,
# you can filter out the "Middle" group:
df_extremes <- df_clustered %>% filter(Cluster != "Middle")

# Optionally, plot the results to visualize the three clusters for each cytokine
ggplot(df_clustered, aes(x = reorder(CellLine, Value), y = Value, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering (k = 3) of Cytokine Values by Cell Line",
       x = "Cell Line",
       y = "Cytokine Value",
       color = "Cluster") +
  facet_wrap(~ Cytokine, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# And you might similarly plot just the extreme groups if needed:
ggplot(df_extremes, aes(x = reorder(CellLine, Value), y = Value, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "Extreme Groups (Low and High) of Cytokine Values",
       x = "Cell Line",
       y = "Cytokine Value",
       color = "Cluster") +
  facet_wrap(~ Cytokine, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Option 1: Pivot wider with both Value and Cluster
df_wide_extremes <- df_extremes %>%
  pivot_wider(
    id_cols = CellLine,
    names_from = Cytokine,
    values_from = c(Value, Cluster)
  )

# View the wide-format data
head(df_wide_extremes)

# Option 2: If you only want the numeric values (without cluster labels), you can do:
df_wide_values <- df_extremes %>%
  select(CellLine, Cytokine, Value) %>%
  pivot_wider(
    id_cols = CellLine,
    names_from = Cytokine,
    values_from = Value
  )

# View the wide-format values data
head(df_wide_values)

# Now, save the desired wide-format data to a CSV file.
# For Option 1 (values and cluster labels):
write.csv(df_wide_extremes, "df_wide_extremes.csv", row.names = FALSE)

# Or, for Option 2 (values only):
write.csv(df_wide_values, "df_wide_values.csv", row.names = FALSE)
