# Read in the CSV file
df <- read.csv("ranks.csv")

# Convert the "Cell.Line" column to character
df$Cell.Line <- as.character(df$Cell.Line)

# Set the number of permutations
set.seed(123)  # for reproducibility
nPerm <- 10000

# Get the unique cytokines and unique cell lines
cytokines <- unique(df$Cytokine)
cell_lines <- unique(df$Cell.Line)

# Create a results data frame to store the observed counts and p-values
results <- data.frame(
  Cell.Line = cell_lines,
  observed_top5 = NA,
  pvalue_top5 = NA,
  observed_bottom5 = NA,
  pvalue_bottom5 = NA,
  stringsAsFactors = FALSE
)

# Create an empty list to store the permutation distribution for top 5 counts
perm_list <- list()

# Create a progress bar for the outer loop (cell lines)
pb <- txtProgressBar(min = 0, max = length(cell_lines), style = 3)

# Loop over each cell line
for (i in seq_along(cell_lines)) {
  
  cl <- cell_lines[i]
  
  # Calculate the observed statistic for this cell line
  obs_top5 <- sum(df$Cell.Line == cl & df$Rank <= 5)
  obs_bottom5 <- sum(df$Cell.Line == cl & df$Rank >= 21)  # bottom 5 if there are 25 ranks
  
  # Vectors to store the test statistic from each permutation
  perm_top5 <- numeric(nPerm)
  perm_bottom5 <- numeric(nPerm)
  
  # Loop over the number of permutations
  for (j in 1:nPerm) {
    # For each cytokine, randomly shuffle the cell line assignments
    permuted_df <- do.call(rbind, lapply(cytokines, function(cyt) {
      sub_df <- subset(df, Cytokine == cyt)
      # Shuffle the cell line labels (keeping the ranking values fixed)
      sub_df$Cell.Line <- sample(sub_df$Cell.Line)
      return(sub_df)
    }))
    
    # Count how many times 'cl' appears in the top 5 and bottom 5 in the permuted data
    perm_top5[j] <- sum(permuted_df$Cell.Line == cl & permuted_df$Rank <= 5)
    perm_bottom5[j] <- sum(permuted_df$Cell.Line == cl & permuted_df$Rank >= 21)
  }
  
  # Calculate one-sided p-values:
  # For top 5: probability that the permuted count is greater than or equal to the observed count.
  p_val_top <- mean(perm_top5 >= obs_top5)
  # For bottom 5:
  p_val_bottom <- mean(perm_bottom5 >= obs_bottom5)
  
  # Save the results for this cell line in the results data frame
  results[results$Cell.Line == cl, "observed_top5"] <- obs_top5
  results[results$Cell.Line == cl, "pvalue_top5"] <- p_val_top
  results[results$Cell.Line == cl, "observed_bottom5"] <- obs_bottom5
  results[results$Cell.Line == cl, "pvalue_bottom5"] <- p_val_bottom
  
  # Save the permutation distribution for this cell line (top 5 counts only)
  perm_list[[cl]] <- perm_top5
  
  # Update the progress bar after processing each cell line
  setTxtProgressBar(pb, i)
}

# Close the progress bar
close(pb)

# Show the results data frame
print(results)

library(dplyr)

# Combine the permutation results into a data frame
all_perm_data <- do.call(rbind, lapply(names(perm_list), function(cl) {
  data.frame(Cell.Line = cl, permCount = perm_list[[cl]])
}))

# Define the selected cell lines
selected_cell_lines <- c("Nike", "CLL1390", "1771", "Yamane", "Vogel")

# Filter the permutation data and the observed results to only the selected cell lines
all_perm_data_subset <- all_perm_data %>% 
  filter(Cell.Line %in% selected_cell_lines)

observed_df_subset <- results %>% 
  filter(Cell.Line %in% selected_cell_lines)

# Define the selected cell lines (if you want to show both, for now we do this for top 5)
selected_cell_lines <- c("Nike", "CLL1390", "1771", "Yamane", "Vogel")

# Filter the results for top 5 permutation data.
# (Assume perm_list_top was created analogously to our previous code,
#  and 'results' contains the observed_top5 and pvalue_top5.)
all_perm_data_top <- do.call(rbind, lapply(selected_cell_lines, function(cl) {
  data.frame(Cell.Line = cl, permCount = perm_list[[cl]])  # perm_list here is the top 5 list
}))

observed_df_top <- results %>% 
  filter(Cell.Line %in% selected_cell_lines)

# Create the top 5 plot:
p_top <- ggplot(all_perm_data_top, aes(x = permCount)) +
  geom_histogram(color = "black", fill = "lightblue", bins = 30) +
  geom_vline(data = observed_df_top, aes(xintercept = observed_top5), 
             color = "red", size = 1.2) +
  facet_wrap(~ Cell.Line, scales = "free_y") +
  labs(title = "Permutation Distribution (Top 5)",
       x = "Count in Top 5",
       y = "Frequency") +
  theme_minimal() +
  geom_text(data = observed_df_top, 
            aes(x = observed_top5, 
                label = paste("p =", signif(pvalue_top5, 3))),
            y = Inf, vjust = 2, color = "red", size = 4)

print(p_top)

# Define the cell line for the top-5 plot (only Nike)
selected_top <- "Nike"

# Extract the permutation data for Nike from the top-5 list.
# (This assumes you have stored it in the list 'perm_list' from your permutation code.)
all_perm_data_top <- data.frame(
  Cell.Line = selected_top, 
  permCount = perm_list[[selected_top]]
)

# Extract Nike’s observed top-5 count and p-value from your results data frame.
observed_df_top <- results %>% 
  filter(Cell.Line == selected_top)

# Create the top-5 plot for Nike
p_top <- ggplot(all_perm_data_top, aes(x = permCount)) +
  geom_histogram(color = "black", fill = "lightblue", bins = 30) +
  geom_vline(data = observed_df_top, 
             aes(xintercept = observed_top5), 
             color = "red", size = 1.2) +
  labs(title = paste("Permutation Distribution (Top 5) for", selected_top),
       x = "Count in Top 5",
       y = "Frequency") +
  theme_minimal() +
  geom_text(data = observed_df_top, 
            aes(x = observed_top5, label = paste("p =", signif(pvalue_top5, 3))),
            y = Inf, vjust = 2, color = "red", size = 4,
            position = position_nudge(x = 0.9))

# Display the plot
print(p_top)
