# Install packages if you haven't already
# install.packages("readxl")
# install.packages("dplyr")

library(readxl)
library(dplyr)

# Read in the two Excel sheets.
# Update the file paths to where your Excel files are located.
z_scores <- read.csv("z_scores.csv")
mutations <- read.csv("mutation.csv")

# Optional: inspect the data
head(z_scores)
head(mutations)

# Check that the join column ("cell.line") exists in both data frames.
# You can also check column names:
names(z_scores)
names(mutations)

# If the join column might have extra spaces or different capitalization, clean it:
z_scores$cell.line <- trimws(z_scores$cell.line)
mutations$cell.line <- trimws(mutations$cell.line)

# Now perform an inner join on the "cell.line" column.
# This will keep only the rows that are present in both datasets.
joined_data <- inner_join(z_scores, mutations, by = "Cell.line")

joined_data$X <- NULL

# View the resulting joined dataset.
print(joined_data)
head(joined_data)

#convert to binary
joined_data <- joined_data %>%
  mutate_at(vars(12:ncol(joined_data)),
            ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes")))

#Handle ERK/AKT
# Verify the unique values are now cleaned:
unique(joined_data$`pAKT`)
unique(joined_data$`pERK`)

# Define the vector of column names
status_columns <- c("pERK", "pAKT")

# Convert the cleaned columns to factors with the specified levels and labels
joined_data[status_columns] <- lapply(joined_data[status_columns], function(x) {
  factor(x, levels = c("l","c","a"), labels = c("Low", "Constitutive", "Activated"))
})

# Check the conversion
str(joined_data[status_columns])

column_classes <- sapply(joined_data, class)
print(column_classes)
str(joined_data)

# Define the columns to convert
percentage_columns <- c("VEGF", "IL.8", "KC.like", "IL.10", "TNF.a")

# Convert each of these columns: remove extra spaces, remove "%" and then convert to numeric
joined_data[percentage_columns] <- lapply(joined_data[percentage_columns], function(x) {
  as.numeric(gsub("%", "", trimws(x)))
})

# Set the first column as row names
rownames(joined_data) <- joined_data[, 1]

# Then remove the first column from the data frame
joined_data <- joined_data[, -1]

library(ggplot2)


write.csv(joined_data, "joined_data.csv", row.names = TRUE)



