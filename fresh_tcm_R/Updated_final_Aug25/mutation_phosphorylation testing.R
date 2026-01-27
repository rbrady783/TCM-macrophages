# =========================================
# Full Script: Mutation-Cytokine Analysis with Proper Data Formatting
# =========================================

# 0) INSTALL / LOAD PACKAGES ----
# install.packages(c("readr","dplyr","purrr","broom","tidyr","stringr","ggplot2","knitr","tibble"))
library(readr)
library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(stringr)
library(ggplot2)
library(knitr)
library(tibble)

# 1) READ THE FORMATTED DATA ----
# Note: Run the data preparation script first to create these files
cyto <- read_csv("means_modz_for_mutations.csv", show_col_types = FALSE) %>%
  rename(treatment = 1) %>%
  mutate(treatment = str_trim(as.character(treatment)))

mut  <- read_csv("mutations_updated.csv", show_col_types = FALSE) %>%
  rename(treatment = 1) %>%
  mutate(treatment = str_trim(as.character(treatment)))

# 2) VERIFY AND MERGE DATA ----
common <- intersect(cyto$treatment, mut$treatment)
print(paste("Common treatments:", length(common)))

cyto   <- filter(cyto,  treatment %in% common)
mut    <- filter(mut,   treatment %in% common)
df     <- inner_join(cyto, mut, by = "treatment")

cytokines <- setdiff(names(cyto),  "treatment")
mutations <- setdiff(names(mut),   "treatment")

print(paste("Cytokines to analyze:", length(cytokines)))
print(paste("Mutations to analyze:", length(mutations)))

# Convert mutation columns to factors (they should already be, but double-check)
df <- df %>% mutate(across(all_of(mutations), ~ as.factor(.x)))

print("\nMutation variable levels:")
for(mut_var in mutations) {
  cat(mut_var, ":", paste(levels(df[[mut_var]]), collapse = ", "), "\n")
}

# 3) RUN NONPARAMETRIC TESTS ----
print("\nRunning statistical tests...")

results <- 
  expand_grid(cytokine = cytokines, mut = mutations) %>%
  mutate(test = map2(cytokine, mut, ~ {
    y <- df[[.x]]  
    x <- df[[.y]]  
    
    # Skip if mutation has no variation
    if (nlevels(x) < 2 | all(table(x) == 0)) {
      return(tibble(method=NA, statistic=NA_real_, p.value=NA_real_))
    } 
    
    # Check if there are enough observations in each group
    group_counts <- table(x)
    if (any(group_counts == 0) || length(group_counts) < 2) {
      return(tibble(method=NA, statistic=NA_real_, p.value=NA_real_))
    }
    
    # Two-group comparison: Wilcoxon rank-sum test
    if (nlevels(x) == 2) {
      tryCatch({
        tidy(wilcox.test(y ~ x, exact = FALSE))
      }, error = function(e) {
        tibble(method="Wilcoxon rank sum test", statistic=NA_real_, p.value=NA_real_)
      })
    } 
    # Multi-group comparison: Kruskal-Wallis test
    else {
      tryCatch({
        tidy(kruskal.test(y ~ x))
      }, error = function(e) {
        tibble(method="Kruskal-Wallis rank sum test", statistic=NA_real_, p.value=NA_real_)
      })
    }
  })) %>%
  unnest(test) %>%
  select(cytokine, mut, method, statistic, p.value) %>%
  # Remove rows where test couldn't be performed
  filter(!is.na(p.value)) %>%
  # Apply multiple testing correction within each mutation
  group_by(mut) %>%
  mutate(p.adj = p.adjust(p.value, "BH")) %>%
  ungroup()

print(paste("Total tests performed:", nrow(results)))

# 4) FILTER EXPLORATORY HITS (unadj p < 0.05) ----
exploratory <- filter(results, p.value < 0.05)
significant <- filter(results, p.adj < 0.05)

print(paste("Exploratory hits (p < 0.05):", nrow(exploratory)))
print(paste("Significant hits (FDR < 0.05):", nrow(significant)))

# 5) DIRECTION FUNCTION (enhanced for multi-level factors) ----
get_direction <- function(df, mut_col, cyto_col) {
  x    <- df[[mut_col]]; y <- df[[cyto_col]]
  lvls <- levels(x)
  
  if (length(lvls) < 2) return(NA_character_)
  
  # Calculate median for each group
  m <- tapply(y, x, median, na.rm = TRUE)
  
  if (all(is.na(m)) || length(m) < 2) return(NA_character_)
  
  # For binary variables (Yes/No)
  if (length(lvls) == 2) {
    if ("Yes" %in% names(m) && "No" %in% names(m)) {
      if (m["Yes"] > m["No"]) {
        return(paste0("higher when ", mut_col, " = Yes"))
      } else {
        return(paste0("higher when ", mut_col, " = No"))
      }
    }
  }
  
  # For multi-level variables
  hi <- names(which.max(m))
  lo <- names(which.min(m))
  
  if (length(hi) == 1 && length(lo) == 1 && hi != lo) {
    return(paste0("highest in ", mut_col, " = ", hi, "; lowest in ", mut_col, " = ", lo))
  } else {
    return(paste0("median values: ", paste(paste(names(m), round(m, 2), sep = "="), collapse = ", ")))
  }
}

# 6) ANNOTATE DIRECTIONS ----
exploratory_dir <- exploratory %>%
  rowwise() %>%
  mutate(direction = get_direction(df, mut, cytokine)) %>%
  ungroup() %>%
  arrange(p.value)

# 7) ENHANCED BOXPLOT FUNCTION ----
plot_mut_cyto <- function(df, mut_col, cyto_col, pval, padj = NULL, out_file) {
  dfp <- df %>%
    select(all_of(c(mut_col, cyto_col))) %>%
    rename(mut = all_of(mut_col), cyto = all_of(cyto_col)) %>%
    mutate(mut = factor(mut))
  
  # Remove empty factor levels
  dfp$mut <- droplevels(dfp$mut)
  
  # Compute vertical position for p-value
  y_max <- max(dfp$cyto, na.rm = TRUE)
  y_min <- min(dfp$cyto, na.rm = TRUE)
  y_range <- y_max - y_min
  y_pos <- y_max + 0.08 * y_range
  
  # Create p-value text
  p_text <- paste0("p = ", sprintf("%.3f", pval))
  if (!is.null(padj)) {
    p_text <- paste0(p_text, "\nFDR = ", sprintf("%.3f", padj))
  }
  
  plt <- ggplot(dfp, aes(x = mut, y = cyto, fill = mut)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.8, size = 2) +
    annotate("text",
             x = Inf, y = y_pos,
             label = p_text,
             hjust = 1.1, vjust = 0,
             size = 3.5, lineheight = 0.9) +
    labs(
      x     = mut_col,
      y     = cyto_col,
      title = paste0(cyto_col, " by ", mut_col)
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1),
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.margin     = margin(5, 5, 20, 5)
    ) +
    coord_cartesian(clip = "off")
  
  ggsave(out_file, plot = plt, width = 5, height = 4, dpi = 300)
  invisible(plt)
}

# 8) GENERATE & SAVE BOXPLOTS FOR EXPLORATORY HITS ----
dir.create("mutation_cytokine_boxplots", showWarnings = FALSE)

print("\nGenerating boxplots...")
pwalk(list(exploratory_dir$mut, exploratory_dir$cytokine, exploratory_dir$p.value, exploratory_dir$p.adj), 
      function(mut_col, cyto_col, pval, padj) {
        fn <- file.path("mutation_cytokine_boxplots",
                        paste0(mut_col, "_vs_", cyto_col, ".png"))
        plot_mut_cyto(df, mut_col, cyto_col, pval, padj, fn)
      })

# 9) POST-HOC PAIRWISE TESTS FOR MULTI-LEVEL MUTATIONS ----
multi_hits <- exploratory_dir %>% 
  filter(str_detect(method, "Kruskal-Wallis")) %>%
  select(mut, cytokine) %>%
  distinct()

if (nrow(multi_hits) > 0) {
  print("\nPerforming post-hoc pairwise tests...")
  
  pairwise_results <- map_dfr(1:nrow(multi_hits), function(i) {
    mut_col <- multi_hits$mut[i]
    cyto_col <- multi_hits$cytokine[i]
    
    # subset data to this pair
    dfp <- df %>% 
      select(all_of(c(mut_col, cyto_col))) %>%
      rename(mut = all_of(mut_col), cyto = all_of(cyto_col)) %>%
      mutate(mut = droplevels(as.factor(mut)))
    
    # Only proceed if we have multiple levels with observations
    if (nlevels(dfp$mut) < 2) return(NULL)
    
    tryCatch({
      # Run pairwise Wilcoxon with BH correction
      pw <- pairwise.wilcox.test(
        x               = dfp$cyto, 
        g               = dfp$mut, 
        p.adjust.method = "BH", 
        exact           = FALSE
      )
      
      # Convert to long format
      if (!is.null(pw$p.value) && !all(is.na(pw$p.value))) {
        pw$p.value %>%
          as.data.frame() %>%
          rownames_to_column("group1") %>%
          pivot_longer(-group1, names_to = "group2", values_to = "p.adj") %>%
          filter(!is.na(p.adj)) %>%
          mutate(
            mutation = mut_col,
            cytokine = cyto_col
          ) %>%
          select(mutation, cytokine, group1, group2, p.adj)
      }
    }, error = function(e) {
      warning(paste("Pairwise test failed for", mut_col, "vs", cyto_col))
      NULL
    })
  })
  
  if (!is.null(pairwise_results) && nrow(pairwise_results) > 0) {
    write_csv(pairwise_results, "pairwise_posthoc_mutation_cytokine.csv")
    print("Post-hoc pairwise results saved to: pairwise_posthoc_mutation_cytokine.csv")
  }
} else {
  print("No multi-level mutations with significant associations found.")
}

# 10) SAVE RESULTS TABLES ----
write_csv(
  exploratory_dir %>% 
    select(mut, cytokine, method, statistic, p.value, p.adj, direction) %>%
    arrange(p.value),
  "exploratory_cytokine_mutation_associations.csv"
)

if (nrow(significant) > 0) {
  significant_dir <- significant %>%
    rowwise() %>%
    mutate(direction = get_direction(df, mut, cytokine)) %>%
    ungroup() %>%
    arrange(p.adj)
  
  write_csv(
    significant_dir %>% 
      select(mut, cytokine, method, statistic, p.value, p.adj, direction),
    "significant_cytokine_mutation_associations.csv"
  )
}

exploratory_table <- read_csv("exploratory_cytokine_mutation_associations.csv")
print(exploratory_table)

# For driver mutations (binary: 0 vs 1)
table(df$`2 or more drivers`)

# For pAKT (three levels: A, C, L)
table(df$pAKT)

Show full exploratory results
exploratory_table <- read_csv("exploratory_cytokine_mutation_associations.csv")
print(exploratory_table, width = Inf)

# Show pairwise results
pairwise <- read_csv("pairwise_posthoc_mutation_cytokine.csv")
print(pairwise)

