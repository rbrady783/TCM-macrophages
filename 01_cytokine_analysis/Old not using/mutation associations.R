# =========================================
# Full Script: Read Data → Run Tests → Extract Exploratory Hits → Boxplots with p‐values
# =========================================

# 0) INSTALL / LOAD PACKAGES ----
# install.packages(c("readr","dplyr","purrr","broom","tidyr","stringr","ggplot2","knitr"))
library(readr)
library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(stringr)
library(ggplot2)
library(knitr)

# 1) READ & MATCH TREATMENTS ----
cyto <- read_csv("heatmap.csv", show_col_types = FALSE) %>%
  rename(treatment = 1) %>%
  mutate(treatment = str_trim(treatment))
mut  <- read_csv("mutations.csv", show_col_types = FALSE) %>%
  rename(treatment = 1) %>%
  mutate(treatment = str_trim(treatment))

common <- intersect(cyto$treatment, mut$treatment)
cyto   <- filter(cyto,  treatment %in% common)
mut    <- filter(mut,   treatment %in% common)

# 2) MERGE & PREPARE ----
df        <- inner_join(cyto, mut, by = "treatment")
cytokines <- setdiff(names(cyto),  "treatment")
mutations <- setdiff(names(mut),   "treatment")
df        <- mutate(df, across(all_of(mutations), ~ factor(.)))

# 3) RUN NONPARAMETRIC TESTS ----
results <- 
  expand_grid(cytokine = cytokines, mut = mutations) %>%
  mutate(test = map2(cytokine, mut, ~ {
    y <- df[[.x]]  
    x <- df[[.y]]  
    if (nlevels(x) < 2) {
      tibble(method=NA, statistic=NA_real_, p.value=NA_real_)
    } else if (nlevels(x) == 2) {
      tidy(wilcox.test(y ~ x, exact = FALSE))
    } else {
      tidy(kruskal.test(y ~ x))
    }
  })) %>%
  unnest(test) %>%
  select(cytokine, mut, method, statistic, p.value) %>%
  group_by(mut) %>%
  mutate(p.adj = p.adjust(p.value, "BH")) %>%
  ungroup()

# 4) FILTER EXPLORATORY HITS (unadj p < 0.05) ----
exploratory <- filter(results, p.value < 0.05)

# 5) DIRECTION FUNCTION ----
get_direction <- function(df, mut_col, cyto_col) {
  x    <- df[[mut_col]]; y <- df[[cyto_col]]
  lvls <- levels(x)
  if (length(lvls) < 2) return(NA_character_)
  m    <- tapply(y, x, median, na.rm = TRUE)
  if (all(is.na(m))) return(NA_character_)
  hi <- names(which.max(m)); lo <- names(which.min(m))
  paste0("highest in ", mut_col, "=", hi, 
         "; lowest in ", mut_col, "=", lo)
}

# 6) ANNOTATE DIRECTIONS ----
exploratory_dir <- exploratory %>%
  rowwise() %>%
  mutate(direction = get_direction(df, mut, cytokine)) %>%
  ungroup()

# 7) BOXPLOT HELPER (with uncluttered p‐value) ----
plot_mut_cyto <- function(df, mut_col, cyto_col, pval, out_file) {
  dfp <- df %>%
    select(all_of(c(mut_col, cyto_col))) %>%
    rename(mut = all_of(mut_col), cyto = all_of(cyto_col)) %>%
    mutate(mut = factor(mut))
  
  # compute vertical position for p‐value
  y_max <- max(dfp$cyto, na.rm = TRUE)
  y_min <- min(dfp$cyto, na.rm = TRUE)
  y_pos <- y_max + 0.05 * (y_max - y_min)
  
  plt <- ggplot(dfp, aes(x = mut, y = cyto, fill = mut)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
    annotate("text",
             x = Inf, y = y_pos,
             label = paste0("p = ", sprintf("%.3f", pval)),
             hjust = 1.1, vjust = 0,
             size = 4) +
    labs(
      x     = mut_col,
      y     = cyto_col,
      title = paste0(cyto_col, " by ", mut_col)
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.margin     = margin(5, 5, 20, 5)
    ) +
    coord_cartesian(clip = "off")
  
  ggsave(out_file, plot = plt, width = 4, height = 4, dpi = 300)
  invisible(plt)
}

# 8) GENERATE & SAVE BOXPLOTS ----
dir.create("mutation_cytokine_boxplots", showWarnings = FALSE)
walk2(
  exploratory_dir$mut,
  exploratory_dir$cytokine,
  ~ {
    mut_col  <- .x
    cyto_col <- .y
    pval     <- exploratory_dir$p.value[
      exploratory_dir$mut == mut_col &
        exploratory_dir$cytokine == cyto_col
    ]
    fn <- file.path("mutation_cytokine_boxplots",
                    paste0(mut_col, "_vs_", cyto_col, ".png"))
    plot_mut_cyto(df, mut_col, cyto_col, pval, fn)
  }
)

# 9) (OPTIONAL) SAVE EXPLORATORY TABLE ----
write_csv(
  exploratory_dir %>% select(mut, cytokine, method, statistic, p.value, direction),
  "exploratory_cytokine_mutation_associations.csv"
)

# ————— POST-HOC PAIRWISE TESTS FOR MULTI-LEVEL MUTATIONS —————

# We’ll use pairwise Wilcoxon with BH correction:
# (you could swap in FSA::dunnTest() if you prefer Dunn’s test)
library(tibble)   # for rownames_to_column()


multi_hits <- exploratory_dir %>% 
  filter(method == "Kruskal-Wallis rank sum test") %>%
  pull(mut) %>% unique()

pairwise_results <- map_dfr(multi_hits, function(mut_col) {
  # find all cytokines tested for this multi‐level mutation
  cyto_list <- exploratory_dir %>% 
    filter(mut == mut_col) %>% 
    pull(cytokine)
  
  map_dfr(cyto_list, function(cyto_col) {
    # subset data to this pair
    dfp <- df %>% select(all_of(c(mut_col, cyto_col))) %>%
      rename(mut = all_of(mut_col), cyto = all_of(cyto_col)) %>%
      mutate(mut = factor(mut))
    
    # run pairwise Wilcoxon with BH
    pw <- pairwise.wilcox.test(
      x              = dfp$cyto, 
      g              = dfp$mut, 
      p.adjust.method= "BH", 
      exact          = FALSE
    )
    
    # tidy up into a long data.frame
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
  })
})

# ————— SAVE PAIRWISE RESULTS —————
write_csv(pairwise_results, "pairwise_posthoc_mutation_cytokine.csv")

# ————— OPTIONAL: VIEW A SAMPLE —————
print(pairwise_results)
