# ————— INSTALL / LOAD PACKAGES —————
# install.packages(c("ggplot2","dplyr","purrr","readr","stringr"))
library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(stringr)

# assume:
# df              : your merged cytokine + mutation data.frame
# exploratory_dir : data.frame with columns mut, cytokine, p.value

# ————— 1) Helper to plot one mutation–cytokine pair —————
plot_mut_cyto <- function(df, mut_col, cyto_col, pval, out_file) {
  dfp <- df %>%
    select(all_of(c(mut_col, cyto_col))) %>%
    rename(mut = all_of(mut_col), cyto = all_of(cyto_col)) %>%
    mutate(mut = factor(mut))
  
  # find the top of your data range
  y_max <- max(dfp$cyto, na.rm = TRUE)
  y_pos <- y_max + 0.05 * (y_max - min(dfp$cyto, na.rm = TRUE))
  
  plt <- ggplot(dfp, aes(x = mut, y = cyto, fill = mut)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
    # place the p-value at y_pos, with a little horizontal padding
    annotate("text", 
             x = Inf, y = y_pos,
             label = paste0("p = ", sprintf("%.3f", pval)),
             hjust = 1.1, vjust = 0,        # hjust>1 pushes left of the right edge
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
      plot.margin     = margin(5, 5, 20, 5)  # add extra space at top
    ) +
    coord_cartesian(clip = "off")           # allow text outside panel
  
  ggsave(out_file, plot = plt, width = 4, height = 4, dpi = 300)
  invisible(plt)
}

  ggsave(out_file, plot = plt, width = 4, height = 4, dpi = 300)
  invisible(plt)
}

# ————— 2) Loop over all exploratory hits and save plots —————
# create an output directory
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
    fn       <- paste0("mutation_cytokine_boxplots/",
                       mut_col, "_vs_", cyto_col, ".png")
    plot_mut_cyto(df, mut_col, cyto_col, pval, fn)
  }
)
