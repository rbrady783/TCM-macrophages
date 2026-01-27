# install.packages(c("readr","dplyr","stringr","purrr","tibble"))  # if needed
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# 1) point to your folder
res_dir <- "C:/Users/brady/OneDrive/Desktop/TCM_RNA/Pairwise code/GSEA results"   # ← adjust this

# 2) list all your .csv files
all_files <- list.files(
  res_dir, 
  pattern    = "_GSEA_.*_preRanked\\.csv$", 
  full.names = TRUE
)

# 3) build a metadata table
file_info <- tibble(path = all_files) %>%
  mutate(
    file       = basename(path),
    cytokine   = str_extract(file, "^[^_]+"),
    collection = str_extract(file, "(?<=_GSEA_)[^_]+(?=_preRanked)")
  )

# 4) read them all in, keeping metadata
all_gsea <- file_info %>%
  mutate(data = map(path, read_csv)) %>%
  dplyr::select(-path, -file)

# 5) unnest into one big tibble
combined <- all_gsea %>%
  unnest(data)

# 6) split back by collection
hallmark_df <- combined %>% filter(collection == "Hallmark")
C2_df       <- combined %>% filter(collection == "C2")
C5_df       <- combined %>% filter(collection == "C5")
C6_df       <- combined %>% filter(collection == "C6")   # added C6
C7_df       <- combined %>% filter(collection == "C7")

# 7) write one file per collection
readr::write_csv(hallmark_df, file.path(res_dir, "Hallmark_GSEA_results.csv"))
readr::write_csv(C2_df,       file.path(res_dir, "C2_GSEA_results.csv"))
readr::write_csv(C5_df,       file.path(res_dir, "C5_GSEA_results.csv"))
readr::write_csv(C6_df,       file.path(res_dir, "C6_GSEA_results.csv"))  # new
readr::write_csv(C7_df,       file.path(res_dir, "C7_GSEA_results.csv"))