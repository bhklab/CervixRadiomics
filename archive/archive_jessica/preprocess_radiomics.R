library(tidyverse)

fnames <- list.files("../data/features/", pattern = ".+bincount.+")
# fnames <- c("features_081x081_all.csv", "features_162x162_all.csv")

dfs <- lapply(fnames, function(f) {
  spacing <- str_split(f, "_")[[1]][2]
  df <- read_csv(paste("../data/features", f, sep = "/")) %>%
    select(-contains("diagnostics"), -Image, -Mask) %>%
    rename_at(vars(-patient_id), list(~paste("spacing", spacing, ., sep = "_")))
  return(df)
})

dfs %>% 
  reduce(left_join, by = "patient_id") %>%
  write_csv("../data/features_all.csv")