# uncomment before first use to install required packages
# install.packages("tidyverse")
# install.packages("survival")
# install.packages("BiocManager)
# BiocManager::install("survcomp")
# install.packages("broom")
# install.packages("caret")
# install.packages(c("devtools", "githubinstall"))
# library(githubinstall)
# githubinstall("bhklab/mRMRe", ref = "bo_upgrade")

library(tidyverse)
library(survival)
library(survcomp)
library(broom)
library(caret)
library(mRMRe)

Sys.setenv(KMP_DUPLICATE_LIB_OK=TRUE) # when installing the new mRMRe from GitHub, libomp is loaded twice,
# crashing R, probably related to my installation of clang-omp.

clinical <- read_csv("../data/clinical_clean.csv")
radiomics <- read_csv("../data/radiomics.csv") %>% 
  select(-contains("diagnostics"), -contains("Image"), -contains("Mask")) %>%
  rename_at(vars(-patient_id), .funs = ~paste("radiomics_", ., sep = ""))

data <- inner_join(clinical, radiomics)
clinical_cols <- 
  list(
    clinical_basic_sex_f  = sym("SEX_F"),
    clinical_basic_age    = sym("AGE"),
    # clinical_basic_tstage = sym("cT"),
    # clinical_basic_nstage = sym("cN"),
    # clinical_basic_chemo  = sym("Concurrent chemo (1=yes, 0=no)"),
    clinical_basic_volume = sym("GTV Size (cc)"),
    # clinical_basic_hiv    = sym("HIV (0=no, 1=yes)"),
    clinical_tstage = sym("cT"),
    clinical_nstage = sym("cN"),
    clinical_chemo  = sym("Concurrent chemo (1=yes, 0=no)"),
    clinical_hiv    = sym("HIV (0=no, 1=yes)"),
    clinical_grade        = sym("Hist Grade"),
    clinical_size         = sym("Size"),
    clinical_dose         = sym("Primary Dose (cGy)"),
    clinical_hb           = sym("Hb"),
    clinical_hb_100       = sym("Hb >100 (1=yes, 0=no)"),
    clinical_wbc          = sym("WBC"),
    clinical_plt          = sym("Plt"),
    clinical_neut         = sym("Neut"),
    clinical_lymph        = sym("Lymph"),
    clinical_plr          = sym("Platelets:Lymphocyte ratio"),
    clinical_nlr          = sym("Neutrophil:Lymphocyte ratio"),
    clinical_nlr2         = sym("NLR>2 (1= >2, 0=<2)"))

data <- data %>%
  rename(!!!clinical_cols) %>%
  drop_na(contains("clinical"))

mrmr_select <- function(data, feature_count = 10, fixed_feature_count = 0) {
  data_df <- as.data.frame(data)
  data_df$surv <- Surv(data_df$time_col, data_df$event_col)
  data_df$time_col <- NULL
  data_df$event_col <- NULL
  mrmr_data <- mRMR.data(data_df)
  mrmr_obj <- mRMR.ensemble(data = mrmr_data,
                            target_indices = c(length(data_df)),
                            solution_count = 1,
                            feature_count = feature_count,
                            fixed_feature_count = fixed_feature_count)
  return(featureNames(mrmr_obj)[solutions(mrmr_obj)[[1]]])
}

run_cox_model <- function(train, test, feature_count = 10, fixed_feature_count = 0) {
  
  selected_features <- mrmr_select(train, feature_count, fixed_feature_count)
  
  model <- train %>%
    select(selected_features, time_col, event_col) %>%
    coxph(Surv(time_col, event_col) ~ ., data=.)
  
  predicted <- predict(model,
                       test %>% select(selected_features),
                       type = "risk")
  
  # print("predicted:")
  # print(predicted)
  # print("time:")
  # print(test$time_col)
  # print("event:")
  # print(test$event_col)
  # print(paste(rep("=", 20)))
  return(tidy(model, exponentiate = TRUE) %>%
           mutate(ci=concordance.index(predicted,
                                       test$time_col,
                                       test$event_col,
                                       outx = FALSE)$c.index))
}

find_optimal_feature_count <- function(data, lower = 5, upper = 20) {
  folds <- createFolds(as.factor(data$event_col), k = 3)
  results <- tibble(n_feats = numeric(), fold = numeric(), ci = numeric())
  for(n_feats in lower:upper) {
    for(i in 1:length(folds))  {
      train_fold <- data[-folds[[i]],]
      test_fold <- data[folds[[i]],]
      ci <- run_cox_model(train_fold, test_fold, n_feats)$ci[[1]]
      results <- results %>% add_row(n_feats = n_feats,
                                     fold = i,
                                     ci = ci)
    }
  }
  top_n_feats <- results %>% group_by(n_feats) %>% summarise(ci=mean(ci)) %>% top_n(1) %>% select(n_feats)
  return(top_n_feats[[1]])
}

run_analysis <- function(data, feature_count = 10 , n_iter = 100, fixed_feature_count = 0) {
  results <- tibble(
    ci = numeric(),
    variable = character(),
    iteration = numeric(),
    model = character())
  
  for(i in 1:n_iter) {
    train_indices <- createDataPartition(as.factor(data$event_col), p = .8, list = FALSE)
    train <- data[train_indices,]
    test <- data[-train_indices,]
    
    if(feature_count == "optim") {
      feature_count_clin <- train %>%
        select(contains("clinical_basic"), time_col, event_col) %>%
        find_optimal_feature_count(upper = ncol(.) - 2)
    } else {
      feature_count_clin <- min(ncol(data %>%
                                       select(contains("clinical_basic"))), 
                                feature_count)
    }
    
    clinical_basic_model <- run_cox_model(train %>% select(contains("clinical_basic"), time_col, event_col),
                                          test,
                                          feature_count_clin,
                                          fixed_feature_count = 0)
    selected_clinical_basic <- gsub("`", "", clinical_basic_model$term)
    
    if(feature_count == "optim") {
      feature_count_rad <- train %>%
        select(contains("radiomics"), time_col, event_col) %>%
        find_optimal_feature_count(upper = ncol(train %>% select(contains("clinical"))))
    } else {
      feature_count_rad <- feature_count
    }
    
    radiomics_model <- run_cox_model(train %>% select(contains("radiomics"), time_col, event_col),
                                     test,
                                     feature_count_rad)
    selected_radiomics <- gsub("`", "", radiomics_model$term)
    
    if(feature_count == "optim") {
      feature_count_comb <- train %>%
        select(selected_clinical_basic, selected_radiomics, time_col, event_col) %>%
        find_optimal_feature_count(upper = ncol(.) - 2)
    } else {
      feature_count_comb <- feature_count
    }
    
    combined_model <- run_cox_model(train %>% select(selected_clinical_basic, selected_radiomics, time_col, event_col),
                                    test,
                                    feature_count = feature_count_comb,
                                    fixed_feature_count = fixed_feature_count)
    selected_combined <- gsub("`", "", combined_model$term)
    
    if(feature_count == "optim") {
      feature_count_clin_adv <- train %>%
        select(contains("clinical"), time_col, event_col) %>%
        find_optimal_feature_count(upper = ncol(.) - 2)
    } else {
      feature_count_clin_adv <- feature_count
    }
    
    clinical_adv_model <- run_cox_model(train %>% select(contains("clinical"),
                                                         time_col, 
                                                         event_col),
                                        test,
                                        feature_count_clin_adv)
    selected_clinical_adv <- gsub("`", "", clinical_adv_model$term)
    
    for(v in selected_clinical_basic) {
      results <- results %>% add_row(
        ci = clinical_basic_model$ci[[1]],
        variable = v,
        iteration = i,
        model = "clinical_basic")
    }
    
    for(v in selected_radiomics) {
      results <- results %>% add_row(
        ci = radiomics_model$ci[[1]],
        variable = v,
        iteration = i,
        model = "radiomics")
    }
    
    for(v in selected_combined) {
      results <- results %>% add_row(
        ci = combined_model$ci[[1]],
        variable = v,
        iteration = i,
        model = "combined")
    }
    
    for(v in selected_clinical_adv) {
      results <- results %>% add_row(
        ci = clinical_adv_model$ci[[1]],
        variable = v,
        iteration = i,
        model = "clinical_adv")
    }
  }
  return(results)
}

plot_variable_frequencies <- function(data) {
  plt <- data %>% 
    group_by(model, variable) %>% 
    summarise(percent_selected = n_distinct(iteration) / 1000) %>%
    top_n(10, wt = percent_selected) %>%
    ungroup() %>% 
    arrange(model, percent_selected) %>%
    mutate(order = row_number())
  
  p <- plt %>%
    ggplot(aes(x = order, y = percent_selected)) +
    geom_col() +
    facet_grid(rows = vars(model), scales = "free_y") + 
    scale_x_continuous(
      breaks = plt$order,
      labels = plt$variable,
      expand = c(0,0)) + 
    coord_flip()
  return(p)
}

results_os <- run_analysis(data %>% filter(patient_id != 994) %>% select(contains("clinical"), contains("radiomics"), time_col = death_time, event_col = event_os), n_iter = 100, feature_count = 10, fixed_feature_count = data %>% select(contains("clinical_basic")) %>% ncol)