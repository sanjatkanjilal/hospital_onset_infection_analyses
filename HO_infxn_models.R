# /usr/bin/R

#####################################################################################################
# Description: Models for colonization pressure / hospital-onset infection analysis
#
# DEPENDENT SCRIPTS
#
# HO_infxn_functions.R - Contains all functions in data processing scripts
# HO_infxn_micro_prep.R - Bins micro results into organism categories
# HO_infxn_input_data_processing.R - Pre-processes micro, antibiotic, encounter, and admit source datasets
# HO_infxn_build_unmatched_cohorts.R - Generates cohorts without features
# HO_infxn_add_features.R - Adds features to unmatched cohorts
# HO_infxn_build_matched_cohorts.R - Matches cases to controls for final analysis
#
# INPUT DATASETS
#
# Matched cohorts with features
#
# OUTPUTS
#
# Conditional logistic regression model results
# XGBoost model results
#
# Authors: Ziming (Alex) Wei, Sanjat Kanjilal
# Last updated: 2025-02-17
#####################################################################################################

##### LIBRARIES ####
library(tidyverse)
library(survival)
library(xgboost)
library(caret)
library(pROC)
library(PRROC)
library(Matrix)
library(conflicted)
conflicts_prefer(base::`%in%`)

#### ENVIRONMENTAL VARIABLES ####

# Set working directory and output directory to be whatever you generate today
# mainDir <- "/data/tide/projects/ho_infxn_ml/"
mainDir <- "/Users/zimingwei/GitHub/nosocomial-acquisition_colonization-pressure/"
setwd(file.path(mainDir))

#### IMPORT DATASETS ####
cc_final <- readr::read_csv("clean_data/HO_infxn_environmental_analysis.csv")

#### SET UP MODELS ####

# Looping variables
matching_rubric <- unique(cc_final$match)
targets <- unique(cc_final$run)

# Holding datasets
clogit_coefficients <- data.frame()

#### CONDITIONAL LOGISTIC REGRESSION ####

# Test dataset using C difficile only
# dat.run <- cc_final %>%
#   dplyr::filter(match == "environmental" & run == "C_diff") %>%
#   as.data.frame(.)

# Loop through all organisms
clr <- do.call(dplyr::bind_rows,
               lapply(levels(factor(matching_rubric)), function(x) 
                 {
                 dat.match <- cc_final %>% dplyr::filter(match == x)
                 
                 if(x == "environmental")
                   {
                   do.call(dplyr::bind_rows,
                           lapply(levels(factor(dat.match$run)), function(y) 
                             {
                             dat.run <- dat.match %>% 
                               dplyr::filter(run == y)
                             
                             prior_path_target <- paste0("prior_", y)
                             print(paste("Environmental match regression for", y))
                             dat.run <- dat.run %>% dplyr::select(group_binary, group_index, elix_index_mortality, CDiff_cp:DR_PsA_cp)
                             dat.run <- dat.run %>% dplyr::select(where(~dplyr::n_distinct(.) > 1))
                             form1 = as.formula(paste("group_binary ~", paste(names(dat.run)[3:ncol(dat.run)], collapse = "+"), paste("+ strata(group_index)")))
                             fit.clr <- clogit(form1, data = dat.run)
                             summary_logistic <- summary(fit.clr)
                             coef_clogit <- dplyr::bind_cols(summary_logistic[["conf.int"]][, c(1,3,4)], 
                                                      summary_logistic[["coefficients"]][,c(3,5)]) %>%
                               dplyr::mutate(match = x,
                                      target = y,
                                      variable = rownames(summary_logistic[["coefficients"]])) %>%
                               dplyr::mutate(variable = ifelse(variable == "dat.run[, prior_path_target]", prior_path_target, variable))
                             rownames(coef_clogit) <- NULL
                             clogit_coefficients <- rbind(clogit_coefficients, coef_clogit)
                             
                             }))
                   } else {
                     do.call(dplyr::bind_rows,
                             lapply(levels(factor(dat.match$run)), function(y) 
                               {
                               dat.run <- dat.match %>% 
                                 dplyr::filter(run == y)
                               print(paste("Patient match regression for", y))
                               prior_path_target <- paste0("prior_", y)
                               dat.run <- dat.run %>% dplyr::select(group_binary, group_index, 
                                                             age, sex, elix_index_mortality, any_surgery, anti_anaerobe_0_60:tetracycline_0_60)
                               dat.run <- dat.run %>% dplyr::select(where(~dplyr::n_distinct(.) > 1))
                               form1 = as.formula(paste("group_binary ~", paste(names(dat.run)[3:ncol(dat.run)], collapse = "+"), paste("+ strata(group_index)")))
                               fit.clr <- clogit(form1, data = dat.run)
                               summary_logistic <- summary(fit.clr)
                               coef_clogit <- dplyr::bind_cols(summary_logistic[["conf.int"]][, c(1,3,4)], 
                                                        summary_logistic[["coefficients"]][,c(3,5)]) %>%
                                 dplyr::mutate(match = x,
                                        target = y,
                                        variable = rownames(summary_logistic[["coefficients"]])) %>%
                                 dplyr::mutate(variable = ifelse(variable == "test1[, prior_path_target]", prior_path_target, variable))
                               rownames(coef_clogit) <- NULL
                               clogit_coefficients <- rbind(clogit_coefficients, coef_clogit)
                             }))
                     }
                 }))
                           
# Format results
clr.final <- clr %>%
  dplyr::select(match:variable, coef = `exp(coef)`, lower_CI = `lower .95`, upper_CI = `upper .95`, SE_coef = `se(coef)`, pval = `Pr(>|z|)`) %>%
  dplyr::mutate(coef_clean = dplyr::case_when(
    coef < 6 ~ as.character(round(coef, 2)),
    coef > 6 & coef < 10 ~ "6 - 10",
    coef > 10  ~ ">10",
    .default = NULL)) %>%
  dplyr::mutate(lower_CI_clean = dplyr::case_when(
    lower_CI < 6 ~ as.character(round(lower_CI, 2)),
    lower_CI > 6 & lower_CI < 10 ~ "6 - 10",
    lower_CI > 10 ~ ">10",
    .default = NULL)) %>%
  dplyr:: mutate(upper_CI_clean = dplyr::case_when(
    upper_CI < 6 ~ as.character(round(upper_CI, 2)),
    upper_CI > 6 & upper_CI < 10 ~ "6 - 10",
    upper_CI > 10 ~ ">10",
    .default = NULL)) %>%
  dplyr::mutate(SE_coef_clean = dplyr::case_when(
    SE_coef < 6 ~ as.character(round(SE_coef, 2)),
    SE_coef > 6 & SE_coef < 10 ~ "6 - 10",
    SE_coef > 10 ~ ">10",
    .default = NULL)) %>%
  dplyr::mutate(pval_clean = dplyr::case_when(
    pval > 0.05 ~ "p > 0.05",
    pval < 0.05 & pval > 0.01 ~ "p < 0.05",
    pval < 0.01 & pval > 0.001 ~ "p < 0.01",
    pval < 0.001 ~ "p < 0.001",
    .default = NULL)) %>%
  dplyr::mutate(sig_flag = ifelse(pval < 0.05, "*", NA)) %>%
  dplyr::mutate(across(ends_with("clean"), ~ ifelse(. %in% c("1", "0"), paste0(., ".00"), .))) %>%
  dplyr::select(match:variable, coef = coef_clean, lower_CI = lower_CI_clean, upper_CI = upper_CI_clean, 
         SE_coef = SE_coef_clean, pval = pval_clean, sig_flag, coef_raw = coef, lower_CI_raw = lower_CI, 
         upper_CI_raw = upper_CI, SE_coef_raw = SE_coef, pval_raw = pval)

# Save file
readr::write_csv(clr.final, file = paste0("results/model_results/add_dept/clogit_coefficients.csv"))

#### XGBOOST MODELS ####

xgboost_metrics <- data.frame()
xgboost_feature_importance <- data.frame()
env_dat_with_pred <- data.frame()

# Running XGBoost environmental matching
run_xgboost_environmental <- function(y){
  
  # Filter for dataset to contain the focused pathogent
  dat.run <- dat.match %>% dplyr::filter(run == y)
  
  # Select Features to include in training the conditional logistic regression model
  prior_path_target <- paste0("prior_", y)
  print(paste("Environmental match regression for", y))
  dat.run <- dat.run %>% select(group_binary, elix_index_mortality, CDiff_cp:DR_PsA_cp)
  dat.run <- dat.run %>% select(where(~n_distinct(.) > 1)) # Remove features with only 1 value among samples
  
  # One-hot encode non-numeric columns
  data_encoded <- dat.run %>%
    mutate(across(where(is.character), as.factor)) %>%
    model.matrix(~ . - 1, data = .)
  
  # 5-fold cross-validation setup
  folds <- createFolds(data_encoded[, "group_binary"], k = 5, list = TRUE, returnTrain = FALSE)
  auc_values <- c()
  auprc_values <- c()
  importance_list <- list()
  shap_value_list <- list()
  
  for(i in 1:5) {
    # Splitting data into training and test sets for each fold
    test_indices <- folds[[i]]
    train_indices <- setdiff(1:nrow(data_encoded), test_indices)
    train_data <- data_encoded[train_indices, ]
    test_data <- data_encoded[test_indices, ]
    test_matrix <- test_data
    
    # Create XGBoost matrices
    xgb_train <- xgb.DMatrix(data = as.matrix(train_data[, -which(colnames(train_data) == "group_binary")]), 
                             label = train_data[, "group_binary"])
    xgb_test <- xgb.DMatrix(data = as.matrix(test_data[, -which(colnames(test_data) == "group_binary")]), 
                            label = test_data[, "group_binary"])
    
    # Calculate the scale_pos_weight for current fold
    pos_count <- sum(train_data[, "group_binary"] == 1)
    neg_count <- sum(train_data[, "group_binary"] == 0)
    scale_pos_weight <- neg_count / pos_count
    
    # Update XGBoost parameters
    params <- list(objective = "binary:logistic",
                   eval_metric = "auc",
                   max_depth = 6,
                   eta = 0.01,
                   subsample = 0.5,
                   colsample_bytree = 0.5,
                   min_child_weight = 1,
                   scale_pos_weight = scale_pos_weight)
    
    # Train XGBoost model
    xgb_model <- xgb.train(params = params,
                           data = xgb_train,
                           nrounds = 150,
                           watchlist = list(train = xgb_train, test = xgb_test),
                           early_stopping_rounds = 20,
                           print_every_n = 10)
    
    # Make predictions on test data
    test_data <- data.frame(test_data)
    test_data$pred <- predict(xgb_model, xgb_test)
    env_dat_with_pred <<- rbind(env_dat_with_pred, test_data %>% select(elix_index_mortality, CDiff_cp:DR_PsA_cp, group_binary, pred) %>% mutate(run = y, match = 'environmental',fold = i))
    
    # Calculate and store AUC for each fold
    roc_obj <- roc(test_data$group_binary, test_data$pred)
    auc_values <- c(auc_values, auc(roc_obj))
    pr_curve <- pr.curve(scores.class0 = test_data$pred[test_data$group_binary == 1],
                         scores.class1 = test_data$pred[test_data$group_binary == 0], 
                         curve = TRUE)
    auprc_values <- c(auprc_values, pr_curve$auc.integral)
    
    # Store feature importance for each fold
    importance_list[[i]] <- xgb.importance(model = xgb_model)
    
    # Save the models
    output_dir = paste0("results/model_results/","20250603_SensitivityAnalysis","/xgb/model_checkpoints/environmental_", y)
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    } 
    xgb.save(xgb_model, paste0(output_dir, "/fold_", i, ".model"))
  }
  
  # Calculate average AUROC
  avg_auc <- mean(auc_values)
  std_auc <- sd(auc_values) / sqrt(length(auc_values))
  t_value <- qt(0.975, df = length(auc_values) - 1)
  margin_of_error <- t_value * std_auc
  auc_lower_CI <- avg_auc - margin_of_error
  auc_upper_CI <- avg_auc + margin_of_error
  cat("Average AUC across folds:", avg_auc, "\n")
  
  # Calculate average AUPRC
  avg_auprc <- mean(auprc_values)
  std_auprc <- sd(auprc_values)
  t_value <- qt(0.975, df = length(auprc_values) - 1)
  margin_of_error <- t_value * std_auprc
  auprc_lower_CI <- avg_auprc - margin_of_error
  auprc_upper_CI <- avg_auprc + margin_of_error
  cat("Average AUPRC across folds:", avg_auprc, "\n")
  
  # Aggregate and calculate average feature importance
  all_importances <- bind_rows(importance_list)
  avg_importance <- all_importances %>% 
    group_by(Feature) %>% 
    summarize(Gain = mean(Gain), Cover = mean(Cover), Frequency = mean(Frequency))
  
  # Calculating confusion matrix of model prediction
  threshold <- 0.5
  xgb_preds_binary <- ifelse(test_data$pred > threshold, 1, 0)
  cm <- confusionMatrix(factor(xgb_preds_binary), 
                        factor(test_data$group_binary),
                        positive = '1')
  
  # For XGBoost metrics
  metrics <- data.frame(
    Target = y,
    Match = 'environmental',
    AUROC = avg_auc,
    AUROC_lower_CI = auc_lower_CI,
    AURPC_upper_CI = auc_upper_CI,
    AUPRC = avg_auprc,
    AUPRC_lower_CI = auprc_lower_CI,
    AUPRC_upper_CI = auprc_upper_CI,
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    BalancedAccuracy = cm$byClass["Balanced Accuracy"],
    PPV = cm$byClass["Pos Pred Value"],
    NPV = cm$byClass["Neg Pred Value"],
    Precision = cm$byClass["Precision"],
    Recall = cm$byClass["Recall"],
    F1 = cm$byClass["F1"]
  )
  
  rownames(metrics) <- NULL
  xgboost_metrics <<- rbind(xgboost_metrics, metrics)
  
  # Feature importance
  avg_importance$Target <- y
  avg_importance$match <-'environmental'
  xgboost_feature_importance <<- rbind(xgboost_feature_importance, avg_importance)
  
}

# Running environmental match
for (rubric in matching_rubric){
  dat.match <- cc_final %>% dplyr::filter(match == rubric)
  if (rubric == 'environmental'){
    remaining_targets = setdiff(targets, env_dat_with_pred$run)
    for (target in remaining_targets){
      cat(paste0("Running environmental match models for target: ", target), "\n\n")
      run_xgboost_environmental(target)
      write.csv(env_dat_with_pred, paste0("results/model_results/","20250603_SensitivityAnalysis","/xgb/environmental_match_predictions.csv"), row.names = FALSE)
      write.csv(xgboost_metrics, paste0("results/model_results/","20250603_SensitivityAnalysis","/xgb/xgboost_metrics.csv"), row.names = FALSE)
    }
  }
}


#### CLEAN UP ####

# Remove global paths / variables
rm(mainDir, matching_rubric, targets)

# Remove datasets
rm(cc_final, dat.run, clogit_coefficients, clr, clr.final)

# Remove functions


