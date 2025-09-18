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

#### CONDITIONAL LOGISTIC REGRESSION ####
# Loop through all organisms

k <- 5  # number of CV folds

clr <- do.call(dplyr::bind_rows,
               lapply(levels(factor(matching_rubric)), function(x) 
               {
                 dat.match <- cc_final %>% dplyr::filter(match == x)
                 
                 if (x == "environmental")
                 {
                   do.call(dplyr::bind_rows,
                           lapply(levels(factor(dat.match$run)), function(y) 
                           {
                             dat.run <- dat.match %>% 
                               dplyr::filter(run == y)
                             
                             prior_path_target <- paste0("prior_", y)
                             print(paste("Environmental match regression for", y))
                             
                             dat.run <- dat.run %>% 
                               dplyr::select(group_binary, group_index, elix_index_mortality, CDiff_cp:DR_PsA_cp) %>%
                               dplyr::select(where(~dplyr::n_distinct(.) > 1))
                             
                             pred_cols <- setdiff(names(dat.run), c("group_binary","group_index"))
                             
                             # set K if not defined elsewhere
                             k <- if (exists("k")) get("k") else 5
                             
                             set.seed(1)
                             uniq_groups <- sort(unique(dat.run$group_index))
                             fold_assign <- sample(rep(1:k, length.out = length(uniq_groups)))
                             grp2fold <- data.frame(group_index = uniq_groups, fold = fold_assign)
                             
                             dat.run <- dat.run %>% dplyr::left_join(grp2fold, by = "group_index")
                             
                             aucs   <- rep(NA_real_, k)
                             auprcs <- rep(NA_real_, k)
                             ppvs   <- rep(NA_real_, k)
                             npvs   <- rep(NA_real_, k)
                             lrp    <- rep(NA_real_, k)  # LR+
                             lrn    <- rep(NA_real_, k)  # LR-
                             
                             for (i in 1:k) {
                               dat.train <- dat.run %>% dplyr::filter(fold != i)
                               dat.test  <- dat.run %>% dplyr::filter(fold == i)
                               
                               form1 <- as.formula(paste(
                                 "group_binary ~",
                                 paste(pred_cols, collapse = " + "),
                                 "+ strata(group_index)"
                               ))
                               
                               fit.clr.cv <- clogit(form1, data = dat.train, x = TRUE, model = TRUE)
                               
                               beta <- coef(fit.clr.cv)
                               tt <- terms(fit.clr.cv)
                               sp <- attr(tt, "specials")$strata
                               tt_nostrata <- if (!is.null(sp)) drop.terms(tt, sp, keep.response = FALSE) else tt
                               mm_test <- model.matrix(tt_nostrata, data = dat.test)
                               
                               common <- intersect(colnames(mm_test), names(beta))
                               score <- if (length(common) > 0) {
                                 as.numeric(mm_test[, common, drop = FALSE] %*% beta[common])
                               } else {
                                 rep(NA_real_, nrow(dat.test))
                               }
                               dat.test$score <- score
                               
                               ok_pos <- sum(dat.test$group_binary == 1, na.rm = TRUE) > 0
                               ok_neg <- sum(dat.test$group_binary == 0, na.rm = TRUE) > 0
                               ok_var <- sd(dat.test$score, na.rm = TRUE) > 0
                               
                               if (ok_pos && ok_neg && ok_var) {
                                 # AUROC
                                 roc_obj <- pROC::roc(dat.test$group_binary, dat.test$score)
                                 aucs[i] <- as.numeric(roc_obj$auc)
                                 
                                 # AUPRC
                                 pos <- dat.test$score[dat.test$group_binary == 1]
                                 neg <- dat.test$score[dat.test$group_binary == 0]
                                 pr  <- PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg)
                                 auprcs[i] <- as.numeric(pr$auc.integral)
                                 
                                 # Threshold by Youden's J, then PPV/NPV + LR+/LR-
                                 th_metrics <- pROC::coords(
                                   roc_obj, x = "best",
                                   best.method = "youden",
                                   ret = c("threshold","sensitivity","specificity","precision","npv"),
                                   transpose = TRUE
                                 )
                                 sens <- as.numeric(th_metrics["sensitivity"])
                                 spec <- as.numeric(th_metrics["specificity"])
                                 ppvs[i] <- as.numeric(th_metrics["precision"]) # PPV
                                 npvs[i] <- as.numeric(th_metrics["npv"])       # NPV
                                 
                                 # LR+ = sens/(1-spec), LR- = (1-sens)/spec (guard edge cases)
                                 lrp[i] <- if (is.finite(1 - spec) && (1 - spec) > 0) sens / (1 - spec) else NA_real_
                                 lrn[i] <- if (is.finite(spec) && spec > 0) (1 - sens) / spec else NA_real_
                               } else {
                                 aucs[i]   <- NA_real_
                                 auprcs[i] <- NA_real_
                                 ppvs[i]   <- NA_real_
                                 npvs[i]   <- NA_real_
                                 lrp[i]    <- NA_real_
                                 lrn[i]    <- NA_real_
                               }
                             }
                             
                             mean_auc   <- mean(aucs,   na.rm = TRUE)
                             sd_auc     <- sd(aucs,     na.rm = TRUE)
                             min_auc    <- suppressWarnings(min(aucs,   na.rm = TRUE))
                             max_auc    <- suppressWarnings(max(aucs,   na.rm = TRUE))
                             
                             mean_auprc <- mean(auprcs, na.rm = TRUE)
                             sd_auprc   <- sd(auprcs,   na.rm = TRUE)
                             min_auprc  <- suppressWarnings(min(auprcs, na.rm = TRUE))
                             max_auprc  <- suppressWarnings(max(auprcs, na.rm = TRUE))
                             
                             mean_ppv   <- mean(ppvs,   na.rm = TRUE)
                             sd_ppv     <- sd(ppvs,     na.rm = TRUE)
                             mean_npv   <- mean(npvs,   na.rm = TRUE)
                             sd_npv     <- sd(npvs,     na.rm = TRUE)
                             mean_lrp   <- mean(lrp,    na.rm = TRUE)
                             sd_lrp     <- sd(lrp,      na.rm = TRUE)
                             mean_lrn   <- mean(lrn,    na.rm = TRUE)
                             sd_lrn     <- sd(lrn,      na.rm = TRUE)
                             
                             k_used <- sum(is.finite(aucs))
                             
                             cv_row <- dplyr::tibble(
                               match = x,
                               target = y,
                               k = k,
                               k_used = k_used,
                               mean_auc = mean_auc,
                               sd_auc = sd_auc,
                               min_auc = ifelse(is.infinite(min_auc), NA_real_, min_auc),
                               max_auc = ifelse(is.infinite(max_auc), NA_real_, max_auc),
                               mean_auprc = mean_auprc,
                               sd_auprc = sd_auprc,
                               min_auprc = ifelse(is.infinite(min_auprc), NA_real_, min_auprc),
                               max_auprc = ifelse(is.infinite(max_auprc), NA_real_, max_auprc),
                               mean_ppv = mean_ppv,
                               sd_ppv = sd_ppv,
                               mean_npv = mean_npv,
                               sd_npv = sd_npv,
                               mean_lr_pos = mean_lrp,
                               sd_lr_pos = sd_lrp,
                               mean_lr_neg = mean_lrn,
                               sd_lr_neg = sd_lrn
                             )
                             
                             cv_row
                           }))
                 }
               }))

# Save file
readr::write_csv(clr, file = paste0("results/model_results/clogit_performance.csv"))

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
  
  # 10-fold CV -> 8:1:1 (train:val:test). For fold i, test = i, val = i+1.
  k <- 10
  folds <- createFolds(data_encoded[, "group_binary"], k = k, list = TRUE, returnTrain = FALSE)
  auc_values <- c()
  auprc_values <- c()
  importance_list <- list()
  shap_value_list <- list()
  all_test_preds <- list()  # collect test predictions across folds
  
  for(i in seq_len(k)) {
    test_indices <- folds[[i]]
    val_fold_idx <- ifelse(i < k, i + 1, 1)
    val_indices  <- folds[[val_fold_idx]]
    train_indices <- setdiff(seq_len(nrow(data_encoded)), union(test_indices, val_indices))
    
    train_data <- data_encoded[train_indices, ]
    val_data   <- data_encoded[val_indices, ]
    test_data  <- data_encoded[test_indices, ]
    test_matrix <- test_data
    
    # Create XGBoost matrices
    xgb_train <- xgb.DMatrix(
      data = as.matrix(train_data[, -which(colnames(train_data) == "group_binary")]),
      label = train_data[, "group_binary"]
    )
    xgb_val <- xgb.DMatrix(
      data = as.matrix(val_data[, -which(colnames(val_data) == "group_binary")]),
      label = val_data[, "group_binary"]
    )
    xgb_test <- xgb.DMatrix(
      data = as.matrix(test_data[, -which(colnames(test_data) == "group_binary")]),
      label = test_data[, "group_binary"]
    )
    
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
    
    # Train XGBoost model (early stop on validation, not test)
    xgb_model <- xgb.train(params = params,
                           data = xgb_train,
                           nrounds = 150,
                           watchlist = list(train = xgb_train, val = xgb_val),
                           early_stopping_rounds = 20,
                           print_every_n = 10)
    
    # Make predictions on test data (evaluation only on test)
    test_df <- data.frame(test_data)
    test_df$pred <- predict(xgb_model, xgb_test)
    all_test_preds[[i]] <- test_df
    
    # Save per-fold preds to global collector (like before)
    env_dat_with_pred <<- rbind(
      env_dat_with_pred,
      test_df %>%
        dplyr::select(elix_index_mortality, CDiff_cp:DR_PsA_cp, group_binary, pred) %>%
        dplyr::mutate(run = y, match = 'environmental', fold = i)
    )
    
    # Calculate and store AUC/AUPRC for this test fold
    roc_obj <- roc(test_df$group_binary, test_df$pred)
    auc_values <- c(auc_values, auc(roc_obj))
    pr_curve <- pr.curve(scores.class0 = test_df$pred[test_df$group_binary == 1],
                         scores.class1 = test_df$pred[test_df$group_binary == 0], 
                         curve = TRUE)
    auprc_values <- c(auprc_values, pr_curve$auc.integral)
    
    # Store feature importance for each fold
    importance_list[[i]] <- xgb.importance(model = xgb_model)
    
    # Save the models
    output_dir <- paste0("results/model_results/xgb/model_checkpoints/environmental_", y)
    if (!dir.exists(output_dir)){
      dir.create(output_dir, recursive = TRUE)
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
  
  # Calculate average AUPRC (kept same style as your code)
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
    summarize(Gain = mean(Gain), Cover = mean(Cover), Frequency = mean(Frequency), .groups = "drop")
  
  # Confusion matrix on concatenated test predictions across folds
  all_test_df <- bind_rows(all_test_preds)
  threshold <- 0.5
  xgb_preds_binary <- ifelse(all_test_df$pred > threshold, 1, 0)
  cm <- confusionMatrix(factor(xgb_preds_binary), 
                        factor(all_test_df$group_binary),
                        positive = '1')
  
  # For XGBoost metrics (kept your column names as-is)
  metrics <- data.frame(
    Target = y,
    Match = 'environmental',
    AUROC = avg_auc,
    AUROC_lower_CI = auc_lower_CI,
    AURPC_upper_CI = auc_upper_CI,  # (kept your original naming)
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
  avg_importance$match <- 'environmental'
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
      write.csv(env_dat_with_pred, paste0("results/model_results/xgb/environmental_match_predictions.csv"), row.names = FALSE)
      write.csv(xgboost_metrics, paste0("results/model_results/xgb/xgboost_metrics.csv"), row.names = FALSE)
    }
  }
}

#### CLEAN UP ####

# Remove global paths / variables
rm(mainDir, matching_rubric, targets)

# Remove datasets
rm(cc_final, dat.run, clogit_coefficients, clr, clr.final)

# Remove functions


