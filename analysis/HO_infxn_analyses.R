#### HEADER ####
# Description: 
#   Script to run LR and XGBoost models on case/control cohorts for HO-infection study
# Dependencies: 
#   final_dataset_for_models.csv
# Output: 
#   results/model_results/clr/environmental_colonization_pressure_with_pred.csv
#   results/model_results/clr/clr_coefficients_removed_prior_pathogen.csv
#   results/model_results/clr/clr_model_performance.csv
#   results/model_results/xgb/model_checkpoints/
#   results/model_results/xgb/environmental_match_predictions.csv
#   results/model_results/xgb/xgboost_metrics.csv
# Authors: Luke Sagers, Ziming (Alex) Wei, Sanjat Kanjilal
# Last updated: 2023-12-16

##### LIBRARIES ####
library(tidyverse)
library(lubridate)
library(mice)
library(skimr)
library(survival)
library(xgboost)
library(caret)
library(pROC)
library(PRROC)
library(Matrix)
library(conflicted)
library(gridExtra)
conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::lag())
conflicts_prefer(lubridate::month)
conflicts_prefer(lubridate::year)
conflicts_prefer(lubridate::week)
conflicts_prefer(lubridate::quarter)
conflicts_prefer(dplyr::first)
conflicts_prefer(reshape2::melt)
conflicts_prefer(reshape2::dcast)

# Set working directory and output directory to be whatever you generate today
mainDir <- "~/ho_infxn_clean_code/"
setwd(file.path(mainDir))

#### IMPORT DATASETS ####
cc_final <- data.table::fread(file = "results/dataset/final_dataset_for_models.csv")

#### SET UP MODELS ####

# Looping variables
matching_rubric <- unique(cc_final$match)
targets <- unique(cc_final$run)

# Holding datasets
clogit_coefficients <- data.frame()
clogit_metrics <- data.frame()
env_dat_with_pred <- data.frame()
patient_dat_with_pred <- data.frame()

#### CONDITIONAL LOGISTIC REGRESSION ####

# Loop through all organisms
clr <- do.call(bind_rows,
               lapply(levels(factor(matching_rubric)), function(x) 
                 {
                 dat.match <- cc_final %>% filter(match == x)
                 
                 if(x == "environmental")
                   {
                       do.call(bind_rows,
                               lapply(levels(factor(dat.match$run)), function(y) 
                               {
                                 dat.run <- dat.match %>% 
                                   filter(run == y) 
                                 folds = createFolds(unique(dat.run$group_index), k = 5, returnTrain = TRUE)
                                 
                                 for (i in seq(1,5)){
                                   dat.run <- dat.match %>% 
                                     filter(run == y) 
                                   dat.run <- dat.run %>% filter(group_index %in% folds[[i]])
                                   prior_path_target <- paste0("prior_", y)
                                   print(paste("Environmental match regression for", y))
                                   dat.run <- dat.run %>% select(group_binary, group_index, CDiff_cp:DR_PsA_cp)
                                   dat.run <- dat.run %>% select(where(~n_distinct(.) > 1))
                                   
                                   form1 = as.formula(paste("group_binary ~", paste(names(dat.run)[3:ncol(dat.run)], collapse = "+"), paste("+ strata(group_index)")))
                                   fit.clr <- clogit(form1, data = dat.run)
                                   
                                   dat.run$pred <- predict(fit.clr,
                                                           newdata = dat.run,
                                                           type = 'lp')
                                   
                                   dat.run$pred <- plogis(dat.run$pred) 
                                   
                                   env_dat_with_pred <<- rbind(env_dat_with_pred, dat.run %>% select(CDiff_cp:DR_PsA_cp, group_binary, pred) %>% mutate(run = y))
                                   pr_curve <- pr.curve(scores.class0 = dat.run$pred[dat.run$group_binary == 1],
                                                        scores.class1 = dat.run$pred[dat.run$group_binary == 0], 
                                                        curve = TRUE)
                                   auprc_value <- pr_curve$auc.integral
                                   
                                   # Extra Metrics
                                   # Calculating confusion matrix of model prediction
                                   threshold <- 0.5
                                   dat.run$pred_binary <- ifelse(dat.run$pred > threshold, 1, 0)
                                   cm <- confusionMatrix(factor(dat.run$pred_binary), 
                                                         factor(dat.run$group_binary),
                                                         positive = '1')
                                   
                                   # For Metrics
                                   metrics <- data.frame(
                                     Target = y,
                                     Match = 'environmental',
                                     AUROC = auc(dat.run$group_binary,
                                                 dat.run$pred),
                                     AUPRC = auprc_value,
                                     LRP = cm$byClass["Sensitivity"] / ( 1 - cm$byClass["Specificity"]),
                                     LRN = ( 1 - cm$byClass["Sensitivity"]) / cm$byClass["Specificity"],
                                     PPV = cm$byClass["Pos Pred Value"],
                                     NPV = cm$byClass["Neg Pred Value"]
                                   )
                                   
                                   rownames(metrics) <- NULL
                                   
                                   clogit_metrics <<- rbind(clogit_metrics, metrics)
                                   summary_logistic <- summary(fit.clr)
                                   coef_clogit <- bind_cols(summary_logistic[["conf.int"]][, c(1,3,4)], 
                                                            summary_logistic[["coefficients"]][,c(3,5)]) %>%
                                     mutate(match = x,
                                            target = y,
                                            variable = rownames(summary_logistic[["coefficients"]])) %>%
                                     mutate(variable = ifelse(variable == "dat.run[, prior_path_target]", prior_path_target, variable))
                                   rownames(coef_clogit) <- NULL
                                   if(i==1){
                                     clogit_coefficients <<- rbind(clogit_coefficients, coef_clogit)}
                                   }
                               }))
                   } 
                 }))

# Calculate mean, standard error and confidence intervals
mean_ci <- function(x) {
  mean_val <- mean(x)
  se <- sd(x) / sqrt(length(x))  # Standard error
  ci_lower <- mean_val - qt(0.975, df = length(x) - 1) * se  # Lower CI
  ci_upper <- mean_val + qt(0.975, df = length(x) - 1) * se  # Upper CI
  return(c(mean = mean_val, lower = ci_lower, upper = ci_upper))
}

# Calculate the metrics for the conditional logistic regression model
clogit_metrics_integrated <- clogit_metrics %>%
  group_by(Target, Match) %>%
  summarize(AUROC_mean = mean(AUROC), 
            AUROC_lower_ci = mean_ci(AUROC)[2], 
            AUROC_upper_ci = mean_ci(AUROC)[3],
            AUPRC_mean = mean(AUPRC), 
            AUPRC_lower_ci = mean_ci(AUPRC)[2], 
            AUPRC_upper_ci = mean_ci(AUPRC)[3],
            PPV_mean = mean(PPV),
            PPV_lower_ci = mean_ci(PPV)[2],
            PPV_upper_ci = mean_ci(PPV)[3],
            NPV_mean = mean(NPV),
            NPV_lower_ci = mean_ci(NPV)[2],
            NPV_upper_ci = mean_ci(NPV)[3],
            LRP_mean = mean(LRP),
            LRP_lower_ci = mean_ci(LRP)[2],
            LRP_upper_ci = mean_ci(LRP)[3],
            LRN_mean = mean(LRN),
            LRN_lower_ci = mean_ci(LRN)[2],
            LRN_upper_ci = mean_ci(LRN)[3],
            )

# Store the coefficients of the conditional logistic regression                           
clr <- clogit_coefficients %>%
  select(match:variable, coef = `exp(coef)`, lower_CI = `lower .95`, upper_CI = `upper .95`, SE_coef = `se(coef)`, pval = `Pr(>|z|)`) %>%
  mutate(sig_flag = ifelse(pval < 0.05, "*", NA))

# Save Output file for Conditional logistic regression models
write.csv(env_dat_with_pred, file = "results/model_results/clr/environmental_colonization_pressure_with_pred.csv", 
          row.names = FALSE, na = "")
write.csv(clr, file = "results/model_results/clr/clr_coefficients.csv", 
          row.names = FALSE, na = "")
write.csv(clogit_metrics_integrated, file = "results/model_results/clr/clr_model_performance.csv",
          row.names = FALSE, na = "")

#### XGBOOST MODELS ####

# XGB holding datasets
glm_coefficients <- data.frame()
xgboost_metrics <- data.frame()
xgboost_feature_importance <- data.frame()
env_dat_with_pred <- data.frame()
patient_dat_with_pred <- data.frame()

# XGBoost with 5-fold cross-validation for environmental analysis
env_dat_with_pred <- data.table::fread(file = "results/model_results/xgb/environmental_match_predictions.csv")
xgboost_metrics <- data.table::fread(file = "results/model_results/xgb/xgboost_metrics.csv")

# Running XGBoost environmental matching
run_xgboost_environmental <- function(y){
  
    # Filter for dataset to contain the focused pathogent
    dat.run <- dat.match %>% filter(run == y)
    
    # Select Features to include in training the conditional logistic regression model
    prior_path_target <- paste0("prior_", y)
    print(paste("Environmental match regression for", y))
    dat.run <- dat.run %>% select(group_binary, CDiff_cp:DR_PsA_cp)
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
      env_dat_with_pred <<- rbind(env_dat_with_pred, test_data %>% select(CDiff_cp:DR_PsA_cp, group_binary, pred) %>% mutate(run = y, match = 'environmental',fold = i))
      
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
      output_dir = paste0("results/model_results/xgb/model_checkpoints/environmental_", y)
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
    
    # # Get SHAP plot
    # contr <- predict(xgb_model, xgb_test, predcontrib = TRUE)
    # shap_sum <- xgb.ggplot.shap.summary(data = as.matrix(test_matrix[, -which(colnames(test_matrix) == "group_binary")]), 
    #                                     shap_contrib = contr, model = xgb_model) +
    #             theme_minimal() +
    #             labs(x = "Shap Value", y = "Feature", color = "Feature Value")
    # 
    # shap_feature_plots <- xgb.plot.shap(data = as.matrix(test_matrix[, -which(colnames(test_matrix) == "group_binary")]), 
    #                                     shap_contrib = contr, model = xgb_model, top_n = 12, n_col = 3)
    # 
    # shap_summary_plot <- xgb.ggplot.shap.summary(as.matrix(test_matrix[, -which(colnames(test_matrix) == "group_binary")]), 
    #                                              model = xgb_model, 
    #                                              target_class = 1, 
    #                                              top_n = 10)
    # 
    # ggtitle(y)
    # 
    # ggsave(filename = paste0("results/xgb/shap_plots/environmental_match/",y,"_xgb_shap.png"),
    #        plot = shap_summary_plot)
}

# Running environmental match
for (rubric in matching_rubric){
  dat.match <- cc_final %>% filter(match == rubric)
  if (rubric == 'environmental'){
    remaining_targets = setdiff(targets, env_dat_with_pred$run)
    for (target in remaining_targets){
      cat(paste0("Running environmental match models for target: ", target), "\n\n")
      run_xgboost_environmental(target)
      write.csv(env_dat_with_pred, "results/model_results/xgb/environmental_match_predictions.csv", row.names = FALSE)
      write.csv(xgboost_metrics, "results/model_results/xgb/xgboost_metrics.csv", row.names = FALSE)
    }
  }
}


