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
# library(caret)
# library(pROC)
# library(Matrix)
# library(conflicted)
conflicts_prefer(base::`%in%`)

#### ENVIRONMENTAL VARIABLES ####

# Set working directory and output directory to be whatever you generate today
mainDir <- "/data/tide/projects/ho_infxn_ml/"
setwd(file.path(mainDir))

#### IMPORT DATASETS ####
cc_final <- readr::read_csv("clean_data/20250217/final_dataset_for_models_20250217.csv")

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
readr::write_csv(clr.final, file = paste0("results/model_results/20250217/clogit_coefficients.csv"))

#### XGBOOST MODELS ####

#### CLEAN UP ####

# Remove global paths / variables
rm(mainDir, matching_rubric, targets)

# Remove datasets
rm(cc_final, dat.run, clogit_coefficients, clr, clr.final)

# Remove functions


