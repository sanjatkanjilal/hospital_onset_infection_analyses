# /usr/bin/R

#####################################################################################################
# Description: Matches cases to controls for colonization pressure / hospital-onset infection analysis
#
# DEPENDENT SCRIPTS
#
# HO_infxn_functions.R - Contains all functions in data processing scripts
# HO_infxn_micro_prep.R - Bins micro results into organism categories
# HO_infxn_input_data_processing.R - Pre-processes micro, antibiotic, encounter, and admit source datasets
# HO_infxn_build_unmatched_cohorts.R - Generates cohorts without features
# HO_infxn_add_features.R - Adds features to cohorts
#
# INPUT DATASETS
#
# Pathogen categories table
# Unmatched cohort + all features
#
# OUTPUT DATASETS
#
# Matched cohort + all features
#
# ENVIRONMENTAL FACTOR ANALYSIS
# Case / control matching criteria: 
#   - Age +/- 3 years
#   - Sex
#   - Prior surgery in the previous 90 days
#   - Length of stay +/- 3 days or 10%
#   - Prior antibiotic courses (by class) in the previous 60 days
#
# PATIENT-LEVEL ANALYSIS
# Case / control matching criteria: 
#   - Ward 
#   - Month / year
#   - Length of stay +/- 3 days or 10%
#
# Authors: Sanjat Kanjilal
# Last updated: 2025-02-17
#####################################################################################################

#### SET ENVIRONMENT VARIABLES ####

# Import the functions for the pipeline
source("/PHShome/zw852/colonization-pressure_HAI/HO_infxn_functions.R")

# Set working directory for imports
mainDir <- "/data/tide/projects/ho_infxn_ml/"
setwd(file.path(mainDir))

# Global paths
path_cat_table_filename <- "clean_data/20250603_SensitivityAnalysis/path_cat_table_matching.csv"
unmatched_features_filename <- "clean_data/20250603_SensitivityAnalysis/unmatched_case_controls_features_org_group_3.csv"

#### BUILD MATCHED DATASETS ####

path_cat_table_matching <- readr::read_csv(path_cat_table_filename)
unmatched_features <- readr::read_csv(unmatched_features_filename)

unmatched_features <- subset(unmatched_features, select = -c(prior_NA))
# Import datasets as needed
# unmatched_features <- read_csv(file = paste0("clean_data/",  today, "/unmatched_case_controls_features_", pathogen_hierarchy, ".csv"))

# Match on patient characteristics to estimate impact of colonization presssure
environment.matched <- do.call(dplyr::bind_rows,
                           lapply(levels(factor(path_cat_table_matching$pathogen_category)), function(y) 
                               
                             {
                               print(paste("Matching cases to controls for environmental-level analysis for", y))
                               
                               pathogen_hierarchy <- unique(path_cat_table_matching$pathogen_hierarchy)
                               
                               unmatched_features.pathogen <- unmatched_features %>%
                                 dplyr::filter(run == y)
                               
                               # Add matching features to case / control data
                               tic()
                               
                               matched_features.pathogen <- environmental.matching_process(pathogen_category = y,
                                                                                     dat = unmatched_features.pathogen)
                               
                               toc()
                               
                               # Compile data across each pathogen category
                               matched_cc <- as.data.frame(cbind(matched_features.pathogen))
                               
                               }))

# Match on environmental characteristics to estimate impact of antibiotic exposure
patient.matched <- do.call(dplyr::bind_rows,
                               lapply(levels(factor(path_cat_table_matching$pathogen_category)), function(y) 
                                 
                               {
                                 print(paste("Matching cases to controls for patient-level analysis for", y))
                                 
                                 pathogen_hierarchy <- unique(path_cat_table_matching$pathogen_hierarchy)
                                 
                                 unmatched_features.pathogen <- unmatched_features %>%
                                   dplyr::filter(run == y)
                                 
                                 # Add matching features to case / control data
                                 tic()
                                 
                                 matched_features.pathogen <- patient.matching_process(pathogen_category = y,
                                                                                             dat = unmatched_features.pathogen)
                                 
                                 toc()
                                 
                                 # Compile data across each pathogen category
                                 matched_cc <- as.data.frame(cbind(matched_features.pathogen))
                                 
                               }))

# Bind environmental and patient-level matched datasets with features
environment.matched <- environment.matched %>%
  dplyr::mutate(match = "environmental")

patient.matched <- patient.matched %>%
  dplyr::mutate(match = "patient")

matched.final <- dplyr::bind_rows(environment.matched, patient.matched) %>%
  dplyr::select(match, run:group, group_index, PatientID:DR_PsA_cp) %>%
  dplyr::distinct()

# Save matched patient-level dataset + features
readr::write_csv(matched.final, file = paste0("clean_data/20250603_SensitivityAnalysis/final_cohort_org_group_3.csv"))

#### FINAL DATASET CLEAN / PREP FOR MODELS ####

# Check sample sizes
sample_sizes <- matched.final %>%
  dplyr::select(match:PatientID) %>%
  dplyr::distinct() %>%
  dplyr::group_by(match, run, group) %>%
  dplyr::summarize(sample_size = dplyr::n()) %>%
  dplyr::ungroup() %>%
  tidyr::spread(group, sample_size)

# Check for missingness
percent_missing_values <- matched.final %>%
  dplyr::select(-group_index) %>%
  dplyr::group_by(match, run, group) %>%
  dplyr::summarise_all(~sum(is.na(.)) / dplyr::n() * 100) %>%
  tidyr:: pivot_longer(cols = c(PatientID:DR_PsA_cp), names_to = "variable", values_to = "percent_missing") %>%
  dplyr::filter(percent_missing > 0) %>%
  dplyr::arrange(variable, match, run, group)

# Check for empty columns
colSums(is.na(matched.final))

# Binarize group
matched.final.models <- matched.final %>%
  dplyr::mutate(group_binary = ifelse(group == "case", 1, 0))

# Select a limited variable pool
abx_of_interest <- c("penicillin_0_60", "extended_spectrum_penicillin_0_60", "cephalosporin_0_60", 
                     "extended_spectrum_cephalosporin_0_60", "carbapenem_0_60", "anti_staph_beta_lactam_0_60", 
                     "fluoroquinolone_0_60", "glycopeptide_0_60", "anti_anaerobe_0_60", "anti_Cdiff_0_60",
                     "tetracycline_0_60", "macrolide_0_60", "sulfonamide_0_60", "lincosamide_0_60", "any_abx_0_60")

variables_to_drop <- names(matched.final.models)[grepl("_0_60$", names(matched.final.models)) & !(names(matched.final.models) %in% abx_of_interest)]

matched.final.models.with_elix <- matched.final.models %>% 
  dplyr::select(match:group_index, group_binary, PatientID, DTS_in_month, DTS_in_year, duration,
                time_to_infxn, matching_duration, age, sex, starts_with('elix'), 
                elix_index_mortality, any_surgery, admit_source_clean, 
                any_abx_0_60:other_abx_0_60, prior_C_diff:DR_PsA_cp) %>%
  dplyr::select(-all_of(variables_to_drop)) %>%
  dplyr::select(-matches("60_plus"))

matched.final.models <- matched.final.models %>% 
  dplyr::select(match:group_index, group_binary, PatientID, DTS_in_month, DTS_in_year, duration,
         time_to_infxn, matching_duration, age, sex, # starts_with('elix'), 
         elix_index_mortality, any_surgery, admit_source_clean, 
         any_abx_0_60:other_abx_0_60, prior_C_diff:DR_PsA_cp) %>%
  dplyr::select(-all_of(variables_to_drop)) %>%
  dplyr::select(-matches("60_plus"))

# Drop remaining NAs from predictor features
matched.final.models <- matched.final.models %>%
  tidyr::drop_na(age, CDiff_cp:DR_PsA_cp)
matched.final.models.with_elix <- matched.final.models.with_elix %>%
  tidyr::drop_na(age, CDiff_cp:DR_PsA_cp)

# Impute the remaining missing values with MICE (Elixhauser scores)
# imp <- mice(matched.final)
# matched.final <- complete(imp)

# Recheck missingness
colSums(is.na(matched.final.models))

# Save file
readr::write_csv(matched.final.models.with_elix, file = paste0("clean_data/20250603_SensitivityAnalysis/final_dataset_for_models_elix_20250603_SensitivityAnalysis.csv"))
readr::write_csv(matched.final.models, file = paste0("clean_data/20250603_SensitivityAnalysis/final_dataset_for_models_20250603_SensitivityAnalysis.csv"))

#### CLEAN UP ####

# Remove global paths / variables
rm(mainDir, path_cat_table_filename, unmatched_features_filename)

# Remove datasets
rm(path_cat_table_matching, unmatched_features, environment.matched, patient.matched, matched.final, sample_sizes,
   percent_missing_values, matched.final.models, abx_of_interest, variables_to_drop)

# Remove functions
rm(grouping_function, micro_dedup, initial_abx_clean, process_ip_abx_features, calculate_abx_ip_durations, 
   calculate_op_abx_durations, calculate_abx_op_last_date, abx_courses, enc_processing, path_cat_hierarchy,
   admt.mapped, first_eligible_rooms, micro_prep, room_dat_eligible_trimmed, adt_micro_join, time_to_infxn, 
   drop_old_micro, remove_treated_patients, drop_prev_enc, adt_micro_hierarchy, potential_cases, after_first_room, 
   flag_prior_infxn, generate_controls, add_demographics, add_abx, final_abx_clean, add_comorbidities, add_cpt, 
   add_admit.source, add_prior_pathogens, calculate_CP_score, calculate_and_add_cp_scores_parallel, 
   pathogen_selection, adt_micro_initial_prep, generate_unmatched_data, generate.unmatched.dataset, 
   environmental.matching_process, patient.matching_process)


