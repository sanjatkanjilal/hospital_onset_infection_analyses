# /usr/bin/R

#####################################################################################################
# Description: Build unmatched cohorts for colonization pressure / hospital-onset infection analysis
#
# DEPENDENT SCRIPTS
#
# HO_infxn_functions.R - Contains all functions in data processing scripts
# HO_infxn_micro_prep.R - Bins micro results into organism categories
# HO_infxn_input_data_processing.R - Pre-processes micro, antibiotic, encounter, and admit source datasets
#
# INPUT DATASETS
#
# Clean ADT dataset
# Raw microbiology dataset
# Clean (de-duplicated) microbiology dataset
# Abx (post-initial clean)
# Abx grouped into courses
# Cleaned encounters data
#
# OUTPUT DATASETS
# 
# Joined ADT / micro dataset for further refining
# Target pathogen table (only those from org_group_3)
# Datasets with cases / controls for each target organism (no features added)
#
# Authors: Sanjat Kanjilal
# Last updated: 2025-02-17
#####################################################################################################

#### SET ENVIRONMENT VARIABLES / IMPORT FUNCTIONS ####

# Import the functions for the pipeline
source("/PHShome/zw852/colonization-pressure_HAI/HO_infxn_functions.R")
conflicted::conflicts_prefer(dplyr::n)
conflicts_prefer(dplyr::filter)

# Set working directory for imports
mainDir <- "/data/tide/projects/ho_infxn_ml/"
setwd(file.path(mainDir))

# Global paths
# ADT_filename <- "clean_data/20250411/ADT_clean.csv"
# micro_dedup_filename <- "clean_data/20250411/micro_dedup.csv"
# abx_courses_filename <- "clean_data/20250411/abx_courses.csv"
# enc_clean_filename <- "clean_data/20250411/enc_clean.csv"
# micro_filename <- "input_data/20250411/micro.ground_truth_20150525-20250131.csv"
# abx_prelim_filename <- "clean_data/20250411/abx_prelim_clean.csv"
ADT_filename <- "clean_data/20250603_SensitivityAnalysis/ADT_clean.csv"
micro_dedup_filename <- "clean_data/20250603_SensitivityAnalysis/micro_dedup.csv"
abx_courses_filename <- "clean_data/20250603_SensitivityAnalysis/abx_courses.csv"
enc_clean_filename <- "clean_data/20250603_SensitivityAnalysis/enc_clean.csv"
micro_filename <- "input_data/20250411/micro.ground_truth_20150525-20250131.csv"
abx_prelim_filename <- "clean_data/20250603_SensitivityAnalysis/abx_prelim_clean.csv"

# Set pathogen hierarchy
pathogen_hierarchy <- "org_group_3"

#### ADT / MICRO DATA JOIN ####

# Join ADT to micro data using threshold of 3 days to define hospital-onset infection
room_dat <- readr::read_csv(ADT_filename)
micro.dedup <- readr::read_csv(micro_dedup_filename)
abx.courses <- readr::read_csv(abx_courses_filename)
enc_clean <- readr::read_csv(enc_clean_filename)

adt.micro.raw <- adt_micro_initial_prep(room_dat = room_dat,
                                        micro.dedup = micro.dedup,
                                        day_threshold = 3, 
                                        abx.courses = abx.courses,
                                        enc_clean = enc_clean)

#### SAVE PRELIMINARY ADT / MICRO JOINED DATASET
readr::write_csv(adt.micro.raw, file = "clean_data/20250603_SensitivityAnalysis/adt_micro_raw.csv")

#### BUILD UNMATCHED CASE / CONTROL DATASETS ####

# Select hierarchy for looping
# org_group_1 for body niche (skin, enteric, environmental)
# org_group_2 for species specific categories (Staph_aureus, etc)
# org_group_3 for species-subtype specific categories (MRSA, MSSA, etc)
micro  <- readr::read_csv(micro_filename)

path_cat_table <- path_cat_hierarchy(micro, hierarchy = "org_group_3") 

readr::write_csv(path_cat_table, file = "clean_data/20250603_SensitivityAnalysis/path_cat_table.csv")

# Identify cases / controls for each organism in the table above
abx_prelim <- readr::read_csv(abx_prelim_filename)

flow_chart <- data.frame(
  pathogen = character(),
  initial_cases = numeric(),
  case_after_drop_subsequent = numeric(),
  case_after_drop_prior_infection = numeric(),
  stringsAsFactors = FALSE
)

cc.unmatched <- do.call(bind_rows, 
                        lapply(levels(factor(path_cat_table$pathogen_category)), function(y) 
                        {
                          
                          tic()
                          
                          print(paste("Running for", y))
                          
                          pathogen_hierarchy <- unique(path_cat_table$pathogen_hierarchy)
                          
                          # Select cases and controls 
                          unmatched_cc <- generate_unmatched_data(
                            pathogen_hierarchy = pathogen_hierarchy,
                            pathogen_category = y,
                            day_threshold = 3,
                            adt.micro.raw = adt.micro.raw,
                            abx_prelim = abx_prelim)
                          
                          # Add target name for spread
                          unmatched_cc <- cbind(run = y, unmatched_cc) %>% 
                            distinct()
                          
                          toc()
                          
                          # Ensure the last expression in the function is the dataset
                          return(unmatched_cc)
                        }))

cc.unmatched <- cc.unmatched %>%
  arrange(hospitalization_id, room_stay_id) %>%
  spread(key = run, value = group)

#### SAVE UNMATCHED COHORT DATASET ####
readr::write_csv(cc.unmatched, file = paste0("clean_data/20250603_SensitivityAnalysis/unmatched_case_controls_no_features_", pathogen_hierarchy, ".csv"))

#### CLEAN UP ####

# Remove global paths / variables
rm(mainDir, ADT_filename, micro_dedup_filename, Abx_courses_filename, enc_clean_filename, micro_filename,abx_prelim_filename, pathogen_hierarchy)

# Remove datasets
rm(room_dat, micro.dedup, abx.courses, enc_clean, adt.micro.raw, micro, path_cat_table, abx_prelim, cc.unmatched)

# Remove functions
rm(grouping_function, micro_dedup, initial_abx_clean, process_ip_abx_features, calculate_abx_ip_durations, calculate_op_abx_durations,
   calculate_abx_op_last_date, abx_courses, enc_processing, path_cat_hierarchy, admt.mapped, first_eligible_rooms, micro_prep, 
   room_dat_eligible_trimmed, adt_micro_join, time_to_infxn, drop_old_micro, remove_treated_patients, drop_prev_enc, 
   adt_micro_hierarchy, potential_cases, after_first_room, flag_prior_infxn, generate_controls, add_demographics, add_abx, 
   final_abx_clean, add_comorbidities, add_cpt, add_admit.source, add_prior_pathogens, calculate_CP_score, 
   calculate_and_add_cp_scores_parallel, adt_micro_initial_prep, generate_unmatched_data, generate.unmatched.dataset, 
   environmental.matching_process, patient.matching_process)
