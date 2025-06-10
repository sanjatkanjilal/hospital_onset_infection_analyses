# /usr/bin/R

#####################################################################################################
# Description: Add features to cohorts in colonization pressure / hospital-onset infection analysis
#
# DEPENDENT SCRIPTS
#
# HO_infxn_functions.R - Contains all functions in data processing scripts
# HO_infxn_micro_prep.R - Bins micro results into organism categories
# HO_infxn_input_data_processing.R - Pre-processes micro, antibiotic, encounter, and admit source datasets
# HO_infxn_build_unmatched_cohorts.R - Generates cohorts without features
#
# INPUT DATASETS
#
# Unmatched cohorts with no features
# Raw ADT / micro joined dataset
# Raw ADT dataset
# Raw microbiology dataset
# Clean (de-duplicated) microbiology dataset
# Abx grouped into courses
# Cleaned encounters data (grouped)
# Raw demographics dataset
# Raw comorbidities dataset (with Elixhauser index)
# Raw procedures dataset
# Department mapping file
# Pathogen categories table
#
# OUTPUT DATASETS
#
# Demographics featurized for cohort
# Antibiotics featurized for cohort
# Elixhauser comorbidities featurized for cohort
# Procedures featurized for cohort
# Admission disposition featurized for cohort
# Prior room occupant featurized for cohort
# Colonization pressure featurized for cohort
# Unmatched cohort + all features
#
# Authors: Luke Sagers, Sanjat Kanjilal
# Last updated: 2025-02-17
#####################################################################################################

#### SET ENVIRONMENT VARIABLES / IMPORT FUNCTIONS ####

# Import the functions for the pipeline
# source("/PHShome/sk726/Scripts/ml-hc-class/prior_patient_infection/HO_infxn_functions.R")
source("~/colonization-pressure_HAI/HO_infxn_functions.R")
# conflicted::conflicts_prefer(tidyr::replace_na)
# conflicted::conflicts_prefer(dplyr::mutate)
# conflicted::conflicts_prefer(dplyr::select)
# conflicted::conflicts_prefer(dplyr::distinct)
# conflicted::conflicts_prefer(dplyr::if_else)

# Set working directory for imports
mainDir <- "/data/tide/projects/ho_infxn_ml/"
setwd(file.path(mainDir))

# Global paths
unmatched_cohort_no_features <- "clean_data/20250603_SensitivityAnalysis/unmatched_case_controls_no_features_org_group_3.csv"
adt_micro_raw <- "clean_data/20250603_SensitivityAnalysis/adt_micro_raw.csv"
dems_filename <- "input_data/20250411/MGB_demographics_20240715.csv"
abx_courses <- "clean_data/20250603_SensitivityAnalysis/abx_courses.csv"
elix_filename <-"input_data/20250411/MGB_elixhauser_20240715_updated.csv"
cpt_filename <- "input_data/20250411/MGB_procedures_20240715.csv"
admt_filename <- "clean_data/20250603_SensitivityAnalysis/admt_clean.csv"
location_map_filename <- "mappings/department_mapping.csv"
enc_filename <- "clean_data/20250603_SensitivityAnalysis/enc_clean.csv"
ADT_filename <- "input_data/20250411/MGB_ADT_20140131-20240713.csv"
micro_filename <- "input_data/20250411/micro.ground_truth_20150525-20250131.csv"
path_cat_table_filename <- "clean_data/20250603_SensitivityAnalysis/path_cat_table.csv"

# Set pathogen hierarchy
pathogen_hierarchy <- "org_group_3"

#### GENERATE FEATURES FOR ADT / MICRO DATA ####

# Import dataset with unmatched cohort and no features
cc.unmatched <- readr::read_csv(unmatched_cohort_no_features)

# Import cleaned ADT / micro dataset for case / control selection
adt.micro.raw <- readr::read_csv(adt_micro_raw)

# Demographics
dems <- readr::read_csv(dems_filename)

dem_features <- add_demographics(adt.micro.raw, dems) 

# Antibiotics
abx.courses <- readr::read_csv(abx_courses)

abx_features.prelim <- add_abx(abx.courses, adt.micro.raw) 

abx_features <- final_abx_clean(abx_features.prelim, id_vars = c("PatientID", "DTS_in", "coll_datetime_UTC"))

# Comorbidities
elix <- readr::read_csv(elix_filename) %>%
  select(PatientEncounterID, elix_index_mortality, elix_AIDS:elix_CBVD) %>%
  mutate(PatientEncounterID = as.numeric(PatientEncounterID))

elix_features <- add_comorbidities(adt.micro.raw, elix)

# Procedures
cpt <- readr::read_csv(cpt_filename) %>%
  select(PatientID, proc_simple, proc_date) %>%
  distinct()

cpt_features <- add_cpt(adt.micro.raw, cpt) 

# Admission disposition
admt_mapped <- readr::read_csv(admt_filename)

admt_features <- add_admit.source(adt.micro.raw, admt_mapped) 

# Log transform of # of days since a prior pathogen noted in the room from a different patient
room_dat <- readr::read_csv(ADT_filename)
micro <- readr::read_csv(micro_filename)

prior_occupant_features <- add_prior_pathogens(room_dat, micro, pathogen_hierarchy, adt.micro.raw) 

# Colonization pressure (sum of log transform of days since ward patients had pathogen) # ~14 hour run time
location_map <- readr::read_csv(paste0(location_map_filename)) 

cp_features.prelim <- calculate_and_add_cp_scores_parallel(adt.micro.raw, micro, room_dat, location_map) 

cp_features <- cp_features.prelim %>%
  dplyr::select(PatientID, hospitalization_id, CDiff_cp, MSSA_cp, MRSA_cp, 
                DS_Entero_cp, ESBL_cp, VSE_cp, VRE_cp, DS_PsA_cp, DR_PsA_cp) %>%
  distinct()

####  SAVE FEATURE DATASETS ####
readr::write_csv(dem_features, file = paste0("clean_data/20250603_SensitivityAnalysis/dems_", pathogen_hierarchy, ".csv"))
readr::write_csv(abx_features, file = paste0("clean_data/20250603_SensitivityAnalysis/abx_", pathogen_hierarchy, ".csv"))
readr::write_csv(elix_features, file = paste0("clean_data/20250603_SensitivityAnalysis/elix_", pathogen_hierarchy, ".csv"))
readr::write_csv(cpt_features, file = paste0("clean_data/20250603_SensitivityAnalysis/cpt_", pathogen_hierarchy, ".csv"))
readr::write_csv(admt_features, file = paste0("clean_data/20250603_SensitivityAnalysis/admt_", pathogen_hierarchy, ".csv"))
readr::write_csv(prior_occupant_features, file = paste0("clean_data/20250603_SensitivityAnalysis/prior_occupant_", pathogen_hierarchy, ".csv"))
readr::write_csv(cp_features, file = paste0("clean_data/20250603_SensitivityAnalysis/col_pressure_", pathogen_hierarchy, ".csv"))

#### CREATE FULLY FEATURIZED COHORT ####

# Import clean featurized datasets as needed
# adt.micro.raw <- readr::read_csv(file = ("clean_data/20250217/adt_micro_raw.csv"))
# cc.unmatched <- readr::read_csv(file = ("clean_data/20250217/unmatched_case_controls_no_features_org_group_3.csv"))
# dem_features <- readr::read_csv(file = ("clean_data/20250217/dems_org_group_3.csv"))
# abx_features <- readr::read_csv(file = ("clean_data/20250217/abx_org_group_3.csv"))
# elix_features <- readr::read_csv(file = ("clean_data/20250217/elix_org_group_3.csv"))
# cpt_features <- readr::read_csv(file = ("clean_data/20250217/cpt_org_group_3.csv"))
# admt_features <- readr::read_csv(file = ("clean_data/20250217/admt_org_group_3.csv"))
# prior_occupant_features <- readr::read_csv(file = ("clean_data/20250217/prior_occupant_org_group_3.csv"))
# cp_features <- readr::read_csv(file = ("clean_data/20250217/col_pressure_org_group_3.csv"))

# Pick organisms with at least X number of cases in unmatched cohort
path_cat_table <- readr::read_csv(path_cat_table_filename)

# Select organisms with n>300 cases in unmatched cohorts
path_cat_table_matching <- pathogen_selection(cc.unmatched, path_cat_table)
readr::write_csv(path_cat_table_matching, file = paste0("clean_data/20250603_SensitivityAnalysis/path_cat_table_matching.csv"))

unmatched_features <- do.call(dplyr::bind_rows, 
                              lapply(levels(factor(path_cat_table_matching$pathogen_category)), function(y) 
                              {
                                
                                tic()
                                
                                print(paste("Adding features for", y))
                                
                                # Select cases and controls 
                                unmatched_cc_features <- generate.unmatched.dataset(
                                  pathogen_category = y,
                                  adt.micro.raw = adt.micro.raw,
                                  cc.unmatched = cc.unmatched,
                                  dem_features = dem_features,
                                  abx_features = abx_features,
                                  elix_features = elix_features,
                                  cpt_features = cpt_features,
                                  admt_features = admt_features,
                                  prior_occupant_features = prior_occupant_features,
                                  cp_features = cp_features)
                                
                                toc()
                                
                                # Ensure the last expression in the function is the flattened dataset
                                
                                return(unmatched_cc_features)
                                
                              }))

# Replace missing values with 0 
unmatched_features <- unmatched_features %>%
  dplyr::mutate(across(.cols = matches("0_60|60_plus|sx|surgery|prior|elix"), ~replace(., is.na(.), 0))) 

# Clean up columns
unmatched_features <- unmatched_features %>%
  dplyr::select(run, group, PatientID, PatientEncounterID, dept_code, hospitalization_id, room_stay_id, 
         duration, time_to_infxn, matching_duration, DTS_in, DTS_in_month, DTS_in_year, 
         age, sex:ethnicity, any_abx_0_60, any_abx_60_plus, aminoglycoside_0_60:other_abx_0_60, 
         elix_index_mortality:admit_source_clean, prior_C_diff:DR_PsA_cp) %>% 
  dplyr::distinct()

#### SAVE UNMATCHED COHORTS WITH FEATURES ADDED ####
readr::write_csv(unmatched_features, file = paste0("clean_data/20250603_SensitivityAnalysis/unmatched_case_controls_features_", pathogen_hierarchy, ".csv"))

#### CLEAN UP ####

# Remove global paths / variables
rm(mainDir, unmatched_cohort_no_features, adt_micro_raw, micro_dedup_filename, dems_filename, abx_courses_filename, 
   abx_prelim_filename, elix_filename, enc_clean_filename, cpt_filename, admt_filename, location_map_filename, 
   enc_filename, ADT_filename, micro_filename, path_cat_table_filename, pathogen_hierarchy)

# Remove datasets
rm(cc.unmatched, adt.micro.raw, micro.dedup, dems, dem_features, abx.courses, abx_features.prelim, abx_features, 
   elix, elix_features, enc_clean, cpt, cpt_features, admt_mapped, admt_features, room_dat, micro, prior_occupant_features, 
   location_map, cp_features.prelim, cp_features, path_cat_table, path_cat_table_matching, unmatched_features)

# Remove functions
rm(grouping_function, micro_dedup, initial_abx_clean, process_ip_abx_features, calculate_abx_ip_durations, 
   calculate_op_abx_durations, calculate_abx_op_last_date, abx_courses, enc_processing, path_cat_hierarchy,
   admt.mapped, first_eligible_rooms, micro_prep, room_dat_eligible_trimmed, adt_micro_join, time_to_infxn, 
   drop_old_micro, remove_treated_patients, drop_prev_enc, adt_micro_hierarchy, potential_cases, after_first_room, 
   flag_prior_infxn, generate_controls, add_demographics, add_abx, final_abx_clean, add_comorbidities, add_cpt, 
   add_admit.source, add_prior_pathogens, calculate_CP_score, calculate_and_add_cp_scores_parallel, 
   pathogen_selection, adt_micro_initial_prep, generate_unmatched_data, generate.unmatched.dataset, 
   environmental.matching_process, patient.matching_process)

