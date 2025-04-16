# /usr/bin/R

#####################################################################################################
# Description: Input data pre-processing for colonization pressure / hospital-onset infection analysis
#
# DEPENDENT SCRIPTS
#
# HO_infxn_functions.R - Contains all functions in data processing scripts
# HO_infxn_micro_prep.R - Bins micro results into organism categories
#
# INPUT DATASETS
#
# Raw ADT dataset
# Raw microbiology dataset
# Raw abx dataset
# Raw encounters data
#
# OUTPUT DATASETS
#
# Microbiology - Drop duplicate cultures within 30d of the first
# ADT - Filter for room stays >48 hours
# Antibiotics - Simplify into abx classes and define start / end dates
# Encounters - Group encounters occurring within 7d of each other into single index
# Admission disposition - Mapped to categories
#
# Authors: Sanjat Kanjilal
# Last updated: 2025-02-17
#####################################################################################################

#### SET ENVIRONMENT VARIABLES / IMPORT FUNCTIONS ####
conflicted::conflicts_prefer(dplyr::row_number)
conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::lag())
conflicts_prefer(lubridate::month)
conflicts_prefer(lubridate::year)
conflicts_prefer(lubridate::week)
conflicts_prefer(lubridate::quarter)
conflicts_prefer(dplyr::first)
conflicts_prefer(reshape2::melt)
conflicts_prefer(reshape2::dcast)
conflicts_prefer(data.table::shift)

# Set working directory for imports
mainDir <- "/data/tide/projects/ho_infxn_ml/input_data/20250411/"
setwd(file.path(mainDir))

# Global paths
micro_filename <- "micro.ground_truth_20150525-20250131.csv"
adt_filename <- "MGB_ADT_20140131-20240713.csv"
abx_filename <- "MGB_abx_20150323-20240715.csv"
abx_mapping_filename <- "mappings/EDW_abx_map_20150525-20240716.csv"
admt_filename <- "MGB_IP_ED_encounters_20150430-20240715.csv"
enc_filename <- "MGB_encounters_20241105.csv"

# Import the functions for the pipeline
# source("/PHShome/sk726/Scripts/ml-hc-class/prior_patient_infection/HO_infxn_functions.R")
source("~/colonization-pressure_HAI/HO_infxn_functions.R")

#### PROCESS INPUT DATA ####

# Micro data
micro <- readr::read_csv(micro_filename) %>%
  dplyr::filter(!grepl("strep", org_group_2, ignore.case = T)) %>%
  # mutate(org_group_1 = ifelse(org_group_2 == "Strep_pneumo", "Respiratory flora", org_group_1)) %>%
  # mutate(org_group_3 = ifelse(!is.na(org_group_2) & is.na(org_group_3), org_group_2, org_group_3)) %>%
  # gather(hierarchy, strata, org_group_1:org_group_3) %>%
  dplyr::filter(!is.na(org_group_1))

micro.dedup <- micro_dedup(micro)

# ADT data
room_dat <- readr::read_csv(adt_filename) %>%
  dplyr::mutate(dept_room = paste(DepartmentDSC, RoomID)) %>%
  dplyr::mutate(duration = round(duration/24, 1)) %>% # Convert duration scale to days
  dplyr::select(PatientID:PatientClassDSC, dept_room, DepartmentDSC:duration)

# Antibiotic data

# Import raw abx data
abx <- readr::read_csv(abx_filename) %>%
  dplyr::mutate(across(where(is.character), ~ dplyr::na_if(.,""))) 

# Map major antibiotic classes
major <- readr::read_csv(abx_mapping_filename) %>% 
  dplyr::filter(major=="x")

# Pre-processing to drop prns, 'ordered but not given', multiple orders given on a same day, then group by regimen dates
abx_prelim <- initial_abx_clean(abx, major)

# Calculate IP abx course durations by abx class 
abx.ip <- process_ip_abx_features(abx_prelim) 
abx_ip_durations <- calculate_abx_ip_durations(abx.ip) 

# Calculate OP abx course durations by abx class 
abx_op_duration_raw <- calculate_op_abx_durations(abx_prelim)
abx_op_durations <- calculate_abx_op_last_date(abx_op_duration_raw) 

# Create final abx courses dataset
abx.courses <- abx_courses(abx_ip_durations, abx_op_durations)

# Encounter data
enc <- readr::read_csv(enc_filename)

enc_clean <- enc_processing(enc)

# Map admission disposition to clean categories
admt <- readr::read_csv(admt_filename)

admt_clean <- admt.mapped(admt)

#### SAVE CLEANED INPUT DATASETS FOR ADT / MICRO JOIN ####
# Set working directory for exports
mainDir <- "/data/tide/projects/ho_infxn_ml/clean_data/20250411/"
setwd(file.path(mainDir))

readr::write_csv(room_dat, file = "ADT_clean.csv")
readr::write_csv(micro.dedup, file = "micro_dedup.csv")
readr::write_csv(abx_prelim, file = "abx_prelim_clean.csv")
readr::write_csv(abx.courses, file = "abx_courses.csv")
readr::write_csv(enc_clean, file = "enc_clean.csv")
readr::write_csv(admt_clean, file = "admt_clean.csv")

#### CLEAN UP ####

# Remove global paths / variables
rm(mainDir, today, micro_filename, adt_filename, abx_filename, abx_mapping_filename,
   dems_filename, elix_filename, cpt_filename, admt_filename, location_map_filename, 
   enc_filename, pathogen_hierarchy)

# Remove datasets
rm(room_dat, micro, micro.dedup, abx, major, abx.ip, abx_ip_durations, 
   abx_op_duration_raw, abx_op_durations, abx_courses, enc, enc_clean, 
   admt, admt_clean)

# Remove functions
rm(grouping_function, micro_dedup, initial_abx_clean, process_ip_abx_features, calculate_abx_ip_durations, calculate_op_abx_durations,
   calculate_abx_op_last_date, abx_courses, enc_processing, path_cat_hierarchy, admt.mapped, first_eligible_rooms, micro_prep, 
   room_dat_eligible_trimmed, adt_micro_join, time_to_infxn, drop_old_micro, remove_treated_patients, drop_prev_enc, 
   adt_micro_hierarchy, potential_cases, after_first_room, flag_prior_infxn, generate_controls, add_demographics, add_abx, 
   final_abx_clean, add_comorbidities, add_cpt, add_admit.source, add_prior_pathogens, calculate_CP_score, 
   calculate_and_add_cp_scores_parallel, adt_micro_initial_prep, generate_unmatched_data, generate.unmatched.dataset, environmental.matching_process, 
   patient.matching_process)