# /usr/bin/R

#####################################################################################################
# Description: Script to build case / control cohorts for hospital-onset infection analysis: 
# Dependencies: HO_infxn_micro_prep.R, HO_infxn_functions.R
# Authors: Luke Sagers, Ziming (Alex) Wei, Sanjat Kanjilal
# Last updated: 2024-12-16
#####################################################################################################

##### LIBRARIES #####
library(tictoc)
library(tidyverse)
library(data.table)
library(reshape2)
library(MatchIt)
library(conflicted)
library(lubridate)
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

#### SET ENVIRONMENT VARIABLES ####

# Set working directory and output directory to be whatever you generate today
mainDir <- "~/ho_infxn_clean_code/"
setwd(file.path(mainDir))

# Import the functions for the pipeline
source("~/ho_infxn_clean_code/data/HO_infxn_functions.R")

# Global paths
micro_filename <- "micro.ground_truth_20150525-20240701.csv"
adt_filename <- "MGB_ADT_20140131-20240713.csv"
abx_filename <- "MGB_abx_20150323-20240715.csv"
abx_mapping_filename <- "EDW_abx_map_20150525-20240716.csv"
dems_filename <- "MGB_demographics_20240715.csv"
# elix_filename <-"MGB_elixhauser_20240715.csv"
elix_filename <- "MGB_elixhauser_20241105.csv"
cpt_filename <- "MGB_procedures_20240715.csv"
admt_filename <- "MGB_IP_ED_encounters_20150430-20240715.csv"
location_map_filename <- "department_mapping.csv"
# enc_filename <- "MGB_encounters_20240715.csv"
enc_filename <- "MGB_encounters_20241105.csv"

#### IMPORT AND PREP DATASETS FOR COHORT BUILDING PIPELINE ####

# Room data
room_dat <- read_csv(paste0("/data/tide/data/edw/ADT/", adt_filename)) %>%
  mutate(dept_room = paste(DepartmentDSC, RoomID)) %>%
  mutate(duration = round(duration/24, 1)) %>% # Convert duration scale to days
  select(PatientID:PatientClassDSC, dept_room, DepartmentDSC:duration)

# Micro data
micro <- read_csv(paste0("input_data/", micro_filename)) %>%
  filter(!grepl("strep", org_group_2, ignore.case = T)) %>%
  # mutate(org_group_1 = ifelse(org_group_2 == "Strep_pneumo", "Respiratory flora", org_group_1)) %>%
  # mutate(org_group_3 = ifelse(!is.na(org_group_2) & is.na(org_group_3), org_group_2, org_group_3)) %>%
  # gather(hierarchy, strata, org_group_1:org_group_3) %>%
  filter(!is.na(org_group_1))

micro.dedup <- micro_dedup(micro)

# Antibiotic data
abx <- read_csv(paste0("/data/tide/data/edw/medications/", abx_filename)) %>%
  mutate(across(where(is.character), ~ na_if(.,""))) 

# Map major antibiotic classes
major <- read_csv(paste0("/data/tide/data/mappings/edw/", abx_mapping_filename)) %>% 
  filter(major=="x")

# Pre-processing to drop prns, 'ordered but not given', multiple orders given on a same day, then group by regimen dates
abx_prelim <- initial_abx_clean(abx, major)

# Calculate IP abx course durations by abx class 
abx.ip <- process_ip_abx_features(abx_prelim) 
abx_ip_durations <- calculate_abx_ip_durations(abx.ip) 

# Calculate OP abx course durations by abx class 
abx_op_duration_raw <- calculate_op_abx_durations(abx_prelim)
abx_op_durations <- calculate_abx_op_last_date(abx_op_duration_raw) 

# Create final abx courses dataset
abx.courses <-  abx_courses(abx_ip_durations, abx_op_durations)

rm(abx, major, abx.ip, abx_ip_durations, abx_op_duration_raw, abx_op_durations)

# Process encounter data
enc <- data.table::fread(paste0("/data/tide/data/edw/encounters/", enc_filename))

enc_clean <- enc_processing(enc)

rm(enc)

# Save clean input datasets for ADT / micro join
write_csv(micro.dedup, file = paste0("clean_data/",  today, "/micro_dedup.csv"))
write_csv(abx.courses, file = paste0("clean_data/",  today, "/abx_courses.csv"))
write_csv(enc_clean, file = paste0("clean_data/",  today, "/enc_clean.csv"))

#### ADT / MICRO DATA JOIN ####
adt.micro.raw <- adt_micro_initial_prep(room_dat = room_dat,
                                        micro.dedup = micro.dedup,
                                        day_threshold = 3, 
                                        abx.courses = abx.courses,
                                        enc_clean = enc_clean)

write_csv(adt.micro.raw, file = paste0("clean_data/",  today, "/adt_micro_raw.csv"))

#### BUILD UNMATCHED CASE / CONTROL DATASETS ####

# Select hierarchy for looping
# org_group_1 for body niche (skin, enteric, environmental)
# org_group_2 for species specific categories (Staph_aureus, etc)
# org_group_3 for species-subtype specific categories (MRSA, MSSA, etc)
path_cat_table <- path_cat_hierarchy(micro, hierarchy = "org_group_3") 

pathogen_hierarchy <- "org_group_3"

write_csv(path_cat_table, file = paste0("clean_data/",  today, "/path_cat_table.csv"))

# Identify cases / controls for each organism in the table above
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
                
# Write to .csv
write_csv(cc.unmatched, file = paste0("clean_data/",  today, "/unmatched_case_controls_no_features_", pathogen_hierarchy, ".csv"))

#### ADD FEATURES TO CLEAN ADT / MICRO DATA ####

# Filter for organisms / groups with >=300 cases
path_cat_table_matching <- cc.unmatched %>%
  gather(run, group, C_diff:VSE_faecium) %>%
  count(run, group) %>%
  arrange(-n) %>%
  filter(group == "case") %>%
  filter(n >= 300) %>%
  select(pathogen_category = run) %>%
  mutate(pathogen_category = factor(pathogen_category)) %>%
  left_join(path_cat_table)

# Import cleaned ADT / micro dataset for case / control selection
adt.micro.raw <- read_csv(file = paste0("clean_data/",  today, "/adt_micro_raw.csv"))

# Demographics
dems <- read_csv(paste0("/data/tide/data/edw/demographics/", dems_filename))

dem_features <- add_demographics(adt.micro.raw, dems) 

write_csv(dem_features, file = paste0("clean_data/",  today, "/dems_", pathogen_hierarchy, ".csv"))

rm(dems)

# Antibiotics
abx_features.prelim <- add_abx(abx.courses, adt.micro.raw) 

abx_features <- final_abx_clean(abx_features.prelim, id_vars = c("PatientID", "DTS_in", "coll_datetime_UTC"))
  
write_csv(abx_features, file = paste0("clean_data/",  today, "/abx_", pathogen_hierarchy, ".csv"))

rm(abx_features.prelim)

# Comorbidities
elix <- read_csv(paste0("/data/tide/data/edw/diagnoses/", elix_filename)) %>%
  rename("indiv_score" = "elix_index_mortality") %>%
  select(PatientEncounterID, indiv_score, elix_AIDS:elix_CBVD) %>%
  mutate(PatientEncounterID = as.numeric(PatientEncounterID))

elix_features <- add_comorbidities(adt.micro.raw, elix)

write_csv(elix_features, file = paste0("clean_data/",  today, "/elix_", pathogen_hierarchy, ".csv"))

rm(elix)

# Procedures
cpt <- data.table::fread(paste0("/data/tide/data/edw/procedures/", cpt_filename)) %>%
  select(PatientID, proc_simple, proc_date) %>%
  distinct()

cpt_features <- add_cpt(adt.micro.raw, cpt) 

write_csv(cpt_features, file = paste0("clean_data/",  today, "/cpt_", pathogen_hierarchy, ".csv"))

rm(cpt)

# Admission disposition
admt <- data.table::fread(paste0("/data/tide/data/edw/encounters/", admt_filename))

admt_mapped <- admt.mapped(admt)

admt_features <- add_admit.source(adt.micro.raw, admt_mapped) 

write_csv(admt_features, file = paste0("clean_data/",  today, "/admt_", pathogen_hierarchy, ".csv"))

rm(admt, admt_mapped)

# Log transform of # of days since a prior pathogen noted in the room from a different patient
prior_occupant_features <- add_prior_pathogens(room_dat, micro, pathogen_hierarchy, adt.micro.raw) 

write_csv(prior_occupant_features, file = paste0("clean_data/",  today, "/prior_occupant_", pathogen_hierarchy, ".csv"))

# Colonization pressure (sum of log transform of days since ward patients had pathogen) # ~14 hour run time
location_map <- read_csv(paste0("mappings/", location_map_filename)) 

cp_features.prelim <- calculate_and_add_cp_scores_parallel(adt.micro.raw, micro, room_dat, location_map) 

cp_features <- cp_features.prelim %>%
  select(PatientID, hospitalization_id, CDiff_cp, MSSA_cp, MRSA_cp, 
         DS_Entero_cp, ESBL_cp, VSE_cp, VRE_cp, DS_PsA_cp, DR_PsA_cp) %>%
  distinct()

write_csv(cp_features, file = paste0("clean_data/",  today, "/col_pressure_", pathogen_hierarchy, ".csv"))

#### CREATE FULLY FEATURIZED DATASET ####

# Import datasets as needed
# micro <- read_csv(file = paste0("clean_data/",  today, "/micro.csv"))
# path_cat_table <- read_csv(file = paste0("clean_data/", today, "/path_cat_table.csv"))
# adt.micro.raw <- read_csv(file = paste0("clean_data/", today, "/adt_micro_raw.csv"))
# cc.unmatched <- read_csv(file = paste0("clean_data/",  today, "/unmatched_case_controls_no_features_org_group_3.csv"))
# dem_features <- read_csv(file = paste0("clean_data/",  today, "/dems_org_group_3.csv"))
# abx_features <- read_csv(file = paste0("clean_data/",  today, "/abx_org_group_3.csv"))
# elix_features <- read_csv(file = paste0("clean_data/",  today, "/elix_org_group_3.csv"))
# cpt_features <- read_csv(file = paste0("clean_data/",  today, "/cpt_org_group_3.csv"))
# admt_features <- read_csv(file = paste0("clean_data/",  today, "/admt_org_group_3.csv"))
# prior_occupant_features <- read_csv(file = paste0("clean_data/20240715/prior_occupant_org_group_3.csv"))
# cp_features <- read_csv(file = paste0("clean_data/",  today, "/col_pressure_org_group_3.csv"))

path_cat_table_matching <- cc.unmatched %>%
  gather(run, group, C_diff:VSE_faecium) %>%
  count(run, group) %>%
  arrange(-n) %>%
  filter(group == "case") %>%
  filter(n >= 300) %>%
  select(pathogen_category = run) %>%
  mutate(pathogen_category = factor(pathogen_category)) %>%
  left_join(path_cat_table)

unmatched_features <- do.call(bind_rows, 
                              lapply(levels(factor(path_cat_table_matching$pathogen_category)), function(y) 
                                {
                                
                                tic()
                                
                                print(paste("Adding features for", y))
                                
                                # Select cases and controls 
                                unmatched_cc_features <- dataset.build(
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
  mutate(across(.cols = matches("0_60|60_plus|sx|surgery|prior|elix"), ~replace(., is.na(.), 0))) 

unmatched_features <- unmatched_features %>%
  select(run, group, PatientID, PatientEncounterID, dept_code, hospitalization_id, room_stay_id, 
         duration, time_to_infxn, matching_duration, DTS_in, DTS_in_month, DTS_in_year, 
         age, sex:ethnicity, any_abx_0_60, any_abx_60_plus, aminoglycoside_0_60:other_abx_0_60, 
         indiv_score:admit_source_clean, prior_C_diff:DR_PsA_cp) %>% 
  distinct()

# Save unmatched dataset + features
write_csv(unmatched_features, file = paste0("clean_data/",  today, "/unmatched_case_controls_features_", pathogen_hierarchy, ".csv"))

#### BUILD MATCHED DATASETS ####

# Import datasets as needed
# unmatched_features <- read_csv(file = paste0("clean_data/",  today, "/unmatched_case_controls_features_", pathogen_hierarchy, ".csv"))

## Filter out patients >90 yo
unmatched_features <- unmatched_features %>% filter(age <= 89) %>% filter(age >= 18)

## Filter out patients with negative time to infection
unmatched_features <- unmatched_features %>% filter(time_to_infxn >= 0)

environment.matched <- do.call(bind_rows,
                           lapply(levels(factor(path_cat_table_matching$pathogen_category)), function(y) 
                               
                             {
                               print(paste("Matching cases to controls for environmental-level analysis for", y))
                               
                               pathogen_hierarchy <- unique(path_cat_table_matching$pathogen_hierarchy)
                               
                               unmatched_features.pathogen <- unmatched_features %>%
                                 filter(run == y)
                               
                               # Add matching features to case / control data
                               tic()
                               
                               matched_features.pathogen <- environmental.matching_process(pathogen_category = y,
                                                                                     dat = unmatched_features.pathogen)
                               
                               toc()
                               
                               # Compile data across each pathogen category
                               matched_cc <- as.data.frame(cbind(matched_features.pathogen))
                               
                               }))

patient.matched <- do.call(bind_rows,
                               lapply(levels(factor(path_cat_table_matching$pathogen_category)), function(y) 
                                 
                               {
                                 print(paste("Matching cases to controls for patient-level analysis for", y))
                                 
                                 pathogen_hierarchy <- unique(path_cat_table_matching$pathogen_hierarchy)
                                 
                                 unmatched_features.pathogen <- unmatched_features %>%
                                   filter(run == y)
                                 
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
  mutate(match = "environmental")

patient.matched <- patient.matched %>%
  mutate(match = "patient")

matched.final <- bind_rows(environment.matched, patient.matched) %>%
  select(match, run:group, group_index, PatientID:DR_PsA_cp) %>%
  distinct()

# Save matched patient-level dataset + features
# write_csv(matched.final, file = paste0("clean_data/",  today, "/final_cohort_", pathogen_hierarchy, ".csv"))
write_csv(matched.final, file = paste0("clean_data/",  today, "/final_cohort_age_filtered_", pathogen_hierarchy, ".csv"))

#### FINAL DATASET CLEAN / PREP FOR MODELS####

# Check sample sizes
sample_sizes <- matched.final %>%
  select(match:PatientID) %>%
  distinct() %>%
  group_by(match, run, group) %>%
  summarize(sample_size = n()) %>%
  ungroup()

# Check for missingness
percent_missing_values <- matched.final %>%
  select(-group_index) %>%
  group_by(match, run, group) %>%
  summarise_all(~sum(is.na(.)) / n() * 100) %>%
  pivot_longer(cols = c(PatientID:DR_PsA_cp), names_to = "variable", values_to = "percent_missing") %>%
  filter(percent_missing > 0) %>%
  arrange(variable, match, run, group)

# Check for empty columns
colSums(is.na(matched.final))

# Binarize group
matched.final.models <- matched.final %>%
  mutate(group_binary = ifelse(group == "case", 1, 0))

# Select a limited variable pool
abx_of_interest <- c("penicillin_0_60", "extended_spectrum_penicillin_0_60", "cephalosporin_0_60", 
                     "extended_spectrum_cephalosporin_0_60", "carbapenem_0_60", "anti_staph_beta_lactam_0_60", 
                     "fluoroquinolone_0_60", "glycopeptide_0_60", "anti_anaerobe_0_60", "anti_Cdiff_0_60",
                     "tetracycline_0_60", "macrolide_0_60", "sulfonamide_0_60", "lincosamide_0_60", "any_abx_0_60")

variables_to_drop <- names(matched.final.models)[grepl("_0_60$", names(matched.final.models)) & !(names(matched.final.models) %in% abx_of_interest)]

matched.final.models <- matched.final.models %>% 
  select(match:group, group_index, group_binary, PatientID, DTS_in_month, DTS_in_year, duration,
         time_to_infxn, age, sex, indiv_score, any_surgery, admit_source_clean, 
         any_abx_0_60:other_abx_0_60, prior_C_diff:DR_PsA_cp) %>%
  select(-all_of(variables_to_drop)) %>%
  select(-matches("60_plus"))

# Drop remaining NAs
matched.final.models <- matched.final.models %>%
  drop_na()

# Impute the remaining missing values with MICE (Elixhauser scores)
# imp <- mice(matched.final)
# matched.final <- complete(imp)

# Recheck missingness
colSums(is.na(matched.final.models))

today = '20241210'

# Save file
write_csv(sample_sizes,file = paste0("clean_data/", today, "/sample_sizes_", today, ".csv"))

# write_csv(matched.final.models, file = paste0("clean_data/", today, "/both_matched_final_dataset_for_models_", today, ".csv"))
write_csv(matched.final.models, file = paste0("clean_data/", today, "/age_filtered_final_dataset_for_models_", today, ".csv"))
write_csv(matched.final.models, file = '~/hospital_onset_personal_folder/patient_below90_final_dataset.csv')

#### CLEAN UP ####

# Remove global paths / variables
rm(mainDir, today, micro_filename, adt_filename, abx_filename, abx_mapping_filename, dems_filename, elix_filename, cpt_filename,
   admt_filename, location_map_filename, enc_filename, pathogen_hierarchy)

# Remove global datasets
rm(path_cat_table, path_cat_table_matching, location_map, abx_prelim, abx.courses, abx_features, 
   room_dat, micro, micro.dedup, admt_features, adt.micro.raw, cc.unmatched, cp_features, 
   cp_features.prelim, cpt_features, dem_features, elix_features, enc_clean, environment.matched, 
   matched.final, patient.matched, prior_occupant_features, unmatched_features, abx_of_interest, 
   variables_to_drop, percent_missing_values, sample_sizes)

# Remove functions
rm(grouping_function, micro_dedup, initial_abx_clean, process_ip_abx_features, calculate_abx_ip_durations, calculate_op_abx_durations,
   calculate_abx_op_last_date, abx_courses, enc_processing, path_cat_hierarchy, admt.mapped, first_eligible_rooms, micro_prep, 
   room_dat_eligible_trimmed, adt_micro_join, time_to_infxn, matching_duration, drop_old_micro, remove_treated_patients, drop_prev_enc, 
   adt_micro_hierarchy, potential_cases, after_first_room, flag_prior_infxn, generate_controls, add_demographics, add_abx, 
   final_abx_clean, add_comorbidities, add_cpt, add_admit.source, add_prior_pathogens, calculate_CP_score, 
   calculate_and_add_cp_scores_parallel, adt_micro_initial_prep, generate_unmatched_data, dataset.build, environmental.matching_process, 
   patient.matching_process)


