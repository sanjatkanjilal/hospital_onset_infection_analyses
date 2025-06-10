# /usr/bin/R

#####################################################################################################
# Description: Functions for pre-processing and cohort building for 
# colonization pressure / hospital-onset infection analysis
#
# DEPENDENT SCRIPTS
# None
#
# INPUT DATASETS
# None
#
# OUTPUT DATASETS
# None
#
# FUNCTION SETS
# SET 1: HELPER FUNCTIONS
# SET 2: INPUT DATASET PRE-PROCESSING
# SET 3: BUILD ADT - MICRO GROUND TRUTH DATASET FOR ALL PATHOGENS
# SET 4: PATHOGEN-SPECIFIC CASE / CONTROL SELECTION
# SET 5: GENERATE FEATURES FOR MATCHING AND MODELING
#
# METAFUNCTIONS
# METAFUNCTION 1: PREP ADT AND MICRO DATA (NOT PATHOGEN-SPECIFIC)
# METAFUNCTION 2: GENERATE UNMATCHED CASES & CONTROLS (PATHOGEN-SPECIFIC)
# METAFUNCTION 3: CREATION OF FINAL DATASET
# METAFUNCTION 4: ENVIRONMENTAL-LEVEL ANALYSIS: MATCHING ALGORITHM
# METAFUNCTION 5: PATIENT-LEVEL ANALYSIS: MATCHING ALGORITHM
#
# Author: Luke Sagers, Sanjat Kanjilal
# Last updated: 2025-02-17
#####################################################################################################

#### LIBRARIES ####
library(data.table)
# library(tidytable)
library(tidyverse)
library(reshape2)
library(MatchIt)
library(tictoc)
library(conflicted)
library(parallel)
conflicts_prefer(base::`%in%`)

#### FUNCTION SET 1: HELPER FUNCTIONS ####

# Function to group room stays within a specific time interval (in days) with each other
# 'type' refers to the user-defined name for the grouping (ie episodes, infections; defined below)
grouping_function <- function(data, id, index, start, end, time_interval) 
{
  
  data <- data %>%  
    dplyr::arrange({{id}}, {{index}}) %>%
    dplyr::mutate(group = ifelse({{index}} == 1, 1, NA)) %>%
    dplyr::select({{id}}, {{index}}, {{start}}, {{end}}, group) %>%
    dplyr::group_by({{id}}) %>%
    dplyr::mutate(tempdate = dplyr::lag({{end}})) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tempdiff = interval(tempdate, {{start}})/days(x = 1)) %>%
    dplyr::group_by({{id}})
  
  print("Calculation of time intervals complete")
  
  group1 <- data %>% 
    dplyr::filter(!is.na(group))
  
  mydata <- tibble()
  
  while(sum(is.na(data$group)) > 0) {
    
    # print(dplyr::filter(data,PatientID=='Z10051770'), n=41)
    # print(dplyr::filter(data,is.na(group)))
    
    if (sum(is.na(data$group)) == 0)  {
      
      break
      
    } else {
      
      group.temp <- data %>% 
        dplyr::mutate(lag_group = dplyr::lag(group)) %>%
        dplyr::filter(is.na(group)) 
      
      group.temp1 <- group.temp %>%
        dplyr::mutate(group = dplyr::if_else(tempdiff[1] <= time_interval, lag_group[1], lag_group[1] + 1)) %>%
        dplyr::mutate(group = dplyr::if_else({{index}} == {{index}}[1], group[1], NA))
      
      group.new <- group.temp1 %>%
        dplyr::filter(!is.na(group))
      
      data <- group.temp1
      
      mydata <- dplyr::bind_rows(mydata, group.new)
      
      print("Iteration complete")
      
    }
  }
  
  output <- dplyr::full_join(mydata, group1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-tempdate, -tempdiff, -lag_group) %>% 
    dplyr::arrange({{id}}, {{start}})
  
  print("Grouping complete")
  
  return(output)
}

#### FUNCTION SET 2:INPUT DATASET PRE-PROCESSING ####

# MICROBIOLOGY

# Drop repeat micro specimens with the same organisms drawn within 30 days of a first positive
micro_dedup <- function(micro)
{
  micro.dedup <- micro %>%
    # head(1000) %>%
    dplyr::arrange(PatientID, org_group_3, coll_datetime_UTC) %>%
    dplyr::group_by(PatientID, org_group_3) %>%
    dplyr::mutate(sample_index = row_number(), 
                  prev_sample = dplyr::lag(coll_datetime_UTC)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(micro_diff = difftime(coll_datetime_UTC, prev_sample, units = "days"),
                  micro_diff = round(micro_diff, digits = 1),
                  flag = ifelse(micro_diff <= 30, "X", NA)) %>%
    dplyr::filter(is.na(flag))
  
  return(micro.dedup)
  
}

# ANTIBIOTICS

# Process IP antibiotic data prior to joining to case / control data
initial_abx_clean <- function(abx, major)
{
  abx_prelim <- abx %>%
    # head(10000) %>%
    dplyr::filter(!grepl("need|prn|on call", order_freq, ignore.case = T, perl = T)) %>%
    dplyr::filter(!(is.na(order_dose) & is.na(order_unit) & is.na(order_freq) & is.na(order_quantity))) %>%
    dplyr::mutate(order_DT = floor_date(order_DTS, unit = "days"),
                  order_DT = ymd(order_DT),
                  taken_DT = floor_date(taken_DTS, unit = "days"),
                  taken_DT = ymd(taken_DT)) %>%
    dplyr::group_by(patientID, order_DT, antibiotic) %>%
    dplyr::mutate(index_order = row_number()) %>% 
    dplyr::filter(index_order == 1) %>%
    dplyr::group_by(patientID, taken_DT, antibiotic) %>%
    dplyr::mutate(index_taken = row_number()) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(drop = ifelse(!is.na(taken_DT) & index_taken > 1, "X", NA)) %>%
    dplyr::filter(is.na(drop)) %>%
    dplyr::select(-drop) %>%
    dplyr::filter(antibiotic %in% major$antibiotic) %>% 
    dplyr::left_join(major %>% dplyr::select(antibiotic, drug_class), by = c('antibiotic' = 'antibiotic')) %>%
    dplyr::mutate(drug_class = ifelse(antibiotic == "VANCOMYCIN" & route_abbr == "PO", "anti_Cdiff", drug_class)) %>%
    dplyr::filter(!is.na(drug_class)) %>%
    dplyr::arrange(patientID, antibiotic, order_DT, taken_DT) %>%
    dplyr::select(patientID, order_DT, taken_DT, antibiotic, drug_class, order_freq, order_quantity, route_abbr) %>%
    dplyr::distinct()
  
  return(abx_prelim)
}

# Antibiotic features
process_ip_abx_features <- function(abx_prelim) 
{
  tic()
  abx.ip.prelim <- as.data.table(abx_prelim)
  abx.ip.prelim <- abx.ip.prelim[!is.na(taken_DT)]
  abx.ip.prelim <- abx.ip.prelim[order(patientID, drug_class, taken_DT)]
  abx.ip.prelim[, id := seq_len(.N), by = c("patientID", "taken_DT", "drug_class")]
  abx.ip.prelim <- abx.ip.prelim[id == 1]
  abx.ip <- unique(abx.ip.prelim[, .(patientID, taken_DT, antibiotic, drug_class)])
  abx.ip <- as.data.frame(abx.ip)
  toc()
  
  return(abx.ip)
}

calculate_abx_ip_durations <- function(abx.ip) 
{
  # Calculate time intervals and group in 7 day courses
  abx.ip.prelim <- as.data.table(abx.ip)
  abx.ip.prelim[, 
                dt_lag := data.table::shift(.SD, 1), 
                by = .(patientID, drug_class), 
                .SDcols="taken_DT"]
  
  abx.ip.prelim <- abx.ip.prelim %>%
    dplyr::mutate(dt_lag = dplyr::if_else(is.na(dt_lag), taken_DT, dt_lag)) %>%
    dplyr::mutate(day_diff = taken_DT - dt_lag) %>%
    dplyr::select(-dt_lag) %>% 
    dplyr::mutate(new_prescription = dplyr::case_when(day_diff == 0 ~ 1, 
                                                      day_diff >= 7 ~ 1)) %>% 
    dplyr::group_by(patientID, drug_class) %>% 
    dplyr::mutate(new_prescription_1 = cumsum(!is.na(new_prescription))) %>%
    dplyr::ungroup()
  
  # Get start and end dates
  abx.ip.prelim <- as.data.table(abx.ip.prelim)
  abx.ip.prelim[, 
                `:=`(start = min(taken_DT), end = max(taken_DT)), 
                by = .(patientID, drug_class, new_prescription_1)]
  
  # Calculate durations and keep only summary row
  abx_ip_durations <- abx.ip.prelim %>% 
    dplyr::mutate(abx_duration_days = as.numeric(end - start + 1)) %>%
    tidytable::slice_tail(n=1, .by = c(patientID, drug_class, new_prescription_1)) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(end_DT = taken_DT) %>% 
    dplyr::select(patientID, drug_class, end_DT, abx_duration_days) %>%
    dplyr::distinct()
  
  abx_ip_durations <- as.data.frame(abx_ip_durations)
  
  return(abx_ip_durations)
}

# Calculate durations for OP antibiotic courses
calculate_op_abx_durations <- function(abx_prelim) 
{
  abx_op_duration_raw <- as.data.table(abx_prelim)
  abx_op_duration_raw <- abx_op_duration_raw[is.na(taken_DT)]
  abx_op_duration_raw <- abx_op_duration_raw[order(patientID, drug_class, order_DT)]
  abx_op_duration_raw[, id := seq_len(.N), by = c("patientID", "order_DT", "drug_class")]
  abx_op_duration_raw <- abx_op_duration_raw[id == 1]  
  abx_op_duration_raw <- unique(abx_op_duration_raw[, .(patientID, order_DT, antibiotic, drug_class, order_freq, order_quantity)])
  
  abx_op_duration_raw <- abx_op_duration_raw %>% 
    dplyr::mutate(numeric_frequency_per_day = dplyr::case_when(
      stringr::str_detect(tolower(order_freq), "(?i)^every\\s+\\d{1,2}\\s+hours") ~ 24/as.numeric(stringr::str_match(tolower(order_freq), "(?i)^every\\s+(\\d{1,2})\\s+hours")[,2]),
      stringr::str_detect(tolower(order_freq), "(?i)once") ~ 1,
      stringr::str_detect(tolower(order_freq), "(?i)\\d+\\s+times\\s+daily") ~ as.numeric(stringr::str_match(tolower(order_freq), "(?i)(\\d+)\\s+times\\s+daily")[,2]),
      stringr::str_detect(tolower(order_freq), "(?i)\\d+\\s+times\\s+weekly") ~ as.numeric(stringr::str_match(tolower(order_freq), "(?i)(\\d+)\\s+times\\s+weekly")[,2])/7,
      stringr::str_detect(tolower(order_freq), "(?i)^every\\s+\\d{1,2}\\s+days") ~ 1/as.numeric(stringr::str_match(tolower(order_freq), "(?i)^every\\s+(\\d{1,2})\\s+days")[,2]),
      tolower(order_freq) %in% c("daily", "daily ") ~ 1,
      tolower(order_freq) %in% c("nightly", "nightly ") ~ 1,
      tolower(order_freq) %in% c("daily with breakfast", "daily with breakfast ", "every morning", "every morning ", "every night", "every night ", "every evening",
                                 "daily before breakfast", "daily with dinner", "daily after dinner", "daily with lunch", "every afternoon", "daily before dinner") ~ 1,
      tolower(order_freq) %in% c("every other day") ~ 0.5,
      tolower(order_freq) %in% c("every week", "every week ", "weekly") ~ 1/7,
      ((is.na(order_freq) | order_freq == "See admin instructions") & antibiotic %in% c("AZITHROMYCIN", "FLUCONAZOLE", "LEVOFLOXACIN", "NITROFURANTOIN", "FOSFOMYCIN")) ~ 1,
      ((is.na(order_freq) | order_freq == "See admin instructions") & antibiotic %in% c("AMOXICILLIN", "TRIMETHOPRIM-SULFAMETHOXAZOLE", "AMOXICILLIN-CLAVULANATE", 
                                                                                        "DOXYCYCLINE", "CIPROFLOXACIN", "CEFDINIR", "CEFPODOXIME", "MINOCYCLINE", "CEFADROXIL", "LINEZOLID")) ~ 2,
      ((is.na(order_freq) | order_freq == "See admin instructions") & antibiotic %in% c("CEPHALEXIN", "CLINDAMYCIN", "METRONIDAZOLE")) ~ 3,
      TRUE ~ NA_real_
    )) %>%
    dplyr::mutate(numeric_frequency_per_day = round(numeric_frequency_per_day, 1)) %>% 
    dplyr::filter(numeric_frequency_per_day > 0) %>% 
    dplyr::mutate(quantity_numeric = dplyr::case_when(
      stringr::str_detect(tolower(order_quantity), "(?i)\\d+\\s+tablet") ~ as.numeric(stringr::str_match(tolower(order_quantity), "(?i)(\\d+)\\s+tablet")[,2]),
      stringr::str_detect(tolower(order_quantity), "(?i)\\d+\\s+capsule") ~ as.numeric(stringr::str_match(tolower(order_quantity), "(?i)(\\d+)\\s+capsule")[,2]),
      stringr::str_detect(tolower(order_quantity), "(?i)\\d+\\s+packet") ~ as.numeric(stringr::str_match(tolower(order_quantity), "(?i)(\\d+)\\s+packet")[,2]),
      stringr::str_detect(tolower(order_quantity), "(?i)\\d+\\s+dose") ~ as.numeric(stringr::str_match(tolower(order_quantity), "(?i)(\\d+)\\s+dose")[,2])
    )) %>% 
    dplyr::mutate(quantity_numeric = ifelse(is.na(quantity_numeric), 5, quantity_numeric)) %>%
    dplyr::mutate(abx_duration_days = round(quantity_numeric / numeric_frequency_per_day, 1)) %>% 
    dplyr::filter(!is.na(abx_duration_days)) %>%
    dplyr::mutate(end_DT = order_DT + abx_duration_days) %>%
    dplyr::select(patientID, drug_class, order_DT, end_DT, abx_duration_days) %>%
    dplyr::distinct()
  
  abx_op_duration_raw <- as.data.frame(abx_op_duration_raw)
  
  return(abx_op_duration_raw)
}

# Group courses within 7 days of each other and calculate last taken date
calculate_abx_op_last_date <- function(abx_op_duration_raw) 
{
  
  tic()
  
  print("Removing duplicate orders")
  abx_op_durations.prelim <- abx_op_duration_raw %>% 
    dplyr::arrange(patientID, drug_class, order_DT, end_DT) %>%
    tidyr::unite(key, patientID, drug_class, remove = F, sep = "-") %>%
    dplyr::group_by(key, drug_class, order_DT) %>%
    dplyr::mutate(course_index_raw = row_number()) %>%
    dplyr::filter(course_index_raw == max(course_index_raw)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-course_index_raw) %>%
    dplyr::distinct() %>%
    dplyr::group_by(key) %>%
    dplyr::mutate(course_index = row_number()) %>%
    dplyr::ungroup()
  
  print("Performing grouping")
  abx_op_episodes <- grouping_function(data = abx_op_durations.prelim, id = key, index = course_index, 
                                       start = order_DT, end = end_DT, time_interval = 7) 
  
  print("Calculating courses after grouping")
  abx_op_episodes <- as.data.table(abx_op_episodes)
  abx_op_episodes <- abx_op_episodes[, `:=`(start_new = min(order_DT), end_new = max(end_DT)), by = .(key, group)]
  abx_op_episodes <- abx_op_episodes %>%
    dplyr::select(-course_index, -order_DT, -end_DT) %>%
    dplyr::distinct() %>%
    dplyr::mutate(abx_duration_days_new = end_new - start_new) 
  
  print("Join back to main abx dataset")
  abx_op_durations <- abx_op_durations.prelim %>%
    dplyr::left_join(abx_op_episodes, relationship = "many-to-many") %>%
    dplyr::select(key, end_DT = end_new, abx_duration_days = abx_duration_days_new) %>%
    dplyr::distinct() %>%
    dplyr::mutate(abx_duration_days = as.integer(abx_duration_days)) %>%
    tidyr::separate(key, c("patientID", "drug_class"), sep = "-")
  
  toc()
  
  return(abx_op_durations)
}

# Combine IP and OP abx courses
abx_courses <- function(abx_ip_durations, abx_op_durations)
{
  abx.courses <- dplyr::bind_rows(abx_ip_durations, abx_op_durations) %>%
    dplyr::arrange(patientID, drug_class, end_DT)
  
  return(abx.courses)
}

# ENCOUNTERS

# Process encounter data
enc_processing <- function(enc, sensitivity_analysis=FALSE)
{
  # Filter to encounters lasting >48 hours and then indexing
  
  tic()
  print("Initial cleaning of encounter data")
  if(sensitivity_analysis==TRUE){
    enc.prelim <- enc %>%
      dplyr::filter(visit_type == "Inpatient") %>%
      dplyr::select(PatientID, visit_type, hospitalAdmitDTS, hospitalDischargeDTS) %>%
      dplyr::distinct() %>%
      dplyr::arrange(PatientID, hospitalAdmitDTS) %>%
      dplyr::mutate(duration = as.numeric(as.difftime(hospitalDischargeDTS - hospitalAdmitDTS), units = "days")) %>%
      dplyr::filter(duration > 0) %>%
      dplyr::group_by(PatientID) %>%
      dplyr::mutate(enc_index = row_number()) %>%
      dplyr::ungroup()
  }
  else{
    enc.prelim <- enc %>%
      dplyr::filter(visit_type == "Inpatient") %>%
      dplyr::select(PatientID, visit_type, hospitalAdmitDTS, hospitalDischargeDTS) %>%
      dplyr::distinct() %>%
      dplyr::arrange(PatientID, hospitalAdmitDTS) %>%
      dplyr::mutate(duration = as.numeric(as.difftime(hospitalDischargeDTS - hospitalAdmitDTS), units = "days")) %>%
      dplyr::filter(duration > 2) %>%
      dplyr::group_by(PatientID) %>%
      dplyr::mutate(enc_index = row_number()) %>%
      dplyr::ungroup()
  }  
  toc()
  
  # Group admit/discharge episodes within 7 days of each other into the same hospitalization episode
  tic()
  print("Group encounters occuring within 7 days of eachother into a single encounter episode")
  enc.episodes <- grouping_function(data = enc.prelim, id = PatientID, index = enc_index, 
                                    start = hospitalAdmitDTS, end = hospitalDischargeDTS, time_interval = 7) 
  toc()
  
  tic()
  print("Re-calculating start / end dates for encounters with new grouping")
  enc.episodes <- as.data.table(enc.episodes)
  enc.episodes <- enc.episodes[, `:=`(start_new = min(hospitalAdmitDTS), end_new = max(hospitalDischargeDTS)), by = .(PatientID, group)]
  enc_clean <- enc.episodes %>%
    dplyr::select(PatientID, IP_enc_endDTS = end_new) %>%
    dplyr::distinct()
  toc()
  
  return(enc_clean)
}

# ADMISSION DISPOSITION

# Map admission dispositions to clean categories
admt.mapped <- function(admt)
{
  admt_mapped <- admt %>% 
    dplyr::select(PatientID, PatientEncounterID, AdmitSourceDSC) %>% 
    dplyr::mutate(PatientEncounterID = as.numeric(PatientEncounterID)) %>%
    dplyr::mutate(admit_source_clean = dplyr::case_when(
      AdmitSourceDSC %in% c("Ambulatory Surgery Center", "Born Outside this Hospital", "Discharge and Readmit", "NSM Emergency",
                            "Outside Health Care Facility", "Outside Hospital", "Psych, Substance Abuse, or Rehab Hospital",
                            "Transfer from another Health Care Facility") ~ "Transfer from outside facility",
      AdmitSourceDSC %in% c("Court/Law Enforcement", "Self Referral") ~ "Self referral",
      AdmitSourceDSC %in% c("Skilled Nursing Facility", "Hospice") ~ "Skilled Nursing facility",
      AdmitSourceDSC %in% c("Physician or Clinic Referral") ~ "Clinic referral",
      AdmitSourceDSC %in% c("", "Information Unavailable") ~ "Unknown",
      AdmitSourceDSC %in% c("Born Inside this Hospital") ~ "Born in hospital")) %>%
    dplyr::select(PatientID, PatientEncounterID, admit_source_clean) %>%
    dplyr::distinct()
  
  return(admt_mapped)
}

# Set target hierarchy for looping through pathogens 
path_cat_hierarchy <- function(micro, hierarchy)
{
  path_cat_table <- micro %>%
    tidyr::gather(pathogen_hierarchy, pathogen_category, org_group_1:org_group_3) %>%
    dplyr::filter(!is.na(pathogen_category)) %>%
    dplyr::select(pathogen_hierarchy, pathogen_category) %>%
    dplyr::distinct()
  
  path_cat_table <- path_cat_table %>%
    dplyr::filter(pathogen_hierarchy == hierarchy) 
  
  return(path_cat_table)
  
}

#### FUNCTION SET 3: BUILD ADT - MICRO GROUND TRUTH DATASET FOR ALL PATHOGENS ####

# Filter room data to identify eligible hospital stays (room stays within 30d of each other) and filter out ER, periop and residential stays

# Ennumerate rooms and for ones that have the same DTS_in, sort and choose the latest one only
first_eligible_rooms <- function(room_dat)
  {
  
  room_dat$dept_room <- paste(room_dat$DepartmentDSC, room_dat$RoomID, sep = " ")
  
  room_dat.eligible.prelim <- room_dat %>%
    # head(10000) %>%
    dplyr::filter(!grepl("emerg|periop|unspecified|home ", DepartmentDSC, ignore.case = T, perl = T)) %>%
    dplyr::arrange(PatientID, DTS_in, DTS_out) %>%
    dplyr::group_by(PatientID, DTS_in) %>%
    dplyr::mutate(room_stay_index_raw = row_number()) %>%
    dplyr::filter(room_stay_index_raw == max(room_stay_index_raw)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-room_stay_index_raw) %>%
    dplyr::distinct() %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(room_stay_index = row_number()) %>%
    dplyr::ungroup()
  
  # Group room stays within 7 days of each other into episodes
  room_dat.episode <- grouping_function(room_dat.eligible.prelim, id = PatientID, index = room_stay_index, 
                                        start = DTS_in, end = DTS_out, time_interval = 7)

  room_dat.eligible.prelim2 <- room_dat.eligible.prelim %>%
    dplyr::left_join(room_dat.episode) %>%
    dplyr::arrange(PatientID, DTS_in) %>%
    dplyr::filter(duration >= 2) %>% # Trim to rooms where patient stayed between 3 & 30 days
    dplyr::mutate(drop = ifelse(duration > 30, "X", NA)) %>% # Drop all room stays that occur after any room stay >30 days in duration
    dplyr::group_by(PatientID, group) %>%
    tidyr::fill(drop, .direction = "downup") %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.na(drop)) %>%
    dplyr::group_by(PatientID, group) %>% # Re-index room stays after having taken out redundant rooms and rooms with durations outside of 48 to 720 hrs
    dplyr::mutate(room_stay_index = row_number(),
           n_rooms_hosp = max(room_stay_index),
           first_room = ifelse(room_stay_index == 1, "X", NA),
           hospitalization_id = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(PatientID, group, room_stay_index) %>% 
    dplyr::mutate(room_stay_id = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::select(PatientID:PatientClassDSC, DepartmentDSC, dept_room, first_room, room_stay_index, n_rooms_hosp, 
           room_stay_id, hospitalization_index = group, hospitalization_id, DTS_in, DTS_out, duration) %>%
    dplyr::distinct()
  
  rm(room_dat.episode)
  
  return(room_dat.eligible.prelim2)
  }

# Filter micro data for patients in ADT dataset
micro_prep <- function(micro.dedup, room_dat.eligible.prelim2)
  {
  micro.eligible <- micro.dedup %>% 
    dplyr::filter(PatientID %in% room_dat.eligible.prelim2$PatientID) %>% 
    dplyr::select(PatientID, accession, coll_datetime_UTC, org_group_1:org_group_3) %>%
    dplyr::distinct()
  
  return(micro.eligible)
  }

# Ensure room data begins after the start of the micro data
room_dat_eligible_trimmed <- function(room_dat.eligible.prelim2, micro.eligible)
  {
  room_dat.eligible.trimmed <- room_dat.eligible.prelim2 %>%
    dplyr::filter(DTS_in > min(micro.eligible$coll_datetime_UTC)) 
  
  return(room_dat.eligible.trimmed)
  }

# Add micro data to ADT / room data
adt_micro_join <- function(room_dat.eligible.trimmed, micro.eligible)
  {
  adt.micro.raw <- room_dat.eligible.trimmed %>%
    dplyr::left_join(micro.eligible, by = c("PatientID"), relationship = "many-to-many") %>%
    dplyr::distinct() %>%
    dplyr::arrange(hospitalization_id, coll_datetime_UTC)
  
  return(adt.micro.raw)
  }

# Add time to infection relative to first room in a hospitalization 
time_to_infxn <- function(adt.micro.raw)
  {
  adt.micro.raw <- adt.micro.raw %>%
    dplyr::mutate(time_to_infxn = difftime(coll_datetime_UTC, DTS_in, units = "days"),
           time_to_infxn = round(time_to_infxn, digits = 1))
  
  return(adt.micro.raw)
  }

# Drop micro data for rooms where sample was drawn >1 year prior to and >30 days after entry (skin flora & Enterococcus)
# or >6 months prior to and >30 days after entry (other)
drop_old_micro <- function(adt.micro.raw)
  {
  adt.micro.raw <- adt.micro.raw %>%
    dplyr::mutate(accession = dplyr::case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        (time_to_infxn < -365 | time_to_infxn > 30) ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        (time_to_infxn < -182.5 | time_to_infxn > 30) ~ NA,
      .default = accession)) %>%
    dplyr::mutate(coll_datetime_UTC = dplyr::case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = coll_datetime_UTC)) %>%
    dplyr::mutate(org_group_1 = dplyr::case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = org_group_1)) %>%
    dplyr::mutate(org_group_2 = dplyr::case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = org_group_2)) %>%
    dplyr::mutate(org_group_3 = dplyr::case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = org_group_3)) %>%
    dplyr::mutate(time_to_infxn = dplyr::case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = time_to_infxn)) %>%
    dplyr::distinct() %>%
  dplyr::group_by(room_stay_id) %>%
  dplyr::mutate(group_size = n(),
         drop = ifelse(group_size > 1 & is.na(accession), "X", NA)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.na(drop)) %>%
  dplyr::select(-group_size, -drop)
  
  return(adt.micro.raw)
  }

# Drop hospitalizations where the patient received antibiotics within 7 days of their stay in their 1st room 
remove_treated_patients <- function(adt.micro.raw, abx.courses, day_threshold)
  {
  abx.courses.filtered <- abx.courses %>%
    dplyr::rename(PatientID = patientID) %>%
    dplyr::filter(abx_duration_days >=2)
  
  adt.micro.raw <- adt.micro.raw %>%
    # head(10000) %>%
    dplyr::left_join(abx.courses.filtered, relationship = "many-to-many") %>%
    dplyr::mutate(diff_days = as.Date(floor_date(end_DT, unit = "days")) - as.Date(floor_date(DTS_in, unit = "days"))) %>%
    dplyr::mutate(drop = ifelse(diff_days >= -7 & diff_days <= day_threshold & !is.na(first_room), "X", NA)) %>%
    dplyr::group_by(hospitalization_id) %>%
    tidyr::fill(drop, .direction = "downup") %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.na(drop)) %>%
    dplyr::select(PatientID:time_to_infxn) %>%
    dplyr::distinct()
  
  return(adt.micro.raw)
}

# Drop hospitalizations where the patient had a previous IP encounter >48 hours in length
drop_prev_enc <- function(adt.micro.raw, enc_clean)
  {
  adt.micro.raw <- adt.micro.raw %>%
    dplyr::left_join(enc_clean, relationship = "many-to-many") %>%
    dplyr::mutate(diff_days  = as.Date(floor_date(IP_enc_endDTS, unit = "days")) - as.Date(floor_date(DTS_in, unit = "days"))) %>%
    dplyr::mutate(drop = ifelse(diff_days < 0 & diff_days > -90, "X", NA)) %>%
    dplyr::group_by(hospitalization_id) %>%
    tidyr::fill(drop, .direction = "downup") %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.na(drop)) %>%
    dplyr::select(PatientID:time_to_infxn) %>%
    dplyr::distinct()
  
  return(adt.micro.raw)
}

# Add feature for matching cases (time to infxn) to control (LOS in room)
# matching_duration <- function(adt.micro.raw)
#   {
#   adt.micro.raw <- adt.micro.raw %>%
#     dplyr::mutate(time_to_infxn = as.numeric(time_to_infxn)) %>%
#     dplyr::mutate(matching_duration = dplyr::case_when(
#        time_to_infxn > 0 ~ time_to_infxn,
#        .default = duration))
#   
#   return(adt.micro.raw)
#   
# }

#### FUNCTION SET 4: PATHOGEN-SPECIFIC CASE / CONTROL SELECTION ####

# Select proper pathogen hierarchy in micro data
adt_micro_hierarchy <- function(adt.micro.final, pathogen_hierarchy)
{
  adt.micro.raw.temp <- adt.micro.raw %>%
    dplyr::select(PatientID:coll_datetime_UTC, time_to_infxn)
  
  adt.micro.raw.hierarchy <- adt.micro.raw %>%
    dplyr::select(org_group_1:org_group_3) %>%
    dplyr::select(matches(pathogen_hierarchy)) %>%
    dplyr::rename(org_group = matches(pathogen_hierarchy))
  
  adt.micro.raw.filtered <- dplyr::bind_cols(adt.micro.raw.temp, adt.micro.raw.hierarchy) %>%
    dplyr::select(PatientID:coll_datetime_UTC, org_group, time_to_infxn)
  
  rm(adt.micro.raw.temp, adt.micro.raw.hierarchy)
  
  return(adt.micro.raw.filtered)
}

# Flag hospitalizations where the pathogen of interest was identified after day X (set by user) and before day 30 in the first eligible room 
potential_cases <- function(adt.micro.raw.filtered, pathogen_category, day_threshold) 
  {
  adt.micro.raw.filtered <- adt.micro.raw.filtered %>% 
    dplyr::mutate(potential_case = ifelse(!is.na(first_room) &
                                     org_group == pathogen_category & 
                                     time_to_infxn >= day_threshold & time_to_infxn <= 30, "X", NA))
  
  return(adt.micro.raw.filtered)
  }

# Drop potential cases where micro timestamp falls during occupation of a later room during that hospitalization episode
after_first_room <- function(adt.micro.raw.filtered, pathogen_category, day_threshold) 
{
  secondary.case.diff.room <- adt.micro.raw.filtered %>%
    dplyr::filter(org_group == pathogen_category) %>%
    dplyr::select(PatientID, hospitalization_id, room_stay_id, org_group, accession, time_to_infxn) %>%
    dplyr::distinct() %>%
    dplyr::filter(time_to_infxn >= day_threshold & time_to_infxn <= 30) %>%
    dplyr::group_by(hospitalization_id, accession) %>%
    dplyr::mutate(num_rooms = n()) %>%
    dplyr::mutate(drop = ifelse(num_rooms > 1, "X", NA)) %>%
    tidyr::fill(drop, .direction = "downup") %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(drop))
  
  adt.micro.raw.filtered <- adt.micro.raw.filtered %>%
    dplyr::filter(!hospitalization_id %in% secondary.case.diff.room$hospitalization_id) 
  
  return(adt.micro.raw.filtered)
}

# Drop potential cases where the patient had the pathogen of interest between day -365 (S. aureus, Enterococcus) / day -182.5 (all others)
# and day +X (set by user), relative to entry into first room (suggesting not nosocomially acquired)
flag_prior_infxn <- function(adt.micro.raw.filtered, pathogen_category, day_threshold) 
  {
  prior_infxn <- adt.micro.raw.filtered %>% 
    dplyr::arrange(PatientID, coll_datetime_UTC) %>%
    dplyr::filter(org_group == pathogen_category) %>%
    dplyr::mutate(prior_infxn = 
             dplyr::case_when(
               !is.na(first_room) &
                 org_group %in% c("Skin flora", "Staph_aureus", "MSSA", "MRSA", "E_faecalis", "E_faecium", 
                                  "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") & 
                 time_to_infxn >= -365 & time_to_infxn <= day_threshold ~ "X",
               !is.na(first_room) &
                 !org_group %in% c("Skin flora", "Staph_aureus", "MSSA", "MRSA", "E_faecalis", "E_faecium", 
                                   "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") & 
                 time_to_infxn >= -182.5 & time_to_infxn <= day_threshold ~ "X",
               .default = NA)) %>%
    dplyr::group_by(hospitalization_id) %>%
    tidyr::fill(prior_infxn, .direction = "downup") %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(prior_infxn))
  
  adt.micro.cases <- adt.micro.raw.filtered %>%
    dplyr::filter(!hospitalization_id %in% prior_infxn$hospitalization_id) 
  
  return(adt.micro.cases)
}

# Identify controls (for matching later)
generate_controls <- function(adt.micro.cases, pathogen_category) 
{
  adt.micro.cc.temp <- adt.micro.cases %>%
    dplyr::mutate(group = ifelse(org_group == pathogen_category & !is.na(potential_case), "case", NA)) %>%
    dplyr::mutate(group = ifelse(org_group != pathogen_category & !is.na(first_room), "control", group)) %>%
    dplyr::filter(!is.na(group))
  
  cc.overlap <- adt.micro.cc.temp %>%
    dplyr::select(hospitalization_id, group) %>%
    dplyr::distinct() %>%
    dplyr::group_by(hospitalization_id) %>%
    dplyr::mutate(overlap = n(), 
           overlap_flag = ifelse(overlap > 1, "X", NA)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(overlap_flag))
  
  adt.micro.cc <- adt.micro.cc.temp %>%
    dplyr::mutate(group = ifelse(hospitalization_id %in% cc.overlap$hospitalization_id & group == "control", NA, group)) %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::select(hospitalization_id, room_stay_id, group) %>% 
    dplyr::distinct()
  
  return(adt.micro.cc)
}

#### FUNCTION SET 5: GENERATE FEATURES FOR MATCHING AND MODELING ####

# Demographics
add_demographics <- function(adt.micro.raw, dems) 
  {
  dem_features <- adt.micro.raw %>%
    dplyr::select(PatientID, room_stay_id, hospitalization_id, DTS_in) %>%
    dplyr::distinct() %>%
    dplyr::left_join(dems) %>%
    dplyr::mutate(age = round(time_length(difftime(DTS_in, DOB), "years"), digits = 1)) %>%
    dplyr::mutate(sex = ifelse(sex %in% c("Male", "Female"), sex, "Other")) %>%
    dplyr::mutate(across(where(is.character), ~na_if(., ""))) %>%
    dplyr::select(-DTS_in, -DOB)
  
  # Impute missing values for admission source and race/ethnicity
  dem_features <- dem_features %>% 
    dplyr::mutate(race = ifelse(is.na(race), "Declined", race)) %>%
    dplyr::mutate(ethnicity = ifelse(is.na(ethnicity), "Prefer not to say/Decline", ethnicity)) 

  return(dem_features)
}

# Antibiotics
add_abx <- function(abx.courses, adt.micro.raw) 
  {
  # Code to create the feature set for inpatient abx history
  abx_features.prelim <- adt.micro.raw %>% 
    dplyr::select(PatientID, room_stay_id, DTS_in, coll_datetime_UTC) %>%
    dplyr::distinct() %>%
    # head(1000) %>%
    dplyr::left_join(abx.courses %>% dplyr::select(PatientID = patientID, drug_class, end_DT, abx_duration_days), 
              relationship = "many-to-many") %>% 
    dplyr::mutate(diff_days = as.Date(floor_date(DTS_in, unit = "days")) - end_DT) %>%
    dplyr::filter(diff_days > 0 & diff_days < 182.5) %>% 
    dplyr::mutate(interval = dplyr::case_when(diff_days <= 60 ~ "0_60",
                                diff_days > 60 ~ "60_plus")) %>%
    dplyr::mutate(duration_binary = dplyr::case_when(abx_duration_days <= 2 ~ "under2days",
                                       abx_duration_days > 2 ~ "over2days")) %>% 
    dplyr::filter(duration_binary == "over2days") %>%
    dplyr::mutate(val = 1) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(PatientID) %>% 
    tidyr::pivot_wider(names_from = c(drug_class, duration_binary, interval), names_sep = "_", values_from = val) %>% 
    dplyr::ungroup() %>%
    dplyr::select(-c(DTS_in, end_DT, abx_duration_days, diff_days)) %>%
    dplyr::distinct() %>%
    dplyr::mutate_if(is.numeric, list(~ tidyr::replace_na(., 0)))
  
  return(abx_features.prelim)
  }

final_abx_clean <- function(abx_features.prelim, id_vars) 
  {
  id_vars = c("PatientID", "room_stay_id", "coll_datetime_UTC")
  
  abx_features.melt <- reshape2::melt(abx_features.prelim, id.vars = id_vars)
  
  abx_features.melt <- abx_features.melt %>%
    dplyr::mutate(variable = gsub("_over2days", "", variable))
  
  abx_features <- reshape2::dcast(abx_features.melt, PatientID + room_stay_id ~ variable, value.var = "value", fun.aggregate = sum)
  
  abx_features <- abx_features %>%
    dplyr::mutate_if(is.numeric, list(~ tidyr::replace_na(., 0)))
  
  abx_of_interest <- c("penicillin_0_60", "extended_spectrum_penicillin_0_60", "cephalosporin_0_60", 
                       "extended_spectrum_cephalosporin_0_60", "carbapenem_0_60", "anti_staph_beta_lactam_0_60", 
                       "fluoroquinolone_0_60", "glycopeptide_0_60", "anti_anaerobe_0_60", "anti_Cdiff_0_60",
                       "tetracycline_0_60", "macrolide_0_60", "sulfonamide_0_60", "lincosamide_0_60")
  
  variables_to_group <- names(abx_features)[grepl("_0_60$", names(abx_features)) & !(names(abx_features) %in% abx_of_interest)]
  
  # Add flag for antibiotics not in existing categories
  abx_features <- abx_features %>%
    dplyr::mutate(other_abx_0_60 = as.integer(rowSums(dplyr::select(., all_of(variables_to_group)), na.rm = TRUE) >= 1))
  
  abx_features <- abx_features %>% 
    dplyr::mutate(any_abx_0_60 = ifelse(rowSums(dplyr::select(., contains("_0_60"))) > 0, 1, 0)) %>% 
    dplyr::mutate(any_abx_60_plus = ifelse(rowSums(dplyr::select(., contains("_60_plus"))) > 0, 1, 0)) 
  
  return(abx_features)
  }

# Comorbidities
add_comorbidities <- function(adt.micro.raw, elix)
{
  elix_features <- adt.micro.raw %>%
    dplyr::select(PatientID, PatientEncounterID, hospitalization_id) %>%
    dplyr::distinct() %>%
    dplyr::left_join(elix)
  
  return(elix_features)
}

# CPT codes
add_cpt <- function(adt.micro.raw, cpt) 
  {
  cpt_features <- adt.micro.raw %>%
    dplyr::select(PatientID, hospitalization_id, DTS_in) %>%
    dplyr::distinct() %>%
    dplyr::left_join(cpt, relationship = "many-to-many") %>% 
    dplyr::filter(proc_date < DTS_in) %>% 
    dplyr::mutate(diff_days = difftime(proc_date, DTS_in, units="days")) %>% 
    dplyr::mutate(surgery_past_90 = ifelse((diff_days < 0) & (diff_days >= -90), proc_simple, NA)) %>% 
    dplyr::filter(!is.na(surgery_past_90)) %>% 
    dplyr::select(PatientID, hospitalization_id, surgery_past_90) %>% 
    dplyr::distinct() %>%
    dplyr::mutate(val = 1) %>% 
    tidyr::pivot_wider(names_from = surgery_past_90, values_from = val) %>%
    dplyr::mutate_if(is.numeric, list(~ replace_na(., 0))) %>%
    dplyr::arrange(hospitalization_id, PatientID) %>%
    dplyr::mutate(any_surgery = 1)
  
  cpt_features <- cpt_features %>%
    dplyr::rename(ortho_sx = `Orthopedic surgery`, abd_sx = `Intra_abdominal surgery`, rp_sx = `Retroperitoneal surgery`,
           pelvic_sx = `Pelvic surgery`, vasc_sx = `Vascular surgery`, card_sx = `Cardiac surgery`,
           thoracic_sx = `Thoracic surgery`, neuro_sx = `Neurosurgery`, ent_sx = `ENT surgery`, 
           spinal_sx = `Spinal surgery`, thyroid_sx = `Thyroid surgery`, transplant_sx = `Transplant surgery`,
           breast_sx = `Breast surgery`, hernia_sx = `Hernia surgery`)
    
  return(cpt_features)
  }

# Add admission source
add_admit.source <- function(adt.micro.raw, admt_mapped) 
{
  admt_features <- adt.micro.raw %>%
    dplyr::select(PatientID, PatientEncounterID, hospitalization_id) %>%
    dplyr::distinct() %>%
    dplyr::left_join(admt_mapped, relationship = "many-to-many") 
  
  # Impute missing values for admission source and race/ethnicity
  admt_features <- admt_features %>% 
    dplyr::mutate(admit_source_clean = dplyr::case_when(
      is.na(admit_source_clean) ~ "Information Unavailable",
      admit_source_clean %in% c("Self referral", "Clinic referral", "Skilled Nursing facility", 
                                 "Transfer from outside facility", "Information Unavailable") ~ admit_source_clean,
      .default = "Other or unknown"))
  
  # Re-level
  admt_features <- admt_features %>%
    dplyr::mutate(admit_source_clean = factor(admit_source_clean)) %>%
    dplyr::mutate(admit_source_clean = fct_relevel(admit_source_clean, "Self referral", after = 0)) %>%
    dplyr::mutate(admit_source_clean = fct_relevel(admit_source_clean, "Information Unavailable", after = 5))
  
  return(admt_features)
}

# Prior room occupant
add_prior_pathogens <- function(room_dat, micro, pathogen_hierarchy, adt.micro.raw) 
  {
  # Prep total room dataset
  
  room_dat$dept_room <- paste(room_dat$DepartmentDSC, room_dat$RoomID, sep = " ")
  room_dat$PatientEncounterID <- as.character(room_dat$PatientEncounterID)
  
  room_dat.simple <- room_dat %>%
    dplyr::select(PatientID, PatientEncounterID, dept_room, DTS_in, DTS_out, duration) %>%
    dplyr::distinct() 
  
  # Prep total micro dataset
  micro.simple <- micro %>%
    dplyr::select(PatientID, PatientEncounterID = CSN, coll_datetime_UTC) 
  
  micro.hierarchy <- micro %>%
    dplyr::select(org_group_1:org_group_3) %>%
    dplyr::select(matches(pathogen_hierarchy)) %>%
    dplyr::rename(org_group = matches(pathogen_hierarchy))
  
  micro.filtered <- dplyr::bind_cols(micro.simple, micro.hierarchy)
  
  # Join total room and total micro datasets
  room_micro <- room_dat.simple %>%
    dplyr::left_join(micro.filtered, relationship = "many-to-many") %>%
    dplyr::filter(coll_datetime_UTC > DTS_in & coll_datetime_UTC < DTS_out) %>%
    dplyr::mutate(prior_colldate = floor_date(coll_datetime_UTC, unit = "days")) %>%
    dplyr::select(dept_room, prior_PatientID = PatientID, prior_colldate, prior_org_group = org_group) %>%
    dplyr::distinct()
    
  # Join total room+micro data to clean ADT+micro dataset and filter for when a different patient had a pathogen in the previous 365 days
  prior_pathogen.prelim <- adt.micro.raw %>%
    # head(1000) %>%
    dplyr::filter(!is.na(first_room)) %>%
    dplyr::mutate(DT_in = floor_date(DTS_in, unit = "days")) %>%
    dplyr::select(PatientID, dept_room, hospitalization_id, DT_in) %>%
    dplyr::distinct() %>%
    dplyr::left_join(room_micro, relationship = "many-to-many") %>%
    dplyr::mutate(diff_days = difftime(DT_in, prior_colldate, units="days")) %>%
    dplyr::filter(diff_days > 0 & diff_days < 365) %>%
    dplyr::filter(PatientID != prior_PatientID)
  
  # Select closest date to DTS_in if pathogen present more than once in the past
  prior_pathogen.prelim.dedup <- as.data.table(prior_pathogen.prelim)
  prior_pathogen.prelim.dedup <- prior_pathogen.prelim.dedup[order(dept_room, prior_org_group, diff_days)]
  prior_pathogen.prelim.dedup[,
                         diff_max := max(diff_days),
                         by = .(dept_room, prior_org_group)]
  prior_pathogen.prelim.dedup <- prior_pathogen.prelim.dedup %>%
    dplyr::filter(diff_days == diff_max) %>%
    dplyr::mutate(diff_days = as.integer(diff_days))
  
  # Negative exponential transform of days since last pathogen
  diff_days_transform <- prior_pathogen.prelim.dedup %>% 
    dplyr::select(dept_room, hospitalization_id, diff_days) %>%
    dplyr::mutate(diff_days_transform = exp(-0.01 * (diff_days))) %>%
    dplyr::mutate(diff_days_transform = round(diff_days_transform, 2)) %>%
    dplyr::mutate(diff_days_transform = diff_days_transform * 5) # Additional up-weighting of prior occupant above that of other ward occupants

  # Spread transformed data
  prior_pathogen.prelim.transform <- prior_pathogen.prelim.dedup %>%
    dplyr::left_join(diff_days_transform, relationship = "many-to-many") %>%
    dplyr::select(PatientID:hospitalization_id, prior_org_group, diff_days_transform) %>%
    dplyr::distinct() %>%
    dplyr::mutate(prior_org_group = paste0("prior_", prior_org_group)) %>% 
    spread(prior_org_group, diff_days_transform) %>%
    dplyr::arrange(hospitalization_id) %>%
    dplyr::mutate_if(is.numeric, list(~ replace_na(., 0)))
  
   # Add back to main ADT / micro dataset 
  prior_occupant_features <- adt.micro.raw %>%
    dplyr::select(PatientID, dept_room, hospitalization_id) %>%
    dplyr::distinct() %>%
    dplyr::left_join(prior_pathogen.prelim.transform) %>%
    dplyr::mutate_if(is.numeric, list(~ replace_na(., 0)))

  return(prior_occupant_features)
  }

# Colonization pressure
calculate_CP_score <- function(date, location_specific, room_dat, numerator_category, micro_df, lambda = 0.01) 
  {
  # Define the date 30 days prior to the patient of interest 
  ward_previous_start_date <- date - days(30)
  
  prev_month_occupants <- room_dat %>% 
    dplyr::filter(DepartmentDSC == location_specific) %>% 
    dplyr::filter((DTS_in >= ward_previous_start_date & DTS_in < date) | (DTS_out > ward_previous_start_date & DTS_in < ward_previous_start_date))
  
  # Filter the microbiology data by date range and hospital
  filtered_data <- micro_df %>%
    dplyr::filter(PatientID %in% prev_month_occupants$PatientID) %>% 
    dplyr::filter(coll_datetime_UTC >= date - days(90), coll_datetime_UTC <= date)
  
  # Filter the data based on the numerator and denominator categories
  numerator_categories <- strsplit(numerator_category, " \\+ ")[[1]]
  # denominator_categories <- strsplit(denominator_category, " \\+ ")[[1]]
  
  filtered_data <- filtered_data %>%
    dplyr::filter(org_group_3 %in% numerator_categories) %>% 
    dplyr::distinct(PatientID, loc_cat, loc_specific, coll_datetime_UTC, org_group_3)
  
  tic()
  
  # Calculate the colonization pressure
  if(nrow(filtered_data) == 0) {
    colonization_pressure <- 0
  } else {
    colonization_pressure <- filtered_data %>%
      dplyr::group_by(PatientID, org_group_3) %>%
      dplyr::filter(coll_datetime_UTC == max(coll_datetime_UTC)) %>% # Retain only the most recent positive sample
      dplyr::summarize(days_since_last_positive = as.numeric(date - coll_datetime_UTC), .groups = "drop") %>%
      dplyr::mutate(score = exp(-lambda * days_since_last_positive)) %>% # Exponential decay
      dplyr::pull(score) %>% 
      sum(na.rm = TRUE)
  }
  
  return(colonization_pressure)
  
  toc()
  }

calculate_and_add_cp_scores_parallel <- function(adt.micro.raw, micro, room_dat, location_map) 
  {
  loc_map_single <- location_map %>% 
    dplyr::distinct(DepartmentDSC, .keep_all = TRUE)
  
  cp_features <- adt.micro.raw %>%
    dplyr::filter(!is.na(first_room)) %>%
    dplyr::select(PatientID, hospitalization_id, DepartmentDSC, DTS_in) %>%
    dplyr::distinct() %>%
    dplyr::left_join(loc_map_single, by = "DepartmentDSC") %>%
    dplyr::mutate(filter_date = DTS_in) %>%
    dplyr::mutate(micro_loc = if_else(is.na(micro_loc), DepartmentDSC, micro_loc))
  
  cp_score_function <- function(date, location_specific, micro_df, room_dat, numerator_category, lambda = 0.01) 
    {
    calculate_CP_score(date, location_specific, room_dat, numerator_category, micro_df, lambda)
    }
  
  # Categories, numerator categories, and corresponding column names
  categories_info <- list(
    C_diff = list(numerator_category = "C_diff", column_name = "CDiff_cp"),
    MSSA = list(numerator_category = "MSSA", column_name = "MSSA_cp"),
    MRSA = list(numerator_category = "MRSA", column_name = "MRSA_cp"),
    DS_Enterobacterales = list(numerator_category = "DS_E_coli + DS_P_mirabilis + DS_K_aerogenes + DS_C_koseri + DS_K_pneumoniae + 
                               DS_K_oxytoca + DS_E_cloacae + DS_C_freundii + DS_S_marcescens + DS_Providencia + DS_M_morganii", 
                               column_name = "DS_Entero_cp"),
    ESBL = list(numerator_category = "ESBL_E_coli + ESBL_K_pneumoniae + ESBL_K_oxytoca + ESBL_K_aerogenes + ESBL_E_cloacae + 
                                     ESBL_C_freundii + ESBL_C_koseri + ESBL_S_marcescens + ESBL_Providencia + ESBL_M_morganii + 
                                     ESBL_P_mirabilis", column_name = "ESBL_cp"),
    VSE = list(numerator_category = "VSE_faecalis + VSE_faecium", column_name = "VSE_cp"),
    VRE = list(numerator_category = "VRE_faecalis + VRE_faecium", column_name = "VRE_cp"),
    DS_P_aeruginosa = list(numerator_category = "DS_P_aeruginosa", column_name = "DS_PsA_cp"),
    DR_P_aeruginosa = list(numerator_category = "DR_P_aeruginosa", column_name = "DR_PsA_cp")
  )
  
  for (category_info in categories_info) 
    {
    
    # Create a list of parameters for each mapply call
    parameters_list <- lapply(1:nrow(cp_features), function(i) 
      {
      list(
        date = cp_features$filter_date[i],
        location_specific = cp_features$DepartmentDSC[i],
        micro_df = micro,
        room_dat = room_dat,
        numerator_category = category_info$numerator_category
        )
      })
    
    # Apply mclapply for parallel processing
    cp_scores <- mclapply(parameters_list, function(x) do.call(cp_score_function, x), mc.cores = 16)
    
    # Assign scores back to cp_features
    cp_features[[category_info$column_name]] <- unlist(cp_scores)
    }
  
  return(cp_features)
  }

# Select pathogens with n>300 cases in unmatched cohort

pathogen_selection <- function(unmatched.data, path_categories)
  {
  path_cat_table_matching <- unmatched.data %>%
    tidyr::gather(run, group, C_diff:VSE_faecium) %>%
    dplyr::count(run, group) %>%
    dplyr::arrange(-n) %>%
    dplyr::filter(group == "case") %>%
    dplyr::filter(n >= 300) %>%
    dplyr::select(pathogen_category = run) %>%
    dplyr::mutate(pathogen_category = factor(pathogen_category)) %>%
    dplyr::left_join(path_categories)
  
  return(path_cat_table_matching)
}

#### METAFUNCTION 1: PREP ADT AND MICRO DATA (NOT PATHOGEN-SPECIFIC) ####

adt_micro_initial_prep <- function(room_dat, 
                                   micro.dedup, 
                                   day_threshold,
                                   abx.courses,
                                   enc_clean)
{
  # Group hospital room stays into episodes 
  print("Group hospital room stays into episodes")
  tic()
  room_dat.eligible <<- first_eligible_rooms(room_dat)
  toc()
  
  # Filter micro data for patients in ADT dataset
  print("Filter micro data for patients in ADT dataset")
  tic()
  micro.eligible <<- micro_prep(micro.dedup, room_dat.eligible)
  toc()
  
  # Ensure room data begins after the start of the micro data
  print("Ensure room data begins after the start of the micro data")
  tic()
  room_dat.eligible.trimmed <<- room_dat_eligible_trimmed(room_dat.eligible, micro.eligible)
  toc()
  
  # Add prepped micro data to prepped ADT data
  print("Add prepped micro data to prepped ADT data")
  tic()
  adt.micro.raw <- adt_micro_join(room_dat.eligible.trimmed, micro.eligible)
  toc()
  
  # Add time interval between entry into first room in a hospitalization and a positive micro specimen
  print("Add time interval between entry into first room in a hospitalization and a positive micro specimen")
  tic()
  adt.micro.raw <- time_to_infxn(adt.micro.raw)
  toc()
  
  # Drop micro data for rooms where sample was drawn >1 year prior to and >30 days after entry (skin flora & Enterococcus)
  # or >6 months prior to and >30 days after entry (other)
  print("Drop old micro data")
  tic()
  adt.micro.raw <- drop_old_micro(adt.micro.raw)
  print(nrow(adt.micro.raw))
  toc()
  
  # Print Total Hospitalization
  cohort_size <- adt.micro.raw %>% group_by(org_group_3) %>% summarise(unique_encounters = n_distinct(hospitalization_id))
  print(cohort_size)
  
  # Remove hospitalization episodes where the patient received inpatient antibiotics prior to their stay in the room based on when their abx course finished
  print("Remove potential cases where the patient received concomitant IP/OP antibiotics")
  tic()
  adt.micro.raw <- remove_treated_patients(adt.micro.raw, abx.courses, day_threshold)
  print(nrow(adt.micro.raw))
  toc()
  
  # Print Total Hospitalization
  cohort_size <- adt.micro.raw %>% group_by(org_group_3) %>% summarise(unique_encounters = n_distinct(hospitalization_id))
  print(cohort_size)
  
  # Remove hospitalization episodes where patients had another IP encounter >2d in duration in the previous 90 days
  print("Drop potential cases where the patient had another hospital encounter in the previous 90 days")
  tic()
  adt.micro.raw <- drop_prev_enc(adt.micro.raw, enc_clean)
  print(nrow(adt.micro.raw))
  toc()
  
  # Print Total Hospitalization
  cohort_size <- adt.micro.raw %>% group_by(org_group_3) %>% summarise(unique_encounters = n_distinct(hospitalization_id))
  print(cohort_size)
  
  # Add feature for matching cases (time to infxn) to control (LOS in room)
  # print("Add feature for matching cases to control on time to infxn (cases) or LOS (controls)")
  # tic()
  # adt.micro.raw <- matching_duration(adt.micro.raw)
  # toc()
  
  return(adt.micro.raw)
  
}

#### METAFUNCTION 2: GENERATE UNMATCHED CASES & CONTROLS (PATHOGEN-SPECIFIC) ####

# Generate cases and controls (pathogen specific)
generate_unmatched_data <- function(pathogen_hierarchy,
                                    pathogen_category, 
                                    day_threshold,
                                    adt.micro.raw, 
                                    abx_prelim) 
{
  
  # Select proper pathogen hierarchy in micro data
  print("Filter ADT micro data for proper micro hierarchy")
  tic()
  adt.micro.raw.filtered <- adt_micro_hierarchy(adt.micro.raw, pathogen_hierarchy)
  toc()
  
  # Flag hospitalizations where the pathogen of interest was identified after day X (set by user) and before day 30 in the first eligible room 
  print(paste0("Identify potential cases of ", pathogen_category))
  tic()
  adt.micro.raw.filtered <- potential_cases(adt.micro.raw.filtered, pathogen_category, day_threshold) 
  initial_case = nrow(adt.micro.raw.filtered %>% filter(potential_case == 'X'))
  toc()
  
  # Drop subsequent potential cases after first case identified in an eligible room (including ones with same timestamp)
  print(paste0("Remove potential cases identified in room stays after the first case in an episode"))
  tic()
  adt.micro.raw.filtered <- after_first_room(adt.micro.raw.filtered, pathogen_category, day_threshold) 
  drop_subsequent = nrow(adt.micro.raw.filtered %>% filter(potential_case == 'X'))
  toc()
  
  # Exclude hospitalizations where the patient had the pathogen of interest between day -365 (S. aureus, Enterococcus) / day -182.5 (all others)
  # and day +X (set by user), relative to entry into first room (suggesting not nosocomially acquired)
  print(paste0("Remove potential cases that have evidence of prior infection with ", pathogen_category))
  tic()
  adt.micro.cases <- flag_prior_infxn(adt.micro.raw.filtered, pathogen_category, day_threshold) 
  drop_prior_infxn = nrow(adt.micro.cases %>% filter(potential_case == 'X'))
  toc()
  
  # Select controls (for matching later)
  print("Select preliminary controls")
  tic()
  adt.micro.cc <- generate_controls(adt.micro.cases, pathogen_category) 
  toc()
  
  row <- data.frame(
    pathogen = pathogen_category,
    initial_cases = initial_case,
    case_after_drop_subsequent = drop_subsequent,
    case_after_drop_prior_infection = drop_prior_infxn
  )
  
  flow_chart <<- rbind(flow_chart, row)
  
  return(adt.micro.cc)
}

#### METAFUNCTION 3: CREATION OF FINAL DATASET ####

generate.unmatched.dataset <- function(pathogen_category,
                                       adt.micro.raw,
                                       cc.unmatched,
                                       dem_features,
                                       abx_features,
                                       elix_features,
                                       cpt_features,
                                       admt_features,
                                       prior_occupant_features,
                                       cp_features)
  
  {
  # Prep ADT / micro dataset
  adt.micro <- adt.micro.raw %>%
    dplyr::mutate(DTS_in_month = floor_date(DTS_in, unit = "month"),
           DTS_in_year = floor_date(DTS_in, unit = "year")) %>%
    dplyr::group_by(DepartmentDSC) %>%
    dplyr::mutate(dept_code = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(org_group_3 = ifelse(org_group_3 == pathogen_category, org_group_3, NA)) %>%
    dplyr::mutate(time_to_infxn = ifelse(org_group_3 == pathogen_category, time_to_infxn, NA)) %>%
    dplyr::group_by(room_stay_id) %>%
    tidyr::fill(time_to_infxn, .direction = "downup") %>%
    tidyr::fill(org_group_3, .direction = "downup") %>%
    dplyr::ungroup() %>%
    dplyr::select(PatientID, PatientEncounterID, dept_code, hospitalization_id, room_stay_id, 
           DTS_in, DTS_in_month, DTS_in_year, duration, time_to_infxn, org_group_3) %>%
    dplyr::distinct() %>%
    dplyr::mutate(matching_duration = dplyr::case_when(
      is.na(org_group_3) ~ duration,
      !is.na(org_group_3) & duration > time_to_infxn ~ time_to_infxn,
      !is.na(org_group_3) & duration <= time_to_infxn ~ duration))
  
  # Prep case/controls
  cc.pathogen.temp <- cc.unmatched %>%
    dplyr::select(hospitalization_id:room_stay_id)
  
  cc.pathogen <- cc.unmatched %>%
    dplyr::select(matches(pathogen_category))
  
  cc <- dplyr::bind_cols(cc.pathogen.temp, cc.pathogen) %>%
    dplyr::rename(group = 3) %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::mutate(run = pathogen_category)
  
  cc <- cc %>%
    dplyr::left_join(adt.micro)
  
  # Add demographics
  cc.d <- cc %>%
    dplyr::left_join(dem_features, relationship = "many-to-many") %>%
    dplyr::filter(!is.na(sex) | !is.na(age))
  
  # Add antibiotics
  cc.da <- cc.d %>%
    dplyr::left_join(abx_features, relationship = "many-to-many")

  # Add comorbidities
  cc.dac <- cc.da %>%
    dplyr::left_join(elix_features, relationship = "many-to-many")
  
  # Add procedures
  cc.dacp <- cc.dac %>%
    dplyr::left_join(cpt_features, relationship = "many-to-many")
  
  # Add admission disposition
  cc.dacpa <- cc.dacp %>%
    dplyr::left_join(admt_features, relationship = "many-to-many")
  
  # Add prior occupant pathogen history
  cc.dacpapo <- cc.dacpa %>%
    dplyr::left_join(prior_occupant_features, relationship = "many-to-many")
  
  # Add colonization pressure
  cc.unmatched.features <- cc.dacpapo %>%
    dplyr::left_join(cp_features, relationship = "many-to-many")
  
  return(cc.unmatched.features)
  }

#### METAFUNCTION 4: ENVIRONMENTAL-LEVEL ANALYSIS: MATCHING ALGORITHM ####

environmental.matching_process <- function(pathogen_category, 
                                           dat) 
  {
  
  # Cohort identifier
  dat <- dat %>% 
    dplyr::mutate(group_binary = ifelse(group == "case", 1, 0))
  
  # Drop NAs in relevant matching columns
  dat <- dat %>% 
    tidyr::drop_na(age, sex, elix_index_mortality, any_surgery)
  
  # Match on age, sex, surgical history, abx exposures, LOS
  matched_out <- match.data(matchit(group_binary ~ age + sex + any_surgery + matching_duration + 
                                      penicillin_0_60 + extended_spectrum_penicillin_0_60 + cephalosporin_0_60 + 
                                      extended_spectrum_cephalosporin_0_60 + carbapenem_0_60 + anti_staph_beta_lactam_0_60 + 
                                      fluoroquinolone_0_60 + glycopeptide_0_60 + anti_anaerobe_0_60 + anti_Cdiff_0_60 + 
                                      tetracycline_0_60 + macrolide_0_60 + sulfonamide_0_60 + lincosamide_0_60,  
                                    data = dat, 
                                    method = "nearest", 
                                    exact = c("sex", "any_surgery", "penicillin_0_60", "extended_spectrum_penicillin_0_60", "cephalosporin_0_60", 
                                              "extended_spectrum_cephalosporin_0_60", "carbapenem_0_60", "anti_staph_beta_lactam_0_60", 
                                              "fluoroquinolone_0_60", "glycopeptide_0_60", "anti_anaerobe_0_60", "anti_Cdiff_0_60",
                                              "tetracycline_0_60", "macrolide_0_60", "sulfonamide_0_60", "lincosamide_0_60"), 
                                    caliper = 1,
                                    ratio = 3))
  
  m_out <- matched_out %>% 
    dplyr::group_by(subclass) %>% 
    dplyr::mutate(age_diff = abs(age[group == "case"] - age)) %>% 
    # dplyr::mutate(elix_diff = abs(elix_index_mortality[group == "case"] - elix_index_mortality)) %>%
    dplyr::mutate(matching_duration_diff = abs(matching_duration[group == "case"] - matching_duration))
  
  # Apply tolerance for some matching
  m_out <- m_out %>% 
    dplyr::filter(group == "case" | group == "control" & (age_diff <= 5)) %>% 
    # dplyr::filter(group == "case" | group == "control" & (elix_diff <= 3)) %>% 
    dplyr::filter(group == "case" | group == "control" & (matching_duration_diff <= 3 | matching_duration <= (0.1 * matching_duration[group == "case"]))) %>% 
    dplyr::filter(any(group == "case") & any(group == "control")) 
  
  # Final clean up
  m_out <- m_out %>%
    dplyr::mutate(group_index = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(group_index, group)
  
  return(m_out)
  
  }

#### METAFUNCTION 5: PATIENT-LEVEL ANALYSIS: MATCHING ALGORITHM ####

patient.matching_process <- function(pathogen_category,
                                     dat) 
  {
  
  # Cohort identifier
  dat <- dat %>% 
    dplyr::mutate(group_binary = ifelse(group == "case", 1, 0))
  
  # Drop NAs in relevant matching columns
  dat <- dat %>% 
    tidyr::drop_na(dept_code, DTS_in_month, DTS_in_year, duration)
  
  # Match on ward, month and length of stay in first room +/- 3 days or 10%
  matched_out <- match.data(matchit(group_binary ~ dept_code + DTS_in_month + DTS_in_year + matching_duration,
                                    data = dat, 
                                    method = "nearest", 
                                    exact = c("dept_code", "DTS_in_month", "DTS_in_year"), 
                                    caliper = 1,
                                    ratio = 3))
  
  m_out <- matched_out %>% 
    dplyr::group_by(subclass) %>% 
    dplyr::mutate(matching_duration_diff = abs(matching_duration[group == "case"] - matching_duration)) 
  
  # Apply tolerance for some matching
  m_out <- m_out %>% 
    dplyr::filter(group == "case" | group == "control" & (matching_duration_diff <= 3 | matching_duration <= (0.1 * matching_duration[group == "case"]))) %>% 
    dplyr::filter(any(group == "case") & any(group == "control")) 
  
  # Final clean up
  m_out <- m_out %>%
    dplyr::mutate(group_index = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(group_index, group)
  
  return(m_out)
}
