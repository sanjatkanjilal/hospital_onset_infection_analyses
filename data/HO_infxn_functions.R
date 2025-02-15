# /usr/bin/R

#####################################################################################################
# Description: Script to define functions for building cohorts for hospital-onset infection analysis: 
#
# ENVIRONMENTAL FACTOR ANALYSIS
# Case / control matching criteria: 
#   - Age +/- 3 years
#   - Sex
#   - Comorbidities +/- 3 (Elixhauser index)
#   - Prior surgery in the previous 90 days
#   - Prior antibiotic exposure in the previous 60 days
#
# PATIENT-LEVEL ANALYSIS
# Case / control matching criteria: 
#   - Ward 
#   - Month / year
#   - Length of stay +/- 3 days or 10%
#
# Author: Luke Sagers, Ziming (Alex) Wei, Sanjat Kanjilal
# Last updated: 2024-05-08
#####################################################################################################

#### LIBRARIES ####
library(tidyverse)
# library(data.table)
library(reshape2)
library(MatchIt)
library(tictoc)
library(conflicted)
library(parallel)
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

#### FUNCTION SET 1: HELPER FUNCTIONS ####

# Function to group room stays within a specific time interval (in days) with each other
# 'type' refers to the user-defined name for the grouping (ie episodes, infections; defined below)
grouping_function <- function(data, id, index, start, end, time_interval) 
  {
  
  data <- data %>%  
    arrange({{id}}, {{index}}) %>%
    mutate(group = ifelse({{index}} == 1, 1, NA)) %>%
    select({{id}}, {{index}}, {{start}}, {{end}}, group) %>%
    group_by({{id}}) %>%
    mutate(tempdate = lag({{end}})) %>%
    ungroup() %>%
    mutate(tempdiff = interval(tempdate, {{start}})/days(x = 1)) %>%
    group_by({{id}})
  
  print("Calculation of time intervals complete")
  
  group1 <- data %>% 
    filter(!is.na(group))
  
  mydata <- tibble()
  
  while(sum(is.na(data$group)) > 0) {
    
    if (sum(is.na(data$group)) == 0)  {
      
      break
      
    } else {
      
      group.temp <- data %>% 
        mutate(lag_group = lag(group)) %>%
        filter(is.na(group)) 
      
      group.temp1 <- group.temp %>%
        mutate(group = if_else(tempdiff[1] <= time_interval, lag_group[1], lag_group[1] + 1)) %>%
        mutate(group = if_else({{index}} == {{index}}[1], group[1], NA))
      
      group.new <- group.temp1 %>%
        filter(!is.na(group))
      
      data <- group.temp1
      
      mydata <- bind_rows(mydata, group.new)

      print("Iteration complete")
      
          }
  }
  
  output <- full_join(mydata, group1) %>%
    ungroup() %>%
    select(-tempdate, -tempdiff, -lag_group) %>% 
    arrange({{id}}, {{start}})
  
  print("Grouping complete")

    return(output)
  }

#### FUNCTION SET 2:GLOBAL DATASET PRE-PROCESSING ####

# Drop repeat micro specimens with the same organisms drawn within 30 days of a first positive
micro_dedup <- function(micro)
{
  micro.dedup <- micro %>%
    # head(1000) %>%
    arrange(PatientID, org_group_3, coll_datetime_UTC) %>%
    group_by(PatientID, org_group_3) %>%
    mutate(sample_index = row_number(), 
           prev_sample = lag(coll_datetime_UTC)) %>%
    ungroup() %>%
    mutate(micro_diff = difftime(coll_datetime_UTC, prev_sample, units = "days"),
           micro_diff = round(micro_diff, digits = 1),
           flag = ifelse(micro_diff <= 30, "X", NA)) %>%
    filter(is.na(flag))
  
  return(micro.dedup)
  
}

# Process IP antibiotic data prior to joining to case / control data
initial_abx_clean <- function(abx, major)
{
  abx_prelim <- abx %>%
    # head(10000) %>%
    filter(!grepl("need|prn|on call", order_freq, ignore.case = T, perl = T)) %>%
    filter(!(is.na(order_dose) & is.na(order_unit) & is.na(order_freq) & is.na(order_quantity))) %>%
    mutate(order_DT = floor_date(order_DTS, unit = "days"),
           order_DT = ymd(order_DT),
           taken_DT = floor_date(taken_DTS, unit = "days"),
           taken_DT = ymd(taken_DT)) %>%
    group_by(patientID, order_DT, antibiotic) %>%
    mutate(index_order = row_number()) %>% 
    filter(index_order == 1) %>%
    group_by(patientID, taken_DT, antibiotic) %>%
    mutate(index_taken = row_number()) %>% 
    ungroup() %>%
    mutate(drop = ifelse(!is.na(taken_DT) & index_taken > 1, "X", NA)) %>%
    filter(is.na(drop)) %>%
    select(-drop) %>%
    filter(antibiotic %in% major$antibiotic) %>% 
    left_join(major %>% select(antibiotic, drug_class), by = c('antibiotic' = 'antibiotic')) %>%
    mutate(drug_class = ifelse(antibiotic == "VANCOMYCIN" & route_abbr == "PO", "anti_Cdiff", drug_class)) %>%
    filter(!is.na(drug_class)) %>%
    arrange(patientID, antibiotic, order_DT, taken_DT) %>%
    select(patientID, order_DT, taken_DT, antibiotic, drug_class, order_freq, order_quantity, route_abbr) %>%
    distinct()
  
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
    mutate(dt_lag = if_else(is.na(dt_lag), taken_DT, dt_lag)) %>%
    mutate(day_diff = taken_DT - dt_lag) %>%
    select(-dt_lag) %>% 
    mutate(new_prescription = case_when(day_diff == 0 ~ 1, 
                                        day_diff >= 7 ~ 1)) %>% 
    group_by(patientID, drug_class) %>% 
    mutate(new_prescription_1 = cumsum(!is.na(new_prescription))) %>%
    ungroup()
  
  # Get start and end dates
  abx.ip.prelim <- as.data.table(abx.ip.prelim)
  abx.ip.prelim[, 
                `:=`(start = min(taken_DT), end = max(taken_DT)), 
                by = .(patientID, drug_class, new_prescription_1)]
  
  # Calculate durations and keep only summary row
  abx_ip_durations <- abx.ip.prelim %>% 
    mutate(abx_duration_days = as.numeric(end - start + 1)) %>%
    tidytable::slice_tail(n=1, .by = c(patientID, drug_class, new_prescription_1)) %>% 
    ungroup() %>% 
    rename(end_DT = taken_DT) %>% 
    select(patientID, drug_class, end_DT, abx_duration_days) %>%
    distinct()
  
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
    dplyr::mutate(numeric_frequency_per_day = case_when(
      str_detect(tolower(order_freq), "(?i)^every\\s+\\d{1,2}\\s+hours") ~ 24/as.numeric(str_match(tolower(order_freq), "(?i)^every\\s+(\\d{1,2})\\s+hours")[,2]),
      str_detect(tolower(order_freq), "(?i)once") ~ 1,
      str_detect(tolower(order_freq), "(?i)\\d+\\s+times\\s+daily") ~ as.numeric(str_match(tolower(order_freq), "(?i)(\\d+)\\s+times\\s+daily")[,2]),
      str_detect(tolower(order_freq), "(?i)\\d+\\s+times\\s+weekly") ~ as.numeric(str_match(tolower(order_freq), "(?i)(\\d+)\\s+times\\s+weekly")[,2])/7,
      str_detect(tolower(order_freq), "(?i)^every\\s+\\d{1,2}\\s+days") ~ 1/as.numeric(str_match(tolower(order_freq), "(?i)^every\\s+(\\d{1,2})\\s+days")[,2]),
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
    dplyr::mutate(quantity_numeric = case_when(
      str_detect(tolower(order_quantity), "(?i)\\d+\\s+tablet") ~ as.numeric(str_match(tolower(order_quantity), "(?i)(\\d+)\\s+tablet")[,2]),
      str_detect(tolower(order_quantity), "(?i)\\d+\\s+capsule") ~ as.numeric(str_match(tolower(order_quantity), "(?i)(\\d+)\\s+capsule")[,2]),
      str_detect(tolower(order_quantity), "(?i)\\d+\\s+packet") ~ as.numeric(str_match(tolower(order_quantity), "(?i)(\\d+)\\s+packet")[,2]),
      str_detect(tolower(order_quantity), "(?i)\\d+\\s+dose") ~ as.numeric(str_match(tolower(order_quantity), "(?i)(\\d+)\\s+dose")[,2])
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
    arrange(patientID, drug_class, order_DT, end_DT) %>%
    unite(key, patientID, drug_class, remove = F, sep = "-") %>%
    group_by(key, drug_class, order_DT) %>%
    mutate(course_index_raw = row_number()) %>%
    filter(course_index_raw == max(course_index_raw)) %>%
    ungroup() %>%
    select(-course_index_raw) %>%
    distinct() %>%
    group_by(key) %>%
    mutate(course_index = row_number()) %>%
    ungroup()

  print("Performing grouping")
  abx_op_episodes <- grouping_function(data = abx_op_durations.prelim, id = key, index = course_index, 
                                       start = order_DT, end = end_DT, time_interval = 7) 

  print("Calculating courses after grouping")
  abx_op_episodes <- as.data.table(abx_op_episodes)
  abx_op_episodes <- abx_op_episodes[, `:=`(start_new = min(order_DT), end_new = max(end_DT)), by = .(key, group)]
  abx_op_episodes <- abx_op_episodes %>%
    select(-course_index, -order_DT, -end_DT) %>%
    distinct() %>%
    mutate(abx_duration_days_new = end_new - start_new) 

  print("Join back to main abx dataset")
  abx_op_durations <- abx_op_durations.prelim %>%
    left_join(abx_op_episodes, relationship = "many-to-many") %>%
    select(key, end_DT = end_new, abx_duration_days = abx_duration_days_new) %>%
    distinct() %>%
    mutate(abx_duration_days = as.integer(abx_duration_days)) %>%
    separate(key, c("patientID", "drug_class"), sep = "-")
  
  toc()
  
  return(abx_op_durations)
}

# Combine IP and OP abx courses
abx_courses <- function(abx_ip_durations, abx_op_durations)
{
  abx.courses <- bind_rows(abx_ip_durations, abx_op_durations) %>%
    arrange(patientID, drug_class, end_DT)
  
  return(abx.courses)
}

# Process encounter data
enc_processing <- function(enc)
{
  # Filter to encounters lasting >48 hours and then indexing

  tic()
  print("Initial cleaning of encounter data")
  enc.prelim <- enc %>%
    filter(visit_type == "Inpatient") %>%
    select(PatientID, visit_type, hospitalAdmitDTS, hospitalDischargeDTS) %>%
    distinct() %>%
    arrange(PatientID, hospitalAdmitDTS) %>%
    mutate(duration = as.numeric(as.difftime(hospitalDischargeDTS - hospitalAdmitDTS), units = "days")) %>%
    filter(duration > 2) %>%
    group_by(PatientID) %>%
    mutate(enc_index = row_number()) %>%
    ungroup()
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
    select(PatientID, IP_enc_endDTS = end_new) %>%
    distinct()
  toc()

  return(enc_clean)
}

# Set target hierarchy for looping through pathogens 
path_cat_hierarchy <- function(micro, hierarchy)
{
  path_cat_table <- micro %>%
    gather(pathogen_hierarchy, pathogen_category, org_group_1:org_group_3) %>%
    filter(!is.na(pathogen_category)) %>%
    select(pathogen_hierarchy, pathogen_category) %>%
    distinct()
  
  path_cat_table <- path_cat_table %>%
    filter(pathogen_hierarchy == hierarchy) 
  
  return(path_cat_table)
  
}

# Map admission dispositions to clean categories
admt.mapped <- function(admt)
{
  admt_mapped <- admt %>% 
    select(PatientID, PatientEncounterID, AdmitSourceDSC) %>% 
    mutate(PatientEncounterID = as.numeric(PatientEncounterID)) %>%
    mutate(admit_source_clean = case_when(
      AdmitSourceDSC %in% c("Ambulatory Surgery Center", "Born Outside this Hospital", "Discharge and Readmit", "NSM Emergency",
                            "Outside Health Care Facility", "Outside Hospital", "Psych, Substance Abuse, or Rehab Hospital",
                            "Transfer from another Health Care Facility") ~ "Transfer from outside facility",
      AdmitSourceDSC %in% c("Court/Law Enforcement", "Self Referral") ~ "Self referral",
      AdmitSourceDSC %in% c("Skilled Nursing Facility", "Hospice") ~ "Skilled Nursing facility",
      AdmitSourceDSC %in% c("Physician or Clinic Referral") ~ "Clinic referral",
      AdmitSourceDSC %in% c("", "Information Unavailable") ~ "Unknown",
      AdmitSourceDSC %in% c("Born Inside this Hospital") ~ "Born in hospital")) %>%
    select(PatientID, PatientEncounterID, admit_source_clean) %>%
    distinct()
  
  return(admt_mapped)
}

#### FUNCTION SET 3: BUILD ADT - MICRO GROUND TRUTH DATASET FOR ALL PATHOGENS ####

# Filter room data to identify eligible hospital stays (room stays within 30d of each other) and filter out ER, periop and residential stays

# Ennumerate rooms and for ones that have the same DTS_in, sort and choose the latest one only
first_eligible_rooms <- function(room_dat)
  {
  room_dat.eligible.prelim <- room_dat %>%
    # head(10000) %>%
    filter(!grepl("emerg|periop|unspecified|home ", DepartmentDSC, ignore.case = T, perl = T)) %>%
    arrange(PatientID, DTS_in, DTS_out) %>%
    group_by(PatientID, DTS_in) %>%
    mutate(room_stay_index_raw = row_number()) %>%
    filter(room_stay_index_raw == max(room_stay_index_raw)) %>%
    ungroup() %>%
    select(-room_stay_index_raw) %>%
    distinct() %>%
    group_by(PatientID) %>%
    mutate(room_stay_index = row_number()) %>%
    ungroup()
  
  # Group room stays within 7 days of each other into episodes
  room_dat.episode <- grouping_function(room_dat.eligible.prelim, id = PatientID, index = room_stay_index, 
                                        start = DTS_in, end = DTS_out, time_interval = 7)

  room_dat.eligible.prelim2 <- room_dat.eligible.prelim %>%
    left_join(room_dat.episode) %>%
    arrange(PatientID, DTS_in) %>%
    filter(duration >= 2) %>% # Trim to rooms where patient stayed between 3 & 30 days
    mutate(drop = ifelse(duration > 30, "X", NA)) %>% # Drop all room stays that occur after any room stay >30 days in duration
    group_by(PatientID, group) %>%
    fill(drop, .direction = "downup") %>%
    ungroup() %>%
    filter(is.na(drop)) %>%
    group_by(PatientID, group) %>% # Re-index room stays after having taken out redundant rooms and rooms with durations outside of 48 to 720 hrs
    mutate(room_stay_index = row_number(),
           n_rooms_hosp = max(room_stay_index),
           first_room = ifelse(room_stay_index == 1, "X", NA),
           hospitalization_id = cur_group_id()) %>%
    ungroup() %>%
    group_by(PatientID, group, room_stay_index) %>% 
    mutate(room_stay_id = cur_group_id()) %>%
    ungroup() %>%
    select(PatientID:PatientClassDSC, DepartmentDSC, dept_room, first_room, room_stay_index, n_rooms_hosp, 
           room_stay_id, hospitalization_index = group, hospitalization_id, DTS_in, DTS_out, duration) %>%
    distinct()
  
  rm(room_dat.episode)
  
  return(room_dat.eligible.prelim2)
  }

# Filter micro data for patients in ADT dataset
micro_prep <- function(micro.dedup, room_dat.eligible.prelim2)
  {
  micro.eligible <- micro.dedup %>% 
    filter(PatientID %in% room_dat.eligible.prelim2$PatientID) %>% 
    select(PatientID, accession, coll_datetime_UTC, org_group_1:org_group_3) %>%
    distinct()
  
  return(micro.eligible)
  }

# Ensure room data begins after the start of the micro data
room_dat_eligible_trimmed <- function(room_dat.eligible.prelim2, micro.eligible)
  {
  room_dat.eligible.trimmed <- room_dat.eligible.prelim2 %>%
    filter(DTS_in > min(micro.eligible$coll_datetime_UTC)) 
  
  return(room_dat.eligible.trimmed)
  }

# Add micro data to ADT / room data
adt_micro_join <- function(room_dat.eligible.trimmed, micro.eligible)
  {
  adt.micro.raw <- room_dat.eligible.trimmed %>%
    left_join(micro.eligible, by = c("PatientID"), relationship = "many-to-many") %>%
    distinct() %>%
    arrange(hospitalization_id, coll_datetime_UTC)
  
  return(adt.micro.raw)
  }

# Add time to infection relative to first room in a hospitalization 
time_to_infxn <- function(adt.micro.raw)
  {
  adt.micro.raw <- adt.micro.raw %>%
    mutate(time_to_infxn = difftime(coll_datetime_UTC, DTS_in, units = "days"),
           time_to_infxn = round(time_to_infxn, digits = 1))
  
  return(adt.micro.raw)
  }

# Drop micro data for rooms where sample was drawn >1 year prior to and >30 days after entry (skin flora & Enterococcus)
# or >6 months prior to and >30 days after entry (other)
drop_old_micro <- function(adt.micro.raw)
  {
  adt.micro.raw <- adt.micro.raw %>%
    mutate(accession = case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        (time_to_infxn < -365 | time_to_infxn > 30) ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        (time_to_infxn < -182.5 | time_to_infxn > 30) ~ NA,
      .default = accession)) %>%
    mutate(coll_datetime_UTC = case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = coll_datetime_UTC)) %>%
    mutate(org_group_1 = case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = org_group_1)) %>%
    mutate(org_group_2 = case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = org_group_2)) %>%
    mutate(org_group_3 = case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = org_group_3)) %>%
    mutate(time_to_infxn = case_when(
      org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -365 | time_to_infxn > 30 ~ NA,
      !org_group_3 %in% c("MSSA", "MRSA", "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") &
        time_to_infxn < -182.5 | time_to_infxn > 30 ~ NA,
      .default = time_to_infxn)) %>%
    distinct() %>%
  group_by(room_stay_id) %>%
  mutate(group_size = n(),
         drop = ifelse(group_size > 1 & is.na(accession), "X", NA)) %>%
  ungroup() %>%
  filter(is.na(drop)) %>%
  select(-group_size, -drop)
  
  return(adt.micro.raw)
  }

# Drop hospitalizations where the patient received antibiotics within 7 days of their stay in their 1st room 
remove_treated_patients <- function(adt.micro.raw, abx.courses, day_threshold)
  {
  abx.courses.filtered <- abx.courses %>%
    rename(PatientID = patientID) %>%
    filter(abx_duration_days >=2)
  
  adt.micro.raw <- adt.micro.raw %>%
    # head(10000) %>%
    left_join(abx.courses.filtered, relationship = "many-to-many") %>%
    mutate(diff_days = as.Date(floor_date(end_DT, unit = "days")) - as.Date(floor_date(DTS_in, unit = "days"))) %>%
    mutate(drop = ifelse(diff_days >= -7 & diff_days <= day_threshold & !is.na(first_room), "X", NA)) %>%
    group_by(hospitalization_id) %>%
    fill(drop, .direction = "downup") %>%
    ungroup() %>%
    filter(is.na(drop)) %>%
    select(PatientID:time_to_infxn) %>%
    distinct()
  
  return(adt.micro.raw)
}

# Drop hospitalizations where the patient had a previous IP encounter >48 hours in length
drop_prev_enc <- function(adt.micro.raw, enc_clean)
  {
  adt.micro.raw <- adt.micro.raw %>%
    left_join(enc_clean, relationship = "many-to-many") %>%
    mutate(diff_days  = as.Date(floor_date(IP_enc_endDTS, unit = "days")) - as.Date(floor_date(DTS_in, unit = "days"))) %>%
    mutate(drop = ifelse(diff_days < 0 & diff_days > -90, "X", NA)) %>%
    group_by(hospitalization_id) %>%
    fill(drop, .direction = "downup") %>%
    ungroup() %>%
    filter(is.na(drop)) %>%
    select(PatientID:time_to_infxn) %>%
    distinct()
  
  return(adt.micro.raw)
}

# Add feature for matching cases (time to infxn) to control (LOS in room)
matching_duration <- function(adt.micro.raw)
  {
  adt.micro.raw <- adt.micro.raw %>%
    mutate(time_to_infxn = as.numeric(time_to_infxn)) %>%
    mutate(matching_duration = case_when(
       time_to_infxn > 0 ~ time_to_infxn,
       .default = duration))
  
  return(adt.micro.raw)
  
}

#### FUNCTION SET 4: PATHOGEN-SPECIFIC CASE / CONTROL SELECTION ####

# Select proper pathogen hierarchy in micro data
adt_micro_hierarchy <- function(adt.micro.final, pathogen_hierarchy)
{
  adt.micro.raw.temp <- adt.micro.raw %>%
    select(PatientID:coll_datetime_UTC, time_to_infxn)
  
  adt.micro.raw.hierarchy <- adt.micro.raw %>%
    select(org_group_1:org_group_3) %>%
    select(matches(pathogen_hierarchy)) %>%
    rename(org_group = matches(pathogen_hierarchy))
  
  adt.micro.raw.filtered <- bind_cols(adt.micro.raw.temp, adt.micro.raw.hierarchy) %>%
    select(PatientID:coll_datetime_UTC, org_group, time_to_infxn)
  
  rm(adt.micro.raw.temp, adt.micro.raw.hierarchy)
  
  return(adt.micro.raw.filtered)
}

# Flag hospitalizations where the pathogen of interest was identified after day X (set by user) and before day 30 in the first eligible room 
potential_cases <- function(adt.micro.raw.filtered, pathogen_category, day_threshold) 
  {
  adt.micro.raw.filtered <- adt.micro.raw.filtered %>% 
    mutate(potential_case = ifelse(!is.na(first_room) &
                                     org_group == pathogen_category & 
                                     time_to_infxn >= day_threshold & time_to_infxn <= 30, "X", NA))
  
  return(adt.micro.raw.filtered)
  }

# Drop potential cases where micro timestamp falls during occupation of a later room during that hospitalization episode
after_first_room <- function(adt.micro.raw.filtered, pathogen_category, day_threshold) 
{
  secondary.case.diff.room <- adt.micro.raw.filtered %>%
    filter(org_group == pathogen_category) %>%
    select(PatientID, hospitalization_id, room_stay_id, org_group, accession, time_to_infxn) %>%
    distinct() %>%
    filter(time_to_infxn >= day_threshold & time_to_infxn <= 30) %>%
    group_by(hospitalization_id, accession) %>%
    mutate(num_rooms = n()) %>%
    mutate(drop = ifelse(num_rooms > 1, "X", NA)) %>%
    fill(drop, .direction = "downup") %>%
    ungroup() %>%
    filter(!is.na(drop))
  
  adt.micro.raw.filtered <- adt.micro.raw.filtered %>%
    filter(!hospitalization_id %in% secondary.case.diff.room$hospitalization_id) 
  
  return(adt.micro.raw.filtered)
}

# Drop potential cases where the patient had the pathogen of interest between day -365 (S. aureus, Enterococcus) / day -182.5 (all others)
# and day +X (set by user), relative to entry into first room (suggesting not nosocomially acquired)
flag_prior_infxn <- function(adt.micro.raw.filtered, pathogen_category, day_threshold) 
  {
  prior_infxn <- adt.micro.raw.filtered %>% 
    arrange(PatientID, coll_datetime_UTC) %>%
    filter(org_group == pathogen_category) %>%
    mutate(prior_infxn = 
             case_when(
               !is.na(first_room) &
                 org_group %in% c("Skin flora", "Staph_aureus", "MSSA", "MRSA", "E_faecalis", "E_faecium", 
                                  "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") & 
                 time_to_infxn >= -365 & time_to_infxn <= day_threshold ~ "X",
               !is.na(first_room) &
                 !org_group %in% c("Skin flora", "Staph_aureus", "MSSA", "MRSA", "E_faecalis", "E_faecium", 
                                   "VSE_faecalis", "VSE_faecium", "VRE_faecalis", "VRE_faecium") & 
                 time_to_infxn >= -182.5 & time_to_infxn <= day_threshold ~ "X",
               .default = NA)) %>%
    group_by(hospitalization_id) %>%
    fill(prior_infxn, .direction = "downup") %>%
    ungroup() %>%
    filter(!is.na(prior_infxn))
  
  adt.micro.cases <- adt.micro.raw.filtered %>%
    filter(!hospitalization_id %in% prior_infxn$hospitalization_id) 
  
  return(adt.micro.cases)
}

# Identify controls (for matching later)
generate_controls <- function(adt.micro.cases, pathogen_category) 
{
  adt.micro.cc.temp <- adt.micro.cases %>%
    mutate(group = ifelse(org_group == pathogen_category & !is.na(potential_case), "case", NA)) %>%
    mutate(group = ifelse(org_group != pathogen_category & !is.na(first_room), "control", group)) %>%
    filter(!is.na(group))
  
  cc.overlap <- adt.micro.cc.temp %>%
    select(hospitalization_id, group) %>%
    distinct() %>%
    group_by(hospitalization_id) %>%
    mutate(overlap = n(), 
           overlap_flag = ifelse(overlap > 1, "X", NA)) %>%
    ungroup() %>%
    filter(!is.na(overlap_flag))
  
  adt.micro.cc <- adt.micro.cc.temp %>%
    mutate(group = ifelse(hospitalization_id %in% cc.overlap$hospitalization_id & group == "control", NA, group)) %>%
    filter(!is.na(group)) %>%
    select(hospitalization_id, room_stay_id, group) %>% 
    distinct()
  
  return(adt.micro.cc)
}

#### FUNCTION SET 5: GENERATE FEATURES FOR MATCHING AND MODELING ####

# Demographics
add_demographics <- function(adt.micro.raw, dems) 
  {
  dem_features <- adt.micro.raw %>%
    select(PatientID, room_stay_id, hospitalization_id, DTS_in) %>%
    distinct() %>%
    left_join(dems) %>%
    mutate(age = round(time_length(difftime(DTS_in, DOB), "years"), digits = 1)) %>%
    mutate(sex = ifelse(sex %in% c("Male", "Female"), sex, "Other")) %>%
    mutate(across(where(is.character), ~na_if(., ""))) %>%
    select(-DTS_in, -DOB)
  
  # Impute missing values for admission source and race/ethnicity
  dem_features <- dem_features %>% 
    mutate(race = ifelse(is.na(race), "Declined", race)) %>%
    mutate(ethnicity = ifelse(is.na(ethnicity), "Prefer not to say/Decline", ethnicity)) 

  return(dem_features)
}

# Antibiotics
add_abx <- function(abx.courses, adt.micro.raw) 
  {
  # Code to create the feature set for inpatient abx history
  abx_features.prelim <- adt.micro.raw %>% 
    select(PatientID, room_stay_id, DTS_in, coll_datetime_UTC) %>%
    distinct() %>%
    # head(1000) %>%
    left_join(abx.courses %>% select(PatientID = patientID, drug_class, end_DT, abx_duration_days), 
              relationship = "many-to-many") %>% 
    mutate(diff_days = as.Date(floor_date(DTS_in, unit = "days")) - end_DT) %>%
    filter(diff_days > 0 & diff_days < 182.5) %>% 
    mutate(interval = case_when(diff_days <= 60 ~ "0_60",
                                diff_days > 60 ~ "60_plus")) %>%
    mutate(duration_binary = case_when(abx_duration_days <= 2 ~ "under2days",
                                       abx_duration_days > 2 ~ "over2days")) %>% 
    filter(duration_binary == "over2days") %>%
    mutate(val = 1) %>% 
    distinct() %>% 
    group_by(PatientID) %>% 
    pivot_wider(names_from = c(drug_class, duration_binary, interval), names_sep = "_", values_from = val) %>% 
    ungroup() %>%
    select(-c(DTS_in, end_DT, abx_duration_days, diff_days)) %>%
    distinct() %>%
    mutate_if(is.numeric, list(~ replace_na(., 0)))
  
  return(abx_features.prelim)
  }

final_abx_clean <- function(abx_features.prelim, id_vars) 
  {
  id_vars = c("PatientID", "room_stay_id", "coll_datetime_UTC")
  
  abx_features.melt <- reshape2::melt(abx_features.prelim, id.vars = id_vars)
  
  abx_features.melt <- abx_features.melt %>%
    mutate(variable = gsub("_over2days", "", variable))
  
  abx_features <- reshape2::dcast(abx_features.melt, PatientID + room_stay_id ~ variable, value.var = "value", fun.aggregate = sum)
  
  abx_features <- abx_features %>%
    mutate_if(is.numeric, list(~ replace_na(., 0)))
  
  abx_of_interest <- c("penicillin_0_60", "extended_spectrum_penicillin_0_60", "cephalosporin_0_60", 
                       "extended_spectrum_cephalosporin_0_60", "carbapenem_0_60", "anti_staph_beta_lactam_0_60", 
                       "fluoroquinolone_0_60", "glycopeptide_0_60", "anti_anaerobe_0_60", "anti_Cdiff_0_60",
                       "tetracycline_0_60", "macrolide_0_60", "sulfonamide_0_60", "lincosamide_0_60")
  
  variables_to_group <- names(abx_features)[grepl("_0_60$", names(abx_features)) & !(names(abx_features) %in% abx_of_interest)]
  
  # Add flag for antibiotics not in existing categories
  abx_features <- abx_features %>%
    mutate(other_abx_0_60 = as.integer(rowSums(select(., all_of(variables_to_group)), na.rm = TRUE) >= 1))
  
  abx_features <- abx_features %>% 
    mutate(any_abx_0_60 = ifelse(rowSums(select(., contains("_0_60"))) > 0, 1, 0)) %>% 
    mutate(any_abx_60_plus = ifelse(rowSums(select(., contains("_60_plus"))) > 0, 1, 0)) 
  
  return(abx_features)
  }

# Comorbidities
add_comorbidities <- function(adt.micro.raw, elix)
{
  elix_features <- adt.micro.raw %>%
    select(PatientID, PatientEncounterID, hospitalization_id) %>%
    distinct() %>%
    left_join(elix)
  
  return(elix_features)
}

# CPT codes
add_cpt <- function(adt.micro.raw, cpt) 
  {
  cpt_features <- adt.micro.raw %>%
    select(PatientID, hospitalization_id, DTS_in) %>%
    distinct() %>%
    left_join(cpt, relationship = "many-to-many") %>% 
    filter(proc_date < DTS_in) %>% 
    mutate(diff_days = difftime(proc_date, DTS_in, units="days")) %>% 
    mutate(surgery_past_90 = ifelse((diff_days < 0) & (diff_days >= -90), proc_simple, NA)) %>% 
    filter(!is.na(surgery_past_90)) %>% 
    select(PatientID, hospitalization_id, surgery_past_90) %>% 
    distinct() %>%
    mutate(val = 1) %>% 
    pivot_wider(names_from = surgery_past_90, values_from = val) %>%
    mutate_if(is.numeric, list(~ replace_na(., 0))) %>%
    arrange(hospitalization_id, PatientID) %>%
    mutate(any_surgery = 1)
  
  cpt_features <- cpt_features %>%
    rename(ortho_sx = `Orthopedic surgery`, abd_sx = `Intra_abdominal surgery`, rp_sx = `Retroperitoneal surgery`,
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
    select(PatientID, PatientEncounterID, hospitalization_id) %>%
    distinct() %>%
    left_join(admt_mapped, relationship = "many-to-many") 
  
  # Impute missing values for admission source and race/ethnicity
  admt_features <- admt_features %>% 
    mutate(admit_source_clean = case_when(
      is.na(admit_source_clean) ~ "Information Unavailable",
      admit_source_clean %in% c("Self referral", "Clinic referral", "Skilled Nursing facility", 
                                 "Transfer from outside facility", "Information Unavailable") ~ admit_source_clean,
      .default = "Other or unknown"))
  
  # Re-level
  admt_features <- admt_features %>%
    mutate(admit_source_clean = factor(admit_source_clean)) %>%
    mutate(admit_source_clean = fct_relevel(admit_source_clean, "Self referral", after = 0)) %>%
    mutate(admit_source_clean = fct_relevel(admit_source_clean, "Information Unavailable", after = 5))
  
  return(admt_features)
}

# Prior room occupant
add_prior_pathogens <- function(room_dat, micro, pathogen_hierarchy, adt.micro.raw) 
  {
  # Prep total room dataset
  room_dat.simple <- room_dat %>%
    select(PatientID, PatientEncounterID, dept_room, DTS_in, DTS_out, duration) %>%
    distinct() 
  
  # Prep total micro dataset
  micro.simple <- micro %>%
    select(PatientID, PatientEncounterID = CSN, coll_datetime_UTC) 
  
  micro.hierarchy <- micro %>%
    select(org_group_1:org_group_3) %>%
    select(matches(pathogen_hierarchy)) %>%
    rename(org_group = matches(pathogen_hierarchy))
  
  micro.filtered <- bind_cols(micro.simple, micro.hierarchy)
  
  # Join total room and total micro datasets
  room_micro <- room_dat.simple %>%
    left_join(micro.filtered, relationship = "many-to-many") %>%
    filter(coll_datetime_UTC > DTS_in & coll_datetime_UTC < DTS_out) %>%
    mutate(prior_colldate = floor_date(coll_datetime_UTC, unit = "days")) %>%
    select(dept_room, prior_PatientID = PatientID, prior_colldate, prior_org_group = org_group) %>%
    distinct()
    
  # Join total room+micro data to clean ADT+micro dataset and filter for when a different patient had a pathogen in the previous 365 days
  prior_pathogen.prelim <- adt.micro.raw %>%
    # head(1000) %>%
    filter(!is.na(first_room)) %>%
    mutate(DT_in = floor_date(DTS_in, unit = "days")) %>%
    select(PatientID, dept_room, hospitalization_id, DT_in) %>%
    distinct() %>%
    left_join(room_micro, relationship = "many-to-many") %>%
    mutate(diff_days = difftime(DT_in, prior_colldate, units="days")) %>%
    filter(diff_days > 0 & diff_days < 365) %>%
    filter(PatientID != prior_PatientID)
  
  # Select closest date to DTS_in if pathogen present more than once in the past
  prior_pathogen.prelim.dedup <- as.data.table(prior_pathogen.prelim)
  prior_pathogen.prelim.dedup <- prior_pathogen.prelim.dedup[order(dept_room, prior_org_group, diff_days)]
  prior_pathogen.prelim.dedup[,
                         diff_max := max(diff_days),
                         by = .(dept_room, prior_org_group)]
  prior_pathogen.prelim.dedup <- prior_pathogen.prelim.dedup %>%
    filter(diff_days == diff_max) %>%
    mutate(diff_days = as.integer(diff_days))
  
  # Negative exponential transform of days since last pathogen
  diff_days_transform <- prior_pathogen.prelim.dedup %>% 
    select(dept_room, hospitalization_id, diff_days) %>%
    mutate(diff_days_transform = exp(-0.01 * (diff_days))) %>%
    mutate(diff_days_transform = round(diff_days_transform, 2)) %>%
    mutate(diff_days_transform = diff_days_transform * 5) # Additional up-weighting of prior occupant above that of other ward occupants

  # Spread transformed data
  prior_pathogen.prelim.transform <- prior_pathogen.prelim.dedup %>%
    left_join(diff_days_transform, relationship = "many-to-many") %>%
    select(PatientID:hospitalization_id, prior_org_group, diff_days_transform) %>%
    distinct() %>%
    mutate(prior_org_group = paste0("prior_", prior_org_group)) %>% 
    spread(prior_org_group, diff_days_transform) %>%
    arrange(hospitalization_id) %>%
    mutate_if(is.numeric, list(~ replace_na(., 0)))
  
   # Add back to main ADT / micro dataset 
  prior_occupant_features <- adt.micro.raw %>%
    select(PatientID, dept_room, hospitalization_id) %>%
    distinct() %>%
    left_join(prior_pathogen.prelim.transform) %>%
    mutate_if(is.numeric, list(~ replace_na(., 0)))

  return(prior_occupant_features)
  }

# Colonization pressure
calculate_CP_score <- function(date, location_specific, room_dat, numerator_category, micro_df, lambda = 0.01) 
  {
  # Define the date 30 days prior to the patient of interest 
  ward_previous_start_date <- date - days(30)
  
  prev_month_occupants <- room_dat %>% 
    filter(DepartmentDSC == location_specific) %>% 
    filter((DTS_in >= ward_previous_start_date & DTS_in < date) | (DTS_out > ward_previous_start_date & DTS_in < ward_previous_start_date))
  
  # Filter the microbiology data by date range and hospital
  filtered_data <- micro_df %>%
    filter(PatientID %in% prev_month_occupants$PatientID) %>% 
    filter(coll_datetime_UTC >= date - days(90), coll_datetime_UTC <= date)
  
  # Filter the data based on the numerator and denominator categories
  numerator_categories <- strsplit(numerator_category, " \\+ ")[[1]]
  # denominator_categories <- strsplit(denominator_category, " \\+ ")[[1]]
  
  filtered_data <- filtered_data %>%
    filter(org_group_3 %in% numerator_categories) %>% 
    distinct(PatientID, loc_cat, loc_specific, coll_datetime_UTC, org_group_3)
  
  tic()
  
  # Calculate the colonization pressure
  if(nrow(filtered_data) == 0) {
    colonization_pressure <- 0
  } else {
    colonization_pressure <- filtered_data %>%
      group_by(PatientID, org_group_3) %>%
      filter(coll_datetime_UTC == max(coll_datetime_UTC)) %>% # Retain only the most recent positive sample
      summarize(days_since_last_positive = as.numeric(date - coll_datetime_UTC), .groups = "drop") %>%
      mutate(score = exp(-lambda * days_since_last_positive)) %>% # Exponential decay
      pull(score) %>% 
      sum(na.rm = TRUE)
  }
  
  return(colonization_pressure)
  
  toc()
  }

calculate_and_add_cp_scores_parallel <- function(adt.micro.raw, micro, room_dat, location_map) 
  {
  loc_map_single <- location_map %>% 
    distinct(DepartmentDSC, .keep_all = TRUE)
  
  cp_features <- adt.micro.raw %>%
    filter(!is.na(first_room)) %>%
    select(PatientID, hospitalization_id, DepartmentDSC, DTS_in) %>%
    distinct() %>%
    left_join(loc_map_single, by = "DepartmentDSC") %>%
    mutate(filter_date = DTS_in) %>%
    mutate(micro_loc = if_else(is.na(micro_loc), DepartmentDSC, micro_loc))
  
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
  room_dat.eligible <- first_eligible_rooms(room_dat)
  toc()
  
  # Filter micro data for patients in ADT dataset
  print("Filter micro data for patients in ADT dataset")
  tic()
  micro.eligible <- micro_prep(micro.dedup, room_dat.eligible)
  toc()
  
  # Ensure room data begins after the start of the micro data
  print("Ensure room data begins after the start of the micro data")
  tic()
  room_dat.eligible.trimmed <- room_dat_eligible_trimmed(room_dat.eligible, micro.eligible)
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
  toc()
  
  # Remove hospitalization episodes where the patient received inpatient antibiotics prior to their stay in the room based on when their abx course finished
  print("Remove potential cases where the patient received concomitant IP/OP antibiotics")
  tic()
  adt.micro.raw <- remove_treated_patients(adt.micro.raw, abx.courses, day_threshold)
  toc()
  
  # Remove hospitalization episodes where patients had another IP encounter >2d in duration in the previous 90 days
  print("Drop potential cases where the patient had another hospital encounter in the previous 90 days")
  tic()
  adt.micro.raw <- drop_prev_enc(adt.micro.raw, enc_clean)
  toc()
  
  # Add feature for matching cases (time to infxn) to control (LOS in room)
  print("Add feature for matching cases to control on time to infxn (cases) or LOS (controls)")
  tic()
  adt.micro.raw <- matching_duration(adt.micro.raw)
  toc()
  
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
  print(length(unique(adt.micro.raw.filtered[['PatientEncounterID']])))
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
  print(drop_prior_infxn)
  toc()
  
  potential_case <- adt.micro.cases %>% filter(potential_case == 'X')
  
  # ## Exclude hospitalizations where have prior antibiotic prescription
  # subgroup_abx_prelim <- abx_prelim %>% filter(patientID %in% potential_case$PatientID)
  # subgroup_abx_prelim <- merge(subgroup_abx_prelim,
  #                              potential_case,
  #                              by.x='patientID',
  #                              by.y='PatientID')
  # subgroup_abx_prelim <- subgroup_abx_prelim %>% select(patientID, hospitalization_id, taken_DT, DTS_in)
  # subgroup_abx_prelim[['abx_bound_lower']] = as.Date(subgroup_abx_prelim[['DTS_in']]) - 7
  # subgroup_abx_prelim[['abx_bound_upper']] = as.Date(subgroup_abx_prelim[['DTS_in']])
  # subgroup_abx_prelim <- subgroup_abx_prelim %>% filter(taken_DT > abx_bound_lower) %>% filter(taken_DT < abx_bound_upper)
  
  # adt.micro.cases <- adt.micro.cases %>% filter(! (hospitalization_id %in% subgroup_abx_prelim$hospitalization_id))
  
  
  drop_abx <- nrow(adt.micro.cases %>% filter(potential_case == 'X'))
  print(drop_abx)
  
  # Select controls (for matching later)
  print("Select preliminary controls")
  tic()
  adt.micro.cc <- generate_controls(adt.micro.cases, pathogen_category) 
  toc()
  
  row <- data.frame(
    pathogen = pathogen_category,
    initial_cases = initial_case,
    case_after_drop_subsequent = drop_subsequent,
    case_after_drop_prior_infection = drop_prior_infxn,
    case_after_drop_abx = drop_abx
  )
  
  flow_chart <<- rbind(flow_chart, row)
  
  return(adt.micro.cc)
}

#### METAFUNCTION 3: CREATION OF FINAL DATASET ####

dataset.build <- function(pathogen_category,
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
    mutate(DTS_in_month = floor_date(DTS_in, unit = "month"),
           DTS_in_year = floor_date(DTS_in, unit = "year")) %>%
    group_by(DepartmentDSC) %>%
    mutate(dept_code = cur_group_id()) %>%
    ungroup() %>%
    select(PatientID, PatientEncounterID, dept_code, hospitalization_id, room_stay_id, 
           DTS_in, DTS_in_month, DTS_in_year, duration, time_to_infxn, matching_duration) %>%
    distinct()
  
  # Prep case/controls
  cc.pathogen.temp <- cc.unmatched %>%
    select(hospitalization_id:room_stay_id)
  
  cc.pathogen <- cc.unmatched %>%
    select(matches(pathogen_category))
  
  cc <- bind_cols(cc.pathogen.temp, cc.pathogen) %>%
    rename(group = 3) %>%
    filter(!is.na(group)) %>%
    mutate(run = pathogen_category)
  
  cc <- cc %>%
    left_join(adt.micro)
  
  # Add demographics
  cc.d <- cc %>%
    left_join(dem_features, relationship = "many-to-many") %>%
    filter(!is.na(sex) | !is.na(age))
  
  # Add antibiotics
  cc.da <- cc.d %>%
    left_join(abx_features, relationship = "many-to-many")

  # Add comorbidities
  cc.dac <- cc.da %>%
    left_join(elix_features, relationship = "many-to-many")
  
  # Add procedures
  cc.dacp <- cc.dac %>%
    left_join(cpt_features, relationship = "many-to-many")
  
  # Add admission disposition
  cc.dacpa <- cc.dacp %>%
    left_join(admt_features, relationship = "many-to-many")
  
  # Add prior occupant pathogen history
  cc.dacpapo <- cc.dacpa %>%
    left_join(prior_occupant_features, relationship = "many-to-many")
  
  # Add colonization pressure
  cc.unmatched.features <- cc.dacpapo %>%
    left_join(cp_features, relationship = "many-to-many")
  
  return(cc.unmatched.features)
  }

#### METAFUNCTION 4: ENVIRONMENTAL-LEVEL ANALYSIS: MATCHING ALGORITHM ####

environmental.matching_process <- function(pathogen_category, 
                                           dat) 
  {
  
  # Cohort identifier
  dat <- dat %>% 
    mutate(group_binary = ifelse(group == "case", 1, 0))
  
  # Drop NAs in relevant matching columns
  dat <- dat %>% 
    drop_na(age, sex, indiv_score, any_surgery)
  
  # Match on age, sex, comorbidities, surgical history, abx exposures
  matched_out <- match.data(matchit(group_binary ~ age + sex + indiv_score + any_surgery + matching_duration + 
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
    group_by(subclass) %>% 
    mutate(age_diff = abs(age[group == "case"] - age)) %>% 
    mutate(elix_diff = abs(indiv_score[group == "case"] - indiv_score)) %>%
    mutate(matching_duration = abs(matching_duration[group == "case"] - matching_duration))
  
  # Apply tolerance for some matching
  m_out <- m_out %>% 
    filter(group == "case" | group == "control" & (age_diff <= 5)) %>% 
    filter(group == "case" | group == "control" & (elix_diff <= 5)) %>% 
    filter(group == "case" | group == "control" & (matching_duration <= 3 | matching_duration <= (0.1 * matching_duration[group == "case"]))) %>% 
    filter(any(group == "case") & any(group == "control")) 
  
  # Final clean up
  m_out <- m_out %>%
    mutate(group_index = cur_group_id()) %>%
    ungroup() %>%
    arrange(group_index, group)
  
  return(m_out)
  
  }

#### METAFUNCTION 5: PATIENT-LEVEL ANALYSIS: MATCHING ALGORITHM ####

patient.matching_process <- function(pathogen_category,
                                     dat) 
  {
  
  # Cohort identifier
  dat <- dat %>% 
    mutate(group_binary = ifelse(group == "case", 1, 0))
  
  # Drop NAs in relevant matching columns
  dat <- dat %>% 
    drop_na(dept_code, DTS_in_month, DTS_in_year, duration)
  dat <- dat %>% 
    drop_na(age, sex, indiv_score, any_surgery)
  
  print(nrow(dat))
  
  # Match on ward, month and length of stay in first room +/- 3 days or 10%
  matched_out <- match.data(matchit(# group_binary ~ dept_code + DTS_in_month + DTS_in_year + matching_duration,
                                    group_binary ~ dept_code + DTS_in_month + age + sex + indiv_score + any_surgery + matching_duration,
                                    data = dat, 
                                    method = "nearest", 
                                    exact = c("sex", "any_surgery","dept_code", "DTS_in_month"), 
                                    caliper = 1,
                                    ratio = 3))
  
  # print(nrow(matched_out))
  
  m_out <- matched_out %>% 
    group_by(subclass) %>% 
    mutate(age_diff = abs(age[group == "case"] - age)) %>% 
    mutate(elix_diff = abs(indiv_score[group == "case"] - indiv_score)) %>%
    mutate(matching_duration = abs(matching_duration[group == "case"] - matching_duration)) 
  
  # Apply tolerance for some matching
  m_out <- m_out %>% 
    filter(group == "case" | group == "control" & (age_diff <= 5)) 
  
  print("After Matching Age")
  # print(nrow(m_out))
  
  m_out <- m_out%>% 
    filter(group == "case" | group == "control" & (elix_diff <= 3)) 
    
  print("After Matching Elix")
  # print(nrow(m_out))
  
  m_out <- m_out %>%
    filter(group == "case" | group == "control" & (matching_duration <= 3 | matching_duration <= (0.1 * matching_duration[group == "case"]))) 
  
  print("After Matching Duration")
  # print(nrow(m_out))
  
  m_out <- m_out %>% 
    filter(any(group == "case") & any(group == "control")) 
  
  print(nrow(m_out))
  
  # Final clean up
  m_out <- m_out %>%
    mutate(group_index = cur_group_id()) %>%
    ungroup() %>%
    arrange(group_index, group)
  
  return(m_out)
}
