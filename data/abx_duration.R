####################################################################################
# Description: Script to pre-process antibiotic data for hHO infection analysis
# Author: Sanjat Kanjilal
# Last updated: 2024-05-02
####################################################################################

#### LIBRARIES ####
library(tidyverse)
library(data.table)
library(reshape2)
library(MatchIt)
library(tictoc)
library(conflicted)
conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::lag())
conflicts_prefer(lubridate::month)
conflicts_prefer(lubridate::year)
conflicts_prefer(lubridate::week)
conflicts_prefer(lubridate::quarter)
conflicts_prefer(dplyr::first)
conflicts_prefer(reshape2::melt)
conflicts_prefer(reshape2::dcast)

#### FUNCTION SET 1: HELPER FUNCTIONS ####

# Function to group events (room stays / drug orders) within a specific time interval (in days) with each other
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

#### FUNCTION SET 2: DATASET PROCESSING ####

# Process raw antibiotic data 
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
    left_join(major %>% select(antibiotic, drug_class), by = c('antibiotic'='antibiotic')) %>%
    arrange(patientID, antibiotic, order_DT, taken_DT) %>%
    select(patientID, order_DT, taken_DT, antibiotic, drug_class, order_freq, order_quantity, route_abbr) %>%
    distinct()
  
  return(abx_prelim)
}

# Pre-process inpatient antibiotic data
process_ip_abx_features <- function(abx_prelim) 
{
  abx.ip <- abx_prelim %>% 
    filter(!is.na(taken_DT)) %>%
    arrange(patientID, drug_class, taken_DT) %>%
    group_by(patientID, taken_DT, drug_class) %>%
    mutate(index = row_number()) %>% 
    ungroup() %>%
    filter(index == 1) %>%
    select(patientID, taken_DT, antibiotic, drug_class) %>%
    distinct()
  
  return(abx.ip)
}

# Calculate inpatient antibiotic duration
calculate_ip_abx_durations <- function(abx.ip)
{
  # Code calculates difference in days and duration
  abx_ip_durations <- abx.ip %>% 
    # head(1000) %>%
    group_by(patientID, drug_class) %>% 
    mutate(day_diff = taken_DT - lag(taken_DT, default = first(taken_DT))) %>% 
    mutate(new_prescription = case_when(day_diff == 0 ~ 1, 
                                        day_diff >= 7 ~ 1)) %>% 
    mutate(new_prescription_1 = cumsum(!is.na(new_prescription))) %>% 
    group_by(patientID, drug_class, new_prescription_1) %>% 
    mutate(abx_duration_days = as.numeric(max(taken_DT) - min(taken_DT) + 1)) %>% 
    slice_tail(n=1) %>% 
    ungroup() %>% 
    rename(end_DT = taken_DT) %>% 
    select(patientID, drug_class, end_DT, abx_duration_days) %>%
    distinct()
  
  return(abx_ip_durations)
}

# Calculate durations for OP antibiotic courses
calculate_op_abx_durations <- function(abx_prelim) 
{
  tic()
  abx_op_durations_raw <- abx_prelim %>% 
    filter(is.na(taken_DT)) %>% 
    # head(1000) %>%
    arrange(patientID, drug_class, order_DT) %>%
    group_by(patientID, order_DT, drug_class) %>%
    mutate(index = row_number()) %>% 
    ungroup() %>%
    filter(index == 1) %>%
    select(patientID, order_DT, antibiotic, drug_class:order_quantity) %>%
    distinct() %>%
    mutate(numeric_frequency_per_day = case_when(
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
    mutate(numeric_frequency_per_day = round(numeric_frequency_per_day, 1)) %>% 
    filter(numeric_frequency_per_day > 0) %>% 
    mutate(quantity_numeric = case_when(
      str_detect(tolower(order_quantity), "(?i)\\d+\\s+tablet") ~ as.numeric(str_match(tolower(order_quantity), "(?i)(\\d+)\\s+tablet")[,2]),
      str_detect(tolower(order_quantity), "(?i)\\d+\\s+capsule") ~ as.numeric(str_match(tolower(order_quantity), "(?i)(\\d+)\\s+capsule")[,2]),
      str_detect(tolower(order_quantity), "(?i)\\d+\\s+packet") ~ as.numeric(str_match(tolower(order_quantity), "(?i)(\\d+)\\s+packet")[,2]),
      str_detect(tolower(order_quantity), "(?i)\\d+\\s+dose") ~ as.numeric(str_match(tolower(order_quantity), "(?i)(\\d+)\\s+dose")[,2])
    )) %>% 
    mutate(quantity_numeric = ifelse(is.na(quantity_numeric), 5, quantity_numeric)) %>%
    mutate(abx_duration_days = round(quantity_numeric / numeric_frequency_per_day, 1)) %>% 
    filter(!is.na(abx_duration_days)) %>%
    mutate(end_DT = order_DT + abx_duration_days) %>%
    select(patientID, drug_class, order_DT, end_DT, abx_duration_days) %>%
    distinct()
  toc()
  
  return(abx_op_durations_raw)
}

# Group outpatient ourses within 7 days of each other and calculate last taken date
calculate_op_abx_last_date <- function(abx_op_durations_raw) 
{
  tic()
  print("Removing duplicate orders")
  # Remove duplicate courses and re-index
  abx_op_durations.prelim <- abx_op_durations_raw %>% 
    # head(10000) %>%
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
  
  toc()
  
  abx_op_episodes <- grouping_function(data = abx_op_durations.prelim, id = key, index = course_index, 
                                       start = order_DT, end = end_DT, time_interval = 7) 
  
  tic()
  print("Calculating courses after grouping")
  abx_op_episodes <- abx_op_episodes %>%
    group_by(key, group) %>%
    mutate(start_new = min(order_DT),
           end_new = max(end_DT)) %>%
    select(-course_index, -order_DT, -end_DT) %>%
    distinct() %>%
    mutate(abx_duration_days_new = end_new - start_new) 
  toc()
  
  tic()
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
  abx_courses <- bind_rows(abx_ip_durations, abx_op_durations) %>%
    arrange(patientID, drug_class, end_DT)
  
  return(abx_courses)
}

#### EXECUTE SCRIPT ####

# Global paths
adt_path <- "/data/clinicalml/data/EDW/ADT/MGB_ADT_20140131-20231126.csv"
abx_path <- "/data/clinicalml/data/EDW/all_MGB/MGB_abx_20150323-20231201.csv"
abx_mapping_path <- "/data/clinicalml/data/EDW/Auxiliary_docs/EDW_abx_map_20150525- 20230125.csv"
output_path <- "/data/clinicalml/skanjilal/Projects/prior_occupant_analysis/"

# Import room data to filter for PatientIDs
patientIDs <- read_csv(adt_path) %>%
  select(PatientID) %>%
  distinct() %>%
  arrange(PatientID)

# Antibiotic data
abx <- read_csv(abx_path) %>%
  mutate(across(where(is.character), ~ na_if(.,""))) %>%
  filter(patientID %in% patientIDs$PatientID)

# Flagging major antibiotic classes
major <- read_csv(abx_mapping_path) %>% 
  filter(major=="x")

# Pre-processing to drop prns, 'ordered but not given', multiple orders given on a same day, then group by regimen dates
abx_prelim <- initial_abx_clean(abx, major)

start <- gsub("-", "", min(abx_prelim$order_DT))
end <- gsub("-", "", max(abx_prelim$order_DT))

# Generate inpatient abx durations by abx class 
abx.ip <- process_ip_abx_features(abx_prelim) 
abx_ip_durations <- calculate_ip_abx_durations(abx.ip) 

write_csv(abx_ip_durations, file = paste0(output_path, "abx_ip_durations_", start, "-", end, ".csv"))

# Generate outpatient abx durations by abx class 
abx_op_durations_raw <- calculate_op_abx_durations(abx_prelim)
abx_op_durations <- calculate_op_abx_last_date(abx_op_durations_raw)

write_csv(abx_ip_durations, file = paste0(output_path, "abx_ip_durations_", start, "-", end, ".csv"))

# Combine IP and OP abx duration
abx_courses <-  abx_courses(abx_ip_durations, abx_op_durations)

write_csv(abx_courses, file = paste0(output_path, "abx_courses_", start, "-", end, ".csv"))

#### CLEAN UP ####

rm(major, abx, abx_prelim, abx.ip, abx_ip_durations, abx_op_durations_raw, abx_op_durations)

