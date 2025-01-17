#### HEADER ####
# Description: Script to deidentify the clean dataset for data sharing
# Dependencies: 
#   results/dataset/final_dataset_for_models.csv
# Output:
#   results/dataset/deidentified_final_dataset_for_models.csv
# Authors: Ziming (Alex) Wei
# Last updated: 2023-12-16

##### LIBRARIES ####
library(tidyverse)
library(lubridate)
library(mice)
library(skimr)
#library(DataExplorer)
library(survival)
library(xgboost)
library(caret)
library(pROC)
library(Matrix)
library(conflicted)
library(hash)
conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::lag())
conflicts_prefer(lubridate::month)
conflicts_prefer(lubridate::year)
conflicts_prefer(lubridate::week)
conflicts_prefer(lubridate::quarter)
conflicts_prefer(dplyr::first)
conflicts_prefer(reshape2::melt)
conflicts_prefer(reshape2::dcast)

mainDir <- "~/ho_infxn_clean_code/"
clean_data <- data.table::fread(file = paste0(mainDir,'results/dataset/final_dataset_for_models.csv'))
clean_data <- clean_data %>% filter(match == 'environmental')

# Remove Admission Source
deidentified_data <- clean_data %>% select(-admit_source_clean)

# Map Patient ID to Random Number
unique_patient_id <- unique(deidentified_data$PatientID)
deidentified_patient_id <-sample(1:length(unique_patient_id), length(unique_patient_id), replace = FALSE)
patient_id_dict <- unlist(set_names(deidentified_patient_id, unique_patient_id))
deidentified_data$deidentified_patient_id <- patient_id_dict[deidentified_data$PatientID]

# Map Date to Random Dates
earliest <- deidentified_data %>%
  group_by(deidentified_patient_id) %>%
  summarise(earliest_time = min(DTS_in_month))

deidentified_data <- deidentified_data %>%
  inner_join(earliest, by = c("deidentified_patient_id"))

# Maintain the date relativity for different encounters of the same patient
deidentified_data$date_diff <- difftime(deidentified_data$DTS_in_month, 
                                        deidentified_data$earliest_time, 
                                        units = c("days")) 

unique_deidentified_patient_id = unique(deidentified_data$deidentified_patient_id)

deidentified_date <-sample(seq(as.Date('2000/01/01'), as.Date('2100/01/01'), by="day"), 
                           length(deidentified_patient_id), 
                           replace = TRUE)

# Create a dictionary and map different time to the original dataset
id_date_pair_dict <- unlist(set_names(deidentified_date, unique_deidentified_patient_id))
deidentified_data$deidentified_month <- round_date(id_date_pair_dict[deidentified_data$deidentified_patient_id] + deidentified_data$date_diff, unit=c('month'))
deidentified_data$deidentified_year <- year(deidentified_data$deidentified_month)

deidentified_data <- deidentified_data %>% select(match:group_binary, deidentified_patient_id, deidentified_month, deidentified_year, duration:DR_PsA_cp)

# Save Data
write_csv(deidentified_data,
          file = paste0(mainDir,'results/dataset/deidentified_final_dataset_for_models.csv'))




