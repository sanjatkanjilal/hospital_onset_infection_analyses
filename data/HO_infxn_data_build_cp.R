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

# Import the functions for the pipeline
source("/PHShome/zw852/ml-hc-class/prior_patient_infection/HO_infxn_functions.R")

# Set working directory and output directory to be whatever you generate today
mainDir <- "/data/tide/projects/ho_infxn_ml/"
today <- 20240911
setwd(file.path(mainDir))

micro_filename <- "micro.ground_truth_20150525-20240701.csv"
adt_filename <- "MGB_ADT_20140131-20240713.csv"
abx_filename <- "MGB_abx_20150323-20240715.csv"
abx_mapping_filename <- "EDW_abx_map_20150525-20240716.csv"
dems_filename <- "MGB_demographics_20240715.csv"
elix_filename <-"MGB_elixhauser_20240715.csv"
cpt_filename <- "MGB_procedures_20240715.csv"
admt_filename <- "MGB_IP_ED_encounters_20150430-20240715.csv"
location_map_filename <- "department_mapping.csv"
enc_filename <- "MGB_encounters_20240715.csv"

adt.micro.raw <- read_csv(file = paste0("clean_data/",  today, "/adt_micro_raw.csv"))
micro <- read_csv(paste0("input_data/", micro_filename)) %>%
  filter(!grepl("strep", org_group_2, ignore.case = T)) %>%
  # mutate(org_group_1 = ifelse(org_group_2 == "Strep_pneumo", "Respiratory flora", org_group_1)) %>%
  # mutate(org_group_3 = ifelse(!is.na(org_group_2) & is.na(org_group_3), org_group_2, org_group_3)) %>%
  # gather(hierarchy, strata, org_group_1:org_group_3) %>%
  filter(!is.na(org_group_1))
room_dat <- read_csv(paste0("/data/tide/data/edw/ADT/", adt_filename)) %>%
  mutate(dept_room = paste(DepartmentDSC, RoomID)) %>%
  mutate(duration = round(duration/24, 1)) %>% # Convert duration scale to days
  select(PatientID:PatientClassDSC, dept_room, DepartmentDSC:duration)

location_map <- read_csv(paste0("mappings/", location_map_filename)) 

cp_features.prelim <- calculate_and_add_cp_scores_parallel(adt.micro.raw, micro, room_dat, location_map) 

cp_features <- cp_features.prelim %>%
  select(PatientID, hospitalization_id, CDiff_cp, MSSA_cp, MRSA_cp, 
         DS_Entero_cp, ESBL_cp, VSE_cp, VRE_cp, DS_PsA_cp, DR_PsA_cp) %>%
  distinct()

pathogen_hierarchy <- "org_group_1"
write_csv(cp_features, file = paste0("clean_data/",  today, "/col_pressure_", pathogen_hierarchy, ".csv"))