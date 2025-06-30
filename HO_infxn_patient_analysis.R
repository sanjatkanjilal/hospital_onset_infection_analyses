# /usr/bin/R

#####################################################################################################
# Description: Data visualizations for colonization pressure / hospital-onset infection analysis
#
# DEPENDENT SCRIPTS
#
# HO_infxn_functions.R - Contains all functions in data processing scripts
# HO_infxn_micro_prep.R - Bins micro results into organism categories
# HO_infxn_input_data_processing.R - Pre-processes micro, antibiotic, encounter, and admit source datasets
# HO_infxn_build_unmatched_cohorts.R - Generates cohorts without features
# HO_infxn_add_features.R - Adds features to unmatched cohorts
# HO_infxn_build_matched_cohorts.R - Matches cases to controls for final analysis
# HO_infxn_models.R - Runs CLR and XGB models
#
# INPUT DATASETS
#
# Matched cohorts with features
# CLR and XGB model outputs
#
# OUTPUTs
#
# Tables (baseline characteristics)
# Figures (CLR results; CP distributions; SHAP values from XGB models)
#
# Authors: Ziming (Alex) Wei, Sanjat Kanjilal
# Last updated: 2025-04-11
#####################################################################################################

##### LIBRARIES ####
library(tidyverse)
library(conflicted)
library(SHAPforxgboost)
library(xgboost)
library(dplyr)
library(scales)
conflicts_prefer(dplyr::filter())
conflicts_prefer(base::`%in%`)

#### ENVIRONMENTAL VARIABLES ####

# Set working directory and output directory to be whatever you generate today
mainDir <- "/data/tide/projects/ho_infxn_ml/"
setwd(file.path(mainDir))

#### IMPORT DATASETS ####
# cc_final <- readr::read_csv("clean_data/20250411/final_dataset_for_models_20250411.csv")
cc_final <- read.csv('/data/tide/projects/ho_infxn_ml/clean_data/20250411/final_dataset_for_models_elix_20250411.csv')
clr.results <- readr::read_csv("results/model_results/20250411/clogit_coefficients.csv")

# Filter for result of colonization pressure analysis and set factor leve.s
clean.data.table <- cc_final %>%
  dplyr::filter(match == "patient") %>%
  dplyr::filter(run %in% c("DS_E_coli", "ESBL_E_coli", "DS_K_pneumoniae", "ESBL_K_pneumoniae", 
                           "C_diff", "VSE_faecalis", "VRE_faecium", 
                           "MSSA", "MRSA", "DS_P_aeruginosa", "DR_P_aeruginosa")) %>%
  dplyr::mutate(run = factor(run, levels = c("DS_E_coli", "ESBL_E_coli", "DS_K_pneumoniae", "ESBL_K_pneumoniae", 
                                             "C_diff", "VSE_faecalis", "VRE_faecium", 
                                             "MSSA", "MRSA", "DS_P_aeruginosa", "DR_P_aeruginosa")))


#### TABLE: BASELINE CHARACTERISTICS ####

# Sample size
sample_size <- clean.data.table %>%
  dplyr::group_by(run, group) %>%
  dplyr::summarise(sample_size = dplyr::n())
sample_size_pooled <- clean.data.table %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(sample_size = dplyr::n())

# Age
age <- clean.data.table %>%
  dplyr::group_by(run, group) %>%
  dplyr::summarise(n = dplyr::n(),
                   age_mean = mean(age),
                   age_sd = sd(age),
                   age_se = age_sd / sqrt(n))

age_pooled <- clean.data.table %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n(),
                   age_mean = mean(age),
                   age_sd = sd(age),
                   age_se = age_sd / sqrt(n))

# Elixhauser index
elix <- clean.data.table %>%
  dplyr::group_by(run, group) %>%
  dplyr::summarise(n = dplyr::n(),
                   elix_mean = mean(elix_index_mortality),
                   elix_sd = sd(elix_index_mortality),
                   elix_se = elix_sd / sqrt(n),
                   elix_median = median(elix_index_mortality),
                   elix_q1 = quantile(elix_index_mortality, 0.25),
                   elix_q3 = quantile(elix_index_mortality, 0.75))

elix_pooled <- clean.data.table %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n(),
                   elix_mean = mean(elix_index_mortality),
                   elix_sd = sd(elix_index_mortality),
                   elix_se = elix_sd / sqrt(n),
                   elix_median = median(elix_index_mortality),
                   elix_q1 = quantile(elix_index_mortality, 0.25),
                   elix_q3 = quantile(elix_index_mortality, 0.75))

elix_indiv <- clean.data.table %>%
  dplyr::group_by(run, group) %>%
  dplyr::summarise(across(starts_with("elix"), ~ mean(. != 0, na.rm = TRUE))) %>% 
  select(-elix_index_mortality)

elix_indiv_pooled <- clean.data.table %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(across(starts_with("elix"), ~ mean(. != 0, na.rm = TRUE))) %>% 
  select(-elix_index_mortality)

write.csv(elix_indiv, "results/model_results/20250411/patient_match/elix_table_1.csv")
write.csv(elix_indiv, "results/model_results/20250411/patient_match/pooled_elix_table_1.csv")


# Time to infection (cases only)
time.to.infxn <- clean.data.table %>%
  dplyr::filter(group == "case") %>%
  dplyr::group_by(run, group) %>%
  dplyr::summarise(n = dplyr::n(),
                   time_to_infxn_mean = mean(time_to_infxn),
                   time_to_infxn_sd = sd(time_to_infxn),
                   time_to_infxn_se = time_to_infxn_sd / sqrt(n),
                   time_to_infxn_median = median(time_to_infxn),
                   time_to_infxn_q1 = quantile(time_to_infxn, 0.25),
                   time_to_infxn_q3 = quantile(time_to_infxn, 0.75))

time.to.infxn_pooled <- clean.data.table %>%
  dplyr::filter(group == "case") %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n(),
                   time_to_infxn_mean = mean(time_to_infxn),
                   time_to_infxn_sd = sd(time_to_infxn),
                   time_to_infxn_se = time_to_infxn_sd / sqrt(n),
                   time_to_infxn_median = median(time_to_infxn),
                   time_to_infxn_q1 = quantile(time_to_infxn, 0.25),
                   time_to_infxn_q3 = quantile(time_to_infxn, 0.75))


# Length of stay (LOS) in 1st room
LOS <- clean.data.table %>%
  dplyr::group_by(run, group) %>%
  dplyr::summarise(n = dplyr::n(),
                   LOS_mean = mean(duration),
                   LOS_sd = sd(duration),
                   LOS_se = LOS_sd / sqrt(n),
                   LOS_median = median(duration),
                   LOS_q1 = quantile(duration, 0.25),
                   LOS_q3 = quantile(duration, 0.75))

LOS_pooled <- clean.data.table %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n(),
                   LOS_mean = mean(duration),
                   LOS_sd = sd(duration),
                   LOS_se = LOS_sd / sqrt(n),
                   LOS_median = median(duration),
                   LOS_q1 = quantile(duration, 0.25),
                   LOS_q3 = quantile(duration, 0.75))

# Matching duration
matching.duration <- clean.data.table %>%
  dplyr::group_by(run, group) %>%
  dplyr::summarise(n = dplyr::n(),
                   matching_duration_mean = mean(matching_duration),
                   matching_duration_sd = sd(matching_duration),
                   matching_duration_se = matching_duration_sd / sqrt(n),
                   matching_duration_median = median(matching_duration),
                   matching_duration_q1 = quantile(matching_duration, 0.25),
                   matching_duration_q3 = quantile(matching_duration, 0.75))

matching.duration_pooled <- clean.data.table %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = dplyr::n(),
                   matching_duration_mean = mean(matching_duration),
                   matching_duration_sd = sd(matching_duration),
                   matching_duration_se = matching_duration_sd / sqrt(n),
                   matching_duration_median = median(matching_duration),
                   matching_duration_q1 = quantile(matching_duration, 0.25),
                   matching_duration_q3 = quantile(matching_duration, 0.75))

# Categorical variables

# Sex
sex <- clean.data.table %>%
  dplyr::select(run, group, group_index, PatientID, sex) %>%
  dplyr::distinct() %>%
  dplyr::count(run, group, sex) %>%
  dplyr::group_by(run, group) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(female = (n / total)*100) %>%
  dplyr::filter(sex != "Male")

sex_pooled <- clean.data.table %>%
  dplyr::select(run, group, group_index, PatientID, sex) %>%
  dplyr::distinct() %>%
  dplyr::count(group, sex) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(female = (n / total)*100) %>%
  dplyr::filter(sex != "Male")

# Prior surgery
surgery <- clean.data.table %>%
  dplyr::select(run, group, group_index, PatientID, any_surgery) %>%
  dplyr::distinct() %>%
  dplyr::count(run, group, any_surgery) %>%
  dplyr::group_by(run, group) %>%
  dplyr::mutate(total =sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent_prior_surgery = (n / total)*100) %>%
  dplyr::filter(any_surgery != 0)

surgery_pooled <- clean.data.table %>%
  dplyr::select(run, group, group_index, PatientID, any_surgery) %>%
  dplyr::distinct() %>%
  dplyr::count(group, any_surgery) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(total =sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent_prior_surgery = (n / total)*100) %>%
  dplyr::filter(any_surgery != 0)

# Prior abx exposure
abx <- clean.data.table %>%
  dplyr::select(run, group, group_index, PatientID, anti_anaerobe_0_60:tetracycline_0_60) %>%
  dplyr::distinct() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total_courses_0_60 = sum(dplyr::c_across(contains("_0_60")), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(anti_anaerobe_0_60:total_courses_0_60, names_to = "abx", values_to = "courses") %>%
  dplyr::filter(courses > 0) %>%
  dplyr::count(run, group, PatientID, abx) %>%
  dplyr::group_by(run, group, abx) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::select(run, group, abx, total) %>%
  dplyr::distinct() %>%
  dplyr::left_join(sample_size) %>%
  dplyr::mutate(courses_group = (total / (sample_size/100))) %>%
  dplyr::mutate(abx = factor(abx, levels = c("total_courses_0_60", "penicillin_0_60", "anti_staph_beta_lactam_0_60", "extended_spectrum_penicillin_0_60", 
                                             "cephalosporin_0_60", "extended_spectrum_cephalosporin_0_60", "carbapenem_0_60",
                                             "glycopeptide_0_60", "fluoroquinolone_0_60", "macrolide_0_60", "lincosamide_0_60", 
                                             "tetracycline_0_60", "sulfonamide_0_60", "anti_anaerobe_0_60", "anti_Cdiff_0_60")))

abx_se <- clean.data.table %>%
  dplyr::select(run, group, group_index, PatientID, anti_anaerobe_0_60:tetracycline_0_60) %>%
  dplyr::distinct() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total_courses_0_60 = sum(dplyr::c_across(contains("_0_60")), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(anti_anaerobe_0_60:total_courses_0_60, names_to = "abx", values_to = "courses") %>%
  dplyr::group_by(run, group, abx) %>%
  dplyr::summarize(
    se_courses = sd(courses, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(se_courses = se_courses * 100)

abx <- dplyr::left_join(abx,
                        abx_se,
                        by = c('run','group','abx'))


abx_pooled <- clean.data.table %>%
  dplyr::select(run, group, group_index, PatientID, anti_anaerobe_0_60:tetracycline_0_60) %>%
  dplyr::distinct() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total_courses_0_60 = sum(dplyr::c_across(contains("_0_60")), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(anti_anaerobe_0_60:total_courses_0_60, names_to = "abx", values_to = "courses") %>%
  dplyr::filter(courses > 0) %>%
  dplyr::count(run, group, PatientID, abx) %>%
  dplyr::group_by(group, abx) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::select(group, abx, total) %>%
  dplyr::distinct() %>%
  dplyr::left_join(sample_size_pooled) %>%
  dplyr::mutate(courses_group = (total / (sample_size/100))) %>%
  dplyr::mutate(abx = factor(abx, levels = c("total_courses_0_60", "penicillin_0_60", "anti_staph_beta_lactam_0_60", "extended_spectrum_penicillin_0_60", 
                                             "cephalosporin_0_60", "extended_spectrum_cephalosporin_0_60", "carbapenem_0_60",
                                             "glycopeptide_0_60", "fluoroquinolone_0_60", "macrolide_0_60", "lincosamide_0_60", 
                                             "tetracycline_0_60", "sulfonamide_0_60", "anti_anaerobe_0_60", "anti_Cdiff_0_60")))

abx_pooled_se <- clean.data.table %>%
  dplyr::select(run, group, group_index, PatientID, anti_anaerobe_0_60:tetracycline_0_60) %>%
  dplyr::distinct() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total_courses_0_60 = sum(dplyr::c_across(contains("_0_60")), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(anti_anaerobe_0_60:total_courses_0_60, names_to = "abx", values_to = "courses") %>%
  dplyr::group_by(group, abx) %>%
  dplyr::summarize(
    se_courses = sd(courses, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(se_courses = se_courses * 100)

abx_pooled <- dplyr::left_join(abx_pooled,
                               abx_se,
                               by = c('group','abx'))

# Combine into a single table
age.temp <- age %>% 
  dplyr::select(-n) %>% 
  dplyr::mutate(across(starts_with("age"), ~round(., 1)))

sex.temp <- sex %>% 
  dplyr::select(run, group, perc_female = female) %>% 
  dplyr::mutate(perc_female = round(perc_female, 1))

surgery.temp <- surgery %>% 
  dplyr::select(run, group, perc_prior_surgery = percent_prior_surgery) %>% 
  dplyr::mutate(perc_prior_surgery = round(perc_prior_surgery, 1))

elix.temp <- elix %>% 
  dplyr::select(-n) %>% 
  dplyr::mutate(across(starts_with("elix"), ~round(., 1)))

LOS.temp <- LOS %>% 
  dplyr::select(-n) %>% 
  dplyr::mutate(across(starts_with("LOS"), ~round(., 1)))

matching.duration.temp <- matching.duration %>% 
  dplyr::select(-n) %>% 
  dplyr::mutate(across(starts_with("matching_duration"), ~round(., 1)))

time.to.infxn.temp <- time.to.infxn %>% 
  dplyr::select(-n) %>% 
  dplyr::mutate(across(starts_with("time_to_infxn"), ~round(., 1)))

abx.temp <- abx %>%
  dplyr::select(run, group, abx, courses_group) %>%
  tidyr::pivot_wider(names_from = abx, values_from = courses_group) %>%
  dplyr::mutate(across(ends_with("0_60"), ~tidyr::replace_na(., 0))) %>%
  dplyr::mutate(across(ends_with("0_60"), ~round(., 1)))

table1 <- sample_size %>%
  dplyr::left_join(age.temp) %>%
  dplyr::left_join(sex.temp) %>%
  dplyr::left_join(surgery.temp) %>%
  dplyr::left_join(elix.temp) %>%
  dplyr::left_join(LOS.temp) %>%
  dplyr::left_join(matching.duration.temp) %>%
  dplyr::left_join(time.to.infxn.temp) %>%
  dplyr::left_join(abx.temp) 

# Save table 1
readr::write_csv(table1, file = "results/model_results/20250411/patient_match/table1.csv")

rm(age.temp, sex.temp, surgery.temp, elix.temp, LOS.temp, 
   matching.duration.temp, time.to.infxn.temp, abx.temp)

#### FIGURE: BASELINE CHARACTERISTICS ####
sample_size %>%
  ggplot(aes(x = run, y = sample_size, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  labs(title = "Sample size", x = "Pathogen", y = "Sample size")

age %>% 
  ggplot(aes(x = run, y = age_mean, fill = group)) +
  geom_point(aes(color = group), 
             position = position_dodge(width = 0.9), 
             size = 3) +
  geom_errorbar(
    aes(ymin = age_mean - age_se, ymax = age_mean + age_se, color = group),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  # ylim(55, 75) +
  labs(title = "Age", x = "Pathogen")

LOS %>% 
  ggplot(aes(x = run, y = LOS_mean, fill = group)) +
  geom_point(aes(color = group), 
             position = position_dodge(width = 0.9), 
             size = 3) +
  geom_errorbar(
    aes(ymin = LOS_mean - LOS_se, ymax = LOS_mean + LOS_se, color = group),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  # ylim(55, 75) +
  labs(title = "LOS", x = "Pathogen")

matching.duration %>% 
  ggplot(aes(x = run, y = matching_duration_mean, fill = group)) +
  geom_point(aes(color = group), 
             position = position_dodge(width = 0.9), 
             size = 3) +
  geom_errorbar(
    aes(ymin = matching_duration_mean - matching_duration_se, ymax = matching_duration_mean + matching_duration_se, color = group),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  # ylim(55, 75) +
  labs(title = "matching_duration", x = "Pathogen")

time.to.infxn %>% 
  ggplot(aes(x = run, y = time_to_infxn_mean, fill = group)) +
  geom_point(aes(color = group), 
             position = position_dodge(width = 0.9), 
             size = 3) +
  geom_errorbar(
    aes(ymin = time_to_infxn_mean - time_to_infxn_se, ymax = time_to_infxn_mean + time_to_infxn_se, color = group),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  # ylim(55, 75) +
  labs(title = "time_to_infxn", x = "Pathogen")

elix %>% 
  ggplot(aes(x = run, y = elix_mean)) +
  geom_point(aes(color = group), 
             position = position_dodge(width = 0.9), size = 3) +
  geom_errorbar(
    aes(ymin = elix_mean - elix_se, ymax = elix_mean + elix_se, color = group),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  labs(title = "Elixhauser index", x = "Pathogen")

elix %>% 
  ggplot(aes(x = run, y = elix_median)) +
  geom_point(aes(color = group), 
             position = position_dodge(width = 0.9), size = 3) +
  geom_errorbar(
    aes(ymin = elix_q1, ymax = elix_q3, color = group),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  labs(title = "Elixhauser index", x = "Pathogen")

sex %>%
  ggplot(aes(x = run, y = female)) +
  geom_point(aes(color = group), position = position_dodge(width = 0.9)) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  ylim(0, 100) +
  labs(title = "Sex (% female)", x = "Pathogen")

surgery %>%
  ggplot(aes(x = run, y = percent_prior_surgery)) +
  geom_point(aes(color = group), position = position_dodge(width = 0.9)) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  ylim(0, 100) +
  labs(title = "Prior surgery", x = "Pathogen")

# Combine sex and prior surgery into single plot
surgery %>%
  dplyr::select(run, group, percent_prior_surgery) %>% 
  dplyr::distinct() %>%
  dplyr::left_join(sex %>% dplyr::select(run, group, female) %>% dplyr::distinct()) %>%
  tidyr::pivot_longer(percent_prior_surgery:female, names_to = "var", values_to = "percent") %>%
  ggplot(aes(x = run, y = percent, shape = var)) +
  geom_point(aes(color = group, 
                 shape = var,
                 size = 3), 
             position = position_dodge(width = 0.9)) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  ylim(0, 100) +
  labs(title = "Sex and prior surgery", x = "Pathogen")

# Plot total abx courses
abx %>%
  dplyr::filter(abx == "total_courses_0_60") %>%
  ggplot(aes(x = run, y = courses_group, fill = group)) +
  geom_point(aes(color = group,
                 size = 3), 
             position = position_dodge(width = 0.9)) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  labs(title = "Number of antibiotic courses", 
       y = "Courses / 100 people", x = "Pathogen")

# Plot Individual antibiotic classes
abx %>%
  dplyr::filter(abx != "total_courses_0_60") %>%
  ggplot(aes(x = run, y = courses_group, 
             fill = group,
             color = group)) +
  geom_point(position = position_dodge(width = 0.9),
  ) +
  scale_x_discrete(drop=FALSE) + 
  facet_wrap(vars(abx), nrow = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    panel.background = element_blank()) +
  labs(
    title = "Antibiotic courses by class",
    x = "Pathogen"
  )

# Clean up
rm(sample_size, age, time.to.infxn, LOS, matching.duration, elix, sex, surgery, abx, table1)

#### FIGURE: HOSPITAL-ONSET INFECTION AND COLONIZATION PRESSURE ####

# Data prep
orgs_to_keep <- c("C_diff", "DS_E_coli", "DS_K_pneumoniae", "ESBL_E_coli", "ESBL_K_pneumoniae",
                  "MRSA", "MSSA", "VRE_faecium", "VSE_faecium", "VSE_faecalis",
                  "DR_P_aeruginosa", "DS_P_aeruginosa")

# Prep CLR for plotting
clr.prepped <- clr.results %>%
  dplyr::filter(target %in% orgs_to_keep) %>%
  dplyr::filter(match == 'patient') %>%
  dplyr::select(match:variable, coef_raw, lower_CI_raw, upper_CI_raw, SE_coef_raw, pval_raw, sig_flag) %>%
  dplyr::distinct() %>%
  dplyr::mutate(perc_coef_rounded = round((coef_raw - 1)*100, 1),
                perc_LCI_rounded = round((lower_CI_raw - 1)*100, 1),
                perc_UCI_rounded = round((upper_CI_raw - 1)*100, 1)) %>%
  dplyr::mutate(direction = ifelse(coef_raw > 1, "pos", "neg")) %>%
  dplyr::mutate(org_group = dplyr::case_when(
    target %in% c("C_diff", "DS_E_cloacae", "DS_E_coli", "DS_K_oxytoca", "DS_K_pneumoniae", 
                  "DS_P_mirabilis", "DS_S_marcescens", "ESBL_E_cloacae", "ESBL_E_coli",
                  "ESBL_K_pneumoniae", "ESBL_P_mirabilis", "VRE_faecium", "VSE_faecalis",
                  "VSE_faecium") ~ "enterics",
    target %in% c("MRSA", "MSSA") ~ "skin_flora",
    target %in% c("DR_P_aeruginosa", "DS_P_aeruginosa") ~ "environmentals"
  )) %>%
  dplyr::mutate(target = factor(target, levels = c("DS_E_coli", "ESBL_E_coli", "DS_K_pneumoniae", "ESBL_K_pneumoniae", 
                                                   "DS_E_cloacae", "ESBL_E_cloacae", "C_diff", "VSE_faecalis", "VRE_faecium", 
                                                   "VSE_faecium", "MRSA", "MSSA", "DR_P_aeruginosa", "DS_P_aeruginosa"
  ))) %>%
  dplyr::mutate(variable = factor(variable, levels = c(
    "age",
    "sexMale",
    "elix_index_mortality",
    "any_surgery",
    "anti_anaerobe_0_60",
    "anti_Cdiff_0_60",
    "anti_staph_beta_lactam_0_60",
    "carbapenem_0_60",
    "cephalosporin_0_60",
    "extended_spectrum_cephalosporin_0_60",
    "extended_spectrum_penicillin_0_60",
    "fluoroquinolone_0_60",
    "glycopeptide_0_60",
    "lincosamide_0_60",
    "macrolide_0_60",
    "penicillin_0_60",
    "sulfonamide_0_60",
    "tetracycline_0_60"
  )))


# Function to filter CLR data for target and CP sets
target.abx.data <- function(data, target_set, comparison_groups) 
{
  
  filtered.data <- data %>%
    
    # Filter for target set
    dplyr::filter(target %in% target_set) %>%
    dplyr::mutate(target = factor(target, levels = target_set)) %>%
    
    # Filter for CP set
    dplyr::mutate(comparison = comparison_groups) %>%
    
    dplyr::arrange(variable, target)
  
  return(filtered.data)
  
}

# Define target sets
enteric.target <- c("DS_E_coli", "ESBL_E_coli", "DS_K_pneumoniae", "ESBL_K_pneumoniae", "C_diff", "VSE_faecalis", "VRE_faecium")
skin.target <- c("MSSA", "MRSA")
environmental.target <-c("DS_P_aeruginosa", "DR_P_aeruginosa")

# Build cognate datasets
enteric_dat <- target.cp.data(clr.prepped, enteric.target, comparison_groups = "enteric")
skin_dat <- target.cp.data(clr.prepped, skin.target, comparison_groups = "skin")
environmental_dat <- target.cp.data(clr.prepped, environmental.target, comparison_groups = "environmental")

# Combine all datasets
clr.cp.results <- dplyr::bind_rows(enteric_dat, skin_dat, environmental_dat) %>%
  dplyr::select(match:variable, comparison, coef = perc_coef_rounded, lower_CI = perc_LCI_rounded, upper_CI = perc_UCI_rounded,
                pval_raw:sig_flag)

# Clean up
rm(orgs_to_keep, target.cp.data, enteric.cp, enteric.target, skin.cp, environmental.cp, 
   skin.target, environmental.target, cognate.enteric, cognate.skin, cognate.environmental, 
   noncognate.enteric.skin, noncognate.enteric.env, noncognate.skin.enteric, 
   noncognate.skin.env, noncognate.env.enteric, noncognate.env.skin)

# Function to plot CLR data for target and CP sets
target.abx.plot <- function(filtered.data, comparison_groups) 
{
  
  filtered.plot <- filtered.data %>%
    dplyr::filter(comparison == comparison_groups) %>%
    dplyr::filter(sig_flag=='*') %>%
    dplyr::mutate(across(where(is.factor), fct_drop)) %>%
    ggplot(aes(x = variable, y = coef, ymin=lower_CI, ymax=upper_CI,
               middle = coef, lower = lower_CI, upper = upper_CI)) +
    geom_boxplot(stat = "identity", aes(fill = target)) +
    geom_hline(yintercept = 0, color = "red", linetype='dotted') +
    scale_x_discrete(drop=FALSE) +
    theme_minimal() +
    theme(panel.border = element_rect(linewidth = 1, fill = NA),
          panel.background = element_blank()) +
    labs(x = "Prior Abx", y = "%_change_in_odds")
  
  return(filtered.plot)
  
}

# Plot cognate datasets
enteric.plot <- target.abx.plot(clr.cp.results, comparison_groups = "enteric")
skin.plot <- target.abx.plot(clr.cp.results, comparison_groups = "skin")
environmental.plot <- target.abx.plot(clr.cp.results, comparison_groups = "environmental")

# Plot CLR Coef Heatmap
ggplot(clr.cp.results, aes(x = target, y = variable, fill = coef)) +
  geom_tile(color = "white") +                     # white grid lines
  scale_fill_gradient2(                            # diverging palette centered at 0
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0, name = "% Change in Odds"
  ) +
  labs(x = "Target", y = "Colonization Pressure") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

# Clean up
rm(clr.results, clr.cp.results, clr.prepped, target.cp.plot, cognate.enteric.plot, cognate.skin.plot, cognate.environmental.plot,
   noncognate.enteric.skin.plot, noncognate.enteric.env.plot, 
   noncognate.skin.enteric.plot, noncognate.skin.env.plot, 
   noncognate.env.enteric.plot, noncognate.env.skin.plot)

#### FIGURE: DISTRIBUTION OF COLONIZATION PRESSURE ####

cp <- clean.data.table %>%
  dplyr::select(run, group, group_index, DS_Entero_cp, ESBL_cp, VSE_cp, VRE_cp, CDiff_cp, MSSA_cp, MRSA_cp, DS_PsA_cp, DR_PsA_cp) %>%
  dplyr::distinct() %>%
  tidyr::pivot_longer(DS_Entero_cp:DR_PsA_cp, names_to = "cp_type", values_to = "cp_val") %>%
  dplyr::mutate(cp_type = factor(cp_type, levels = c("DS_Entero_cp", "ESBL_cp", "VSE_cp", 
                                                     "VRE_cp", "CDiff_cp", "MSSA_cp", "MRSA_cp", "DS_PsA_cp",
                                                     "DR_PsA_cp"))) %>%
  dplyr::mutate(percentile = percent_rank(cp_val))

# Compare CP between drug-susceptible and drug-resistant organism pairs (t-tests)

# Drug-susceptible and ESBL Enterobacterales
enterobacterales.cp <- cp %>% 
  dplyr::filter(cp_type %in% c("DS_Entero_cp")) %>%
  dplyr::select(cp_type, cp_val) %>%
  dplyr::distinct() 

esbl.cp <- cp %>% 
  dplyr::filter(cp_type %in% c("ESBL_cp")) %>%
  dplyr::select(cp_type, cp_val) %>%
  dplyr::distinct()

t.test(enterobacterales.cp$cp_val, esbl.cp$cp_val)

# VSE and VRE
vse.cp <- cp %>% 
  dplyr::filter(cp_type %in% c("VSE_cp")) %>%
  dplyr::select(cp_type, cp_val) %>%
  dplyr::distinct() 

vre.cp <- cp %>% 
  dplyr::filter(cp_type %in% c("VRE_cp")) %>%
  dplyr::select(cp_type, cp_val) %>%
  dplyr::distinct()

t.test(vse.cp$cp_val, vre.cp$cp_val)

# MSSA annd MRSA
mssa.cp <- cp %>% 
  dplyr::filter(cp_type %in% c("MSSA_cp")) %>%
  dplyr::select(cp_type, cp_val) %>%
  dplyr::distinct() 

mrsa.cp <- cp %>% 
  dplyr::filter(cp_type %in% c("MRSA_cp")) %>%
  dplyr::select(cp_type, cp_val) %>%
  dplyr::distinct()

t.test(mssa.cp$cp_val, mrsa.cp$cp_val)

# Drug-susceptible and drug-resistant P. aeruginosa
psa.ds.cp <- cp %>% 
  dplyr::filter(cp_type %in% c("DS_PsA_cp")) %>%
  dplyr::select(cp_type, cp_val) %>%
  dplyr::distinct() 

psa.dr.cp <- cp %>% 
  dplyr::filter(cp_type %in% c("DR_PsA_cp")) %>%
  dplyr::select(cp_type, cp_val) %>%
  dplyr::distinct()

t.test(psa.ds.cp$cp_val, psa.dr.cp$cp_val)

# Get percentiles
cp_percentile <- cp %>%
  dplyr::group_by(cp_type) %>%
  dplyr::mutate(overall_mean = mean(cp_val), 
                overall_median = median(cp_val),
                overall_se = (sd(cp_val) / dplyr::n())) %>%
  dplyr::mutate(mean_percentile = mean(cp_val <= overall_mean),
                median_percentile = 0.5) %>%
  dplyr::mutate(org = dplyr::case_when(
    cp_type %in% c("DS_Entero_cp", "ESBL_cp") ~ "Enterobacterales",
    cp_type %in% c("VSE_cp", "VRE_cp") ~ "Enterococcus",
    cp_type %in% c("CDiff_cp") ~ "C_difficile",
    cp_type %in% c("MSSA_cp", "MRSA_cp") ~ "S_aureus",
    cp_type %in% c("DS_PsA_cp", "DR_PsA_cp") ~ "P_aeruginosa")) %>%
  dplyr::select(org, cp_type, overall_mean, mean_percentile, overall_median, median_percentile) %>%
  dplyr::distinct()

# Join back to main dataset
cp.final <- cp %>%
  dplyr::left_join(cp_percentile)

# Plot distribution of CP by organism pairs
cp.final %>%
  # dplyr::filter(cp_type %in% c("DS_Entero_cp", "ESBL_cp")) %>%
  # dplyr::filter(cp_type %in% c("DS_Entero_cp", "ESBL_cp", "VSE_cp", "VRE_cp", "CDiff_cp")) %>%
  # dplyr::filter(cp_type %in% c("MSSA_cp", "MRSA_cp")) %>%
  # dplyr::filter(cp_type %in% c("DS_PsA_cp", "DR_PsA_cp")) %>%
  ggplot(aes(x = cp_val, color = cp_type)) +
  geom_point(aes(x = overall_median, y = median_percentile)) +
  geom_segment(aes(x = overall_median, y = 0, xend = overall_median, yend = 0.5)) +
  stat_ecdf() + 
  facet_grid(org ~ .) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  labs(x = "Colonization pressure value", y = "Percentiile")

# Clean up
rm(cp, enterobacterales.cp, esbl.cp, vse.cp, vre.cp, mssa.cp,
   mrsa.cp, psa.ds.cp, psa.dr.cp, cp_percentile, cp.final)

#### SHAP value plots ####
all_runs = unique(cc_final$run)

cp_colnames = cc_final %>% select(elix_index_mortality, CDiff_cp:DR_PsA_cp) %>% colnames()

shap_value_matrix = data.frame()
elix_shap_matrix = data.frame()
for (run_name in all_runs){
  for (fold in c(1,2,3,4,5)){
    model = xgb.load(paste0('results/model_results/20250411/patient_match/xgb/model_checkpoints/environmental_', run_name ,'/fold_',fold,'.model'))
    
    dat.run <- cc_final %>% filter(run == run_name & match == 'environmental')
    dat.run <- dat.run %>% select(elix_index_mortality, CDiff_cp:DR_PsA_cp)
    dat.run <- dat.run %>% select(where(~n_distinct(.) > 1)) # Remove features with only 1 value among samples
    dat.run <- as.matrix(dat.run)
    
    shap_values <- shap.values(xgb_model = model, X_train = dat.run)
    shap_values <- shap_values$shap_score
    
    shap_long <- shap.prep(xgb_model = model, X_train = dat.run)
    shap_long <- shap.prep(shap_contrib = shap_values, X_train = dat.run)
    
    shap_long <- shap_long[order(shap_long$ID, factor(shap_long$variable, levels = cp_colnames)),]
    
    if (fold == 1) {
      shap_list <- shap_long
    }else{
      shap_list[,3:6] = shap_list[,3:6] + shap_long[,3:6,]
    }
    
  }
  
  # Print the result
  # shap_list[,3:6] = shap_list[,3:6] / 5
  
  elix_shap = shap_long %>% filter(variable == 'elix_index_mortality')
  elix_shap$variable = run_name
  if (nrow(elix_shap_matrix) == 0){
    elix_shap_matrix = elix_shap
  }else{
    elix_shap_matrix <- rbind(elix_shap_matrix, elix_shap)
  }
  
  shap_summary <- shap_long %>%
    group_by(variable) %>%
    summarize(mean_abs_value = mean(abs(value), na.rm = TRUE),
              se = sd(abs(value), na.rm = TRUE) / sqrt(sum(!is.na(value))))
  
  shap_summary <- shap_summary %>%
    mutate(mean_se = paste0(round(mean_abs_value, 4), " ± ", round(se, 4))) %>%
    select(variable, mean_se)
  
  if (nrow(shap_value_matrix) == 0){
    shap_value_matrix = shap_summary
    colnames(shap_value_matrix)[ncol(shap_value_matrix)] <- run_name
  }else{
    shap_value_matrix <- merge(shap_value_matrix, shap_summary, by=c('variable'))
    colnames(shap_value_matrix)[ncol(shap_value_matrix)] <- run_name
  }
  
  plot <- shap.plot.summary(shap_long)
  
  ggsave(filename = paste0('results/model_results/20250411/patient_match/xgb/shap_plots/environmental_', run_name,'_shap_summary.pdf'), 
         plot = plot, 
         width = 8, 
         height = 6)
  
}

elix_shap_matrix <- elix_shap_matrix %>% filter(variable %in% all_runs)
# elix_shap_matrix <- elix_shap_matrix[order(elix_shap_matrix$ID, factor(elix_shap_matrix$variable)),]
elix_shap_matrix$variable = factor(elix_shap_matrix$variable)
elix_plot <- shap.plot.summary(elix_shap_matrix)

order_list <- c("DS_Entero_cp","ESBL_cp","VSE_cp","VRE_cp","CDiff_cp",
                "MSSA_cp","MRSA_cp","DS_PsA_cp","DR_PsA_cp","elix_index_mortality")
shap_value_matrix$variable <- factor(shap_value_matrix$variable, levels = order_list)
shap_value_matrix_ordered <- shap_value_matrix[order(shap_value_matrix$variable), ]

write.csv(x = shap_value_matrix, file = 'results/model_results/20250411/patient_match/xgb/shap_plots/shap_values.csv',row.names = FALSE)
ggsave(filename = paste0('results/model_results/20250411/patient_match/xgb/shap_plots/elixhauser_shap_summary.pdf'), 
       plot = elix_plot, 
       width = 8, 
       height = 6)

# Plot CP Correlation
corr_plot_cp = c("DS_Entero_cp","ESBL_cp","VSE_cp","VRE_cp","CDiff_cp")
corr_plot_cohort = c("DS_E_coli","ESBL_K_pneumoniae","VSE_faecalis","VRE_faecium")

for (r in corr_plot_cohort) {
  for (cp_col in corr_plot_cp) {
    plot_data <- cc_final %>%
      filter(run == r) %>%
      select(group, !!sym(cp_col)) %>%
      rename(corr_cp = !!sym(cp_col)) %>%
      group_by(group)
    
    # Plot
    cdf_data <- plot_data %>%
      group_modify(~ {
        d <- .x %>%
          arrange(corr_cp) %>%
          mutate(
            y = cumsum(rep(1, n())) / n(),  # proportion ≤ value
            x = corr_cp
          )
        return(d)
      })
    
    p <- ggplot(cdf_data, aes(x = x, y = y, color = group)) +
      geom_line(size = 1.2) +
      labs(
        title = paste("Run:", r, "| CP:", cp_col),
        x = "corr_cp (ascending)",
        y = "Cumulative distribution (≤ cp)"
      ) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      scale_color_manual(values = c("case" = "#D4A64A", "control" = "#6A7B9B")) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(filename = paste0("/data/tide/projects/ho_infxn_ml/results/figures_tables/20250411/patient_match/cp_correlation_plots/", r, "_", cp_col, ".pdf"), plot = p, width = 6, height = 4)
    
  }
  
}

cdf_data <- cc_final %>%
  filter(run %in% corr_plot_cohort) %>%
  pivot_longer(
    cols = all_of(corr_plot_cp),
    names_to = "cp_col",
    values_to = "corr_cp"
  ) %>%
  group_by(run, cp_col, group) %>%
  arrange(corr_cp, .by_group = TRUE) %>%
  mutate(
    y = cumsum(rep(1, n())) / n(),  # empirical CDF
    x = corr_cp
  ) %>%
  ungroup()

# Plot with facet_grid
p <- ggplot(cdf_data, aes(x = x, y = y, color = group)) +
  geom_line(size = 1.2) +
  facet_grid(run ~ cp_col, scales = "free_x") +
  labs(
    x = "corr_cp (ascending)",
    y = "Cumulative distribution (≤ cp)",
    title = "Cumulative Distribution by Run and CP"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = c("case" = "#D4A64A", "control" = "#6A7B9B")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10)
  )

# Save full grid plot (adjust size as needed)
ggsave(
  filename = "/data/tide/projects/ho_infxn_ml/results/figures_tables/20250411/patient_match/cp_correlation_plots/all_facet_grid_plot.pdf",
  plot = p,
  width = 16, height = 10
)

#### GLOBAL CLEAN UP ####
rm(mainDir, cc_final, clr.results)

