#### HEADER ####
# Description: Figures for HO infxn paper
# Dependencies: None
# Authors: Sanjat Kanjilal
# Last updated: 2023-12-16

##### LIBRARIES #####
library(tidyverse)
library(plotly)
library(conflicted)
conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::lag())
conflicts_prefer(lubridate::month)
conflicts_prefer(lubridate::year)
conflicts_prefer(lubridate::week)
conflicts_prefer(lubridate::quarter)
conflicts_prefer(dplyr::first)

#### IMPORT ####
main_dir = '~/ho_infxn_clean_code/'
clr <- read_csv(file = paste0(main_dir,"clr_coefficients.csv"))
clrnew<- read_csv(file = paste0(main_dir,"clr_coefficients_removed_prior_pathogen.csv"))
table1 <- read_csv(file = paste0(main_dir,"Table1_Baseline_cohort_characteristics.csv"))
sample.size <- read_csv(file = paste0(main_dir,"new_sample_sizes.csv"))
clean.data <- read_csv(file = paste0(main_dir,"age_filtered_final_dataset_for_models_20241207.csv"))

#### TABLE 1 ####

# Data prep
clean.data.table1 <- clean.data %>%
  filter(match == "environmental") %>%
  filter(!(group == "case" & time_to_infxn < 0)) %>%
  mutate(run = factor(run, levels = c("DS_E_coli", "ESBL_E_coli", "DS_K_pneumoniae", "ESBL_K_pneumoniae", 
                                      "C_diff", "VSE_faecalis", "VRE_faecium", 
                                      "MSSA", "MRSA", "DS_P_aeruginosa", "DR_P_aeruginosa"
  ))) %>%
  filter(!is.na(run))

# Sample size
sample_size <- clean.data.table1 %>%
  group_by(run, group) %>%
  summarise(sample_size = n())

sample_size %>%
  ggplot(aes(x = run, y = sample_size, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  labs(title = "Sample size", x = "Pathogen", y = "Sample size")

# Age
age <- clean.data.table1 %>%
  group_by(run, group) %>%
  summarise(n = n(),
            age_mean = mean(age),
            age_sd = sd(age),
            age_se = age_sd / sqrt(n))

age %>% 
  ggplot(aes(x = run, y = age_mean, fill = group)) +
  geom_point(aes(color = group), 
             position = position_dodge(width = 0.9), 
             size = 3) +
  # geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
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

# Elixhauser index
elix <- clean.data.table1 %>%
  group_by(run, group) %>%
  summarise(n = n(),
            elix_mean = mean(indiv_score),
            elix_sd = sd(indiv_score),
            elix_se = elix_sd / sqrt(n))

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

# Categorical variables
# Sex
sex <- clean.data.table1 %>%
  select(run, group, group_index, PatientID, sex) %>%
  distinct() %>%
  count(run, group, sex) %>%
  group_by(run, group) %>%
  mutate(total =sum(n)) %>%
  ungroup() %>%
  mutate(female = (n / total)*100) %>%
  filter(sex != "Male")

sex %>%
  ggplot(aes(x = run, y = female)) +
  geom_point(aes(color = group), position = position_dodge(width = 0.9)) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  ylim(0, 100) +
  labs(title = "Sex (% female)", x = "Pathogen")

# Prior surgery
surgery <- clean.data.table1 %>%
  select(run, group, group_index, PatientID, any_surgery) %>%
  distinct() %>%
  count(run, group, any_surgery) %>%
  group_by(run, group) %>%
  mutate(total =sum(n)) %>%
  ungroup() %>%
  mutate(percent_prior_surgery = (n / total)*100) %>%
  filter(any_surgery != 0)

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
  select(run, group, percent_prior_surgery) %>% 
  distinct() %>%
  left_join(sex %>% select(run, group, female) %>% distinct()) %>%
  pivot_longer(percent_prior_surgery:female, names_to = "var", values_to = "percent") %>%
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
  
  
# Prior abx exposure
abx <- clean.data.table1 %>%
  select(run, group, group_index, PatientID, anti_anaerobe_0_60:tetracycline_0_60) %>%
  distinct() %>%
  rowwise() %>%
  mutate(total_courses_0_60 = sum(c_across(contains("_0_60")), na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_longer(anti_anaerobe_0_60:total_courses_0_60, names_to = "abx", values_to = "courses") %>%
  filter(courses > 0) %>%
  count(run, group, PatientID, abx) %>%
  group_by(run, group, abx) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  select(run, group, abx, total) %>%
  distinct() %>%
  left_join(sample_size) %>%
  mutate(courses_group = (total / sample_size)*100) %>%
  mutate(abx = factor(abx, levels = c("total_courses_0_60", "penicillin_0_60", "anti_staph_beta_lactam_0_60", "extended_spectrum_penicillin_0_60", 
                                      "cephalosporin_0_60", "extended_spectrum_cephalosporin_0_60", "carbapenem_0_60",
                                      "glycopeptide_0_60", "fluoroquinolone_0_60", "macrolide_0_60", "lincosamide_0_60", 
                                      "tetracycline_0_60", "sulfonamide_0_60", "anti_anaerobe_0_60", "anti_Cdiff_0_60")))

# Plot sex, prior surgery and total 
# Plot total abx courses
abx %>%
  filter(abx == "total_courses_0_60") %>%
  ggplot(aes(x = run, y = courses_group, fill = group)) +
  geom_point(aes(color = group,
                 size = 3), 
             position = position_dodge(width = 0.9)) +
  theme_minimal() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_blank()) +
  labs(title = "Number of antibiotic courses", x = "Pathogen")

# Plot Individual antibiotic classes
abx %>%
  filter(abx != "total_courses_0_60") %>%
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
    panel.background = element_blank()
  ) +
  labs(
    title = "Antibiotic courses by class",
    x = "Pathogen"
  )

# Clean up
rm(age, elix, sample_size, sex, surgery, abx)

#### FIGURE: HOSPITAL-ONSET INFECTION AND  COLONIZATION PRESSURE PLOTS ####

# Data prep
orgs_to_keep <- c("C_diff", "DS_E_coli", "DS_K_pneumoniae", "ESBL_E_coli", "ESBL_K_pneumoniae",
                  "MRSA", "MSSA", "VRE_faecium", "VSE_faecium", "VSE_faecalis",
                  "DR_P_aeruginosa", "DS_P_aeruginosa")

clr.prepped <- clrnew %>%
  filter(target %in% orgs_to_keep) %>%
  mutate(direction = ifelse(coef > 1, "pos", "neg")) %>%
  mutate(org_group = case_when(
    target %in% c("C_diff", "DS_E_cloacae", "DS_E_coli", "DS_K_oxytoca", "DS_K_pneumoniae", 
                  "DS_P_mirabilis", "DS_S_marcescens", "ESBL_E_cloacae", "ESBL_E_coli",
                  "ESBL_K_pneumoniae", "ESBL_P_mirabilis", "VRE_faecium", "VSE_faecalis",
                  "VSE_faecium") ~ "enterics",
    target %in% c("MRSA", "MSSA") ~ "skin_flora",
    target %in% c("DR_P_aeruginosa", "DS_P_aeruginosa") ~ "environmentals"
  )) %>%
  mutate(target = factor(target, levels = c("DS_E_coli", "ESBL_E_coli", "DS_K_pneumoniae", "ESBL_K_pneumoniae", 
                                            "DS_E_cloacae", "ESBL_E_cloacae", "C_diff", "VSE_faecalis", "VRE_faecium", 
                                            "VSE_faecium", "MRSA", "MSSA", "DR_P_aeruginosa", "DS_P_aeruginosa"
  ))) %>%
  mutate(variable = factor(variable, levels = c("DS_Entero_cp", "ESBL_cp", "CDiff_cp", "VSE_cp", "VRE_cp", 
                                                "MSSA_cp", "MRSA_cp", "DS_PsA_cp", "DR_PsA_cp"))) 

# clr.preppednew <- clr.preppednew %>%
#   select(match:variable, coefnew = coef, lower_CInew = lower_CI, upper_CInew = upper_CI, sig_flagnew = sig_flag) %>%
#   distinct()
# 
# test <- full_join(clr.prepped, clr.preppednew) %>%
#   select(match:variable, coefnew, coef, lower_CInew, lower_CI, upper_CInew, upper_CI, sig_flagnew, sig_flag)
  
# Function to filter CLR data for target and CP sets
target.cp.data <- function(data, target_set, cp_set) 
  {
  
  filtered.data <- data %>%
    
    # Filter for target set
    filter(target %in% target_set) %>%
    mutate(target = factor(target, levels = target_set)) %>%
    
    # Filter for CP set
    filter(variable %in% cp_set) %>%
    mutate(variable = factor(variable, levels = cp_set)) %>%
    
    arrange(variable, target)
    
    return(filtered.data)
  
}

# Function to plot CLR data for target and CP sets
target.cp.plot <- function(filtered.data) 
  {
  filtered.plot <- filtered.data %>%
    ggplot(aes(x = variable, y = coef, ymin=lower_CI, ymax=upper_CI,
               middle = coef, lower = lower_CI, upper = upper_CI)) +
    geom_boxplot(stat = "identity", aes(fill = target)) +
    geom_hline(yintercept = 1, color = "red", linetype='dotted') +
    scale_x_discrete(drop=FALSE) +
    theme_minimal() +
    theme(panel.border = element_rect(linewidth = 1, fill = NA),
          panel.background = element_blank()) +
    ylim(0.7, 1.51) +
    labs(x = "Colonization pressure", y = "Coefficient")
  
  return(filtered.plot)
  
}

# Define target sets
enteric.target <- c("DS_E_coli", "ESBL_E_coli", "DS_K_pneumoniae", "ESBL_K_pneumoniae", "C_diff", "VSE_faecalis", "VRE_faecium")
skin.target <- c("MSSA", "MRSA")
environmental.target <-c("DS_P_aeruginosa", "DR_P_aeruginosa")

# Define CP sets
enteric.cp <- c("DS_Entero_cp", "ESBL_cp", "CDiff_cp", "VSE_cp", "VRE_cp")
skin.cp <-c("MSSA_cp", "MRSA_cp")
environmental.cp <- c("DS_PsA_cp", "DR_PsA_cp")

# Build cognate datasets
cognate.enteric <- target.cp.data(clr.prepped, enteric.target, enteric.cp)
cognate.skin <- target.cp.data(clr.prepped, skin.target, skin.cp)
cognate.environmental <- target.cp.data(clr.prepped, environmental.target, environmental.cp)

# Build non-cognate datasets
noncognate.enteric.skin <- target.cp.data(clr.prepped, enteric.target, skin.cp)
noncognate.enteric.env <- target.cp.data(clr.prepped, enteric.target, environmental.cp)

noncognate.skin.enteric <- target.cp.data(clr.prepped, skin.target, enteric.cp)
noncognate.skin.env <- target.cp.data(clr.prepped, skin.target, environmental.cp)

noncognate.env.enteric <- target.cp.data(clr.prepped, environmental.target, enteric.cp)
noncognate.env.skin <- target.cp.data(clr.prepped, environmental.target, skin.cp)

# Plot cognate datasets
cognate.enteric.plot <- target.cp.plot(cognate.enteric)
cognate.skin.plot <- target.cp.plot(cognate.skin)
cognate.environmental.plot <- target.cp.plot(cognate.environmental)

# Plot non-cognate datastes
noncognate.enteric.skin.plot <- target.cp.plot(noncognate.enteric.skin)
noncognate.enteric.env.plot <- target.cp.plot(noncognate.enteric.env)

noncognate.skin.enteric.plot <- target.cp.plot(noncognate.skin.enteric)
noncognate.skin.env.plot <- target.cp.plot(noncognate.skin.env)

noncognate.env.enteric.plot <- target.cp.plot(noncognate.env.enteric)
noncognate.env.skin.plot <- target.cp.plot(noncognate.env.skin)

# Clean up
rm(target.cp.data, target.cp.plot, enteric.target, skin.target, environmental.target,
   cognate.enteric, enteric.cp, skin.cp, cognate.skin, environmental.cp, cognate.environmental, 
   noncognate.enteric.skin, noncognate.enteric.env, noncognate.skin.enteric, 
   noncognate.skin.env, noncognate.env.enteric, noncognate.env.skin, noncognate.enteric.skin.plot, 
   noncognate.enteric.env.plot, noncognate.skin.enteric.plot, noncognate.skin.env.plot, 
   noncognate.env.enteric.plot, noncognate.env.skin.plot)

#### POST-HOC ANALYSES ####

# Abx pressure

# Clean up
rm(clr, clr.prepped, orgs_to_keep, sample.size, table1)

