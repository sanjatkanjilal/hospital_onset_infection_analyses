#### HEADER ####
# Description: Script to process micro data for HO-infection study
# Authors: Ted Pak, Sanjat Kanjilal
# Last updated: 2024-07-12

#### LIBRARIES ####
library(tidyverse)
library(data.table)
library(conflicted)
conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::lag())
conflicts_prefer(lubridate::month)
conflicts_prefer(lubridate::year)
conflicts_prefer(lubridate::week)
conflicts_prefer(lubridate::quarter)
conflicts_prefer(dplyr::first)

#### FUNCTIONS ####

# Categorize organisms
categorizeOrgUsingAST <- function(SD, org_clean, org_name_raw1, org_name_raw2) {
  neg_cx <- SD$neg_cx[1] 
  if (isTRUE(neg_cx) || (!is.na(neg_cx) && neg_cx == 'X')) return(as.character(NA)); 
  
  ast <- SD$pheno_merged[!is.na(SD$pheno_merged)]
  names(ast) <- SD$AST_code[!is.na(SD$pheno_merged)]
  org <- toupper(org_clean)
  org_name_raw <- toupper(paste0(org_name_raw1, "\n", org_name_raw2))
  acc <- SD$accession[1]
  
  # NOTE: We use \b within many of the regexes to ensure the presence of a word boundary
  # Staph aureus and CoNS
  if (grepl("\\bMRSA\\b", org) && !grepl("\\bNOT ", org) && !grepl("\\bNEGATIVE ", org)) return("MRSA")
  if (grepl("\\bSTAPH", org)) {
    if (grepl("\\bAUREUS\\b", org)) {
      if (grepl("\\bMETHICILLIN RESISTANT\\b", org, ignore.case=TRUE)) return("MRSA")
      if (is.na(ast["OXA"]) && is.na(ast["FOX"])) return ("M?SA") 
      if (!is.na(ast["OXA"])) return(ifelse(ast[["OXA"]] == "Susceptible", "MSSA", "MRSA"))
      if (!is.na(ast["FOX"])) return(ifelse(ast[["FOX"]] %in% c("Susceptible", "Negative"), "MSSA", "MRSA")) ## TODO sometimes positive/negative?
      return("MRSA")
    } else if (grepl("\\b(COAGULASE POSITIVE|LUGDUNENSIS|SCHLEIFERI)\\b", org, ignore.case=TRUE)) {
      if (is.na(ast["OXA"])) return ("M?-CoPS")
      if (ast[["OXA"]] == "Susceptible") return("MS-CoPS")
      return("MR-CoPS")
    } else {
      if (is.na(ast["OXA"])) return ("M?-CoNS")
      if (ast[["OXA"]] == "Susceptible") return("MS-CoNS")
      return("MR-CoNS")
    }
  }
  
  # VSE and VRE
  if (grepl("\\bVRE\\b", org) && !grepl("\\bNEGATIVE ", org)) return("VRE")
  if (grepl("\\bENTEROCOCCUS", org) && !grepl("\\bNOT +ENTEROCOCCUS", org)) {
    if (grepl("(\\b|-)VRE|VANCOMYCIN RESISTANT", org_name_raw, ignore.case=TRUE)) return("VRE")
    if (is.na(ast["VAN"])) return("V?E")
    if (ast[["VAN"]] == "Susceptible") return("VSE")
    return("VRE")
  }
  
  # Non-enterococcal Streptococci
  if (grepl("\\bSTREP", org)) return("Streptococcus")
  
  # Pseudomonas
  if (grepl("\\bPSEUDOMONAS\\b", org)) {
    if (!is.na(ast["IPM"]) && ast[["IPM"]] != "Susceptible") return("CR-Pseudomonas")
    if (!is.na(ast["MPM"]) && ast[["MPM"]] != "Susceptible") return("CR-Pseudomonas")
    if (!is.na(ast["CAZ"]) && ast[["CAZ"]] != "Susceptible") return("DR-Pseudomonas")
    if (!is.na(ast["FEP"]) && ast[["FEP"]] != "Susceptible") return("DR-Pseudomonas")
    if (!is.na(ast["TZP"]) && ast[["TZP"]] != "Susceptible") return("DR-Pseudomonas")
    if (!any(is.na(ast[c("CAZ", "FEP", "TZP")])) && all(ast[c("CAZ", "FEP", "TZP")] == "Susceptible"))
      return("DS-Pseudomonas")
    return("D?-Pseudomonas")
  }
  # Acinetobacter and Stenotrophomonas
  if (grepl("\\bACINETOBACTER\\b", org)) {
    if (!is.na(ast["SAM"])) return(ifelse(ast[["SAM"]] == "Susceptible", "DS-Acinetobacter", "DR-Acinetobacter"))
    return("D?-Acinetobacter")
  }
  if (grepl("\\bSTENOTROPHOMONAS\\b", org)) {
    if (!is.na(ast["SXT"]) && ast[["SXT"]] != "Susceptible") return("DR-Stenotrophomonas")
    if (!is.na(ast["MIN"]) && ast[["MIN"]] != "Susceptible") return("DR-Stenotrophomonas")
    if (!is.na(ast["SXT"]) && any(!is.na(ast[c("MIN", "TET")])) && 
        all(ast[c("SXT", "MIN", "TET")] == "Susceptible", na.rm = TRUE))
      return("DS-Stenotrophomonas")
    return("D?-Stenotrophomonas")
  }
  
  # AmpC Enterobacterales
  # Defined per IDSA 2022 guidelines - https://www.idsociety.org/practice-guideline/amr-guidance-2.0/
  # Of note, we define AmpCs EXCLUSIVE to the ESBL phenotypes below;
  #   it has to be sensitive to >=3 of 4 antibiotics that typically define an ESBL
  #   while harboring a clinically significant chance of *inducible* resistance.
  # Also note: although per IDSA, FEP MIC should be checked (and <=4) to permit assignment
  #   to an AmpC category that is deemed OK for FEP treatment, but we are ignoring this right now.
  if (grepl("\\bENTEROBACTER CLOACAE\\b|\\bAEROGENES\\b|\\bCITROBACTER FREUNDII\\b", org)) {
    if (!is.na(ast["IPM"]) && ast[["IPM"]] != "Susceptible") return("CR-Enterobacterales")
    if (!is.na(ast["MPM"]) && ast[["MPM"]] != "Susceptible") return("CR-Enterobacterales")
    if (!is.na(ast["CRO"]) && ast[["CRO"]] != "Susceptible") return("ESBL-Enterobacterales")
    if (!is.na(ast["FEP"]) && ast[["FEP"]] != "Susceptible") return("ESBL-Enterobacterales")
    if (!is.na(ast["CAZ"]) && ast[["CAZ"]] != "Susceptible") return("ESBL-Enterobacterales")
    if (!is.na(ast["TZP"]) && ast[["TZP"]] != "Susceptible") return("ESBL-Enterobacterales")
    if (sum(!is.na(ast[c("CRO", "CAZ", "FEP", "TZP")])) >= 3 && 
        all(ast[c("CRO", "CAZ", "FEP", "TZP")] == "Susceptible", na.rm = TRUE))
      return("AmpC-Enterobacterales")
    return("?ESBL-Enterobacterales")
  }
  
  # Other Enterobacterales
  # Regex is based on the list of genera in https://en.wikipedia.org/wiki/Enterobacterales
  # Per discussion w SJK, if *any* of the ESBL-defining abx are present and all of those present are S,
  #    ok to call it as non-ESBL 
  if (grepl(paste0("\\b((CITRO|ENTERO|PLURALI|CRONO)BACTER|ESCHERICHIA|KLUYVERA|PANTOEA|CEDECEA|LELLIOTTIA",
                   "|PLESIOMONAS|PROTEUS|PROVIDENCIA|SERRATIA|HAFNIA|LECLERCIA",
                   "|YERSINIA|(EDWARDSI|LEMINOR|KLEBSI|MORGAN|SHIG|RAOULT|SALMON)ELLA)\\b"), org)) {
    if (!is.na(ast["IPM"]) && ast[["IPM"]] != "Susceptible") return("CR-Enterobacterales")
    if (!is.na(ast["MPM"]) && ast[["MPM"]] != "Susceptible") return("CR-Enterobacterales")
    if (grepl("\\b(ESBL|EXTENDED SPECTRUM BETA LACTAMASE)\\b", org, ignore.case=TRUE))
      return ("ESBL-Enterobacterales")
    if (!is.na(ast["CRO"]) && ast[["CRO"]] != "Susceptible") return("ESBL-Enterobacterales")
    if (!is.na(ast["FEP"]) && ast[["FEP"]] != "Susceptible") return("ESBL-Enterobacterales")
    if (!is.na(ast["CAZ"]) && ast[["CAZ"]] != "Susceptible") return("ESBL-Enterobacterales")
    if (!is.na(ast["TZP"]) && ast[["TZP"]] != "Susceptible") return("ESBL-Enterobacterales")
    if (any(!is.na(ast[c("CRO", "CAZ", "FEP", "TZP")])) && 
        all(ast[c("CRO", "CAZ", "FEP", "TZP")] == "Susceptible", na.rm = TRUE))
      return("Non-ESBL-non-AmpC-Enterobacterales")
    return("?ESBL-Enterobacterales")
  }
  
  # Other Gram positives
  if (grepl(paste0("\\b(GRAM POSITIVE|(LACTO)?BACILLUS|(ARCANO|CUTI|PROPIONI|BIFIDO|EU|BREVI)BACTERIUM",
                   "|(GEM|GARDNER|EGGERTH|GLOBICAT|TRUEPER)ELLA|NOCARDIA|ROTHIA|KOCURIA|(ARTHRO|DERMA)BACTER",
                   "|ACTINOMYCES|PARVIMONAS|FINEGOLDIA|(DERMA|PEPTOSTREPTO|AERO|MICRO|RHODO|LACTO)COCCUS",
                   "|CELLULOMONAS|CLOSTRIDIUM|LEUCONOSTOC|LISTERIA|ERYSIPELOTHRIX|CORYNE\\w+",
                   "|LEIFSONIA|POS COCCI)\\b"), org))
    return("Other-gram-positive")
  
  # Other Gram negatives
  if (grepl(paste0("\\b(GRAM NEGATIVE|GNRS?|HAEMOPHILUS|(AGGREGATI|ACHROMO|CAMPYLO|HELICO)BACTER",
                   "|(CARDIO|FLAVO|CHRYSEO|FUSO|SPHINGO)BACTERIUM|BURKHOLDERIA|DELFTIA|NEISSERIA",
                   "|(AERO|BREVUNDI|CHRYSEO|FLAVI|ROSEO|SPHINGO|COMA)MONAS|RHIZOBIUM|VIBRIO",
                   "|OCHROBACTRUM|ALCALIGENES|MYROIDES|ELIZABETHKINGIA|CAPNOCYTOPHAGA|BACTEROIDES",
                   "|(VEILLON|KING|SHEWAN|EIKEN|MORAX|BORDET|PREVOT|OLIG|PASTEUR|EWING|WEEKS)ELLA",
                   "|CUPRIAVIDUS)\\b"), org))
    return("Other-gram-negative")
  
  # Candida albicans
  if (grepl(paste0("\\bCANDIDA ALBICANS\\b"), org))
    return("C_albicans")

  # Other Candida species
  if (grepl(paste0("\\b(AURIS|DUBLINIENSIS|GLABRATA|GUILLIERMONDII|HAEMULONII|KRUSEI|LUSITANIAE|METAPSILOSIS",
                   "|ORTHOPSILOSIS|PARAPSILOSIS|CANDIDA SPECIES|TROPICALIS)\\b"), org))
    return("Candida_non_albicans")
  
  return("Unknown")
}

# Create classification categories
createOrgCategoryColumn <- function(dt.micro) {
  dt.micro[
    , 
    org_category := categorizeOrgUsingAST(.SD, org_clean, org_name_raw1, org_name_raw2),
    by = c("accession", "org_clean", "org_name_raw1", "org_name_raw2"),
    .SDcols = c("neg_cx", "prelim_AST", "accession", "AST_code", "pheno_merged")
  ]
}

# Create treatment categories
createOrgTreatmentCategoryColumn <- function(dt.micro) {
  pathogenTreatmentGroups <- list(
    c("MRSA", "MR-CoNS", "M?-CoNS", "M?SA", "MR-CoPS", "M?-CoPS", "Other-gram-positive"),
    c("MSSA", "MS-CoNS", "MS-CoPS"),
    c("VRE", "V?E"),
    c("ESBL-Enterobacterales", "?ESBL-Enterobacterales", "D?-Pseudomonas", "DR-Pseudomonas"),
    c("DR-Acinetobacter", "D?-Acinetobacter"),
    c("DR-Stenotrophomonas", "D?-Stenotrophomonas"),
    c("CR-Pseudomonas", "CR-Enterobacterales"),
    c("C_albicans", "Candida non_albicans")
  )
  hashmap <- new.env(hash = TRUE)
  sapply(unique(dt.micro$org_category), function(o_cat) { hashmap[[o_cat]] <- o_cat })
  for (group in pathogenTreatmentGroups) {
    for (o_cat in group) { hashmap[[o_cat]] <- paste(group, collapse = "/") }
  }
  dt.micro[
    !is.na(org_category),
    org_treatment_category := sapply(org_category, function(o_cat) { hashmap[[o_cat]] })
  ]
}

#### IMPORT DATA ####

project_path <- "/data/tide/projects/ho_infxn_ml/input_data/"

micro.raw <- data.table::fread(paste0(project_path, "AllMGB_micro_20150525-20240702.csv"))
  
#### DATA CLEANING ####

# Add organism / treatment categories
createOrgCategoryColumn(micro.raw)

createOrgTreatmentCategoryColumn(micro.raw)

write_csv(micro.raw, file = paste0(project_path, "micro_org_cat_raw_v20240702.csv"))
# micro.raw <- data.table::fread(paste0(project_path, "micro_org_cat_raw_v20240702.csv"))

# Add EDW PatientIDs
mrns <- read_csv("/data/tide/data/edw/MRN_map/EDW_MRN_map_20240713.csv") %>%
  mutate(MRN = as.character(MRN))

mrn.IDs.only <- mrns %>%
  select(PatientID, hospital, MRN) %>%
  distinct() %>%
  mutate(MRN = as.character(MRN))

micro.raw <- micro.raw %>%
  # slice(1:100000) %>%
  left_join(mrn.IDs.only)

# Fill in missing PatientIDs
missingIDs <- micro.raw %>%
  filter(is.na(PatientID)) %>%
  select(name, DOB) %>% 
  distinct()

missingIDs.micro <- mrns %>%
  select(PatientNM, DOB, PatientID) %>% 
  distinct() %>%
  right_join(missingIDs, by = join_by(PatientNM == name, DOB == DOB)) %>%
  arrange(PatientNM, PatientID) %>%
  filter(!is.na(PatientID)) %>%
  rename(name = PatientNM, PatientID.missing = PatientID)

micro.raw <- micro.raw %>%
  rename(PatientID.old = PatientID) %>%
  left_join(missingIDs.micro) %>%
  mutate(PatientID = ifelse(!is.na(PatientID.old), PatientID.old, PatientID.missing)) %>%
  select(-PatientID.old, -PatientID.missing) %>%
  distinct()

# Clean up duplicate names (take the longest name)
micro.names.clean <- micro.raw %>%
  select(hospital, MRN, name) %>%
  distinct() %>%
  group_by(hospital, MRN) %>%
  mutate(count = n()) %>%
  mutate(name_length = str_length(name),
         name_length_max = max(name_length)) %>%
  mutate(name.clean = ifelse(name_length == name_length_max, name, NA)) %>%
  fill(name.clean, .direction = "updown") %>%
  ungroup() %>%
  select(hospital, MRN, name, name.clean) %>%
  distinct()

micro.raw <- micro.raw %>%
  left_join(micro.names.clean) %>%
  select(hospital, name = name.clean, PatientID, MRN:org_treatment_category) %>%
  distinct()

rm(mrns, mrn.IDs.only, missingIDs, missingIDs.micro, micro.names.clean)

# Limit to necessary fields only and filter for relevant cultures
micro <- micro.raw %>%
  # head(10000) %>%
  select(hospital, MRN, PatientID, CSN, accession, coll_datetime, loc_cat:loc_specific, specimen_code,
         site_clean, neg_cx, org_clean, org_treatment_category, org_name_raw1:org_name_raw2, 
         prelim_AST, AST_code, AST_panel, pheno_merged) %>%
  distinct()

micro <- micro %>%
  filter(is.na(prelim_AST)) %>%
  filter(is.na(neg_cx)) %>%
  filter(org_treatment_category != "Unknown") %>%
  filter(org_clean != "") %>%
  select(-prelim_AST, -neg_cx)

# Flag organisms of interest
micro <- micro %>%
  mutate(keep = case_when(
    org_clean %in% c("ACINETOBACTER BAUMANNII COMPLEX", 
                     "BURKHOLDERIA SPECIES",
                     "CITROBACTER FREUNDII COMPLEX", "CITROBACTER KOSERI (DIVERSUS)",
                     "ENTEROBACTER CLOACAE COMPLEX",
                     "ENTEROCOCCUS FAECALIS", "ENTEROCOCCUS FAECIUM",
                     "ESCHERICHIA COLI",
                     "KLEBSIELLA AEROGENES", "KLEBSIELLA OXYTOCA", "KLEBSIELLA PNEUMONIAE",
                     "MORGANELLA MORGANII",
                     "PROTEUS MIRABILIS", "PROTEUS VULGARIS",
                     "PROVIDENCIA SPECIES",
                     "PSEUDOMONAS AERUGINOSA", "PSEUDOMONAS AERUGINOSA (MUCOID)",
                     "SALMONELLA SPECIES",
                     "SERRATIA MARCESCENS",
                     "SHIGELLA SPECIES",
                     "STAPHYLOCOCCUS AUREUS",
                     "STENOTROPHOMONAS MALTOPHILIA",
                     "STREPTOCOCCUS AGALACTIAE (GROUP B)", "STREPTOCOCCUS ANGINOSUS GROUP", "STREPTOCOCCUS MITIS/ORALIS GROUP", 
                     "STREPTOCOCCUS PNEUMONIAE", "STREPTOCOCCUS PYOGENES (GROUP A)",
                     "CANDIDA AURIS", "CANDIDA DUBLINIENSIS", "CANDIDA GLABRATA", "CANDIDA GUILLIERMONDII", 
                     "CANDIDA HAEMULONII", "CANDIDA KRUSEI", "CANDIDA LUSITANIAE", "CANDIDA METAPSILOSIS", 
                     "CANDIDA ORTHOPSILOSIS", "CANDIDA PARAPSILOSIS", "CANDIDA SPECIES", "CANDIDA TROPICALIS") ~ "X")) %>%
  filter(!is.na(keep)) %>%
  select(-keep) %>%
  distinct()

# Set date/timestamps to UTC to align with EDW data
micro <- micro %>%
  rename(coll_datetime_UTC = coll_datetime) %>%
  mutate(coll_datetime = with_tz(coll_datetime_UTC, "America/New_York")) %>%
  select(hospital:coll_datetime_UTC, coll_datetime, loc_cat:pheno_merged)
  
# Create unique ID for each organism in an accession
test <- micro %>%
  group_by(hospital, MRN, accession, specimen_code) %>%
  mutate(organism_index = cur_group_id()) %>%
  ungroup() %>%
  select(hospital:accession, organism_index, coll_datetime_UTC:pheno_merged)

# Generate middle level of groupings (species-specific)
micro <- micro %>%
  # head(1000) %>%
  mutate(org_group_temp_2 = case_when(

    # Staphylococcus aureus
    
    grepl("aureus", org_clean, ignore.case = T) ~ "Staph_aureus",

    # Streptococcus species
    
    grepl("agalactiae", org_clean, ignore.case = T) ~ "Strep_agalactiae",
    grepl("pyogenes", org_clean, ignore.case = T) ~ "Strep_pyogenes",
    grepl("streptococcus pneumoniae", org_clean, ignore.case = T) ~ "Strep_pneumo",
    grepl("anginosus", org_clean, ignore.case = T) ~ "Strep_anginosus",
    grepl("mitis", org_clean, ignore.case = T) ~ "Strep_mitis",

    # Enterococcus species
    
    grepl("faecalis", org_clean, ignore.case = T) ~ "E_faecalis",
    grepl("faecium", org_clean, ignore.case = T) ~ "E_faecium",
    
    # Enterobacterales (non-AmpC) family
    
    grepl("escherichia coli", org_clean, ignore.case = T) ~ "E_coli",
    grepl("klebsiella pneumo", org_clean, ignore.case = T) ~ "K_pneumoniae",
    grepl("oxytoca", org_clean, ignore.case = T) ~ "K_oxytoca",
    grepl("proteus mirab", org_clean, ignore.case = T) ~ "P_mirabilis",
    grepl("serratia", org_clean, ignore.case = T) ~ "S_marcescens",
    grepl("providencia", org_clean, ignore.case = T) ~ "Providencia",
    grepl("koseri", org_clean, ignore.case = T) ~ "C_koseri",
    grepl("morganella", org_clean, ignore.case = T) ~ "M_morganii",
    
    # Enterobacterales (AmpC) family
    
    grepl("cloacae", org_clean, ignore.case = T) ~ "E_cloacae",
    grepl("aerogenes", org_clean, ignore.case = T) ~ "K_aerogenes",
    grepl("citrobacter freundii", org_clean, ignore.case = T) ~ "C_freundii",
    
    # Gram negative non-lactose fermenters
    
    grepl("aeruginosa", org_clean, ignore.case = T) ~ "P_aeruginosa",
    grepl("maltophilia", org_clean, ignore.case = T) ~ "S_maltophilia",
    grepl("baumannii", org_clean, ignore.case = T) ~ "A_baumannii",
    
    # Candida
    grepl("albicans", org_clean, ignore.case = T) ~ "C_albicans",
    grepl("glabrata", org_clean, ignore.case = T) ~ "C_glabrata",
    grepl("parapsilosis", org_clean, ignore.case = T) ~ "C_parapsilosis"
    
  ))

# Generate most granular level of groupings (AST-profile specific)
micro <- micro %>%
  # head(1000) %>%
  mutate(org_group_temp_3 = case_when(
    
    # MRSA / MSSA
    
    org_group_temp_2 == "Staph_aureus" & org_treatment_category == "MRSA/MR-CoNS/M?-CoNS/M?SA/MR-CoPS/M?-CoPS/Other-gram-positive" ~ "MRSA",
    org_group_temp_2 == "Staph_aureus" & org_treatment_category == "MSSA/MS-CoNS/MS-CoPS" ~ "MSSA",
    
    # VRE / VSE
    
    org_group_temp_2 == "E_faecalis" & org_treatment_category == "VSE" ~ "VSE_faecalis",
    org_group_temp_2 == "E_faecium" & org_treatment_category == "VSE" ~ "VSE_faecium",
    org_group_temp_2 == "E_faecalis" & org_treatment_category == "VRE/V?E" ~ "VRE_faecalis",
    org_group_temp_2 == "E_faecium" & org_treatment_category == "VRE/V?E" ~ "VRE_faecium",
    
    
    # Non_ESBL vs ESBL+CRE Enterobacterales (non-AmpC) family 
    
    org_group_temp_2 == "E_coli" & org_treatment_category == "Non-ESBL-non-AmpC-Enterobacterales" ~ "DS_E_coli",
    org_group_temp_2 == "E_coli" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_E_coli",
    org_group_temp_2 == "K_pneumoniae" & org_treatment_category == "Non-ESBL-non-AmpC-Enterobacterales" ~ "DS_K_pneumoniae",
    org_group_temp_2 == "K_pneumoniae" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_K_pneumoniae",
    org_group_temp_2 == "K_oxytoca" & org_treatment_category == "Non-ESBL-non-AmpC-Enterobacterales" ~ "DS_K_oxytoca",
    org_group_temp_2 == "K_oxytoca" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_K_oxytoca",
    org_group_temp_2 == "P_mirabilis" & org_treatment_category == "Non-ESBL-non-AmpC-Enterobacterales" ~ "DS_P_mirabilis",
    org_group_temp_2 == "P_mirabilis" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_P_mirabilis",
    org_group_temp_2 == "S_marcescens" & org_treatment_category == "Non-ESBL-non-AmpC-Enterobacterales" ~ "DS_S_marcescens",
    org_group_temp_2 == "S_marcescens" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_S_marcescens",
    org_group_temp_2 == "Providencia" & org_treatment_category == "Non-ESBL-non-AmpC-Enterobacterales" ~ "DS_Providencia",
    org_group_temp_2 == "Providencia" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_Providencia",
    org_group_temp_2 == "C_koseri" & org_treatment_category == "Non-ESBL-non-AmpC-Enterobacterales" ~ "DS_C_koseri",
    org_group_temp_2 == "C_koseri" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_C_koseri",
    org_group_temp_2 == "M_morganii" & org_treatment_category == "Non-ESBL-non-AmpC-Enterobacterales" ~ "DS_M_morganii",
    org_group_temp_2 == "M_morganii" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_M_morganii",
    
    # Non_ESBL vs ESBL+CRE Enterobacterales (AmpC) family 
    
    org_group_temp_2 == "E_cloacae" & org_treatment_category == "AmpC-Enterobacterales" ~ "DS_E_cloacae",
    org_group_temp_2 == "E_cloacae" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_E_cloacae",
    org_group_temp_2 == "K_aerogenes" & org_treatment_category == "AmpC-Enterobacterales" ~ "DS_K_aerogenes",
    org_group_temp_2 == "K_aerogenes" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_K_aerogenes",
    org_group_temp_2 == "C_freundii" & org_treatment_category == "AmpC-Enterobacterales" ~ "DS_C_freundii",
    org_group_temp_2 == "C_freundii" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "ESBL_C_freundii",
    
    # 'Wild-type' vs drug-resistant Gram negative lactose non-fermenters
    
    org_group_temp_2 == "P_aeruginosa" & org_treatment_category == "DS-Pseudomonas" ~ "DS_P_aeruginosa",
    org_group_temp_2 == "P_aeruginosa" & 
      (org_treatment_category == "ESBL-Enterobacterales/?ESBL-Enterobacterales/D?-Pseudomonas/DR-Pseudomonas" | 
         org_treatment_category == "CR-Pseudomonas/CR-Enterobacterales") ~ "DR_P_aeruginosa",
    org_group_temp_2 == "S_maltophilia" & org_treatment_category == "DS-Stenotrophomonas" ~ "DS_S_maltophilia",
    org_group_temp_2 == "S_maltophilia" & org_treatment_category == "DR-Stenotrophomonas/D?-Stenotrophomonas" ~ "DR_S_maltophilia",
    org_group_temp_2 == "A_baumannii" & org_treatment_category == "DS-Acinetobacter" ~ "DS_A_baumannii",
    org_group_temp_2 == "A_baumannii" & org_treatment_category == "DR-Acinetobacter/D?-Acinetobacter" ~ "DR_A_baumannii"))

# Generate highest level of groupings (body-niche)
micro <- micro %>%
  mutate(org_group_temp_1 = case_when(
    
    # Skin flora
    
    org_group_temp_2 %in% c("Staph_aureus", "Strep_agalactiae", "Strep_pyogenes") ~ "Skin flora",
    
    # Enteric flora
    
    org_group_temp_2 %in% c("C_freundii", "C_koseri", "E_cloacae", "E_coli", "E_faecalis", "E_faecium",
                       "K_aerogenes", "K_oxytoca", "K_pneumoniae", "M_morganii", "P_mirabilis",
                       "Providencia", "S_marcescens", "Strep_anginosus", "Strep_mitis", "C_albicans", 
                       "Candida_non_albicans") ~ "Enteric flora",
    # Environmental flora
    
    org_group_temp_2 %in% c("A_baumannii", "P_aeruginosa", "S_maltophilia") ~ "Environmental flora"))

# Keep only necessary columns
micro <- micro %>%
  select(hospital, PatientID:accession, CSN, loc_cat:loc_specific, site_clean, 
         coll_datetime_UTC, org_clean, org_group_1 = org_group_temp_1, org_group_2 = org_group_temp_2, 
         org_group_3 = org_group_temp_3) %>%
  distinct()

# Save temp dataset
end.date <- gsub("-", "", as.Date(range(micro$coll_datetime_UTC))[[2]])

write_csv(micro, file = paste0(project_path, "micro.ground_truth_nocdiff_20150525-", end.date, ".csv"), na = "")

#### ADD C DIFF DATA ####
cdiff <- read_csv(file = paste0(project_path, "MGB_Cdiff_20240712.csv"))

# Harmonize column names
cdiff <- cdiff %>%
  filter(result == "Positive") %>%
  select(PatientID, CSN = encounter_id, accession = order_id, coll_date = coll_datetime) %>%
  distinct() %>%
  mutate(accession = as.character(accession),
         coll_date = as.character(coll_date),
         coll_date = paste0(coll_date, " 00:00:00"),
         coll_datetime_UTC = ymd_hms(coll_date),
         org_group_1 = "Enteric flora",
         org_group_2 = "C_diff",
         org_group_3 = "C_diff") %>%
  select(-coll_date)

# Bind to micro data
micro.cdiff <- bind_rows(micro, cdiff) %>%
  arrange(hospital, PatientID)

# Save final file
end.date <- gsub("-", "", as.Date(range(micro$coll_datetime_UTC))[[2]])

write_csv(micro.cdiff, file = paste0(project_path, "micro.ground_truth_20150525-", end.date, ".csv"), na = "")

# Clean up
rm(categorizeOrgUsingAST, createOrgCategoryColumn, createOrgTreatmentCategoryColumn, 
   project_path, micro.raw, cdiff, micro, end.date)

