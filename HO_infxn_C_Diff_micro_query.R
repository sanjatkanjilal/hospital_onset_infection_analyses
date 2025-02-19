#### HEADER ####
# !/usr/bin/R

# Script to pull any patients who had a C.difficile test performed at MGB 
# Author: Sanjat Kanjilal 
# Last modified: 2023-05-24

#### LIBRARIES ####
library(RODBC)
library(tidyverse)
library(lubridate)
library(tictoc)
library(conflicted)
conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::lag())

#### OPEN DATABASE #####

edwDB <- odbcConnect(dsn = "phsedw", uid = "Partners\\sk726", pwd = .rs.api.askForPassword("password"))

#### C.DIFF EDW QUERY ####

edw.cdiff.query <-
  paste0("
    SELECT 
    --TOP (1000) 
      r.[PatientID]
    , r.[PatientEncounterID]
    , r.[OrderProcedureID]
    , r.[OrderDateRealNBR] as orderdate
    --, time(r.[OrderdateRealNBR]) as ordertime
    , r.[ResultDTS]
    , c.[ComponentID]
    , c.[ComponentNM]
    , c.[ExternalNM]
    , c.[BaseNM]
    , c.[ComponentCommonNM]
    , r.[ResultTXT]
    , r.[ComponentCommentTXT]
    
  FROM [Epic].[Reference].[Component] as c 
  INNER JOIN [Epic].[Orders].[Result_PHS] as r ON c.ComponentID = r.ComponentID
  
  WHERE
  
    c.[ComponentCommonNM] like '%diffi%' AND
    r.[ResultStatusCD] = 3 /* Final results only */ 
    -- AND r.[ResultDTS] >= convert(datetime, '2015-05-25')"
  )

tic()

edw.cdiff.raw <- sqlQuery(edwDB, edw.cdiff.query) %>% 
  as.data.frame()

toc()

rm(edw.cdiff.query)

#### DATA CLEAN ####

# Rename
edw.cdiff <- edw.cdiff.raw %>%
  select(PatientID, encounter_id = PatientEncounterID, order_id = OrderProcedureID, orderdate, result_datetime = ResultDTS,
         test_id = ComponentID, test_code = BaseNM, test_name = ComponentNM, external_name = ExternalNM, 
         common_name = ComponentCommonNM, result_raw = ResultTXT, comment = ComponentCommentTXT) %>%
  mutate(result_raw = str_trim(result_raw)) %>%
  mutate(result_raw = ifelse(result_raw %in% c("", " "), NA, result_raw))

# Date parse
edw.cdiff <- edw.cdiff %>%
  mutate(coll_datetime = as_date(orderdate, origin = "1840-12-31")) %>%
  select(PatientID:orderdate, coll_datetime, result_datetime:comment)
  
# Map clean results
cdiff.results.map <- read_csv(file = "/data/tide/data/mappings/edw/EDW_Cdiff_results_map_20240713.csv") %>%
  select(-n) %>%
  distinct()
  
edw.cdiff <- edw.cdiff %>%
  left_join(cdiff.results.map)

# Store results for mapping
cdiff.results.map <- count(cdiff.results.map, result_raw, result, drop) %>% arrange(-n)

write_csv(cdiff.results.map, file = "/data/tide/data/mappings/edw/EDW_Cdiff_results_map_temp.csv", na = "")

# Filter out specimens without a result
edw.cdiff <- edw.cdiff %>%
  filter(is.na(drop)) %>%
  select(-drop) %>% 
  distinct()

#### SAVE FILE ####
end.date <- gsub("-", "", range(edw.cdiff$coll_datetime)[[2]])

write_csv(edw.cdiff, file = paste0("/data/tide/projects/ho_infxn_ml/input_data/MGB_Cdiff_", end.date, ".csv"), na = "")

#### CLEAN UP ####
close(edwDB)

rm(edwDB, query_date, edw.cdiff.raw, cdiff.results.map)
