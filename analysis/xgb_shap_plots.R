library(SHAPforxgboost)
library(xgboost)
library(dplyr)
conflicts_prefer(dplyr::filter())

data <- read.csv('/data/tide/projects/ho_infxn_ml/clean_data/20241210/age_filtered_final_dataset_for_models_20241210.csv')

all_runs = unique(data$run)

for (run_name in all_runs){
  model = xgb.load(paste0('/data/tide/projects/ho_infxn_ml/results/model_results/20241210/xgb/model_checkpoints/environmental_', run_name ,'/fold_1.model'))
  
  dat.run <- data %>% filter(run == run_name & match == 'environmental')
  dat.run <- dat.run %>% select(CDiff_cp:DR_PsA_cp)
  dat.run <- dat.run %>% select(where(~n_distinct(.) > 1)) # Remove features with only 1 value among samples
  dat.run <- as.matrix(dat.run)
  
  shap_values <- shap.values(xgb_model = model, X_train = dat.run)
  shap_values <- shap_values$shap_score
  
  shap_long <- shap.prep(xgb_model = model, X_train = dat.run)
  shap_long <- shap.prep(shap_contrib = shap_values, X_train = dat.run)
  
  plot <- shap.plot.summary(shap_long)
  
  ggsave(filename = paste0('/data/tide/projects/ho_infxn_ml/results/model_results/20241210/xgb/shap_plots/environmental_', run_name,'_shap_summary.pdf'), 
         plot = plot, 
         width = 8, 
         height = 6)
  
  ggsave(filename = paste0('~/ho_infxn_clean_code/results/model_results/xgb/shap_plots/environmental_', run_name,'_shap_summary.pdf'), 
         plot = plot, 
         width = 8, 
         height = 6)
}


