# fit a drc dose-response model for a treatment on a target
get_drm <- function(data, trt, tgt){
  data_subset <- data %>%
    filter(treatment == trt, target == tgt)
  # 4-param logistic model on pre-log-transformed data
  return(drm(activity~conc_logM, data = data_subset, fct = L.4()))
}

# fit a drc dose-response model for a treatment on a target
get_drm_pct <- function(data, trt, tgt){
  data_subset <- data %>%
    filter(treatment == trt, target == tgt)
  # 4-param logistic model on pre-log-transformed data
  return(drm(pct_inhibition~conc_logM, data = data_subset, fct = L.4()))
}