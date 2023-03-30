# get the drc dose-response model for a combination of treatment and target
get_drm <- function(data, trt, tgt){
  data_subset <- data %>%
    filter(treatment == trt, target == tgt)
  # 4-param logistic model on pre-log-transformed data
  return(drm(activity~log.conc, data = data_subset, fct = L.4()))
}