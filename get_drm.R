# get the drc dose-response model for a combination of treatment and target
get_drm <- function(data, trt, tgt){
  cpd_data <- data %>%
    filter(treatment == trt, target == tgt)
  # 4-param logistic model on pre-log-transformed data
  return(drm(activity~log.conc, data=cpd_data, fct=L.4()))
}