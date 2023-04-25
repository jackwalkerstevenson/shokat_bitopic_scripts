get_hill_slope <- function(data, trt, tgt){
  cpd_data <- data %>%
    filter(treatment == trt, target == tgt)
  slope <- coef(drm(activity~conc_logM, data=cpd_data, fct=L.4()))[1]
  return(slope)
}