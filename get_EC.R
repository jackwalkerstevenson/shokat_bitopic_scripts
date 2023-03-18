get_EC_nM <- function(data, cpd, tgt, EC_threshold){
  cpd_data <- data %>%
    filter(compound == cpd, target == tgt)
  # fit 4-param logistic model, then extract EC value
  EC_log <- ED(drm(activity~log.conc, data=cpd_data, fct=L.4()),
           EC_threshold, display = FALSE)[1,1]
  return(10^EC_log*1e9) # 1e9 for nM
}