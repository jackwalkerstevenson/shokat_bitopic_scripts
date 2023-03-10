get_EC <- function(data, cpd, tgt, EC_threshold){
  cpd_data <- data %>%
    filter(compound == cpd, target == tgt)
  # fit 4-param logistic model, then extract EC value
  EC <- ED(drm(activity~log.conc, data=cpd_data, fct=L.4()), EC_threshold, display = FALSE)[1,1]
  return(EC)
}