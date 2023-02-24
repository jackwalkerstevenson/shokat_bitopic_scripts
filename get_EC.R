get_EC <- function(data, cpd, tgt, EC_threshold){
  cpd_data <- data %>%
    filter(compound == cpd, target == tgt)
  EC <- ED(drm(activity~log.conc, data=cpd_data, fct=L.4()), EC_threshold)[1,1]
  return(EC)
}