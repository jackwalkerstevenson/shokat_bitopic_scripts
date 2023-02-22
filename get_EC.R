get_EC <- function(data, cpd, tgt, EC_threshold){
  cpd_data <- data %>%
    filter(compound == cpd, target == tgt)
  EC <- ED(drm(read_norm~log.conc, data=cpd_data, fct=L.4()), EC_threshold)[1,1]
  return(EC)
}