get_hill_slope <- function(data, cpd, tgt){
  cpd_data <- data %>%
    filter(compound == cpd, target == tgt)
  slope <- coef(drm(read_norm~log.conc, data=cpd_data, fct=L.4()))[1]
  return(slope)
}