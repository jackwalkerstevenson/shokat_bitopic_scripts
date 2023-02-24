get_hill_slope <- function(data, cpd, tgt){
  cpd_data <- data %>%
    filter(compound == cpd, target == tgt)
  slope <- coef(drm(activity~log.conc, data=cpd_data, fct=L.4()))[1]
  return(slope)
}