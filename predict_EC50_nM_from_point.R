#' ---
#'title: "predict EC50"
#'author: "Jack Stevenson"
#'date: "2023-03-31"
#' ---
#'function to predict the EC50 of a dose-response curve given fixed other curve parameters
#'parameters:
#'- logmolar_dose: dose (single or vector) in logmolar units (1 nM = -9 logmolar)
#'- pct_inhibition: percent inhibition (single or vector) corresponding to doses
predict_EC50_nM_from_point <- function(logmolar_dose, pct_inhibition, fixed_params){
  # assertthat::assert_that(length(logmolar_dose) == length(pct_inhibition))
  # convert pct inhibition to activity, since prefit model is activity~log.conc
  act <- 100 - pct_inhibition
  print(str_glue("predicting EC50 from {pct_inhibition}% inhibition at logmolar conc {logmolar_dose}"))
  x <- logmolar_dose
  y <- act
  b <- fixed_params["hill_slope"]
  c <- fixed_params["max_effect"]
  d <- fixed_params["min_effect"]
  # rearranged 4-param logistic to solve for EC50
  # -b hacked in for negative hill slope
  EC50_log <- unname(x/(((d-y)/(y-c))^(1/-b)))
  EC50_nM <- 10^EC50_log*1e9
  return(EC50_nM) # convert EC50 to nM
}