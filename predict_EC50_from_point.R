#' ---
#'title: "predict EC50"
#'author: "Jack Stevenson"
#'date: "2023-03-31"
#' ---
#'function to predict the EC50 of a dose-response curve given fixed other curve parameters
#'parameters:
#'- logmolar_dose: dose (single or vector) in logmolar units (1 nM = -9 logmolar)
#'- pct_inhibition: percent inhibition (single or vector) corresponding to doses
predict_EC50_from_point <- function(logmolar_dose, pct_inhibition, fixed_params){
  # assertthat::assert_that(length(logmolar_dose) == length(pct_inhibition))
  # convert pct inhibition to activity, since prefit model is activity~log.conc
  act <- 100 - pct_inhibition
  x <- logmolar_dose
  y <- act
  b <- fixed_params["hill_slope"]
  c <- fixed_params["max_effect"]
  d <- fixed_params["min_effect"]
  # rearranged 4-param logistic to solve for EC50
  log_EC50 <- x/(((d-y)/(y-c))^(1/b))
  print(str_glue("found log EC50 {log_EC50}"))
  return(10^log_EC50*1e9) # convert EC50 to nM
}