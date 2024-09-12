# predict the EC50 of a treatment in selectscreen assuming model parameters
# created 2023-05-22
library(dplyr)
library(doseplotr)
# import selectscreen results, both titrations and single-point measurements
SS_curve_filename <- "ZLYTE_compiled_results_complete.csv"
SS_single_pt_filename <- "ZLYTE_single_point.csv"
SS_data <- import_selectscreen(SS_curve_filename)
single_pt_data <- import_selectscreen(SS_single_pt_filename)
# fit models for PonatiLink-2-7-12 on ABL1 (normal-shaped) and FGFR1 (shallower)
P12_ABL1_SS_model <- SS_data |>
  filter_trt_tgt("PonatiLink-2-7-12", "ABL1") |> 
  get_drda()
P12_FGFR1_SS_model <- SS_data |>
  filter_trt_tgt("PonatiLink-2-7-12", "FGFR1") |> 
  get_drda()
#' Estimate SelectScreen EC50 from one point measurement and an existing model
#'
#' This function uses the parameters of a model fit by drda to estimate the EC50
#' from a single point dose-response measurement.
#' @param conc_nM Concentration at which single point was measured
#' @param pct_inhibition Percent inhibition observed at single point
#'
#' @return The estimated EC50 in nanomolar units.
est_SS_EC50 <- function(model, conc_nM, pct_inhibition){
  EC50_nM_from_point_model(model,
                           conc_nM,
                           100 - pct_inhibition)
}
# wrapper to estimate EC50s using the model fit to PL-2-7-12 on ABL1
est_SS_EC50_P12_ABL1 <- function(conc_nM, pct_inhibition){
  est_SS_EC50(P12_ABL1_SS_model, conc_nM, pct_inhibition)
}
# wrapper to estimate EC50s using the model fit to PL-2-7-12 on FGFR1
est_SS_EC50_P12_FGFR1 <- function(conc_nM, pct_inhibition){
  est_SS_EC50(P12_FGFR1_SS_model, conc_nM, pct_inhibition)
}
# predict a treatment's EC50s for all single-pt-tested kinases with each model
predict_trt_EC50s_from_points <- function(trt){
  single_pt_data |> 
    filter(treatment == trt) |> 
    group_by(target, Compound.Conc) |> # predict separately for each conc
    summarize(conc_nM = mean(Compound.Conc),
              mean_inhib = mean(pct_inhibition)) |>
    mutate(EC50_pred_P12_ABL1 = est_SS_EC50_P12_ABL1(conc_nM, mean_inhib),
          EC50_pred_P12_FGFR1 = est_SS_EC50_P12_FGFR1(conc_nM, mean_inhib),
          input_too_extreme = mean_inhib < 5 | mean_inhib > 95)
}

pt_EC50_predictions_P10 <- predict_trt_EC50s_from_points("PonatiLink-2-7-10")
pt_EC50_predictions_PA <- predict_trt_EC50s_from_points("ponatinib + asciminib")
pt_EC50_predictions_P24 <- predict_trt_EC50s_from_points("PonatiLink-1-24")
