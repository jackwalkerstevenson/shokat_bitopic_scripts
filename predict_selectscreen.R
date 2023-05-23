# predict the EC50 of a treatment in selectscreen assuming model parameters
# created 2023-05-22
library(dplyr)
library(doseplotr)
SS_curve_filename <- "ZLYTE_compiled_results_complete.csv"
SS_single_pt_filename <- "ZLYTE_compiled_results_single_point.csv"
SS_data <- import_selectscreen(SS_curve_filename)
pt_data <- import_selectscreen(SS_single_pt_filename)
P12_ABL1_SS_model <- SS_data |>
  filter_trt_tgt("PonatiLink-2-7-12", "ABL1") |> 
  get_drda()
P12_FGFR1_SS_model <- SS_data |>
  filter_trt_tgt("PonatiLink-2-7-12", "FGFR1") |> 
  get_drda()
#' Estimate SelectScreen EC50 from one point measurement from pona/asc mode
#'
#' This function uses the parameters of a model fit by drda for the activity of
#' ponatinib + asciminib on ABL1 wt to estimate the SelectScreen EC50 from a
#' single point dose-response measurement.
#' @param conc_nM Concentration at which single point was measured
#' @param pct_inhibition Percent inhibition observed at single point
#'
#' @return The estimated EC50 in nanomolar units.
est_SS_EC50 <- function(model, conc_nM, pct_inhibition){
  EC50_nM_from_point_model(model,
                           conc_nM,
                           100 - pct_inhibition)
}
# wrapper to estimate EC50s using PL-2-7-12 on ABL1, which is normal-shaped
est_SS_EC50_P12_ABL1 <- function(conc_nM, pct_inhibition){
  est_SS_EC50(P12_ABL1_SS_model, conc_nM, pct_inhibition)
}
# wrapper to estimate EC50s using PL-2-7-12 on FGFR1, which is shallower
est_SS_EC50_P12_FGFR1 <- function(conc_nM, pct_inhibition){
  est_SS_EC50(P12_FGFR1_SS_model, conc_nM, pct_inhibition)
}
# predict EC50s for PL-2-7-10 on all single-pt kinases using both models
predict_all_EC50s_from_points <- function(trt){
  pt_data |> 
    filter(treatment == trt) |> 
    group_by(target) |>
    summarize(conc_nM = mean(Compound.Conc), # should only be one conc
              mean_inhib = mean(pct_inhibition)) |> 
    mutate(EC50_pred_P12_ABL1 = est_SS_EC50_P12_ABL1(conc_nM, mean_inhib),
           EC50_pred_P12_FGFR1 = est_SS_EC50_P12_FGFR1(conc_nM, mean_inhib),
           disregard = mean_inhib < 5 | mean_inhib > 95)
}

pt_EC50_predictions_P10 <- predict_all_EC50s_from_points("PonatiLink-2-7-10")
pt_EC50_predictions_PA <- predict_all_EC50s_from_points("ponatinib + asciminib")
pt_EC50_predictions_P24 <- predict_all_EC50s_from_points("PonatiLink-1-24")
