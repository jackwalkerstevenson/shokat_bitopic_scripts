#' ---
#'title: "single point comparison"
#'author: "Jack Stevenson"
#'date: "2023-03-30"
#' ---
#'Predict the activity of a set of SelectScreen treatments at the same single dose
#'- import a SelectScreen dose-response dataset for the treatments of interest
#'- fit a dose-response model for each treatment
#'- use the fitted models to predict the activity of each treatment at the given comparison dose
#'- output a CSV with the results
# set up------------------------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
options(dplyr.summarise.inform = FALSE)
# import helper functions
source("import_selectscreen.R")
source("get_drm.R")
source("predict_activity.R")
dir.create("output/", showWarnings = FALSE) # silently create output directory
# set global variables---------------------------------------------------------
source("parameters/treatments.R")
dose_nM <-  1 # dose for which to predict activity of treatments
tgt <- "ABL1" # target for which to predict activity of treatments
# vectorized function to predict activity from targets inside mutate()
predict_targets <- Vectorize(predict_activity, vectorize.args = "trt")
# import data------------------------------------------------------------------
input_filename <- "ZLYTE_compiled_results_complete.csv"
all_data <- import_selectscreen(input_filename, treatments)
# predict activity and write output--------------------------------------------
output <- tibble(treatment = treatments) |>
  mutate(prediction_dose_nM = dose_nM,
         target = tgt,
         predicted_pct_inhibition = 100 - predict_targets(all_data,
                                                          trt = treatment,
                                                          tgt = tgt,
                                                          dose_nM = dose_nM))
write_csv(output, "output/predicted_treatment_inhibition.csv")