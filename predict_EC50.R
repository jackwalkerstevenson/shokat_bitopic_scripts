#' ---
#'title: "predict EC50s"
#'author: "Jack Stevenson"
#'date: "2023-03-31"
#' ---
#'Predict the EC50 of a treatment on a target given a single point measurement
#'Assume curve shape is similar to PonatiLink compounds on Abl
#'- import a SelectScreen single-point dataset for the treatments
#'- fit a dose-response model in which only the EC50 is allowed to vary
#'- output a CSV with the resulting EC50s
# set up------------------------------------------------------------------------
library(tidyverse) # for tidy data handling
library(drc)  # for dose response curves
# import helper functions
source("import_selectscreen.R")
source("nM_to_logmolar.R")
source("logmolar_to_nM.R")
source("get_EC.R")
source("predict_EC50_nM_from_point.R")
dir.create("output/", showWarnings = FALSE) # silently create output directory
# set global variables---------------------------------------------------------
source("parameters/treatments.R")
# 4-param logistic params for PonatiLink-2-7-10 on ABL1, activity~log.conc
P10_params <- c(hill_slope=2.856, max_effect=-2.762, min_effect=98.773, EC50=-9.244)
P10_params_pct <- c(hill_slope=-2.856, max_effect=1.227, min_effect=102.762, EC50=-9.244)
# import data------------------------------------------------------------------
input_filename <- "ZLYTE_compiled_results_single_point.csv"
point_data <- import_selectscreen(input_filename, treatments)
# predict some EC50s from fake data to test------------------------------------
# inhibition EC50 should give
P10_midpt <- 100- mean(c(P10_params["max_effect"], P10_params["min_effect"]))
P10_midpt_pct <- 100- mean(c(P10_params_pct["max_effect"], P10_params_pct["min_effect"]))
test_doses_nM <- c(.1, .5, 1, 2, 5)
test_doses_logmolar <- nM_to_logmolar(test_doses_nM)
test_predicted_EC50 <- predict_EC50_nM_from_point(test_doses_logmolar, 70, fixed_params = P10_params_pct)
test_data <- tibble(test_doses_logmolar, test_predicted_EC50)

# test_pcts <- c(10, 20, 30, 40, P10_midpt_pct, 50, 60)
# test_predicted_EC50 <- predict_EC50_nM_from_point(P10_params["EC50"], test_pcts, fixed_params = P10_params)
# test_data <- tibble(test_pcts, test_predicted_EC50)

# ggplot(test_data, aes(x = test_pcts, y = test_predicted_EC50)) +
ggplot(test_data, aes(x = test_doses_logmolar, y = test_predicted_EC50)) +
  # scale_y_log10() +
  # scale_x_log10() +
  geom_point()
#write_csv(output, "output/predicted_treatment_inhibition.csv")