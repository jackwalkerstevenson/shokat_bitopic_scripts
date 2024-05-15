# generate a list of KinMap style directives from SelectScreen single-point data
# Jack Stevenson started 2024-05
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr)
# set aesthetic parameters for KinMap output-----------------------------------
kinmap_global_shape <- 0
kinmap_ontarget_fill_color <- "green"
kinmap_offtarget_fill_color <- "red"
kinmap_uninhibited_fill_color <- "gray"
kinmap_global_stroke_color <- "black"
kinmap_global_stroke_width <- 1.5
on_targets <- c("ABL1")
# import selectscreen results--------------------------------------------------
SS_single_pt_filename <-
  "input/selectscreen combined results 67309 2024-05-14.csv"
single_pt_data <- import_selectscreen(SS_single_pt_filename)
# add KinMap directives that are constant--------------------------------------
single_pt_data <- single_pt_data |> 
  # color on vs off target
  mutate(kinmap_fill_color = case_when(
    pct_inhibition < 10 ~ kinmap_uninhibited_fill_color,
    target %in% on_targets ~ kinmap_ontarget_fill_color,
    .default = kinmap_offtarget_fill_color),
    # set global aesthetic parameters
    kinmap_stroke_color = kinmap_global_stroke_color,
    kinmap_stroke_width = kinmap_global_stroke_width,
    # scale size by pct_inhibition with floor
    kinmap_size = ifelse(pct_inhibition < 10, 10, pct_inhibition)
         )