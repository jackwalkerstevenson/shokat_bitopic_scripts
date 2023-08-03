#' ---
#'title: "plot_selectscreen"
#'author: "Jack Stevenson"
#' ---
#'script for analyzing Thermo SelectScreen data
#'copied updates from plot_CTG.R 2023-08-03
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(doseplotr) # you bet
# set global parameters---------------------------------------------------------
# the order of the treatment list is the order they will be plotted
source("parameters/treatments.R") # import list of treatments to include in plots
source("parameters/targets.R") # import list of targets to include in plots
input_directory <- "input/" # path to directory containing input files
output_directory <- "output/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # # font size for plots. 14 is theme_prism default
pt_size = 3 # point size for plots
no_legend <- FALSE # whether all plots should have no legend
global_x_lim <- TRUE # whether all plots should use the same x limits
rigid <- FALSE # whether to use rigid low-dose asymptote
grid <- TRUE # whether to plot a background grid
input_filename <- "ZLYTE_compiled_results_complete.csv"
plot_type <- "pdf" # file type for saved output plots
# import and tidy data---------------------------------------------------------
# create input and output directories, since git doesn't track empty directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
plot_data <- import_selectscreen(input_filename) |>
  filter_trt_tgt(treatments, targets)
plot_data <- plot_data |> 
  mutate(treatment = fct_relevel(treatment, treatments))
# temporary script before replacing with function-------------------------------
if(exists("treatments")){
  plot_data <- plot_data |> 
    filter(treatment %in% treatments) # take only specified treatments
  # assert that all treatments listed are actually present in imported data
  imported_treatments <- unique(plot_data$treatment)
  for(treatment in treatments){
    assert_that(treatment %in% imported_treatments,
                msg = str_glue("treatment '{treatment}' from the list of ",
                               "treatments to plot was not found in imported data"))
  }
  plot_data <- plot_data |> 
    mutate(treatment = fct_relevel(treatment, # relevel treatments by input list
                                   treatments))
}
if(exists("targets")){
  plot_data <- plot_data |> 
    filter(target %in% targets) # take only specified targets
  # assert that all targets listed are actually present in imported data
  imported_targets <- unique(plot_data$target)
  for(target in targets){
    assert_that(target %in% imported_targets,
                msg = str_glue("target '{target}' from the list of ",
                               "targets to plot was not found in imported data"))
  }
  plot_data <- plot_data |> 
    mutate(target = fct_relevel(target, # relevel targets by input list
                                targets))
}
# generate data-dependent global plot parameters--------------------------------
if (!exists("treatments")){ # if treatments not specified, use all treatments
  treatments <- as.vector(unique(plot_data$treatment))}
if (!exists("targets")){ # if targets not specified, use all targets
  targets <- as.vector(unique(plot_data$target))}
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plot_data$log_dose))
x_max <- ceiling(max(plot_data$log_dose))
x_limits <- c(x_min, x_max)
# x_limits <- c(-12,-4) # manual x limit backup
# fit models to output EC values------------------------------------------------
EC_summary <- summarize_models(plot_data, response_col = "response")
write_csv(EC_summary, str_glue("output/EC_summary_selectscreen_{get_timestamp()}.csv"))
# plot data for each treatment separately----------------------------------------
for (trt in treatments){
  trt_targets <- as.vector(unique((plot_data |> filter_trt_tgt(trt = trt))$target))
  plot_treatment(plot_data, trt, rigid = rigid, grid = grid,
                 if(global_x_lim){x_limits = x_limits},
                 response_col = "response") |>
    save_plot(
      str_glue("output/plate_treatment_{trt}_{get_timestamp()}.{plot_type}"),
      legend_len = longest(trt_targets))
}  
# plot data for each target separately------------------------------------------
for (tgt in targets){ 
  tgt_treatments <- as.vector(unique((plot_data |> filter_trt_tgt(tgt = tgt))$treatment))
  # tgt_treatments_test <- unique(plot_data$treatment)
  plot_target(plot_data, tgt, rigid = rigid, grid = grid,
              if(global_x_lim){x_limits = x_limits},
              response_col = "response") |>
    save_plot(
      str_glue("output/plate_target_{tgt}_{get_timestamp()}.{plot_type}"),
      # legend_len = longest(treatments))
      legend_len = longest(tgt_treatments))
}
