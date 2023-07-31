#' ---
#'title: "plateplotr"
#'author: "Jack Stevenson"
#'date: "2023"
#' ---
#'plateplotr generates dose-response plots from plate reader data.
#' 
#'Input: platemap spreadsheets, Excel or CSV, formatted as per the plater package
#'
#'An import template is provided to make it easy to make platemaps.
#' 
#'Variables expected in input platemaps:
#' 
#'- 'treatment': name of treatment used
#'
#'note: if a dilution series includes a vehicle (zero-concentration) well, it should be labeled as the same treatment
#'  
#'- 'conc_nM' or 'conc_uM': concentration of treatment used in µM or nM
#'- 'target': target of treatment, e.g. cell line or purified protein
#'- 'read_norm': normalized plate reader data (calculated in the platemap)
#'- 'replicate': replicate of condition (included for possible future QC)
#'
#'How to use plateplotr:
#'
#'1. Fill in "treatments.R" with the list of treatments you want to plot
#'2. Make a copy of the import platemap for each plate of data you want to import
#'3. Fill out each platemap with treatment(s), target(s) (e.g. cell line) and concentrations used
#'4. Copy the corresponding raw plate reader data into each platemap
#'5. Copy input platemaps into a directory called "input" in the same directory as this script
#'6. Run plot.R
#'7. Output plots will be created in an "output" folder in the working directory

# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(readxl) # for excel file handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(scales) # for axis transforms
library(plater)  # for tidy importing of plate data
library(viridis) # for color schemes
library(patchwork) # for plot organization
library(doseplotr) # you bet

# prepare global variables---------------------------------
# the order of the treatment list is the order they will be plotted
# source("parameters/treatments.R") # import list of treatments to include in plots
source("parameters/targets.R") # import list of targets to include in plots
input_directory <- "input/" # path to directory containing input files
output_directory <- "output/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # # font size for plots. 14 is theme_prism default
pt_size = 3 # point size for plots
no_legend <- FALSE # whether all plots should have no legend
global_x_lim <- TRUE # whether all plots should use the same x limits
rigid <- FALSE # whether to use rigid low-dose asymptote
# filename to use if importing data from a single file instead of a directory
input_filename <- "input/2023-07-03 Ivan raw data names edited.csv"
# import and preprocess data----------------------------------------------------
# create input and output directories, since git doesn't track empty directories
dir.create(input_directory, showWarnings = FALSE) # do nothing if directory already exists
dir.create(output_directory, showWarnings = FALSE)
if(exists("input_filename")){
  plate_data <- readr::read_csv(input_filename) |>
    preprocess_plate_data()
} else {plate_data <- import_plates(input_directory)}
# filter and validate imported data---------------------------------------------
filter_validate_reorder <- function(data, colname, values){
  data <- data |> 
    dplyr::filter({{colname}} %in% values)
  data_values <- unique(data[[colname]])
  for(value in values){
    assertthat::assert_that(value %in% data_values, msg = glue::glue(
      "value {value} was not found in the data"))}
  data |>
    dplyr::mutate(colname =
                    forcats::fct_relevel(.data[[colname]], values))
}
# temporary script before replacing with function-------------------------------
if(exists("treatments")){
  plate_data <- plate_data |> 
  filter(treatment %in% treatments) # take only specified treatments
    # assert that all treatments listed are actually present in imported data
  imported_treatments <- unique(plate_data$treatment)
  for(treatment in treatments){
    assert_that(treatment %in% imported_treatments,
                msg = str_glue("treatment '{treatment}' from the list of ",
                "treatments to plot was not found in imported data"))
  }
  plate_data <- plate_data |> 
    mutate(treatment = fct_relevel(treatment, # relevel treatments by input list
                                   treatments))
}
if(exists("targets")){
  plate_data <- plate_data |> 
    filter(target %in% targets) # take only specified targets
  # assert that all targets listed are actually present in imported data
  imported_targets <- unique(plate_data$target)
  for(target in targets){
    assert_that(target %in% imported_targets,
                msg = str_glue("target '{target}' from the list of ",
                               "targets to plot was not found in imported data"))
  }
  plate_data <- plate_data |> 
    mutate(target = fct_relevel(target, # relevel targets by input list
                                targets))
}
plot_data <- plate_data |> filter(log_dose != -Inf)
# generate data-dependent global plot parameters--------------------------------
if (!exists("treatments")){ # if treatments not specified, use all treatments
  treatments <- as.vector(unique(plot_data$treatment))}
if (!exists("targets")){ # if targets not specified, use all targets
  targets <- as.vector(unique(plot_data$target))}
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plot_data$log_dose))
x_max <- ceiling(max(plot_data$log_dose))
x_limits <- c(x_min, x_max)
# x_limits <- c(-11,-5) # manual x limit backup
# create logistic minor breaks for all treatments
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# fit models and report model parameters----------------------------------------
model_summary <- summarize_models(plot_data,
                                  response_col = "response_norm",
                                  rigid = rigid) |> # use global rigid parameter
  select(-model) |>  # remove actual model from report
  mutate(across(where(is.numeric), \(x){signif(x, digits = 4)}))
write_csv(model_summary,
          str_glue("output/plate_model_summary_{get_timestamp()}.csv"))
# plot untreated data by target for QC------------------------------------------
# set color parameters for treatments
vr <- viridis_range(length(treatments))
vr_begin <- vr[[1]]
vr_end <- vr[[2]]
vr_option <- vr[[3]] # color scale
untreated_data_by_tgt <-  plate_data |>
  filter(log_dose == -Inf) |> # untreated data only
  group_by(target) |>
  summarize_response(response_col = "response")
untreated_data_by_tgt_and_trt <- plate_data |> 
  filter(log_dose == -Inf) |> # untreated data only
  group_by(target, treatment) |>
  summarize_response(response_col = "response")
p <- untreated_data_by_tgt_and_trt |> 
  ggplot(aes(y = target, x = mean_response, color = treatment)) +
  geom_point(size = pt_size, alpha = 0.8) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(labels = \(x) format(x, scientific = TRUE)) +
  scale_color_viridis(option = vr_option, discrete = TRUE,
                      begin = vr_begin, end = vr_end) +
  geom_errorbar(aes(xmax = mean_response+sem, xmin = mean_response-sem,
                    width = w)) +
  theme_prism() +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "CTG luminescense",
       title = "Growth of untreated control wells") +
  geom_point(data = untreated_data_by_tgt,
           size = 6, alpha = 0.3, color = "black") +
  geom_text(data = untreated_data_by_tgt, # label average points with value
            aes(label = mean_response |>
                  signif(digits = 3) |>
                  format(scientific = TRUE)),
            show.legend = FALSE, # no legend for the text labels
            position = position_nudge(y = .3)) +
  guides(color = guide_legend(override.aes = list(size = pt_size)))
save_plot(p, str_glue("output/plate_QC_untreated_{get_timestamp()}.{plot_type}"),
          width = 14, height = 8)
# plot data for each treatment separately----------------------------------------
for (trt in treatments){
  trt_targets <- as.vector(unique((plot_data |> filter_trt_tgt(trt = trt))$target))
  plot_treatment(plot_data, trt, rigid = rigid,
                 if(global_x_lim){x_limits = x_limits},
                 response_col = "response_norm") |>
    save_plot(
      str_glue("output/plate_treatment_{trt}_{get_timestamp()}.{plot_type}"),
      legend_len = longest(trt_targets))
}  
# plot data for each target separately------------------------------------------
for (tgt in targets){
  tgt_treatments <- as.vector(unique((plot_data |> filter_trt_tgt(tgt = tgt))$treatment))
  # tgt_treatments_test <- unique(plot_data$treatment)
  plot_target(plot_data, tgt, rigid = rigid,
                 if(global_x_lim){x_limits = x_limits},
                 response_col = "response_norm") |>
    save_plot(
      str_glue("output/plate_target_{tgt}_{get_timestamp()}.{plot_type}"),
      # legend_len = longest(treatments))
      legend_len = longest(tgt_treatments))
}
