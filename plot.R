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
library(plater)  # for tidy importing of plate data
library(viridis) # for color schemes
library(patchwork) # for plot organization
library(doseplotr) # you bet

# import data---------------------------------
# the order of the treatment list is the order they will be plotted
source("parameters/treatments.R") # import list of treatments to include in plots
# create input and output directories, since git doesn't track empty directories
dir.create("input/", showWarnings = FALSE) # do nothing if directory already exists
dir.create("output/", showWarnings = FALSE)
input_directory <- "input/"
plot_type <- "pdf"
# get excel filenames that start with an alphanumeric character (no hidden files)
excel_filenames <- c(list.files(input_directory, pattern = "^[[:alnum:]].*.xls"))
# convert any excel files to csv so plater can import them
excel_to_csv <- function(filename){
  excel_path <- paste0(input_directory, filename) # full file path
  excel_data <- read_excel(excel_path) # get data from excel file
  # create a new file path that replaces the excel extension with csv
  csv_path <- paste0(input_directory, tools::file_path_sans_ext(filename), ".csv")
  write_csv(excel_data, file = csv_path) # write data to new csv file
}
map(excel_filenames, excel_to_csv)
# gather all CSVs in directory and import them
plate_filenames <- c(list.files(input_directory, pattern = "*.csv")) # get file names
plate_paths <- paste0(input_directory, plate_filenames) # full file paths
plate_names <- seq(1,length(plate_filenames))  # create sequential plate IDs
plate_data <- read_plates(plate_paths, plate_names) %>% # import with plater
  rename(activity = read_norm) %>% # seems clearer for now to have "read_norm" in import template
  # backward compatibility for "compound" in imports. rename if present
  rename(any_of(c(treatment = "compound"))) %>%
  filter(treatment != "N/A") %>% # drop empty wells
  filter(treatment %in% treatments) %>% # take only specified treatments
  make_log_conc %>% # convert conc_nM or conc_uM to log molar concentration
  # drop 0 concs before plotting and curve fitting
  # note this is only OK because normalization happens before import
  filter(conc_logM != -Inf)
# assert that all treatments listed are actually present in imported data
imported_treatments <- distinct(plate_data["treatment"])$treatment
for(treatment in treatments){
  assert_that(treatment %in% imported_treatments,
  msg = glue::glue("treatment '{treatment}' from the list of treatments to plot was not found in imported data"))
}
# generate global parameters for all plots------------------------------------------
targets <- distinct(plate_data["target"])$target
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plate_data$conc_logM))
x_max <- ceiling(max(plate_data$conc_logM))
x_limits <- c(x_min, x_max)
# create logistic minor breaks for all treatments
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# set factors so targets get plotted and colored in input order
target_factors <- distinct(plate_data, target)$target # targets in order of appearance in data
plate_data <- plate_data %>% 
  mutate(target = fct_relevel(target, target_factors)) %>%
  mutate(treatment = fct_relevel(treatment, treatments)) # treatments in order of input list
# set default font size for plots
font_base_size <- 14 # 14 is theme_prism default
# set default point size for plots
pt_size = 3
# import helper function for summarizing plate data
source("summarize_activity.R")
# helper function for saving plots----------------------------------------------
scale_facet <- 4.5 # plot width per col/height per row
legend_pad <- 0.3 # extra width for legend icon
text_factor <- font_base_size / 120 # approx width per character of longest legend text
save_plot <- function(filename, legend_len = 0, nrow = 1, ncol = 1, width = 0, height = 0, ...){
  # if width is not provided, calculate width from length of legend text
  if(width == 0){
    width <- ncol * scale_facet + legend_pad + legend_len * text_factor
  }
  if(height == 0){
    height <- nrow * scale_facet
  }
  ggsave(filename, bg = "transparent", width = width, height = height, ...)}
# helper function for getting length of longest string in a list---------------
longest <- function(strings){
  lengths <- map(strings, nchar)
  lengths[[which.max(lengths)]]
}
# helper function to add ggplot objects common to all plots--------------------
plot_global <- function(plot){
  plot +
    scale_x_continuous(guide = "prism_offset_minor", # end at last tick
                       breaks = scales::breaks_width(1),
                       minor_breaks = minor_x) + # manual minor ticks
    scale_y_continuous(guide = "prism_offset",  # end at last tick
                       breaks = c(0,25,50,75,100)) + # manual y axis ticks
    coord_cartesian(xlim = x_limits, # set x axis zoom from global values
                    ylim = c(0,NA)) + # set y axis zoom locally
    theme_prism(base_size = font_base_size) + # make it look fancy like prism
    theme(plot.background = element_blank()) + # need for transparent background
    labs(x = "log [compound] (M)",
         y = "relative cell viability (%)")
}
# fit models to output EC values------------------------------------------------
# seems like you should be able to just pipe group_by into drm(), but nope, so doing this instead
# helper function for getting EC for one treatment, target, and EC threshold
source("get_EC.R")
source("get_hill_slope.R")
EC_summary <- plate_data %>%
  group_by(treatment, target) %>%
  summarize(
    EC50_nM = get_EC_nM(plate_data, treatment, target, 50), # convert M to nM
    # for negative-response data like this, the EC75 is the drop to 25%
    EC75_nM = get_EC_nM(plate_data, treatment, target, 25), # convert M to nM
    hill_slope = get_hill_slope(plate_data, treatment, target)
  )
write_csv(EC_summary, "output/EC_summary.csv")
# set parameters for treatment plots--------------------------------------------
vr <- viridis_range(length(targets))
viridis_begin <- vr[1]
viridis_end <- vr[2]
# helper function to plot one treatment----------------------------------------
plot_treatment <- function(trt){
  plate_summary <- plate_data %>%
    filter(treatment == trt) %>% # get data from one treatment to work with
    group_by(target, conc_logM) %>%  # get set of replicates for each condition
    summarize_activity()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = conc_logM, y = mean_read, color = target)) +
      geom_point(aes(shape = target), size = pt_size) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w)) +
      # second error bars for 95% CI
      # geom_errorbar(aes(ymin = ymin_activity, ymax = ymax_activity, width = w), alpha = 0.4) +
      # use drm method from drc package to fit dose response curve
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1)} %>%
    plot_global() +
    # theme(aspect.ratio = 1) +
    scale_color_viridis(discrete = TRUE, begin = viridis_begin, end = viridis_end) +
    #scale_color_manual(values = c("black","darkred")) +
    labs(title = trt)
}
# plot data for each treatment separately----------------------------------------
for (trt in treatments){
  plot_treatment(trt)
  # save plot with manually optimized aspect ratio
  save_plot(str_glue("output/{trt}.{plot_type}"), legend_len = longest(targets))
}  
# plot data for all treatments in facets----------------------------------
# treatment_plots = list()
# for (trt in treatments){
#   treatment_plots <- append(treatment_plots, list(plot_treatment(trt)))
# }
# 
# plot_mar <- 15 # margin between wrapped plots, in points
# cols = 4
# rows = 2
# wrap_plots(treatment_plots, guides = "collect", ncol = cols, nrow = rows) &
#   theme(plot.margin = unit(c(plot_mar,plot_mar,plot_mar,plot_mar), "pt"),
#         plot.background = element_blank(),
#         legend.text= element_text(face = "bold", size = 16))
# save_plot(str_glue("output/treatment_facets.{plot_type}"), ncol = cols, nrow = rows, legend_len = longest(targets))
# set color parameters for target plots--------------------------------------
alpha_val <- 1
color_scale <- "viridis"
vr <- viridis_range(length(treatments))
viridis_begin <- vr[1]
viridis_end <- vr[2]
grey_start <- 0.7
grey_end <- 0
# plot data for each target separately------------------------------------------
for (t in targets){
  plate_summary <- plate_data %>%
    filter(target == t) %>%
    group_by(treatment, conc_logM) %>% # group into replicates for each condition
    summarize_activity()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = conc_logM, y = mean_read, color = treatment)) +
      geom_point(aes(shape = treatment), size = pt_size) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
      # second error bars for 95% CI
      # geom_errorbar(aes(ymin = ymin_activity, ymax = ymax_activity, width = w), alpha = 0.4) +
      # use drm method from drc package to fit dose response curve
      geom_line(#aes(linetype = treatment),  # linetype for better grayscale
                stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1, alpha = alpha_val)} %>%
    plot_global() +
    #scale_color_grey(start = grey_start, end = grey_end) +
    scale_color_viridis(option = color_scale, discrete = TRUE, begin = viridis_begin, end = viridis_end) +
    labs(title = t)
  save_plot(str_glue("output/{t}.{plot_type}"), legend_len = longest(treatments))
}
# plot data for all targets at once-----------------------------------------
plate_summary <- plate_data %>%
  group_by(target, treatment, conc_logM) %>% # group into replicates for each condition
  summarize_activity()
{ggplot(plate_summary,aes(x = conc_logM, y = mean_read, color = treatment)) +
    geom_point(aes(shape = treatment), size = pt_size) +
    # error bars = mean plus or minus standard error
    geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
    # use drm method from drc package to fit dose response curve
    geom_line(aes(linetype = target), stat = "smooth", method = "drm", method.args = list(fct = L.4()),
              se = FALSE, linewidth = 1, alpha = alpha_val)} %>%
  plot_global() +
  scale_color_viridis(option = color_scale, discrete = TRUE, begin = viridis_begin, end = viridis_end) +
  labs(title = "All data")
save_plot(str_glue("output/all_data.{plot_type}"), legend_len = longest(append(targets, treatments)))
