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
#'- 'compound': name of compound used
#'
#'note: if a dilution series includes a vehicle (zero-concentration) well, it should be labeled as the same compound
#'
#'- 'conc_nM' or 'conc_uM': concentration of compound used in µM or nM
#'- 'target': target of treatment, e.g. cell line or purified protein
#'- 'read_norm': normalized plate reader data (calculated in the platemap)
#'- 'replicate': replicate of condition (included for possible future QC)
#'
#'How to use plateplotr:
#'
#'1. Fill in "compounds.R" with the list of compounds you want to plot
#'2. Make a copy of the import platemap for each plate of data you want to import
#'3. Fill out each platemap with compound(s), target(s) (e.g. cell line) and concentrations used
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

# import data---------------------------------
# helper function for reading concentration from nM or uM input data
make_log_conc <- function(df){
  tryCatch({ # try to convert from conc_uM
    df %>% mutate(log.conc = log10(conc_uM/1e6))},
    error = function(e){ # if no conc_uM, try to convert from conc_nM
      df %>% mutate(log.conc = log10(conc_nM/1e9))})}
# the order of the compound list is the order they will be plotted
source("compounds.R") # import list of compounds to include in plots
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
  rename(activity = read_norm) %>%
  filter(compound != "N/A") %>% # drop empty wells
  filter(compound %in% compounds) %>% # take only specified compounds
  make_log_conc %>% # convert conc_nM or conc_uM to log molar concentration
  # drop 0 concs before plotting and curve fitting
  # note this is only OK because normalization happens before import
  filter(log.conc != -Inf)
# assert that all compounds listed are actually present in imported data
imported_compounds <- distinct(plate_data["compound"])$compound
for(compound in compounds){
  assert_that(compound %in% imported_compounds,
  msg = glue::glue("compound '{compound}' from the list of compounds to plot was not found in imported data"))
}
# generate global parameters for all plots------------------------------------------
targets <- distinct(plate_data["target"])$target
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plate_data$log.conc))
x_max <- ceiling(max(plate_data$log.conc))
x_limits <- c(x_min, x_max)
# create logistic minor breaks for all compounds
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# set factors so targets get plotted and colored in input order
target_factors <- distinct(plate_data, target)$target # targets in order of appearance in data
plate_data <- plate_data %>% 
  mutate(target = fct_relevel(target, target_factors)) %>%
  mutate(compound = fct_relevel(compound, compounds)) # compounds in order of input list
# set default font size for plots
font_base_size <- 14 # 14 is theme_prism default
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
# helper function for summarizing replicate data for plotting------------------
plate_summarize <- function(x){
  summarize(x,
            # standard error for error bars = standard deviation / square root of n
            sem = sd(activity, na.rm = TRUE)/sqrt(n()),
            # get mean normalized readout value for plotting
            mean_read = mean(activity),
            w = 0.06 * n() # necessary for consistent error bar widths across plots
  )
}
# helper function to add ggplot objects common to all plots--------------------
plot_global <- function(plot){
  plot +
    scale_x_continuous(guide = "prism_offset_minor", # end at last tick
                       breaks = breaks_width(1),
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
# helper function for getting EC for one compound, target, and EC threshold
source("get_EC.R")
source("get_hill_slope.R")
EC_summary <- plate_data %>%
  group_by(compound, target) %>%
  summarize(
    EC50_nM = 10^get_EC(plate_data, compound, target, 50) * 1e9, # convert M to nM
    # for negative-response data like this, the EC75 is the drop to 25%
    EC75_nM = 10^get_EC(plate_data, compound, target, 25) * 1e9,
    hill_slope = get_hill_slope(plate_data, compound, target)
  )
write_csv(EC_summary, "output/EC_summary.csv")
# helper function to plot one compound----------------------------------------
plot_compound <- function(cpd){
  plate_summary <- plate_data %>%
    filter(compound == cpd) %>% # get data from one compound to work with
    group_by(target, log.conc) %>%  # get set of replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = target)) +
      geom_point() +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w)) +
      # use drm method from drc package to fit dose response curve
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1)} %>%
    plot_global() +
    theme(aspect.ratio = 1) +
    scale_color_viridis(discrete = TRUE) +
    #scale_color_manual(values = c("black","darkred")) +
    labs(title = cpd)
}
# plot data for each compound separately----------------------------------------
for (cpd in compounds){
  plot_compound(cpd)
  # save plot with manually optimized aspect ratio
  save_plot(str_glue("output/{cpd}.{plot_type}"), legend_len = longest(targets))
}  
# plot data for all compounds in facets----------------------------------
# compound_plots = list()
# for (cpd in compounds){
#   compound_plots <- append(compound_plots, list(plot_compound(cpd)))
# }
# 
# plot_mar <- 15 # margin between wrapped plots, in points
# cols = 4
# rows = 2
# wrap_plots(compound_plots, guides = "collect", ncol = cols, nrow = rows) &
#   theme(plot.margin = unit(c(plot_mar,plot_mar,plot_mar,plot_mar), "pt"),
#         plot.background = element_blank(),
#         legend.text= element_text(face = "bold", size = 16))
# save_plot(str_glue("output/compound_facets.{plot_type}"), ncol = cols, nrow = rows, legend_len = longest(targets))
# set color parameters for overlaid plots--------------------------------------
alpha_val <- 1
color_scale <- "viridis"
viridis_start <- 1
viridis_end <- .1
grey_start <- 0.7
grey_end <- 0
# plot data for each target separately-------------------------------------------------
for (t in targets){
  plate_summary <- plate_data %>%
    filter(target == t) %>%
    group_by(compound, log.conc) %>% # group into replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = compound)) +
      geom_point(aes(shape = compound), size = 3) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
      # use drm method from drc package to fit dose response curve
      geom_line(#aes(linetype = compound),  # linetype for better grayscale
                stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1, alpha = alpha_val)} %>%
    plot_global() +
    #scale_color_grey(start = grey_start, end = grey_end) +
    scale_color_viridis(option = color_scale, discrete = TRUE, begin = viridis_start, end = viridis_end) +
    labs(title = t)
  save_plot(str_glue("output/{t}.{plot_type}"), legend_len = longest(compounds))
}
# plot data for all targets at once-----------------------------------------
plate_summary <- plate_data %>%
  group_by(target, compound, log.conc) %>% # group into replicates for each condition
  plate_summarize()
{ggplot(plate_summary,aes(x = log.conc, y = mean_read, color = compound)) +
    geom_point(aes(shape = compound), size = 3) +
    # error bars = mean plus or minus standard error
    geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
    # use drm method from drc package to fit dose response curve
    geom_line(aes(linetype = target), stat = "smooth", method = "drm", method.args = list(fct = L.4()),
              se = FALSE, linewidth = 1, alpha = alpha_val)} %>%
  plot_global() +
  scale_color_viridis(option = color_scale, discrete = TRUE, begin = viridis_start, end = viridis_end) +
  labs(title = "All data")
save_plot(str_glue("output/all_data.{plot_type}"), legend_len = longest(append(targets, compounds)))
