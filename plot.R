#' ---
#'title: "plateplotr"
#'author: "Jack Stevenson"
#'date: "2022"
#' ---
#'plateplotr generates dose-response plots from plate reader data.
#' 
#'Input: CSVs formatted as per the plater package
#' 
#'Variables expected in import files (import template is already set up for this):
#' 
#'- 'compound': name of compound used (including vehicle wells in dilution series)
#'- 'conc_uM': concentration of compound used, in µM
#'- 'target': target of treatment, e.g. cell line or purified enzyme
#'- 'read_norm': normalized plate reader data (calculated in the spreadsheet)
#'- 'replicate': replicate of curve (included for possible future QC)
#'
#'How to use plateplotr:
#'
#'1. Make one copy of the import template for each plate of data you want to import
#'2. Fill out each sheet with the compound + concs used (check that the layout is correct for your experiment)
#'3. Paste corresponding raw plate reader data into each sheet
#'4. Save input sheets as CSVs
#'5. Copy input sheets into a directory called "Input CSVs" in the same directory as this script
#'(or set the working directory to the directory that contains "Input CSVs/")
#'6. Run plot.R (this file)
#'7. Output plots will be created in a "Plots Output" folder in the working directory

# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(plater)  # for tidy importing of plate data
library(viridis) # for color schemes
library(patchwork) # for plot organization

# specify names of input files and import data---------------------------------
# note the order compounds are imported is the order they will be plotted
input_directory <- "Input CSVs/"
plot_type <- "pdf"
compounds <- c(
  "ponatinib",
  #"PonatiLink-1",
  #"ponatinib + asciminib",
  #"dasatinib",
  #"dasatinib + asciminib",
  #"asciminib",
  #"DasatiLink-1",
  #"DasatiLink-2",
  #"DasatiLink-3",
  #"DasatiLink-4"
  #"PonatiLink-1-12",
  #"PonatiLink-1-16",
  "PonatiLink-1-20",
  "PonatiLink-1-24",
  "PonatiLink-1-28",
  "PonatiLink-2-7-8"
)
plate_filenames <- c(list.files(input_directory, pattern = "*.csv")) #gathers all .csv in directory
plate_paths <- paste0(input_directory, plate_filenames)
plate_names <- seq(1,length(plate_filenames))  # create plate IDs
plate_data <- read_plates(plate_paths, plate_names) %>% # import with plater
  filter(compound != "N/A") %>% # drop empty wells
  filter(compound %in% compounds) %>%
  # drop 0 values before plotting and curve fitting
  # note this is only OK because normalization happens before import
  filter(conc_uM != 0) %>%
  mutate(log.conc = log10(conc_uM/1e6))  # convert conc µM to M and log transform
# generate global parameters for all plots------------------------------------------
targets <- distinct(plate_data["target"])$target
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plate_data$log.conc))
x_max <- ceiling(max(plate_data$log.conc))
x_limits <- c(x_min, x_max)
# create logistic minor breaks for all compounds
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# set factors so targets get plotted and colored in input order
target_factors <- distinct(plate_data, target)$target
compound_factors <- compounds # manual factor order for compounds
plate_data <- plate_data %>% 
  mutate(target = fct_relevel(target, target_factors)) %>%
  mutate(compound = fct_relevel(compound, compound_factors))
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
            sem = sd(read_norm, na.rm = TRUE)/sqrt(n()),
            # get mean normalized readout value for plotting
            mean_read = mean(read_norm),
            w = 0.1 * n() # necessary for consistent error bar widths across plots
  )
}
# helper function to add ggplot objects common to all plots--------------------
plot_global <- function(plot){
  plot +
    scale_x_continuous(guide = "prism_offset_minor", # end at last tick
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
get_EC <- function(cpd, t, EC_threshold){
    cpd_data <- plate_data %>%
      filter(compound == cpd, target == t)
    EC <- ED(drm(read_norm~log.conc, data=cpd_data, fct=L.4()), EC_threshold)[1,1]
    return(EC)
}
EC_summary <- plate_data %>%
  group_by(compound, target) %>%
  summarize(
    EC50_nM = 10^get_EC(compound, target, 50) * 1e9, # convert M to nM
    # for negative-response data like this, the EC75 is the drop to 25%
    EC75_nM = 10^get_EC(compound, target, 25) * 1e9
  )
write_csv(EC_summary, "plots output/EC_summary.csv")
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
                se = FALSE, size = 1)} %>%
    plot_global() +
    theme(aspect.ratio = 1) +
    scale_color_manual(values = c("black","darkred")) +
    labs(title = cpd)
}
# plot data for each compound separately----------------------------------------
for (cpd in compounds){
  plot_compound(cpd)
  # save plot with manually optimized aspect ratio
  save_plot(str_glue("plots output/{cpd}.{plot_type}"), legend_len = longest(targets))
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
# save_plot(str_glue("plots output/compound_facets.{plot_type}"), ncol = cols, nrow = rows, legend_len = longest(targets))
# set color parameters for overlaid plots--------------------------------------
alpha_val <- 1
color_scale <- "viridis"
color_start <- .95
color_end <- 0
# plot data for each target separately-------------------------------------------------
for (t in targets){
  plate_summary <- plate_data %>%
    filter(target == t) %>%
    group_by(compound, log.conc) %>% # group into replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = compound)) +
      geom_point(aes(shape = compound), size = 4) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
      # use drm method from drc package to fit dose response curve
      geom_line(aes(linetype = compound),
                stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, size = 1, alpha = alpha_val)} %>%
    plot_global() +
    scale_color_viridis(option = color_scale, discrete = TRUE, begin = color_start, end = color_end) +
    labs(title = t)
  save_plot(str_glue("plots output/{t}.{plot_type}"), legend_len = longest(compounds))
}
# plot data for all targets at once-----------------------------------------
plate_summary <- plate_data %>%
  group_by(target, compound, log.conc) %>% # group into replicates for each condition
  plate_summarize()
{ggplot(plate_summary,aes(x = log.conc, y = mean_read, color = compound)) +
    geom_point() +
    # error bars = mean plus or minus standard error
    geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
    # use drm method from drc package to fit dose response curve
    geom_line(aes(linetype = target), stat = "smooth", method = "drm", method.args = list(fct = L.4()),
              se = FALSE, size = 1, alpha = alpha_val)} %>%
  plot_global() +
  scale_color_viridis(option = color_scale, discrete = TRUE, begin = color_start, end = color_end) +
  labs(title = "All data")
save_plot(str_glue("Plots Output/all_data.{plot_type}"), legend_len = longest(append(targets, compounds)))
