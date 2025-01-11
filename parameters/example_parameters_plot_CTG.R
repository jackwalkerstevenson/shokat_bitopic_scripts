# parameters for plot_CTG
source("parameters/manual_scales.R")
manually_recolor_treatments <- TRUE # whether to use manual colors for treatments
manually_recolor_targets <- TRUE # whether to use manual colors for targets
manually_reshape_treatments <- TRUE # whether to use manual shapes for treatments
manually_reshape_targets <- TRUE # whether to use manual shapes for targets
manually_relabel_treatments <- TRUE # whether to override treatment legend labels
manually_relabel_targets <- TRUE # whether to override target legend labels
override_treatment_title <- FALSE
override_target_title <- TRUE
target_title <- "cell line"
input_directory <- "input/" # path to directory containing input files
output_directory <- "output/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # font size for plots. 14 is theme_prism default
pt_size = 3 # point size for plots
no_legend <- FALSE # whether all plots should have no legend
global_x_lim <- TRUE # whether all plots should use the same x limits
override_x_lim <- FALSE # whether to apply manual x limits
rigid <- FALSE # whether to use rigid low-dose asymptote
grid <- FALSE # whether to plot a background grid
# filename to use if importing data from a single file instead of a directory
input_filename <- "input/example_data_plot_CTG.csv"
# the order of treatments and targets is the order they will be plotted
treatments <- c(
  "ponatinib",
  "asciminib",
  "ponatinib + asciminib"
)
targets <- c(
  "K562 wt",
  "K562 T315I"
)
