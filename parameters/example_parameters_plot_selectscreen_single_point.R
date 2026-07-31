# parameters for plot_selectscreen_single_point
source("parameters/manual_scales.R")
manually_recolor_treatments <- TRUE # whether to use manual colors for treatments
manually_recolor_targets <- TRUE # whether to use manual colors for targets
manually_reshape_treatments <- TRUE # whether to use manual shapes for treatments
manually_reshape_targets <- TRUE # whether to use manual shapes for targets
manually_relabel_treatments <- TRUE # whether to override treatment legend labels
manually_relabel_targets <- FALSE # whether to override target legend labels
input_directory <- "example_input_dir/" # path to directory containing input files
output_directory <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # font size for plots. 14 is theme_prism default
pt_size = 3 # point size for plots
no_legend <- FALSE # whether all plots should have no legend
input_filename <- "example_input_filename.csv"
# the order of treatments and targets is the order they will be plotted
treatments <- c(
  "example_treatment_1",
  "example_treatment_2"
)
targets <- c(
  "example_target_1",
  "example_target_2"
)
