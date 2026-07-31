# parameters for plot_fold_change
source("parameters/manual_scales.R") # manual color and shape scales
manually_recolor_treatments <- TRUE # whether to use manual colors for treatments
manually_recolor_targets <- TRUE # whether to use manual colors for targets
manually_relabel_treatments <- TRUE # whether to override treatment legend labels
manually_relabel_targets <- TRUE # whether to override target legend labels
input_directory <- "example_input_dir/" # path to directory containing input files
output_directory <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # font size for plots. 14 is theme_prism default
grid <- FALSE # whether to plot a background grid
input_filename <- "example_input_filename.csv"
# wt and control must be among the targets below
wt_target_name <- "example_target_1"
control_target_name <- "example_target_2"
fold_change_axis_title <- "your fold change axis title here"
target_axis_title <- "your target axis title here"
# treatments and targets to include in plots
# the order of these vectors is the order they will be plotted
treatments <- c(
  "example_treatment_1",
  "example_treatment_2"
)
targets <- c(
  "example_target_1",
  "example_target_2"
)
