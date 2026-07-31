# parameters for plot_nanobret_KD
source("parameters/manual_scales.R")
input_directory <- "example_input_dir/" # path to directory containing input files
output_directory <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
input_filename <- "example_input_filename.csv"
treatments <- c(
  "example_treatment_1",
  "example_treatment_2"
)
targets <- c(
  "example_target_1",
  "example_target_2"
)
x_limits <- c(30,300)
tracer_name <- "your tracer name here"
