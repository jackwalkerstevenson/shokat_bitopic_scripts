# parameters for plot_MD_RMSF
source("parameters/manual_scales.R")
input_dir <- "example_input_dir/" # path to directory containing input files
output_dir <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
key_filename <- "example_molecule_key.xlsx"
compounds <- c(
  "example_compound_1",
  "example_compound_2"
)
y_limits <- c(0,12)
