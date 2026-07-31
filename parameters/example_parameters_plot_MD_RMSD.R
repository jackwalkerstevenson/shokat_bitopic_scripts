# parameters for plot_MD_RMSD
source("parameters/manual_scales.R")
input_dir <- "example_input_dir/" # path to directory containing input files
output_dir <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
key_filename <- "example_molecule_key.xlsx"
treatments <- c(
  "example_treatment_1",
  "example_treatment_2"
)
y_limits <- c(1, 4.5)
