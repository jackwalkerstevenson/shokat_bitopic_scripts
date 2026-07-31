# parameters for plot_CFU
source("parameters/manual_scales.R")
input_dir <- "example_input_dir/" # path to directory containing input files
output_dir <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
input_filename <- "example_input_filename.xlsx"
treatments <- c(
  "example_treatment_1",
  "example_treatment_2"
)

display_names_cfu <- c(
  "bfu_e" = "BFU-E",
  "cfu_g" = "CFU-G"
)

color_map_cfu <- c(
  "bfu_e" = "red3",
  "cfu_g" = "skyblue"
)
