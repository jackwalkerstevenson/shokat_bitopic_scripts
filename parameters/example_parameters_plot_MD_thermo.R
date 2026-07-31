# parameters for plot_MD_RMSD
source("parameters/manual_scales.R")
input_dir <- "example_input_dir/" # path to directory containing input files
output_dir <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
thermo_filename <- "example_entropy_data.xlsx"
key_filename <- "example_molecule_key.xlsx"
RMSF_filename <- "example_RMSF_data.csv"
IC50_filename <- "example_IC50s.xlsx"
compounds <- c(
  "example_compound_1",
  "example_compound_2"
)
