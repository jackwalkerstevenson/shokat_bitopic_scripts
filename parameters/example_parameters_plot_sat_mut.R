# parameters for plot_sat_mut
source("parameters/manual_scales_sat_mut.R")
input_directory <- "example_input_dir/" # path to directory containing input files
output_directory <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
input_path <- "example_input_dir"
input_filename <- "example_input_filename.csv"
# single_treatment <- "example_treatment_1"

# dataframe of manual annotations for sequence position plots-------------------
manual_highlight_residues <- tibble(
  protein_start = c(255, 315, 468),
  residue_name = c("E255", "T315", "V468")
)
