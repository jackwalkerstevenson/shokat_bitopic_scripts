# parameters for plot_PRISM
scales_path <- "parameters/manual_scales.R"
source(scales_path) # get manual color, shape and name keys

# input_dir <- "example_input_dir/" # path to directory containing input files
output_dir <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # font size for plots. 14 is theme_prism default
pt_size = 3 # point size for plots

input_path <- "example_input_dir/example_body_weight_data.xlsx"
key_path <- "example_input_dir/example_animal_key.xlsx"

study_name <- "example study name"
