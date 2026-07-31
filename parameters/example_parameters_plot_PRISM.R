# parameters for plot_PRISM
scales_path <- "parameters/manual_scales.R"
source(scales_path) # get manual color, shape and name keys

random_seed <- 16

#input_dir <- "example_input_dir/" # path to directory containing input files
output_dir <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # font size for plots. 14 is theme_prism default
pt_size = 3 # point size for plots

input_filename <- "example_input_dir/example_input_filename.csv"
fusions_filename <- "example_input_dir/example_fusions_filename.csv"
# the order of treatments is the order they will be plotted
treatments <- c(
  "example_treatment_1",
  "example_treatment_2"
)

subtypes_to_label <- c(
  "Chronic Myeloid Leukemia, BCR-ABL1+" = "CML",
  "B-Lymphoblastic Leukemia/Lymphoma with t(9;22)(q34.1;q11.2);BCR-ABL1" = "ALL",
  "T-Lymphoblastic Leukemia/Lymphoma" = "TLL",
  "other"
)
