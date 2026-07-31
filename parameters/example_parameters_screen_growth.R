source("parameters/manual_scales_screen_growth.R")
input_directory <- "example_input_dir/" # path to directory containing input files
output_directory <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # font size for plots. 14 is theme_prism default
pt_size = 3 # point size for plots
input_filename <- "example_input_dir/example_input_filename.csv"
concs_IC25 <- c(
  "ponatinib Shokat" = 0.1,
  # "ponatinib Pritchard" = 0.1,
  "asciminib" = 5,
  "ponatinib + asciminib" = 0.1,
  "PonatiLink-2-7-10" = 1
)
concs_IC50 <- c(
  "ponatinib Shokat" = 0.4,
  # "ponatinib Pritchard" = 0.4,
  "asciminib" = 10,
  "ponatinib + asciminib" = 0.4,
  "PonatiLink-2-7-10" = 5
)
concs_IC75 <- c(
  "ponatinib Shokat" = 1,
  # "ponatinib Pritchard" = 1,
  "asciminib" = 30,
  "ponatinib + asciminib" = 1,
  "PonatiLink-2-7-10" = 10
)
concs_IC90 <- c(
  "ponatinib Shokat" = 5,
  # "ponatinib Pritchard" = 5,
  "asciminib" = 100,
  "ponatinib + asciminib" = 5,
  "PonatiLink-2-7-10" = 30
)
manual_color_four_each <- c(
  "pink", "salmon", "red2", "red3",
  "plum1", "plum", "plum4", "purple4",
  "palegoldenrod", "yellow", "gold", "gold3",
  "palegreen", "green2", "palegreen4", "darkgreen",
  "lightblue", "cadetblue", "blue", "darkblue"

)
manual_id <- c(
  "pon_s_0.1",
  "pon_s_0.4",
  "pon_s_1",
  "pon_s_5",
  "pon_p_0.1",
  "pon_p_0.4",
  "pon_p_1",
  "pon_p_5",
  "asc_5",
  "asc_10",
  "asc_30",
  "asc_100",
  "pa_0.1",
  "pa_0.4",
  "pa_1",
  "pa_5",
  "pl_1",
  "pl_5",
  "pl_10",
  "pl_30"
)
