# parameters for plot_growth_rate
scales_path <- "parameters/manual_scales.R"
source(scales_path) # get manual color, shape and name keys

#random_seed <- 1

input_dir <- "example_input_dir/" # path to directory containing input files
output_dir <- "example_output_dir/" # path to directory in which to write output files
plot_type <- "pdf" # file type of saved plot images
font_base_size <- 14 # font size for plots. 14 is theme_prism default
pt_size = 3 # point size for plots

start_live_cell_count <- 150 # number of live cells all treatments started with

# maximum cell count to plot
count_plot_limit <- 20000

# the order of treatments is the order they will be plotted
# names are display names, values are treatment names as they appear in the data
treatments <- c(
  "Example Treatment 1" = "example_treatment_1",
  "Example Treatment 2" = "example_treatment_2"
)

asc_dose_key_nM <- c(
  "KRA" = 0,
  "KRA-AL" = 100,
  "KRA-AM" = 500,
  "KRA-AH" = 1000
)

display_names <- c(
  "example_treatment_1 10 nM",
  "example_treatment_2 100 nM"
)

color_map_display_names <- c(
  # "ponatinib" = "#A40000", # apparent darkred from chimera
  # "asciminib" = "#FF5300", # apparent orangered from chimera
  # "ponatinib + asciminib" = "gold2",
  # "PonatiLink-2-7-10" = "skyblue3",
  "ponatinib 1 nM" = "#A40000", # apparent darkred from chimera
  "ponatinib 10 nM" = "#A40000", # apparent darkred from chimera
  "ponatinib 30 nM" = "#A40000", # apparent darkred from chimera
  "asciminib 100 nM" = "#FF5300", # apparent orangered from chimera
  "asciminib 500 nM" = "#FF5300", # apparent orangered from chimera
  "asciminib 1000 nM" = "#FF5300", # apparent orangered from chimera
  "ponatinib + asciminib 1 nM + 500 nM" = "gold2",
  "ponatinib + asciminib 10 nM + 500 nM" = "gold2",
  "ponatinib + asciminib 30 nM + 500 nM" = "gold2",
  "PonatiLink-2 1 nM" = "skyblue3",
  "PonatiLink-2 10 nM" = "skyblue3",
  "PonatiLink-2 100 nM" = "skyblue3",
  "PonatiLink-2 1000 nM" = "skyblue3"
)

shape_map_display_names <- c(
  "ponatinib 1 nM" = "circle",
  "ponatinib 10 nM" = "circle",
  "ponatinib 30 nM" = "circle",
  "asciminib 100 nM" = "triangle",
  "asciminib 500 nM" = "triangle",
  "asciminib 1000 nM" = "triangle",
  "ponatinib + asciminib 1 nM + 500 nM" = "square",
  "ponatinib + asciminib 10 nM + 500 nM" = "square",
  "ponatinib + asciminib 30 nM + 500 nM" = "square",
  "PonatiLink-2 1 nM" = "square open",
  "PonatiLink-2 10 nM" = "square open",
  "PonatiLink-2 100 nM" = "square open",
  "PonatiLink-2 1000 nM" = "square open"
)

linetype_map_display_names <- c(
  "ponatinib 1 nM" = "dotted",
  "ponatinib 10 nM" = "dashed",
  "ponatinib 30 nM" = "solid",
  "asciminib 100 nM" = "dotted",
  "asciminib 500 nM" = "dashed",
  "asciminib 1000 nM" = "solid",
  "ponatinib + asciminib 1 nM + 500 nM" = "dotted",
  "ponatinib + asciminib 10 nM + 500 nM" = "dashed",
  "ponatinib + asciminib 30 nM + 500 nM" = "solid",
  "PonatiLink-2 1 nM" = "dotted",
  "PonatiLink-2 10 nM" = "dotdash",
  "PonatiLink-2 100 nM" = "dashed",
  "PonatiLink-2 1000 nM" = "solid"
)
