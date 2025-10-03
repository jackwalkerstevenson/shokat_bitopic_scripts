# plot thermochemistry data from bitopic MD experiments
# depends on output from plot_MD_RMSF
# Jack Stevenson started 2025-07-31
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
library(viridis) # for color palette
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_MD_thermo.R"
source(params_path)
# set up input and output directories
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
doseplotr::file_copy_to_dir("plot_MD/plot_MD_thermo.R", output_dir) # write timestamped code
doseplotr::file_copy_to_dir(params_path, output_dir) # write timestamped params
thermo_path <- str_glue("{input_dir}{thermo_filename}")
doseplotr::file_copy_to_dir(thermo_path, output_dir) # write timestamped thermo
RMSF_path <- str_glue("{input_dir}{RMSF_filename}")
doseplotr::file_copy_to_dir(RMSF_path, output_dir) # write timestamped RMSF file
IC50_path <- str_glue("{input_dir}{IC50_filename}")
doseplotr::file_copy_to_dir(IC50_path, output_dir) # write timestamped IC50 file
key_path <- str_glue("{input_dir}{key_filename}")
doseplotr::file_copy_to_dir(key_path, output_dir) # write timestamped key file
# import, preprocess and report data-----------------------------------------------
runwise_thermo_data <- readxl::read_excel(thermo_path)
# import runwise mean RMSF data from other experiments
runwise_mean_RMSF_data <- readr::read_csv(RMSF_path)
runwise_data <- readxl::read_excel(key_path) |> 
  dplyr::left_join(runwise_thermo_data, by = join_by(kenneth_id)) |> 
  dplyr::left_join(runwise_mean_RMSF_data, by = join_by(compound_name_full, run)) |> 
  doseplotr::filter_validate_reorder("compound_name_full", compounds) |> 
  dplyr::mutate(entropy_per_mass = entropy / molecular_weight,
                entropy_per_linker_atom = entropy / linker_length_atoms,
                entropy_per_atom = entropy / atom_count)

# import IC50 data from other experiments
IC50_ABL1_data <- readxl::read_excel(IC50_path) |>
  dplyr::mutate(
    IC50_nM = as.numeric(IC50_nM),
    linker_length_PEG = as.numeric(linker_length_PEG)) |> 
  dplyr::filter( # filter for just ABL1 wt and desired assays
    variant == "ABL1 wt",
    assay %in% c("SelectScreen", "Kinomescan")
  )

write_csv(runwise_data, # report processed data
          fs::path(output_dir, str_glue("thermo_RMSF_runwise_data_{get_timestamp()}.csv")))

# join runwise to IC50 data
thermo_potency_data <- runwise_data |> 
  dplyr::left_join(IC50_ABL1_data,
                   by = join_by("compound_name_full" == "treatment",
                                "linker_length_PEG" == "linker_length_PEG"),
                   relationship = "many-to-many")
# test correlation between entropy per atom and mean RMSF-----------------------------------------------
cor_spearman <- cor.test(runwise_data$mean_linker_rmsf, runwise_data$entropy_per_atom, method = "spearman")
cor_pearson <- cor.test(runwise_data$mean_linker_rmsf, runwise_data$entropy_per_atom, method = "pearson")
# set up plot constants-----------------------------------------------
scatter_plot_width <- 6
scatter_plot_height <- 4.5
small_scatter_plot_width <- 6
small_scatter_plot_height <- 3.75
point_size_mean <- 3
point_size_individual <- 1
# plot entropy per atom vs. linker length by series------------------------------
runwise_data |> 
  dplyr::group_by(compound_name_full) |> 
  ggplot(aes(x = linker_length_atoms, y = entropy_per_atom, color = series)) +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = point_size_mean) +
  geom_point(size = point_size_individual) +
  scale_color_manual(values = color_map_series) +
  theme_prism() +
  theme(plot.background = element_blank(), # transparent
        legend.title = element_text()) +
  labs(x = "linker length (atoms)",
       y = "cpd. entropy / atom (cal/mol-K-atom)",
       title ="Entropy per atom vs. linker length")
ggsave(str_glue(
  "{output_dir}/entropy_per_atom_vs_linker_atoms_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = scatter_plot_width, height = scatter_plot_height)
# plot entropy per atom vs. SelectScreen potency by series------------------------------
thermo_potency_data |> 
  dplyr::filter(assay == "SelectScreen") |> 
  dplyr::group_by(compound_name_full) |> 
  ggplot(aes(x = entropy_per_atom, y = IC50_nM, color = series)) +
  stat_summary(
    fun = "mean",
    orientation = "y", # for variation on x axis
    geom = "point",
    size = point_size_mean) +
  geom_point(size = point_size_individual) +
  scale_y_continuous(transform = "log10",
                     guide = guide_axis_logticks(long = 1, mid = 0.5, short = 0.5)
                     ) +
  scale_color_manual(values = color_map_series) +
  theme_prism() +
  theme(plot.background = element_blank(), # transparent
        legend.title = element_text()) +
  labs(x = "cpd. entropy / atom (cal/mol-K-atom)",
       y = "ABL1 IC50 (nM)",
       title ="Entropy per atom vs. ABL1 inhibition")
ggsave(str_glue(
  "{output_dir}/entropy_per_atom_vs_selectscreen_IC50_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = small_scatter_plot_width, height = small_scatter_plot_height)
# plot entropy per atom vs. mean linker RMSF by series w/R^2------------------------------
annotation_df <- data.frame(
  x = -Inf,
  y = Inf,
  label = sprintf("R² = %.3f\np < %.3g", 
                  cor_pearson$estimate^2, 
                  cor_pearson$p.value)
)
runwise_data |> 
  ggplot(aes(x = entropy_per_atom, y = mean_linker_rmsf, color = series)) +
  geom_point(size = point_size_mean) +
  scale_color_manual(values = color_map_series) +
  # geom_text(data = annotation_df,
  #           aes(x = x, y = y, label = label),
  #           hjust = -0.1, vjust = 1.5,
  #           size = 15 / .pt,
  #           inherit.aes = FALSE
  #           ) +
  annotate("text",
           x = -Inf, y = Inf,
           hjust = -0.2, vjust = 1.5,
           label = sprintf("R² = %.3f",#\np < %.3g",
                           cor_pearson$estimate^2
                           # cor_pearson$p.value
                           ),
           size = 15 / .pt,
           fontface = "bold") +
  theme_prism() +
  theme(plot.background = element_blank(), # transparent
        legend.title = element_text()) +
  labs(x = "cpd. entropy / atom (cal/mol-K-atom)",
       y = "mean linker atom fluctuation (Å)",
       title ="Entropy per atom vs. linker fluctuation")
ggsave(str_glue(
  "{output_dir}/entropy_per_atom_vs_mean_linker_RMSF_{doseplotr::get_timestamp()}.{plot_type}"),
  bg = "transparent",
  width = scatter_plot_width, height = scatter_plot_height)