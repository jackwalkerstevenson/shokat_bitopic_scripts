# visualizing growth rates from saturation mutagenesis
# Jack Stevenson, started 2023-11-15
# # load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
#library(readxl) # for excel file handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(scales) # for axis transforms
library(viridis) # for color schemes
library(doseplotr) # you bet
# import global parameters and clear environment--------------------------------
rm(list = ls()) # clear environment
source("parameters/parameters_sat_mut.R")
dir.create(output_directory, showWarnings = FALSE)
# import data-----------------------------------------------
all_data <- readr::read_csv(input_filename)
# dataframe of manual annotations for sequence position plots-------------------
manual_highlight_residues <- tibble(
  protein_start = c(255, 315, 468),
  residue_name = c("E255", "T315", "V468")
)
# ggplot chunk for aesthetics of sequence plots---------------------------------
seq_plot <- function(){
  list(
    labs(x = "sequence position"),
    scale_x_continuous(breaks = scales::breaks_width(25),
                       minor_breaks = scales::breaks_width(5),
                       expand = expansion(mult = 0)),
    theme_prism(),
    theme(panel.grid = element_line(color = "black", linewidth = 0.5),
          panel.grid.minor = element_line(color = "black",
                                          linewidth = 0.1,
                                          linetype = "dotted"))
  )
}
# plot number of rows per position---------------------------------------------
(all_data |>
  group_by(protein_start) |> 
  summarize(n = n()) |> 
  ggplot(aes(x = protein_start, y = n)) +
  geom_bar(stat = "identity") +
  labs(title = "number of rows per position") +
  scale_y_continuous(breaks = scales::breaks_width(5),
                     minor_breaks = scales::breaks_width(1)) +
  seq_plot()) |>
save_plot(str_glue("output/num_rows_{get_timestamp()}.{plot_type}"),
          width = 12, height = 6)
# plot number of duplicate mutations per position-------------------------------
(all_data |> 
  group_by(protein_start, alt) |> # group by desired mutation
  summarize(.groups = "keep",
            n = n()) |> 
  ggplot(aes(x = protein_start, y= n)) +
  geom_bar(stat = "identity") +
  labs(title = "number of duplicate mutations per position") +
  seq_plot()) |> 
  save_plot(str_glue("output/num_duplicates_{get_timestamp()}.{plot_type}"),
            width = 12, height = 6)
# function to plot all data by sequence position--------------------------------
plot_entire_sequence2 <- function(data, condition, viridis_option = "turbo"){
  ggplot(data, aes(x = protein_start, y = get(condition))) +
    geom_point(aes(color = alt_aa)) + # color by mutant residue
    geom_vline(xintercept = c(255, 315, 468),
               linetype = "dashed") +
    geom_label(data = manual_highlight_residues,
               mapping = aes(x = protein_start, label = residue_name),
               y = 0) +
    scale_color_viridis(discrete = TRUE, option = viridis_option) +
    labs(title = "Growth rates by sequence",
         x = "sequence position",
         y = get_if(display_names_treatments[condition], condition %in% names(display_names_treatments), condition)
         ) +
    seq_plot()
}
# todo: plot all residues but each one is mean±SEM for multiple conditions
# todo: function to plot all data by sequence, duplicates averaged--------------------