# visualizing growth rates from saturation mutagenesis
# Jack Stevenson, started 2023-11-15
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
#library(readxl) # for excel file handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(scales) # for axis transforms
library(viridis) # for color schemes
library(pals) # color palettes
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
    geom_vline(xintercept = c(321, 394, 465),
               linetype = "dashed"),
    labs(x = "sequence position",
         caption = "Dashed lines represent pool region boundaries"),
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
   count(protein_start) |> 
   ggplot(aes(x = protein_start, y = n)) +
   geom_bar(stat = "identity") +
   labs(title = "dataset rows per position",
        y = "number of rows") +
   scale_y_continuous(breaks = scales::breaks_width(5),
                      minor_breaks = scales::breaks_width(1)) +
   seq_plot()) |>
  save_plot(str_glue("output/num_rows_{get_timestamp()}.{plot_type}"),
            width = 12, height = 6)
# plot number of variants per position------------------------------------------
(all_data |>
   count(protein_start, alt_aa, name = "num_rows_per_variant") |> 
   count(protein_start, name = "num_variants_per_position") |> 
   ggplot(aes(x = protein_start, y = num_variants_per_position)) +
   geom_bar(stat = "identity") +
   geom_hline(yintercept = 19, linetype = "dotted") +
   labs(title = "variant coverage",
        y = "number of unique variants") +
   scale_y_continuous(breaks = scales::breaks_width(5),
                      minor_breaks = scales::breaks_width(1)) +
   seq_plot()) |>
  save_plot(str_glue("output/num_variants_{get_timestamp()}.{plot_type}"),
            width = 12, height = 6)
# plot number of repeat variants per position, stacked residues-------------
(all_data |>
    count(protein_start, alt_aa, name = "num_residue_repeats") |> 
    filter(num_residue_repeats > 1) |>
    ggplot(aes(x = protein_start, y = num_residue_repeats)) +
    geom_bar(aes(fill = alt_aa), stat = "identity", position = "stack") +
    scale_fill_manual(values = unname(pals::cols25(n = 20))) +
    scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
    labs(title = "repeated variants per position",
        y = "number of repeated variants") +
    seq_plot()) |>
  save_plot(str_glue("output/num_repeats_position_stack_{get_timestamp()}",
                     ".{plot_type}"),
            width = 12, height = 6)
# plot number of repeat variants per position total-------------------------
(all_data |>
   # group_by(protein_start, alt_aa) |> # group by desired variant
   # summarize(num_residue_repeats = n()) |> 
   count(protein_start, alt_aa, name = "num_residue_repeats") |> 
   filter(num_residue_repeats > 1) |> 
   group_by(protein_start) |>
   summarize(num_position_repeats = sum(num_residue_repeats)) |> 
   ggplot(aes(x = protein_start, y = num_position_repeats)) +
   geom_bar(stat = "identity") +
   scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
   labs(title = "total repeated variants per position",
        y = "number of repeated variants") +
   seq_plot()) |>  
  save_plot(str_glue("output/num_repeats_position_{get_timestamp()}",
                     ".{plot_type}"),
            width = 12, height = 6)
# function to plot all data by sequence position--------------------------------
plot_entire_sequence2 <- function(data, condition, viridis_option = "turbo"){
  ggplot(data, aes(x = protein_start, y = get(condition))) +
    geom_point(aes(color = alt_aa)) + # color by mutant residue
    # geom_vline(xintercept = c(255, 315, 468),
               # linetype = "dotted") +
    # geom_label(data = manual_highlight_residues,
               # mapping = aes(x = protein_start, label = residue_name),
               # y = 0) +
    scale_color_viridis(discrete = TRUE, option = viridis_option) +
    labs(title = "Growth rates by sequence",
         x = "sequence position",
         y = get_if(display_names_treatments[condition], condition %in% names(display_names_treatments), condition)
         ) +
    seq_plot()
}
# todo: plot all residues but each one is mean±SEM for multiple conditions
# determine how many unique substitutions are in the dataset
test_grouped <- all_data |> 
  count(protein_start, alt_aa) #|> 
  # nrow()
# plot each treatment-----------------------------------------------
# treatments <- unique(all_data$)
# likely todo: function to plot all data by sequence, repeats averaged--------------------