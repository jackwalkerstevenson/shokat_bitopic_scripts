#' ---
#'title: "plot Eide"
#'author: "Jack Stevenson"
#'date: "2023-03-19"
#' ---
#'started 2023-03-19
#'this script compares potency values for the same targets from multiple experiments
#'it was designed to compare asciminib IC50s reported for different vectors in Eide et al. 2019
#'
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(scales) # for fancy plotting scales
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(Cairo) # for special characters in legend label
# set global variables---------------------------------------------------------
input_filename <- "Eide asciminib potency.xlsx"
plot_type <- "pdf" # file type for saved output plots
dir.create("output/", showWarnings = FALSE) # silently create output directory
# import data------------------------------------------------------------------
all_data <- readxl::read_excel(input_filename, .name_repair = "universal") |> 
  # get greater than from e.g. ">2500", then strip values to numbers
  mutate(asc_IC50_greater_than = if_else(startsWith(asciminib_IC50, ">"), TRUE, FALSE)) |> 
  mutate(asciminib_IC50 = str_replace_all(asciminib_IC50, ">", "")) |> 
  mutate(asciminib_IC50 = as.numeric(asciminib_IC50))
# reorder target and variant factors so they plot in input order
targets <- distinct(all_data["target"])$target
variants <- distinct(all_data["variant"])$variant
all_data <- all_data |>
  mutate(target = fct_relevel(target, targets)) |>
  mutate(variant = fct_relevel(variant, variants))
# plot shared targets----------------------------------------------------------
# get targets that are shared between all variants
# filter for the targets of each variant, then get intersection
shared_targets <- map(variants, \(v) filter(all_data, variant == v)$target) |>
  reduce(intersect)
shared_data <- all_data |> filter(target %in% shared_targets)
asc_IC50_greater <- shared_data |> filter(asc_IC50_greater_than) # select data to add 'greater than' note
shared_data |>
  ggplot(aes(y = target)) +
  scale_color_viridis(discrete = TRUE , begin = 0.7, end = 0.2) +
  # scale_color_manual(values = c("green3", "royalblue3")) +
  scale_y_discrete(limits = rev) + # order top to bottom
  scale_x_continuous(trans = c("log10", "reverse")) + # negative log10 for more potency higher
                     # expand = expansion(mult = c(0, 0.05))) + # remove space from bottom of bar
  geom_point(aes(x = asciminib_IC50, color = variant), size = 2) +
  geom_errorbar(aes(x = asciminib_IC50,
                    xmin = asciminib_IC50 - asciminib_IC50_SEM,
                    xmax = asciminib_IC50 + asciminib_IC50_SEM,
                    color = variant),
                width = 0.3,
                show.legend = FALSE) +
  # label appropriate subset of data with 'greater than' symbol
  geom_text(data = asc_IC50_greater,
            aes(x = asciminib_IC50, color = variant),
            nudge_y = 0.4, label = ">", # label above point
            show.legend = FALSE) +
  theme_prism() +
  labs(
    title = "Ba/F3 asciminib potency from Eide et al. 2019",
    x = "asciminib IC50 (nM)",
    y = "Ba/F3 cell line",
    color = "Ba/F3 vector") +
  theme(plot.background = element_blank(), # need for transparent background
        legend.title = element_text(size = 12),  # reinstate legend title bc theme_prism removes)
        axis.text.y = element_text(hjust = 0)) # left-justify axis labels to line up

ggsave(str_glue("output/eide.{plot_type}"),
       bg = "transparent", width = 7, height = 4,
       device = cairo_pdf)
