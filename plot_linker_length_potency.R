#' ---
#'title: "plot_IC50"
#'author: "Jack Stevenson"
#'date: "2023-02"
#' ---
#'This script plots IC50s against linker lengths of linked compounds
#'
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(scales)
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(lemon) # for fancy facet wrapping
library(doseplotr)
# import parameter files and data---------------------------------
# import global parameters and clear environment--------------------------------
rm(list = ls()) # clear environment
source("parameters/parameters_plot_IC50.R")
# note the order treatments are imported is the order they will be plotted
dir.create("output/", showWarnings = FALSE)
plot_type <- "pdf"
IC_data <- readxl::read_excel(input_filename) |>
  mutate(linker_length = as.numeric(linker_length)) |>
  mutate(IC50_nM = as.numeric(IC50_nM)) |>
  # filter for desired treatments, targets and target variants
  filter(treatment %in% treatments) |>
  filter(target %in% targets) |>
  filter(variant %in% variants) |>
  mutate(treatment = fct_relevel(treatment, treatments)) |> # order treatments by list
  mutate(target = fct_relevel(target, targets)) |> # order targets by list
  mutate(variant = fct_relevel(variant, variants)) |> # order variants by list
  mutate(neglog10IC50_nM = -log10(IC50_nM))
# plot ICs in points--------- ------------------------------------------------------
IC_data |>
  ggplot(aes(x = treatment, y = IC50_nM)) +
  geom_point(aes(shape = assay, color = variant), size = 4, alpha = 1) +
  scale_shape(labels = assay_labels,
              guide = guide_legend(order = 1)) + # force to top of legend
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  # reverse and log10 transform to show most potent on top
  scale_y_continuous(trans = c("log10","reverse")) +
  scale_color_manual(values = color_map_variants) +
  theme_prism() + # make it look fancy like prism
  theme(legend.text.align = 0) +
  labs(x = "treatment",
       y = "IC50 (nM)",
       title = (glue::glue("{group_name} series potency"))) +
  theme(plot.background = element_blank()) # need for transparent background
ggsave(str_glue("output/IC50_points_{get_timestamp()}.{plot_type}"),
          bg = "transparent",
          width = 8,
          height = 6)
# calculate linker length scale parameters--------------------------------------
linkers <- IC_data |>
  filter(linker_length > 0) |>
  distinct(linker_length) |>
  pull(linker_length)
linker_min <- min(linkers)
linker_max <- max(linkers)
linker_seq <- seq(linker_min, linker_max, 2)
# plot ECs by linker length-----------------------------------------------------
IC_data |>
  filter(linker_length > 0) |> # only plot treatments with linkers
  ggplot(aes(x = linker_length, y = IC50_nM,
             shape = assay, color = variant)) +
  scale_shape(labels = assay_labels,
              guide = guide_legend(order = 1)) + # force to top of legend
  theme_prism() + # make it look fancy like prism
  scale_x_continuous(
                     guide = "prism_offset_minor", # end at last tick
                     breaks = linker_seq) + # manual x ticks
  scale_y_continuous(trans = c("log10", "reverse"),
                     guide = "prism_offset_minor") +
  scale_color_manual(values = color_map_variants) +
  # remove placeholder shape from color legend with 32, the nonshape
  guides(color=guide_legend(override.aes=list(shape=32))) +
  geom_point(size = 4) +
  geom_line() +
  theme(plot.background = element_blank()) + # need for transparent background +
  theme(legend.text.align = 0) +
  labs(x = "linker length (PEG units)",
       y = "IC50 (nM)",
       title = glue::glue("{group_name} series potency"))
ggsave(str_glue("output/IC50_linker_{get_timestamp()}.{plot_type}"),
        bg = "transparent",
       width = 8,
       height = 4)
# plot ECs across assays--------------------------------------------------------
IC_data |>
  # pivot so IC50s from both assays associate with each linker/variant combo
  pivot_wider(names_from = assay,
              values_from = IC50_nM,
              id_cols = c(linker_length, variant)) |>
  filter(linker_length > 0) |> # only plot treatments with linkers
  ggplot(aes(y = CTG, x = SelectScreen)) +
  theme_prism() + # make it look fancy like prism
  lemon::facet_rep_wrap(vars(variant), repeat.tick.labels = "left") + # repeat y axis
  theme(panel.spacing = unit(.5, "inches")) + # spacing between facets
  # log10 and reverse both axes for intuitive potency direction
  scale_x_continuous(trans = c("log10", "reverse"),
                     expand = expansion(mult = .1),
                     # breaks = c(1.5, 1, .5)
                     ) +
  scale_y_continuous(trans = c("log10", "reverse"),
                     expand = expansion(mult = .1),
                     # limits = c(3, 1),
                     # breaks = c(3, 2, 1)
                     ) +
  coord_fixed() + # even coordinate spacing on both axes
  # geom_path(color = "gray", linewidth = 1) + # path between points
  geom_point(aes(size = linker_length, color = linker_length), alpha = 0.9) +
  # manual limits, breaks and labels for a unified linker length legend
  scale_size_continuous(range = c(3,7),
                        limits = c(linker_min, linker_max),
                        breaks = linkers,
                        labels = linkers,
                        guide = guide_legend(reverse = FALSE,
                                             direction = "horizontal",
                                             title.position = "top",
                                             label.position = "bottom",
                                             nrow = 1,
                                             keywidth = .5,
                                             title = "linker length (PEG units)"),) +
  scale_color_viridis(option = "viridis", # default viridis color scale
                      begin = 1, end = 0,
                      limits = c(linker_min, linker_max),
                      breaks = linkers,
                      labels = linkers,
                      guide = guide_legend(reverse = FALSE,
                                           direction = "horizontal",
                                           title.position = "top",
                                           label.position = "bottom",
                                           nrow = 1,
                                           keywidth = .5,
                                           title = "linker length (PEG units)")) +
  theme(plot.background = element_blank(), # need for transparent background
        legend.title = element_text(size = 10), # reinstate legend title bc theme_prism removes
        # legend.position = "bottom",
        strip.text.x = element_text(size = 14)) + # size facet labels
  labs(y = "CellTiter-Glo IC50 (nM)",
       x = "SelectScreen IC50 (nM)",
       title = str_wrap(glue::glue("{group_name} series potency"), width = 70))
ggsave(str_glue("output/IC50_assays_{get_timestamp()}.{plot_type}"),
       bg = "transparent",
       width = 10,
       height = 5)
