#' ---
#'title: "plateplotr"
#'author: "Jack Stevenson"
#'date: "2022-09-29"
#' ---
#'This is a baby spinoff from plateplotr for plotting EC50s
#'
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(scales)
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes

# specify names of input files and import data---------------------------------
# note the order compounds are imported is the order they will be plotted
input_filename <- "EC50s.csv"
dir.create("output/", showWarnings = FALSE)
plot_type <- "pdf"
# choose and order compounds and targets to plot
source("compounds.R")
source("targets.R")
EC_data <- read_csv(input_filename) %>%
  mutate(linker_length = as.numeric(linker_length)) %>%
  mutate(neglog10EC50_nM = -log10(EC50_nM)) %>%
  # filter for desired compounds
  filter(compound %in% compounds) %>%
  # filter for desired compounds
  filter(target %in% targets) %>%
  mutate(compound = fct_relevel(compound, compounds)) %>% # order compounds by list
  mutate(target = fct_relevel(target, targets)) # order targets by list
# set factors so things get plotted and colored in input order
Abl_factors <- distinct(EC_data, Abl)$Abl
EC_data <- EC_data %>%
  mutate(Abl = fct_relevel(Abl, Abl_factors))
# plot ECs in points---------------------------------------------------------------
EC_data %>%
  ggplot(aes(x = compound, y = EC50_nM)) +
  geom_point(aes(shape = assay, color = Abl), size = 4) +
  scale_shape(labels = c("CellTiter-Glo", "SelectScreen")) +
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  # -log10 transform to show most potent on top
  scale_y_continuous(trans = c("log10","reverse")) +
  scale_color_manual(values = c("black","red3"),
                     labels = c("Abl wt", "Abl T315I")) +
  #scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.8) +
  theme_prism() + # make it look fancy like prism
  guides(shape = guide_legend(order = 1)) + # force shape to top of legend
  labs(x = "compound",
       y = "EC50 (nM)",
       title = ("PonatiLink-2 cell-based vs. biochemical potency")) +
  theme(plot.background = element_blank()) # need for transparent background
ggsave(str_glue("output/EC50_points.{plot_type}"),
          bg = "transparent",
          width = 8,
          height = 6)
# plot ECs by linker length-----------------------------------------------------
EC_data %>%
  filter(linker_length > 0) %>% # only plot compounds with linkers
  ggplot(aes(x = linker_length, y = EC50_nM,
             shape = assay, color = Abl)) +
  theme_prism() + # make it look fancy like prism
  scale_x_continuous(
                     guide = "prism_offset_minor", # end at last tick
                     breaks = seq(11,23,2)) + # manual x ticks
  scale_y_continuous(trans = c("log10", "reverse"),
                     guide = "prism_offset_minor") +
  scale_color_manual(values = c("black","red3")) +
  #scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.8) +
  # remove placeholder shape from color legend with 32, the nonshape
  guides(color=guide_legend(override.aes=list(shape=32))) +
  geom_point(size = 3) +
  geom_line() +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "linker length (PEG units)",
       y = "EC50 (nM)",
       title = "PonatiLink-2 cell-based vs. biochemical potency")
ggsave(str_glue("output/EC50_linker.{plot_type}"),
        bg = "transparent",
       width = 7,
       height = 3.5)
# plot ECs comparing assays-----------------------------------------------------
EC_data %>%
  pivot_wider(names_from = assay, values_from = EC50_nM, id_cols = c(linker_length, Abl)) %>%
  filter(linker_length > 0) %>% # only plot compounds with linkers
  ggplot(aes(x = CTG, y = SelectScreen)) +
  facet_wrap(vars(Abl)) +
  scale_size(range = c(2,8)) +
  geom_point(aes(size = linker_length, color = linker_length)) +
  scale_x_continuous(trans = c("log10", "reverse")) +
  scale_y_continuous(trans = c("log10", "reverse")) +
  coord_fixed() + # even coordinate spacing on both axes
  scale_color_viridis(begin = .95, end = 0) +
  #guides(color=guide_legend(override.aes=list(shape=32))) +
  theme_prism() + # make it look fancy like prism
  theme(panel.spacing = unit(.5, "inches")) +
  # theme(panel.background = element_rect(fill = NA, color = "black")) + # box facets
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "CTG EC50",
       y = "SelectScreen EC50",
       title = str_wrap("Potency of PonatiLink-2 series in cell-based vs biochemical assays", width = 50))
ggsave(str_glue("output/EC50_assays.{plot_type}"),
       bg = "transparent",
       width = 8,
       height = 6)

