#' ---
#'title: "plateplotr"
#'author: "Jack Stevenson"
#'date: "2022-09-29"
#' ---
#'This is a baby spinoff from plateplotr for plotting EC50s
#'
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes

# specify names of input files and import data---------------------------------
# note the order compounds are imported is the order they will be plotted
input_filename <- "EC50s.csv"
plot_type <- "pdf"
# choose and order compounds to plot
source("compounds.R")
EC_data <- read_csv(input_filename) %>%
  mutate(log10EC50_nM = log10(EC50_nM)) %>%
  # filter for desired compounds
  filter(compound %in% compounds) %>%
  mutate(compound = fct_relevel(compound, compounds)) # order compounds by list
# set factors so experiments and targets get plotted and colored in input order
experiment_factors <- distinct(EC_data, experiment)$experiment
Abl_factors <- distinct(EC_data, Abl)$Abl
EC_data <- EC_data %>% 
  mutate(experiment = fct_relevel(experiment, experiment_factors)) %>%
  mutate(Abl = fct_relevel(Abl, Abl_factors))
# plot ECs in points---------------------------------------------------------------
EC_data %>%
  ggplot(aes(x = compound, y = log10EC50_nM)) +
  geom_point(aes(shape = target, color = experiment), size = 3) +
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  theme_prism() + # make it look fancy like prism
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.8) +
  labs(x = "compound",
       y = "log10(EC50) (nM)",
       title = "Relative potency of compounds") +
  theme(plot.background = element_blank()) # need for transparent background
  ggsave(str_glue("plots output/EC50_points.{plot_type}"),
            bg = "transparent",
            width = 7,
            height = 5)

# plot ECs by linker length-----------------------------------------------------
EC_data %>%
  filter(linker_length > 0) %>%
  #group_by(experiment) %>%
  ggplot(aes(x = linker_length, y = EC50_nM,
             color = experiment, shape = Abl)) +
  scale_x_continuous(guide = "prism_offset_minor", # end at last tick
                     breaks = c(12, 16, 20, 24)) + # manual x ticks
  scale_y_log10(guide = "prism_offset") +
  #scale_y_continuous(guide = "prism_offset") +
  scale_color_manual(values = c("blue2","green3"), labels = c("CellTiter-Glo", "SelectScreen")) +
  #scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.8) +
  guides(color=guide_legend(override.aes=list(shape=15))) +
  geom_point(size = 3) +
  geom_line() +
  theme_prism() + # make it look fancy like prism
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "linker length (PEG units)",
       y = "EC50 (nM)",
       title = "Potency of PonatiLink-1 series by linker length") +
  ggsave(str_glue("plots output/EC50_linker.{plot_type}"),
         bg = "transparent",
         width = 7,
         height = 3.5)
