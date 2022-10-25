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
compounds <- c("ponatinib",
               #"dasatinib",
               "asciminib",
               #"ponatinib + asciminib",
               "PonatiLink-1-12",
               "PonatiLink-1-16",
               "PonatiLink-1-20",
               "PonatiLink-1-24")
EC_data <- read_csv(input_filename) %>%
  mutate(pEC50_nM = -log10(EC50_nM)) %>%
  # filter for desired compounds
  filter(compound %in% compounds) %>%
  mutate(compound = fct_relevel(compound, compounds)) # order compounds by list
# set factors so experiments and targets get plotted and colored in input order
experiment_factors <- distinct(EC_data, experiment)$experiment
target_factors <- distinct(EC_data, target)$target
EC_data <- EC_data %>% 
  mutate(experiment = fct_relevel(experiment, experiment_factors)) %>%
  mutate(target = fct_relevel(target, target_factors))
# plot ECs in points---------------------------------------------------------------
EC_data %>%
  ggplot(aes(x = compound, y = pEC50_nM)) +
  geom_point(aes(shape = target, color = experiment), size = 3) +
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  theme_prism() + # make it look fancy like prism
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.8) +
  theme(plot.background = element_blank()) # need for transparent background
  ggsave(str_glue("plots output/EC50_points.{plot_type}"),
            bg = "transparent",
            width = 7,
            height = 5)

# plot ECs by linker length-----------------------------------------------------
EC_data %>%
  filter(linker_length > 0) %>%
  #filter(experiment == "CTG") %>%
  group_by(experiment) %>%
  #filter(target == "wt") %>%
  ggplot(aes(x = linker_length, y = pEC50_nM,
             color = experiment, shape = target)) +
  scale_x_continuous(guide = "prism_offset_minor", # end at last tick
                     breaks = c(12, 16, 20, 24)) + # manual x ticks
  geom_point(size = 3) +
  geom_line() +
  theme_prism() + # make it look fancy like prism
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.8) +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "linker length",
       y = "pEC50 (nM)")
  ggsave(str_glue("plots output/EC50_linker.{plot_type}"),
         bg = "transparent",
         width = 7,
         height = 3.5)
