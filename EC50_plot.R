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
compounds <- c("ponatinib",
               #"dasatinib",
               #"asciminib",
               #"ponatinib + asciminib",
               "PonatiLink-1-12",
               "PonatiLink-1-16",
               "PonatiLink-1-20",
               "PonatiLink-1-24")
EC_data <- read_csv(input_filename) %>%
  # filter for desired compounds
  filter(compound %in% compounds) %>%
  mutate(compound = fct_relevel(compound, compounds)) # order compounds by list
# set factors so cell lines get plotted and colored in input order
experiment_factors <- distinct(EC_data, experiment)$experiment
EC_data <- EC_data %>% 
  mutate(experiment = fct_relevel(experiment, experiment_factors))
# helper function for saving plots----------------------------------------------
scale_facet <- 4 # plot width per col/height per row
#todo: allow overriding width and height if provided
save_plot <- function(plot, nrow = 1, ncol = 1, ...){
  ggsave(plot, bg = "transparent",
         width = ncol*scale_facet + 2,
         height = nrow*scale_facet, ...)}
# plot ECs in bar---------------------------------------------------------------
EC_data %>%
  filter(experiment == "CTG") %>%
  ggplot(aes(x = compound, y = EC50_nM, fill = target)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_prism() + # make it look fancy like prism
  theme(plot.background = element_blank()) # need for transparent background

# plot ECs by linker length-----------------------------------------------------
EC_data %>%
  filter(linker_length > 0) %>%
  #filter(experiment == "CTG") %>%
  group_by(experiment) %>%
  #filter(target == "wt") %>%
  ggplot(aes(x = linker_length, y = EC50_nM, color = target, shape = experiment)) +
  geom_point() +
  geom_line() +
  theme_prism() + # make it look fancy like prism
  theme(plot.background = element_blank()) # need for transparent background



# helper function to add ggplot objects common to all plots--------------------
plot_global <- function(plot){
  plot +
    scale_x_continuous(guide = "prism_offset_minor", # end at last tick
                       minor_breaks = minor_x) + # manual minor ticks
    scale_y_continuous(guide = "prism_offset",  # end at last tick
                       breaks = c(0,25,50,75,100)) + # manual y axis ticks
    coord_cartesian(xlim = x_limits, # set x axis zoom from global values
                    ylim = c(0,NA)) + # set y axis zoom locally
    theme_prism() + # make it look fancy like prism
    theme(plot.background = element_blank()) + # need for transparent background
    labs(x = "log [compound] (M)",
         y = "relative cell viability (%)")
}
# helper function to plot one compound----------------------------------------
plot_compound <- function(cpd){
  plate_summary <- plate_data %>%
    filter(compound == cpd) %>% # get data from one compound to work with
    group_by(cell_line, log.conc) %>%  # get set of replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = cell_line)) +
      geom_point() +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w)) +
      # use drm method from drc package to fit dose response curve
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, size = 1)} %>%
    plot_global() +
    theme(aspect.ratio = 1) +
    scale_color_manual(values = c("black","darkred")) +
    labs(title = cpd)
}