#' ---
#'title: "plateplotr selectscreen"
#'author: "Jack Stevenson"
#'date: "2022-09-09"
#' ---
#'edited version of plateplotr for dealing with SelectScreen data
#'copied from main plot.R 2022-09-09

# load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(patchwork) # for plot organization
library(assertthat) # for QC assertions
options(dplyr.summarise.inform = FALSE)
# set global variables---------------------------------------------------------
input_filename <- "ZLYTE_compiled_results_complete.csv"
plot_type <- "pdf" # file type for saved output plots
dir.create("output/", showWarnings = FALSE) # silently create output directory
source("parameters/compounds.R")
source("parameters/targets.R")
source("get_EC.R")
source("get_hill_slope.R")
# import and tidy data---------------------------------
source("import_selectscreen.R")
plate_data <- import_selectscreen(input_filename, compounds, targets)
# fit models to output EC values------------------------------------------------
EC_summary <- plate_data %>%
  group_by(compound, target) %>%
  summarize(
    EC50_nM = get_EC_nM(plate_data, compound, target, 50), # convert M to nM
    # for negative-response data like this, the EC75 is the drop to 25%
    EC75_nM = get_EC_nM(plate_data, compound, target, 25), # convert M to nM
    hill_slope = get_hill_slope(plate_data, compound, target)
  )
write_csv(EC_summary, "output/EC_summary_selectscreen.csv")
# generate global parameters for all plots------------------------------------------
pt_size = 3 # size for all geom_point
all_compounds <- distinct(plate_data["compound"])$compound
all_targets <- distinct(plate_data["target"])$target
# find x-axis min/max values for consistent zoom window between all plots
x_min <- floor(min(plate_data$log.conc))
# x_min <- -10
x_max <- ceiling(max(plate_data$log.conc))
x_limits <- c(x_min, x_max)
# create logistic minor breaks for all conc plots
minor_x <- log10(rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9)))
# helper function for saving plots----------------------------------------------
scale_facet <- 4 # plot width per col/height per row
#todo: allow overriding width and height if provided
save_plot <- function(plot, nrow = 1, ncol = 1, ...){
  ggsave(plot, bg = "transparent",
         width = ncol*scale_facet + 2,
         height = nrow*scale_facet, ...)}
# helper function for summarizing replicate data for plotting------------------
# source("activity_summarize.R")
plate_summarize <- function(x){
  summarize(x,
            # standard error for error bars = standard deviation / square root of n
            confidence_intervals = list(mean_cl_normal(activity) |>
                                          rename(mean_activity = y, ymin_activity = ymin, ymax_activity = ymax)),
            sem = sd(activity, na.rm = TRUE)/sqrt(n()),
            # get mean normalized readout value for plotting
            mean_read = mean(activity),
            w = 0.12 * n() # necessary for consistent error bar widths across plots
  ) |>
    tidyr::unnest(cols = confidence_intervals)
}
# helper function to add ggplot objects common to all plots--------------------
source("dose_response_global.R")
# helper function to plot one compound----------------------------------------
plot_compound <- function(cpd, viridis_begin = 1, viridis_end = 0){
  plate_summary <- plate_data %>%
    filter(compound == cpd) %>% # get data from one compound to work with
    group_by(target, log.conc) %>%  # get set of replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = target)) +
      geom_point(aes(shape = target), size = pt_size) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w)) +
      # second error bars for 95% CI
      # geom_errorbar(aes(ymin = ymin_activity, ymax = ymax_activity, width = w), alpha = 0.4) +
      # use drm method from drc package to fit dose response curve
      #L.4 method is 4-param logistic curve. The more common LL.4() works on nonlog
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1)} %>%
    dose_response_global(x_limits) +
    theme(aspect.ratio = 1) +
    scale_color_viridis(discrete = TRUE, begin = viridis_begin, end = viridis_end) +
    # scale_color_manual(values = c("black","darkred")) +
    labs(title = cpd,
         y = "kinase activity (%)")
}
# plot data for each compound separately----------------------------------------
source("viridis_range.R")
vr <- viridis_range(length(targets))
viridis_begin <- vr[1]
viridis_end <- vr[2]
for (cpd in all_compounds){
  plot_compound(cpd, viridis_begin, viridis_end)
  # save plot with manually optimized aspect ratio
  save_plot(str_glue("output/{cpd}.{plot_type}"))
}  
# plot data for all compounds in facets----------------------------------
compound_plots = list()
for (cpd in all_compounds){
  compound_plots <- append(compound_plots, list(plot_compound(cpd)))
}

plot_mar <- 15 # margin between wrapped plots, in points
cols = ceiling(sqrt(length(compounds)))
rows = ceiling(length(compounds)/cols)
wrap_plots(compound_plots, guides = "collect", ncol = cols, nrow = rows) &
  theme(plot.margin = unit(c(plot_mar,plot_mar,plot_mar,plot_mar), "pt"),
        plot.background = element_blank(),
        legend.text= element_text(face = "bold", size = 12))
save_plot(str_glue("output/compound_facets.{plot_type}"), ncol = cols, nrow = rows)
# set color parameters for target plots--------------------------------------
alpha_val <- 1
color_scale <- "viridis"
vr <- viridis_range(length(compounds))
viridis_begin <- vr[1]
viridis_end <- vr[2]
# plot data for each target separately------------------------------------------
for (tgt in all_targets){
  plate_summary <- plate_data %>%
    filter(target == tgt) %>%
    group_by(compound, log.conc) %>% # group into replicates for each condition
    plate_summarize()
  # bracket ggplot so it can be piped to helper function
  {ggplot(plate_summary, aes(x = log.conc, y = mean_read, color = compound)) +
      geom_point(aes(shape = compound), size = pt_size) +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
      # second error bars for 95% CI
      # geom_errorbar(aes(ymin = ymin_activity, ymax = ymax_activity, width = w), alpha = 0.4) +
      # use drm method from drc package to fit dose response curve
      geom_line(stat = "smooth", method = "drm", method.args = list(fct = L.4()),
                se = FALSE, linewidth = 1, alpha = alpha_val)} %>%
    dose_response_global(x_limits) +
    scale_color_viridis(option = color_scale, discrete = TRUE, begin = viridis_begin, end = viridis_end) +
    labs(title = tgt,
         y = "kinase activity (%)")
  save_plot(str_glue("output/{tgt}.{plot_type}"))
}
# plot data for all targets at once-----------------------------------------
plate_summary <- plate_data %>%
  group_by(target, compound, log.conc) %>% # group into replicates for each condition
  plate_summarize()
{ggplot(plate_summary,aes(x = log.conc, y = mean_read, color = compound)) +
    geom_point(aes(shape = compound), size = pt_size) +
    # error bars = mean plus or minus standard error
    geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem, width = w), alpha = alpha_val) +
    # second error bars for 95% CI
    # geom_errorbar(aes(ymin = ymin_activity, ymax = ymax_activity, width = w), alpha = 0.4) +
    # use drm method from drc package to fit dose response curve
    geom_line(aes(linetype = target), stat = "smooth", method = "drm", method.args = list(fct = L.4()),
              se = FALSE, linewidth = 1, alpha = alpha_val)} %>%
  dose_response_global(x_limits) +
  scale_color_viridis(option = color_scale, discrete = TRUE, begin = viridis_begin, end = viridis_end) +
  labs(title = "All data",
       y = "kinase activity (%)")
save_plot(str_glue("output/all_data.{plot_type}"))