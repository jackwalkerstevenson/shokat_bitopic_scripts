# plot apparent KD values from a model summary output by plot_nanobret
# load libraries----------------------------------------------------------------
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(scales) # for axis labels
library(doseplotr) # you bet
# import precalculated EC50 table-----------------------------------------------
rm(list = ls()) # clear environment
params_path <- "parameters/parameters_plot_nanobret_KD.R"
scales_path <- "parameters/manual_scales.R"
source(params_path)
source(scales_path)
input_path <- str_glue("{input_directory}{input_filename}")
data <- read_csv(input_path)
# process data---------------------------------------------------------------
if(exists("treatments")){
  data <- filter_validate_reorder(data, "treatment", treatments)
}
if(exists("targets")){
  data <- filter_validate_reorder(data, "target", targets)
}

# write input data, parameters and report of fold changes-----------------------------------
doseplotr::file_copy_to_dir(input_path, output_directory)
doseplotr::file_copy_to_dir(params_path, output_directory)
doseplotr::file_copy_to_dir(scales_path, output_directory)
# write timestamped code to output
doseplotr::file_copy_to_dir("nanobret/plot_nanobret_KD.R", output_directory)
# strip plot of raw EC50s-----------------------------------------------
# # x_min <- floor(min(log10(data$IC50_nM)))
# x_min <- data$IC50_nM |> na.omit() |> log10() |> min() |> floor()
# # x_max <- ceiling(max(log10(data$IC50_nM)))
# x_max <- data$IC50_nM |> na.omit() |> log10() |> max() |> ceiling()
# breaks_x <- 10^rep(x_min : x_max)
# minor_x <- rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9))
p <- data |> 
  dplyr::mutate(target = fct_reorder(target, EC50_nM)) |> 
  ggplot(aes(y = target,
             x = EC50_nM,
             color = target,
             shape = target)) +
  geom_point(size = 5, alpha = 1, stroke = 1) +
  scale_x_continuous(transform = "log10",
                     # guide = "prism_minor", # end at last tick
                     guide = guide_axis_logticks(long = 1, mid = 0.5, short = 0.5),
                     limits = x_limits) +
                     # breaks = breaks_width(50),
                     # labels = label_comma(accuracy = 1, big.mark = ""),
                     # minor_breaks = minor_x,
                     # expand = expansion(mult = .1)) +
  scale_y_discrete(limits = rev, labels = display_names_targets,
                   name = "competitor") +
  scale_color_manual(values = color_map_targets,
                     labels = display_names_targets,
                     name = "competitor") +
  scale_shape_manual(values = shape_map_targets,
                     labels = display_names_targets,
                     name = "competitor") +
  theme_prism() +
  theme(#legend.title = element_text(),
        legend.position = "none", # remove legend entirely
        panel.grid = element_line(color = "black", linewidth = 0.1),
        panel.grid.minor = element_line(color = "black",
                                        linewidth = 0.1,
                                        linetype = "dotted")) +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = bquote(bold(K[D]^apparent~"(nM) of" ~ .(tracer_name))),
       y = "treatment")
save_plot(p, str_glue("output/nanobret_EC50_dot_{get_timestamp()}.{plot_type}"),
          width = 9, height = .5*length(targets) + 0.75)
