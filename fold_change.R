# plot fold changes between targets relative to wt for each treatment
# load libraries----------------------------------------------------------------
library(tidyverse) # for tidy data handling
library(assertthat) # for testing
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(doseplotr) # you bet
# import precalculated IC50 table-----------------------------------------------
input_filename <- "input/model_summary_2023-07-07T164437.csv"
source("parameters/treatments.R") # import list of treatments to include in plots
source("parameters/targets.R") # import list of targets to include in plots
data <- read_csv(input_filename) |> 
  filter(treatment %in% treatments) |> # take only specified treatments
  filter(target %in% targets) |> # take only specified targets
  # factor order treatments and targets for plotting
  mutate(target = fct_relevel(target, targets)) |> 
  mutate(treatment = fct_relevel(treatment, treatments))
# extract wt IC50 and calculate fold change of other targets
wt_IC50s <- data |>
  filter(target == "K562 pUltra BCR-ABL1 wt") |> 
  group_by(treatment) |> 
  summarize(wt_IC50_nM = IC50_nM)
data <- data |> 
  left_join(wt_IC50s, by = "treatment") |> 
  mutate(fold_vs_wt_IC50 = IC50_nM/wt_IC50_nM)
# bar plot of fold changes-----------------------------------------------
plot_type = "pdf"
color_scale <- "viridis"
vr <- viridis_range(length(treatments))
viridis_begin <- vr[1]
viridis_end <- vr[2]
x_min <- floor(min(log10(data$fold_vs_wt_IC50)))
x_max <- ceiling(max(log10(data$fold_vs_wt_IC50)))
breaks_x <- 10^rep(x_min : x_max)
minor_x <- rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9))
p <- data |> 
  filter(!target %in% c("K562 pUltra BCR-ABL1 wt", "K562 pUltra control")) |> 
  ggplot(aes(y = target, x = fold_vs_wt_IC50,
             fill = treatment,
             label = signif(fold_vs_wt_IC50, digits = 2))) +
  geom_bar(stat = "identity",
           position = position_dodge2(reverse = TRUE, padding = 0)) +
  geom_text(position = position_dodge2(width = 1, reverse = TRUE),
            hjust = -0.1,
            parse = FALSE) +
  scale_x_continuous(trans = "log10",
                     expand = expansion(mult = .1),
                     guide = "prism_offset_minor", # end at last tick
                     breaks = breaks_x,
                     labels = label_comma(accuracy = 1, big.mark = ""),
                     minor_breaks = minor_x,) +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis(option = color_scale,
                      discrete = TRUE,
                      begin = viridis_begin, end = viridis_end) +
  theme_prism() +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "fold change in IC50 vs wt")
save_plot(p, str_glue("output/fold_change_bar_{get_timestamp()}.{plot_type}"),
          width = 14, height = 8)
# strip plot of raw IC50s-----------------------------------------------
x_min <- floor(min(log10(data$IC50_nM)))
x_max <- ceiling(max(log10(data$IC50_nM)))
breaks_x <- 10^rep(x_min : x_max)
minor_x <- rep(1:9, x_max - x_min)*(10^rep(x_min:(x_max - 1), each = 9))
p <- data |> 
  ggplot(aes(y = target, x = IC50_nM,
             color = treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_x_continuous(trans = c("log10", "reverse"),
                     guide = "prism_offset_minor", # end at last tick
                     breaks = breaks_x,
                     labels = label_comma(accuracy = 1, big.mark = ""),
                     minor_breaks = minor_x,
                     expand = expansion(mult = .1)) +
  scale_y_discrete(limits = rev) +
  scale_color_viridis(option = color_scale,
                     discrete = TRUE,
                     begin = viridis_begin, end = viridis_end) +
  theme_prism() +
  theme(plot.background = element_blank()) + # need for transparent background
  labs(x = "IC50 (nM)")
save_plot(p, str_glue("output/IC50_bar_{get_timestamp()}.{plot_type}"),
          width = 14, height = 8)
