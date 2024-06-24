# plot large-scale SelectScreen single-point data using KinMap and ggplot2
# Jack Stevenson started 2024-05
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(scales) # for nice breaks
library(ggprism) # for prism theme
library(ggrepel) # for labeling points
# set up input and output directories------------------------------------------
input_dir <- "input"
output_dir <- "output"
dir.create(input_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)
# set aesthetic parameters for KinMap output-----------------------------------
kinmap_global_shape <- 0
kinmap_global_stroke_color <- "black"
kinmap_global_stroke_width <- 1.5
kinmap_ontarget_fill_color <- "green"
kinmap_offtarget_fill_color <- "red"
kinmap_uninhibited_fill_color <- "gray"
# set parameters specific to this dataset--------------------------------------
on_targets <- c("ABL1")
input_SS_single_pt_filename <-
  str_glue("{input_dir}/selectscreen combined results 67309 ZLYTE duplicates removed 2024-05-15.csv")
top_hit_threshold <- 75
# import and process selectscreen results--------------------------------------
raw_single_pt_data <- import_selectscreen(input_SS_single_pt_filename)
# average replicates for Kinmap plotting
avg_single_pt_data <- raw_single_pt_data |> 
  group_by(treatment, target) |> 
  summarize(pct_inhibition = mean(pct_inhibition))

kinase_plot_data <- raw_single_pt_data |> 
  select(treatment, target, pct_inhibition) |> 
  summarize(pct_inhibition = mean(pct_inhibition), .by = c(treatment, target)) |> 
  filter(pct_inhibition > -25) |> # remove large negative outliers
  # pivot to one row per kinase for ggplot plotting
  tidyr::pivot_wider(names_from = treatment,
  names_prefix = "pct_inhibition_",
  values_from = pct_inhibition) |> 
  tidyr::drop_na() |> 
  mutate(top_hit = if_all(starts_with("pct_inhibition"),
                          \(x) x > top_hit_threshold))
# add KinMap directive parameters----------------------------------------------
avg_single_pt_data <- avg_single_pt_data |> 
  dplyr::mutate(
    # add kinmap parameters
    # global aesthetic parameter
    shape = kinmap_global_shape,
    # scale size by pct_inhibition with floor
    size = ifelse(pct_inhibition < 10, 10, pct_inhibition),
    # color on vs off vs uninhibited target
    fill = case_when(
      pct_inhibition < 10 ~ kinmap_uninhibited_fill_color,
      target %in% on_targets ~ kinmap_ontarget_fill_color,
      .default = kinmap_offtarget_fill_color),
    # more global aesthetic parameters
    stroke = kinmap_global_stroke_color,
    strokeWidth = kinmap_global_stroke_width)
# write KinMap directive output table------------------------------------------
write_output_treatment <- function(data, trt){
  data |>
    dplyr::filter(treatment == trt) |> 
    write_csv(file = str_glue("{output_dir}/kinmap_output_{treatment}_{get_timestamp()}.csv"))
}
treatments <- unique(avg_single_pt_data$treatment)
for(treatment in treatments){
  write_output_treatment(avg_single_pt_data, treatment)
}
# plot 2 treatments against each other-----------------------------------------
kinase_plot_data |>
  ggplot(aes(
    y = `pct_inhibition_ponatinib + asciminib`,
    x = `pct_inhibition_PonatiLink-2-7-10`
  )) +
  geom_point(aes(color = top_hit), size = 2, alpha = 0.6) +
  geom_label_repel(
    data = kinase_plot_data |>
      filter(if_any(starts_with("pct_inhibition"), \(x) x >93)),
    aes(label = target),
    # box.padding = 0.1,
    label.padding = 0.15,
    segment.size = .3,
    force = 5,
    min.segment.length = 0,
    max.overlaps = 10) +
  labs(title = "SelectScreen wild-type kinase inhibition",
       caption = "at IC90 for ABL1 T315I",
       x = "% inhibition by PonatiLink-2",
       y = "% inhibition by ponatinib + asciminib") +
  scale_x_continuous(breaks = breaks_width(25)) +
  scale_y_continuous(
    breaks = breaks_width(25),
    expand = expansion(mult = c(0.02, 0.10))) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  theme_prism() +
  theme(panel.grid.major = element_line(linewidth = 0.2, linetype = "dashed"),
        legend.position = "none")
ggsave(str_glue("{output_dir}/scatter_plot_{get_timestamp()}.pdf"),
       width = 7, height = 7)