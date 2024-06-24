# plot large-scale SelectScreen single-point data using KinMap and ggplot2
# Jack Stevenson started 2024-05
# load required libraries------------------------------------------------------
library(tidyverse)
library(purrr)
library(doseplotr)
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
# import and process selectscreen results--------------------------------------
raw_single_pt_data <- import_selectscreen(input_SS_single_pt_filename)
# average replicates for Kinmap plotting
avg_single_pt_data <- raw_single_pt_data |> 
  group_by(treatment, target) |> 
  summarize(pct_inhibition = mean(pct_inhibition))
# pivot to one row per kinase for plotting
kinase_data <- raw_single_pt_data |> 
  select(treatment, target, pct_inhibition) |> 
  tidyr::pivot_wider(names_from = treatment,
                     names_prefix = "pct_inhibition_",
                     values_from = pct_inhibition,
                     values_fn = list) |> 
  # for all pct_inhibition columns, make mean and sem columns from them
  dplyr::mutate(dplyr::across(matches("pct_inhibition")))
  #dplyr::mutate(mean_pct_inhibition_ponatinib =
  #                purrr::map_dbl(pct_inhibition_ponatinib, mean))
  #dplyr::mutate(dplyr::across(matches("pct_inhibition"), ~ mean(.x)))
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
ggplot(kinase_data, aes(x = pct_inhibition_ponatinib,
                        y = `pct_inhibition_PonatiLink-2-7-10`)) +
  geom_point()