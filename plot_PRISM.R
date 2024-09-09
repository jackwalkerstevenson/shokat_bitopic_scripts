#' ---
#'title: "plot PRISM"
#'author: "Jack Stevenson"
#'date: "2024-09-04"
#' ---
#'started 2024-09-04
#'this script plots results from the PRISM cell line screen
#'
# load required libraries------------------------------------------------------
library(tidyverse) # for tidy data handling
library(scales) # for fancy plotting scales
library(ggprism)  # for pretty prism-like plots
library(viridis) # for color schemes
library(doseplotr) # you bet
# set up-----------------------------------------------
rm(list = ls()) # clear environment
params_path <- "parameters/parameters_plot_PRISM.R"
source(params_path)
dir.create("output/", showWarnings = FALSE) # silently create output directory
# import and preprocess data-----------------------------------------------
# import fusion data: DepMap list of cell lines with ABL1 fusions
fusion_data <- readr::read_csv(fusions_filename,
                               name_repair = "universal") |> 
  # drop unneeded columns
  dplyr::select(all_of(c("Cell.Line", "Left.Gene", "Right.Gene"))) |> 
  # boolean for fusions that include BCR
  # note BCR sometimes left, sometimes right. seems inconsistent, so taking all
  dplyr::mutate(BCR = Left.Gene == "BCR" | Right.Gene == "BCR")

# import DRC data: data for each treatment/target curve fit
DRC_data <- readr::read_csv(input_filename) |> 
  dplyr::rename(treatment = varied_iname) |> # rename treatment column
  filter_validate_reorder("treatment", treatments) |> # check desired treatments
  # get simple cell line name from ccle_name
  dplyr::mutate(cell_line = stringr::str_split_i(ccle_name, "_", 1)) |> 
  # join with fusion data to annotate which cell lines have fusions
  dplyr::left_join(
    fusion_data,
    by = c("cell_line" = "Cell.Line"),
    relationship = "many-to-many") |>
  # fill out BCR column and make it a factor
  dplyr::mutate(BCR = tidyr::replace_na(BCR, FALSE) |>
                  factor(levels = c("TRUE", "FALSE"))) |> 
  # annotate all ABL1 fusions
  dplyr::mutate(fusion = Left.Gene == "ABL1" | Right.Gene == "ABL1") |>
  # annotate fusion types
  dplyr::mutate(fusion_type = case_when(
    fusion == TRUE & BCR == TRUE ~ "BCR::ABL1",
    fusion == TRUE & BCR == FALSE ~ "other ABL1 fusion",
    .default = "no ABL1 fusion") |> 
      factor(levels = c("BCR::ABL1", "other ABL1 fusion", "no ABL1 fusion"))
  )

# compare dose-response AUC to Riemann AUC for each treatment-------------------
for (trt in treatments){
  DRC_data |> 
    dplyr::filter(treatment == trt) |> 
    dplyr::mutate(missing_AUC = is.na(auc)) |> # annotate missing AUC
    dplyr::arrange(missing_AUC) |> 
    dplyr::mutate(auc = tidyr::replace_na(auc, 1)) |>  # replace missing AUC with 1
    ggplot(aes(x = auc, y = auc_riemann, color = missing_AUC)) +
    geom_point() +
    theme_prism() +
    theme(legend.title = element_text()) +
    coord_fixed() +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    labs(
      x = "AUC from curve fit",
      y = "Riemann AUC",
      title = "comparing AUC metrics",
      caption = trt,
      color = "missing curve-fit AUC")
  ggsave(str_glue("{output_dir}/AUC_comparison_{trt}_{doseplotr::get_timestamp()}.{plot_type}"),
         width = 8, height = 6)
}
# parameters for all jitter plots-----------------------------------------------
jitter_height = 0.3
alpha_emphasis = 1
alpha_background = 0.6
# jitter plot of cell line sensitivity by BCR or other fusion------------------
set.seed(random_seed)
DRC_data |> 
  ggplot(aes(x = auc,
             y = treatment,
             color = fusion_type,
             shape = fusion_type,
             alpha = fusion_type)) +
  geom_jitter(width = 0, height = jitter_height) + # only vertical jitter
  theme_prism() +
  scale_x_continuous(limits = c(0,1)) +
  scale_color_manual(values = c("BCR::ABL1" = "red",
                                "other ABL1 fusion" = "pink",
                                "no ABL1 fusion" = "black")) +
  scale_shape_manual(values = c("BCR::ABL1" = "triangle",
                                "other ABL1 fusion" = "triangle open",
                                "no ABL1 fusion" = "circle")) +
  scale_alpha_manual(values = c("BCR::ABL1" = alpha_emphasis,
                                "other ABL1 fusion" = alpha_emphasis,
                                "no ABL1 fusion" = alpha_background)) +
  labs(
    x = "AUC",
    y = "treatment",
    title = "Sensitivity of PRISM cell lines to treatments")
ggsave(str_glue("{output_dir}/jitter_fusion_seed_{random_seed}_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 10, height = 3)

# jitter plot of cell line sensitivity by BCR fusion only------------------
set.seed(random_seed)
DRC_data |> 
  ggplot(aes(x = auc, y = treatment, color = BCR, shape = BCR, alpha = BCR)) +
  geom_jitter(width = 0, height = jitter_height) + # only vertical jitter
  theme_prism() +
  scale_x_continuous(limits = c(0,1)) +
  scale_alpha_manual(values = c("TRUE" = alpha_emphasis,
                                "FALSE" = alpha_background),
                     labels = c("TRUE" = "BCR fusion",
                                "FALSE" = "no BCR fusion")) +
  scale_color_manual(values = c("TRUE" = "red",
                                "FALSE" = "black"),
                     labels = c("TRUE" = "BCR fusion",
                                "FALSE" = "no BCR fusion")) +
  scale_shape_manual(values = c("TRUE" = "triangle",
                                "FALSE" = "circle"),
                     labels = c("TRUE" = "BCR fusion",
                                "FALSE" = "no BCR fusion")) +
  labs(
    x = "AUC",
    y = "treatment",
    title = "Sensitivity of PRISM cell lines to treatments")
ggsave(str_glue("{output_dir}/jitter_BCR_seed_{random_seed}_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 10, height = 3)

# waterfall plot of cell line sensitivity by curve fit AUC--------------------
DRC_data |> 
  # only cell lines present for all treatments
  dplyr::group_by(cell_line) |> 
  dplyr::filter(all(!is.na(auc))) |> 
  dplyr::ungroup() |> 
  # rank order within treatment from highest to lowest AUC
  dplyr::group_by(treatment) |>
  dplyr::mutate(plot_rank = dplyr::row_number(-1 * auc)) |>
  ggplot(aes(x = plot_rank, y = auc, color = treatment)) +
  geom_point() +
  theme_prism() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme(legend.title = element_text()) +
  labs(x = "rank order",
       y = "AUC",
       title = "Sensitivity of PRISM cell lines to treatments")
ggsave(str_glue(
    "{output_dir}/waterfall_auc_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 8, height = 4)

# waterfall plot of cell line sensitivity by Riemann AUC----------------------
DRC_data |> 
  # rank order within treatment from highest to lowest Riemann AUC
  dplyr::group_by(treatment) |>
  dplyr::mutate(plot_rank = dplyr::row_number(-1 * auc_riemann)) |>
  ggplot(aes(x = plot_rank, y = auc_riemann, color = treatment)) +
  geom_point() +
  theme_prism() +
  scale_color_manual(values = color_map_treatments,
                     labels = display_names_treatments) +
  theme(legend.title = element_text()) +
  labs(x = "rank order",
       y = "Riemann AUC",
       title = "Sensitivity of PRISM cell lines to treatments")
ggsave(str_glue("{output_dir}/waterfall_auc_riemann_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 8, height = 4)
