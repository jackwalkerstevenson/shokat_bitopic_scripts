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
library(ggprism)  # for pretty prism-like plots
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
  dplyr::mutate(BCR_fusion = Left.Gene == "BCR" | Right.Gene == "BCR") |> 
  # deduplicate cell lines: remove all but one fusion, keeping BCR if present
  dplyr::group_by(Cell.Line) |> 
  dplyr::arrange(Cell.Line, desc(BCR_fusion)) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

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
  # annotate all ABL1 fusions
  dplyr::mutate(ABL1_fusion = Left.Gene == "ABL1" | Right.Gene == "ABL1") |>
  dplyr::mutate(NUP214_fusion = Left.Gene == "NUP214" | Right.Gene == "NUP214") |>
  # annotate BCR::ABL1 or NUP214::ABL1
  dplyr::mutate(fusion_type = case_when(
    ABL1_fusion == TRUE & BCR_fusion == TRUE ~ "BCR::ABL1",
    ABL1_fusion == TRUE &  NUP214_fusion == TRUE ~ "NUP214::ABL1",
    .default = "other") |> 
      factor(levels = c("BCR::ABL1", "NUP214::ABL1", "other"))) |> 
  # old annotations follow
  dplyr::mutate(BCR_or_NUP214 = dplyr::if_any(c(Left.Gene, Right.Gene),
                                              ~ . %in% c("BCR", "NUP214")) |> 
                  factor(levels = c("TRUE", "FALSE"))) |> 
  # fill out BCR column and make it a factor
  dplyr::mutate(BCR_fusion = tidyr::replace_na(BCR_fusion, FALSE) |>
                  factor(levels = c("TRUE", "FALSE")))
  # # annotate fusion types
  # dplyr::mutate(fusion_type = case_when(
  #   ABL1_fusion == TRUE & BCR_fusion == TRUE ~ "BCR::ABL1",
  #   ABL1_fusion == TRUE & BCR_fusion == FALSE ~ "other ABL1 fusion",
  #   .default = "no ABL1 fusion") |> 
  #     factor(levels = c("BCR::ABL1", "other ABL1 fusion", "no ABL1 fusion"))
# report input, raw data and parameters-----------------------------------------
write_csv(DRC_data, fs::path(output_dir,
                               str_glue("PRISM_raw_data_{get_timestamp()}.csv")))
doseplotr::file_copy_to_dir("plot_PRISM.R", output_dir)
doseplotr::file_copy_to_dir(params_path, output_dir)
doseplotr::file_copy_to_dir(scales_path, output_dir)
doseplotr::file_copy_to_dir(input_filename, output_dir)
doseplotr::file_copy_to_dir(fusions_filename, output_dir)
# find fusion cell lines that are in PRISM-----------------------------------------------
PRISM_fusions <- fusion_data |> 
  dplyr::left_join(DRC_data, by = c("Cell.Line" = "cell_line")) |> 
  dplyr::filter(!is.na(max_dose)) |> 
  dplyr::slice(1, .by = Cell.Line)

# plot dose-response AUC vs Riemann AUC for each treatment-------------------
for (trt in treatments){
  DRC_data |> 
    dplyr::filter(treatment == trt) |> # pick out specific treatment
    dplyr::mutate(missing_AUC = is.na(auc)) |> # annotate missing AUC
    dplyr::arrange(missing_AUC) |> # put missing AUC first
    dplyr::mutate(auc = tidyr::replace_na(auc, 1)) |>  # replace missing AUC with 1
    ggplot(aes(x = auc, y = auc_riemann, color = missing_AUC)) +
    geom_point() +
    theme_prism() + # make it look like prism
    theme(legend.title = element_text()) + # reinstate legend title
    coord_fixed() + # same x and y axis scale
    scale_x_continuous(limits = c(0,1)) + # set bounds of axes
    scale_y_continuous(limits = c(0,1)) +
    # set colors for color variable (missing_AUC)
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    labs( # label plot
      x = "AUC from curve fit",
      y = "Riemann AUC",
      title = "comparing AUC metrics",
      caption = trt,
      color = "missing curve-fit AUC")
  # save plot to a file
  ggsave(str_glue("{output_dir}/AUC_comparison_{trt}_{doseplotr::get_timestamp()}.{plot_type}"),
         width = 8, height = 6)
}
# set parameters for jitter plots-----------------------------------------------
jitter_height = 0.3
alpha_emphasis = 1
alpha_background = 0.6
# jitter plot of cell line sensitivity by BCR/NUP214::ABL1----------------------
set.seed(random_seed) # set a manual random seed so output always the same
DRC_data |> 
  dplyr::arrange(desc(fusion_type)) |> # fusion type determines order
  ggplot(aes(x = auc,
             y = treatment,
             color = fusion_type,
             shape = fusion_type,
             alpha = fusion_type)) +
  geom_jitter(width = 0, height = jitter_height) + # only vertical jitter
  theme_prism() +
  scale_x_continuous(limits = c(1,0), transform = "reverse") + # flip axis
  scale_y_discrete(labels = display_names_treatments) +
  # have to set all the scales so the legend will combine into one
  scale_color_manual(values = c("BCR::ABL1" = "red",
                                "NUP214::ABL1" = "pink",
                                "other" = "black")) +
  scale_shape_manual(values = c("BCR::ABL1" = "triangle",
                                "NUP214::ABL1" = "triangle open",
                                "other" = "circle")) +
  scale_alpha_manual(values = c("BCR::ABL1" = alpha_emphasis,
                                "NUP214::ABL1" = alpha_emphasis,
                                "other" = alpha_background)) +
  labs(
    x = "AUC",
    y = "treatment",
    title = "Sensitivity of PRISM cell lines to treatments")
ggsave(str_glue("{output_dir}/jitter_BCR_NUP214_seed_{random_seed}_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 10, height = 3)

# jitter plot of cell line sensitivity by BCR fusion only-----------------------
set.seed(random_seed)
DRC_data |> 
  dplyr::arrange(desc(BCR_fusion)) |>
  ggplot(aes(x = auc,
             y = treatment,
             color = BCR_fusion,
             shape = BCR_fusion,
             alpha = BCR_fusion)) +
  geom_jitter(width = 0, height = jitter_height) + # only vertical jitter
  theme_prism() +
  scale_x_continuous(limits = c(1,0), transform = "reverse") +
  scale_y_discrete(labels = display_names_treatments) +
  scale_alpha_manual(values = c("TRUE" = alpha_emphasis,
                                "FALSE" = alpha_background),
                     labels = c("TRUE" = "BCR::ABL1",
                                "FALSE" = "not BCR::ABL1")) +
  scale_color_manual(values = c("TRUE" = "red",
                                "FALSE" = "black"),
                     labels = c("TRUE" = "BCR::ABL1",
                                "FALSE" = "not BCR::ABL1")) +
  scale_shape_manual(values = c("TRUE" = "triangle",
                                "FALSE" = "circle"),
                     labels = c("TRUE" = "BCR::ABL1",
                                "FALSE" = "not BCR::ABL1")) +
  labs(
    x = "AUC",
    y = "treatment",
    title = "Sensitivity of PRISM cell lines to treatments")
ggsave(str_glue("{output_dir}/jitter_BCR_seed_{random_seed}_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 10, height = 3)

# jitter plot of cell line sensitivity by BCR::ABL1 or NUP214::ABL1-------------
set.seed(random_seed)
DRC_data |> 
  dplyr::arrange(desc(BCR_or_NUP214)) |>
  ggplot(aes(x = auc,
             y = treatment,
             color = BCR_or_NUP214,
             shape = BCR_or_NUP214,
             alpha = BCR_or_NUP214)) +
  geom_point(position = position_jitter(width = 0, # only vertical jitter
                                        height = jitter_height,
                                        seed = random_seed)) + 
  theme_prism() +
  scale_x_continuous(limits = c(1,0), transform = "reverse") +
  scale_y_discrete(labels = display_names_treatments) +
  scale_alpha_manual(values = c("TRUE" = alpha_emphasis,
                                "FALSE" = alpha_background),
                     labels = c("TRUE" = "BCR::ABL1 or NUP214::ABL1",
                                "FALSE" = "other")) +
  scale_color_manual(values = c("TRUE" = "red",
                                "FALSE" = "black"),
                     labels = c("TRUE" = "BCR::ABL1 or NUP214::ABL1",
                                "FALSE" = "other")) +
  scale_shape_manual(values = c("TRUE" = "triangle",
                                "FALSE" = "circle"),
                     labels = c("TRUE" = "BCR::ABL1 or NUP214::ABL1",
                                "FALSE" = "other")) +
  labs(
    x = "AUC",
    y = "treatment",
    title = "Sensitivity of PRISM cell lines to treatments")
ggsave(str_glue("{output_dir}/jitter_BCR_or_NUP214_seed_{random_seed}_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 11, height = 3)

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
  theme(legend.title = element_text()) + # reinstate legend label
  labs(x = "rank order of sensitivity",
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
  theme(legend.title = element_text()) + # reinstate legend label
  labs(x = "rank order of sensitivity",
       y = "AUC (Riemann)",
       title = "Sensitivity of PRISM cell lines to treatments")
ggsave(str_glue("{output_dir}/waterfall_auc_riemann_{doseplotr::get_timestamp()}.{plot_type}"),
       width = 8, height = 4)
