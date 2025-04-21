# plot PK from tumor growth inhibition study
# Jack Stevenson started 2024-10-09
# load required libraries------------------------------------------------------
library(tidyverse)
library(doseplotr) # you bet
library(ggprism) # for prism theme
# set up-----------------------------------------------
rm(list = ls()) # clear environment
# read parameter file
params_path <- "parameters/parameters_plot_PK_TGI.R"
source(params_path)
input_path <- str_glue("{input_directory}{input_filename}")
# write timestamped params to output
doseplotr::file_copy_to_dir(params_path, output_directory)
# write timestamped code to output
doseplotr::file_copy_to_dir("plot_in_vivo/plot_PK_TGI.R", output_directory)
# write timestamped input files to output
doseplotr::file_copy_to_dir(input_path, output_directory)
# set up input and output directories
dir.create(input_directory, showWarnings = FALSE)
dir.create(output_directory, showWarnings = FALSE)
# import, preprocess and report data-----------------------------------------------
PK_data <- readxl::read_excel(
  str_glue("{input_directory}/{input_filename}"),
) |> 
  dplyr::mutate(measurement = fct_inorder(as.factor(measurement)),
                treatment = fct_inorder(as.factor(treatment)))
# plot PK data-----------------------------------------------------------------
PK_data |> 
  group_by(treatment, measurement) |> 
  dplyr::summarize(
    mean_conc_nM = mean(conc_nM),
    # standard error for error bars = standard deviation / square root of n
    sem = sd(conc_nM, na.rm = TRUE)/sqrt(n()),
  ) |> 
  ggplot(aes(x = treatment,
             y = mean_conc_nM,
             fill = treatment,
             group = measurement)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  geom_errorbar(aes(ymin = mean_conc_nM - sem,
                    ymax = mean_conc_nM + sem),
                position = "dodge",
                linewidth = 0.3) +
  geom_label(aes(label = measurement), 
             position = position_dodge(width = 0.9), 
             vjust = -0.5,
             show.legend = FALSE) + 
  scale_fill_manual(values = color_map_treatments) +
  theme_prism() +
  labs(x = "",
       y = "plasma concentration (nM)",
       title = "Peak and trough compound concentrations",
       caption = "after 14 d QD weekday dosing")
ggsave(
  str_glue(
    "{output_directory}/plot_TGI_PK_{get_timestamp()}.{plot_type}",
    bg = "transparent"
  ),
  width = 9, height = 6,
  bg = "transparent")