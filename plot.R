#'plateplotr generates plots with dose-response curves from 96-well plate data.
#'input: data in plater format for tidy import
#'copy raw plate reader data into JWS Excel template and fill in input cells
#'save input file as csv for import
library(drc)  # for dose response curves
library(tidyverse)
library(ggprism)  # for pretty prism-like plots
library(viridis)
library(plater)  # for plate data import
plater_data <- read_plate("JS-B2-79 plater test.csv") %>%
  filter(drug != "N/A") %>%
  # drop 0 values for plotting and curve fitting
  # note this is only OK because normalization is happening in Excel
  filter(concentration != 0) %>%
  mutate(log.conc = log10(concentration/1e6))  # convert conc from µM to M
plater_data.summary <- plater_data %>%
    group_by(cell_line, concentration) %>%
    summarize(
      sem = sd(CTG_normalized, na.rm = TRUE)/sqrt(n()),
      CTG_normalized = mean(CTG_normalized),
      log.conc = log.conc,
    )

plater_data.summary %>% ggplot(aes(x = log.conc, y = CTG_normalized, color = cell_line))+
    geom_point()+
    geom_errorbar(aes(ymax = CTG_normalized+sem, ymin = CTG_normalized-sem))+
    geom_smooth(method = "drm", method.args = list(fct = L.4()), se = FALSE)+
    theme_prism()+
    scale_color_manual(values = c("darkred","black"))+
    labs(x = "Log [compound] (M)",
         y = "Relative cell viability (%)",
         title = "Dasatinib")
