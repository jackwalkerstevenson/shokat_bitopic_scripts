#'plateplotr generates plots with dose-response curves from 96-well plate data.
#'input: data in plater format for tidy import
#'copy raw plate reader data into JWS Excel template and fill in input cells
#'save input file as csv for import
library(drc)  # for dose response curves
library(tidyverse)
library(ggprism)  # for pretty prism-like plots
library(viridis)
library(plater)  # for plate data import
plate_files = c("JS-B2-79 plater 1.csv",
                "JS-B2-79 plater 2.csv",
                "JS-B2-79 plater 3.csv",
                "JS-B2-79 plater 4.csv",
                "JS-B2-79 plater 5.csv")
plate_names = c("Dasatinib",
                "DasatiLink-1",
                "Ponatinib",
                "PonatiLink-1",
                "Asciminib")
plater_data <- read_plates(plate_files, plate_names) %>%
  filter(drug != "N/A") %>%
  # drop 0 values for plotting and curve fitting
  # note this is only OK because normalization is happening in Excel
  filter(concentration != 0) %>%
  mutate(log.conc = log10(concentration/1e6))  # convert conc from µM to M
for(p in distinct(plater_data["Plate"]))
  {print(str_glue("working on plate {p}"))
  plate.summary <- plater_data %>%
    group_by(Plate, cell_line, concentration) %>%
    summarize(
      sem = sd(CTG_normalized, na.rm = TRUE)/sqrt(n()),
      CTG_normalized = mean(CTG_normalized),
      log.conc = log.conc
    )
  
  }



plate.summary %>% ggplot(aes(x = log.conc, y = CTG_normalized, color = cell_line))+
    geom_point()+
    geom_errorbar(aes(ymax = CTG_normalized+sem, ymin = CTG_normalized-sem))+
    geom_smooth(method = "drm", method.args = list(fct = L.4()), se = FALSE)+
    theme_prism()+
    scale_color_manual(values = c("darkred","black"))+
    labs(x = "Log [compound] (M)",
         y = "Relative cell viability (%)",
         title = "Dasatinib")
#ggsave("plots/xxx.pdf", width = 5, height = 4)
