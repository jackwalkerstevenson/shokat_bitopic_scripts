#' ---
#'title: "plateplotr"
#'author: "Jack Stevenson"
#'date: "2022-07-01"
#' ---
#'plateplotr generates dose-response plots from plate reader data.
#' 
#'Input: data formatted as per the plater package. Use the JWS Excel plater template to make this easy.
#' 
#'Variables expected in the import file (JWS Excel template is already set up for this):
#' 
#'- 'compound': name of compound used (including DMSO wells in dilution series)
#'- 'conc': concentration of compound used
#'- 'cell_line': cell line treated
#'- 'readout': raw plate reader data (this is only used to calculate read_norm)
#'- 'read_norm': normalized plate reader data (calculated in the spreadsheet)
#'- 'replicate': replicate of curve (not actually used)
#'
#'How to use plateplotr:
#'
#'1. Make a copy of the JWS Excel plater template for each plate you want to import
#'2. Fill out each copy with the compound + concs used and paste in raw plate reader data
#'3. Save input sheets as csv and put them in your working directory (or move your working directory to where they are)
#'  - alternatively, manually create plater-formatted csv sheets and put them in the working directory
#'4. Fill out input filenames in plate_files below
#'5. Run plot.R (this file)
#'6. Output plots will be created in a "plots" folder in the working directory

# Load required libraries------------------------------------------------------
library(drc)  # for dose response curves
library(tidyverse) # for tidy data handling
library(ggprism)  # for pretty prism-like plots
library(plater)  # for tidy importing of plate data

# specify names of input files and import data---------------------------------
plate_files = c("CSVs/SP-B1-01 Asciminib Read 2.csv",
                "CSVs/SP-B1-01 Ponatinib Read 2.csv",
                "CSVs/SP-B1-01 PonatiLink-1-20 Read 2.csv")
plate_names = seq(1,length(plate_files))  # create plate IDs
plate_data <- read_plates(plate_files, plate_names) %>%
  # drop empty wells
  filter(compound != "N/A") %>%
  # drop 0 values before plotting and curve fitting
  # note this is only OK because normalization happens before import
  filter(conc != 0) %>%
  mutate(log.conc = log10(conc/1e6))  # convert conc from µM to M

# analyze and plot data for each compound--------------------------------------
for (cpd in distinct(plate_data["compound"])$compound){
  print(str_glue("working on compound: {cpd}"))
  # get cell lines in import order
  cell_line_factors <- unique(plate_data[["cell_line"]])
  plate.summary <- plate_data %>%
    # set factors so cell lines get plotted and colored in input order
    mutate(cell_line = fct_relevel(cell_line, cell_line_factors)) %>%
    filter(compound == cpd) %>% # get data from one compound to work with
    group_by(cell_line, conc) %>%  # get set of replicates for each condition
    # create summary table for plotting
    summarize(
      # standard error for error bars = standard deviation / square root of n
      sem = sd(read_norm, na.rm = TRUE)/sqrt(n()),
      # get mean normalized readout value for plotting
      mean_read = mean(read_norm),
      # carry concentration through for plotting
      log.conc = log.conc
      )

  # plot data and fit dose response curve
  plate.summary %>%
    ggplot(aes(x = log.conc, y = mean_read, color = cell_line)) +
      geom_point() +
      # error bars = mean plus or minus standard error
      geom_errorbar(aes(ymax = mean_read+sem, ymin = mean_read-sem)) +
      # use drm method from drc package to fit dose response curve
      geom_smooth(method = "drm", method.args = list(fct = L.4()), se = FALSE) +
      scale_color_manual(values = c("black","darkred")) +
      scale_x_continuous(breaks = ) +
      scale_y_continuous(breaks = c(0,25,50,75,100)) +
      coord_cartesian(xlim = c(-11,-5),
                      ylim = c(0,NA)) +
      theme_prism() + # make it look like prism
      theme(plot.background = element_blank()) +
      labs(x = "Log [compound] (M)",
           y = "Relative cell viability (%)",
           title = cpd)
  # save plot with manually optimized aspect ratio
  ggsave(str_glue("Plot Output/{cpd}.png"), width = 5, height = 4, bg = "transparent")
  print(str_glue("Done plotting compound {cpd}"))
  }
