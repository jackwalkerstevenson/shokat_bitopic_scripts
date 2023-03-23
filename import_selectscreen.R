import_selectscreen <- function(input_filename, treatments, targets, filter = !(missing(treatments) & missing(targets))){
  data <- read_csv(input_filename) %>%
    # rename relevant columns to common names
    rename(target = Kinase) %>%
    # tidy by pivoting duplicates to one measurement per row
    pivot_longer(cols = c(pct_inhibition_1, pct_inhibition_2), names_to = NULL, values_to = "pct_inhibition") %>%
    # wrangle: convert conc to log molar and convert percent inhibition to activity
    mutate(log.conc = log10(Compound_Conc_nM/1e9)) %>%
    mutate(activity = 100 - pct_inhibition)
  # filter for desired treatments and targets, if given, and order by given lists
  if(filter){
    if(!missing(treatments)){
      print("treatments not missing")
      data <- data %>%
        filter(treatment %in% treatments) %>%
        mutate(treatment = fct_relevel(treatment, treatments))
    }
    if(!missing(targets)){
      print("targets not missing")
      data <- data %>%
        filter(target %in% targets) %>%
        mutate(target = fct_relevel(target, targets))
    }
  }
  return(data)
}