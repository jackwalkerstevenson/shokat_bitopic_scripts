import_selectscreen <- function(input_filename, compounds, targets){
  read_csv(input_filename) %>%
    # rename relevant columns to common names
    rename(target = Kinase) %>%
    # filter for desired compounds and targets and order by input lists
    filter(compound %in% compounds) %>%
    mutate(compound = fct_relevel(compound, compounds)) %>%
    filter(target %in% targets) %>%
    mutate(target = fct_relevel(target, targets)) %>%
    # tidy by pivoting duplicates to one measurement per row
    pivot_longer(cols = c(pct_Inhibition_1, pct_Inhibition_2), names_to = NULL, values_to = "pct_inhibition") %>%
    # wrangle: convert conc to log molar and convert percent inhibition to activity
    mutate(log.conc = log10(Compound_Conc_nM/1e9)) %>%
    mutate(activity = 100 - pct_inhibition)
}