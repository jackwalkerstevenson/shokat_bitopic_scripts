# helper function for summarizing replicate data for plotting------------------
summarize_activity <- function(x){
  summarize(x,
            # standard error for error bars = standard deviation / square root of n
            confidence_intervals = list(mean_cl_normal(activity) |>
                        rename(mean_activity = y, ymin_activity = ymin, ymax_activity = ymax)),
            sem = sd(activity, na.rm = TRUE)/sqrt(n()),
            # get mean normalized readout value for plotting
            mean_read = mean(activity),
            w = 0.06 * n() # necessary for consistent error bar widths across plots
  ) |>
    tidyr::unnest(cols = confidence_intervals)
}