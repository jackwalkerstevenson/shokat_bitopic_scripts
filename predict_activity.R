# predict the activity of a treatment on a target at a given dose
predict_activity <- function(data, trt, tgt, dose_nM){
  # models are fitted on log molar conc, so have to convert the desired nM conc
  log_dose <- log10(dose_nM/1e9)
  model <- get_drm(data = data, trt = trt, tgt = tgt)
  # yucky syntax to get predict.drc to predict activity at given value of log.conc
  newdata <- data.frame(log.conc = log_dose)
  return(predict(model, newdata = newdata))
}