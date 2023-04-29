# predict the activity of a treatment on a target at a given dose
predict_activity_drm <- function(data, trt, tgt, dose_nM){
  # models are fitted on log molar conc, so have to convert the desired nM conc
  log_dose <- nM_to_logM(dose_nM)
  model <- get_drm_plateplotr(data = data, trt = trt, tgt = tgt)
  print(model)
  # yucky syntax to get predict.drc to predict activity at given value of log.conc
  newdata <- data.frame(log.conc = log_dose)
  return(predict(model, newdata = newdata))
}