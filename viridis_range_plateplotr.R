# pick a visually pleasing range of the default viridis color palette for discrete values
# returns c(start, end) values
# 2 values should be 0.4ish
viridis_range_plateplotr <- function(n, reverse = TRUE){
  range <- ifelse(n < 6, 0.2*n, 1)
  high <- 0.5 + range / 2
  low <- 0.5 - range / 2
  if(reverse){
    return(c(high, low))
  }
  else{
    return(c(low, high))
  }
}