flag.complex.overlap.loci <- function(data=NULL, minimum=0, far.dist=999999) {
  if (is.null(data)) stop("Mandatory argument is missing")

  # NA values in the distance_to_low_complexity_1 are far from the complex region
  # set NA values of distance_to_low_complexity_1 to far.dist
  data[is.na(data$distance_to_low_complexity_1),]$distance_to_low_complexity_1 <- far.dist
  data$dist_low_complexity <- 0
  data[abs(data$distance_to_low_complexity_1) <= minimum,]$dist_low_complexity <- 'FLAGGED'

  return(data)
  }
