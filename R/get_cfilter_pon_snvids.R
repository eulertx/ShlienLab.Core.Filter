get_cfilter_pon_snvids <- function(data=NULL, pon_max = 1, max_number_of_flagged=2) {
  if (is.null(data)) stop("Mandatory argument data is missing")
  data <- data %>%
    filter(pon_count <= pon_max) %>%
    filter(flagged_count <= max_number_of_flagged)
  return(data$snvid)
  }
