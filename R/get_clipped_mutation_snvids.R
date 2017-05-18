get_clipped_mutation_snvids <- function(data=NULL, hardclip = FALSE, softclip = FALSE) {
  if (is.null(data)) stop("Mandatory argument data is missing")
  if (hardclip) {
    data <- data %>% filter(in_hard_clips == FALSE)
    }
  if (softclip) {
    data <- data %>% filter(in_soft_clips == FALSE)
    }
  return(data$snvid)
  }
