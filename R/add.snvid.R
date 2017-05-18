add.snvid <- function(data=NULL) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  data <- data %>% dplyr::mutate(snvid = paste(sep='_', annovar_chr, annovar_start, annovar_end, annovar_ref, annovar_alt))
  return(data)
  }
