get.recurrent.snvids <- function(data=NULL, id='snvid') {
  if (is.null(data)) stop("Mandatory argument data is missing")
  if (is.null(id)) stop("Mandatory argument id is missing")

  data.counts <- data %>% group_by_(id) %>% summarise(pon_count=n())

  return(data.counts)
  }
