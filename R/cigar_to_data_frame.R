cigar_to_data_frame <- function(cigar = NULL, position = NULL)
{
  if(is.null(cigar)) stop("Mandatory argument cigar is missing")
  if(is.null(position)) stop("Mandatory argument position is missing")

  # get seperate lists of numbers (lengths) and characters (CIGAR operations)
  alpha <- unlist(strsplit(cigar, "[0-9]+"))
  alpha <- tail(alpha, length(alpha)-1)
  num <- as.numeric(unlist(strsplit(cigar, "[MIDNSHP=X]+")))

  # calculate each start and end position based on whole position and previous string
  start <- position
  if (length(alpha) > 1) {
    for (x in 2:length(alpha)) {
      start[x] <- start[x-1] + num[x-1]
    }
  }
  end <- start + num - 1

  # return entire data frame
  data <- as.data.frame(cbind(type=alpha, size=num, start=start, end=end))
  data$type <- sapply(data$type, as.character)
  return(data)
}
