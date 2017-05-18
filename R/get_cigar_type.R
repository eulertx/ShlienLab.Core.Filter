get_cigar_type <- function(operation = NULL) {
  if(is.null(operation)) stop("Mandatory argument 'operation' is missing")

  if (operation=="M") return("alignment_matches")
  if (operation=="I") return("insertions")
  if (operation=="D") return("deletions")
  if (operation=="N") return("skipped_regions")
  if (operation=="S") return("soft_clips")
  if (operation=="H") return("hard_clips")
  if (operation=="P") return("padding")
  if (operation=="=") return("sequence_matches")
  if (operation=="X") return("sequence_mismatches")
}
