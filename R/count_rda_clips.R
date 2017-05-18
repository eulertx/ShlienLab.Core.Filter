count_rda_clips <- function(bam = NULL, rda = NULL, operation = "S", offsets = 10) {
  if(is.null(bam)) stop("Mandatory argument bam is missing")
  if(is.null(rda)) stop("Mandatory argument rda is missing")

  if (is.data.frame(rda)) {
    muts <- rda
  } else {
    snv <- load(rda)
    muts <- get(snv)
  }

  clips_df <- NULL
  for (row in 1:nrow(muts)) {
    new_row <- get_cigar_mutations(bam_file=bam, mut_file=muts[row,], chr=muts[row,"annovar_chr"],
                                   position=muts[row,"annovar_start"], operation=operation, offsets=offsets)
    clips_df <- rbind(clips_df, new_row)
  }
  return(clips_df)
}
