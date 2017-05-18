get_position_reads <- function(bam_file=NULL, chr=NULL, position=NULL, operation="S", offsets=NULL) {

  if(is.null(bam_file)) stop("Mandatory argument bam_file is missing")
  if(is.null(chr)) stop("Mandatory argument chr is missing")
  if(is.null(position)) stop("Mandatory argument position is missing")
  if(is.null(offsets)) stop("Mandatory argument offsets is missing")

  s_cigars <- load_bam_reads(bam_file=bam_file, chr=chr, position=position, operation=operation, offsets=offsets)

  s_clips <- NULL

  # calculate whether the given operation on the read overlaps the given position
  if (nrow(s_cigars) > 0) {
    for (x in 1:nrow(s_cigars)) {
      s_frame <- cigar_to_data_frame(cigar=s_cigars$cigar[x], position=s_cigars$pos[x])

      # if the operation is S and a read starts with S, update position
      if (operation=="S" & s_frame[1,"type"]=="S") {
        s_frame <- cigar_to_data_frame(cigar=s_cigars$cigar[x],
                                       position=s_cigars$pos[x] - as.numeric(as.character(s_frame[1,"size"])))
      }

      s_frame <- s_frame %>% filter(type==operation)
      s_frame$start <- as.numeric(as.character(s_frame$start))
      s_frame$end <- as.numeric(as.character(s_frame$end))

      # if the operation is hard clipping, mark the entire cigar read as hard clipped
      if (operation=="H") {
        s_clips <- rbind(s_clips, s_cigars[x,])
      }

      # if only one side is soft clipped, add if the operation overlaps the position
      else if (nrow(s_frame) == 1) {
        if (s_frame$start <= position & position <= s_frame$end) {
          s_clips <- rbind(s_clips, s_cigars[x,])
        }
      }

      # add if soft clipped on both sides and one side of the soft clips overlaps the position
      else if ((s_frame[1,"start"] <= position & position <= s_frame[1,"end"]) |
               (s_frame[2,"start"] <= position & position <= s_frame[2,"end"])) {
        s_clips <- rbind(s_clips, s_cigars[x,])
      }
    }
  }

  return(s_clips)
}
