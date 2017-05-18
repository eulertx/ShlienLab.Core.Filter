md_snv_position <- function(md = NULL, position = NULL, cigar = NULL, seq = NULL, operation = "S") {
  if(is.null(md)) stop("Mandatory argument md is missing")
  if(is.null(position)) stop("Mandatory argument position is missing")
  if(is.null(cigar)) stop("Mandatory argument cigar is missing")
  if(is.null(seq)) stop("Mandatory argument seq is missing")

  start_pos <- position
  pos_vector <- NULL
  alts <- NULL
  cigar_df <- cigar_to_data_frame(cigar = cigar, position = position)

  if (operation != "S") {
    # does md contain mutation, or is it a matched read (not including insertions)
    if (!grepl("[TCAG]", md)) {
      return("unmutated")
    }

    # does the md read contain ^ (check for deletions)
    if (grepl("[[:punct:]]", md)) {
      return("indel")
    }

    #do something to make sure there are not 2 mutations after one another
    alpha <- unlist(strsplit(md, "[0-9]+"))
    num <- unlist(strsplit(md, "[TCAG]+"))

    doubles <- alpha[nchar(alpha) > 1]
    for(d in doubles) {
      d_vector <- unlist(strsplit(d, ""))
      newd <- paste(d_vector, 0, sep="", collapse="")
      index <- grep(d, alpha)
      alpha[index] <- newd
    }

    while (length(alpha) > length(num)) {
      num <- c(num, "")
    }

    while (length(alpha) < length(num)) {
      alpha <- c(alpha, "")
    }

    if(alpha[1]=="") {
      newmd <- paste0(alpha, num)
    } else {
      newmd <- paste0(num, alpha)
    }
    md <- paste(newmd, collapse="")
    alpha <- unlist(strsplit(md, "[0-9]+"))
    num <- unlist(strsplit(md, "[TCAG]+"))

    # does the cigar string have the same length as the md read (check for insertion)
    clipped_df <- cigar_df[cigar_df$type != "H",]
    clipped_df <- clipped_df[clipped_df$type != "S",]
    cigar_length <- sum(as.numeric(as.character(clipped_df$size)))
    md_length <- sum(as.numeric(num), na.rm = TRUE) + length(alpha[alpha != ""])
    if (md_length != cigar_length) {
      return("indel")
    }

    # shave off the unnecessary front of alpha and num, and update starting position
    if (alpha[1] == "") {
      position <- position + as.numeric(num[1])
      num <- tail(num, length(num)-1)
      alpha <- tail(alpha, length(alpha)-1)
    } else {
      num <- tail(num, length(num)-1)
    }

    # calculate the position of each of the mutations
    for (mut in alpha) {
      pos_vector <- c(pos_vector, position)
      alts <- c(alts, substr(seq, position-start_pos+1, position-start_pos+1))
      position <- position + as.numeric(num[1]) + 1
      num <- tail(num, length(num)-1)
    }
  }

  # for soft clips, add every mutation within the soft clip
  else {
    # what happens with indels?

    # get all bases in soft clip at the beginning of read
    if (cigar_df[1,"type"]=="S") {
      clip_len <- cigar_df[1,"size"]
      for (base_pos in 1:as.numeric(as.character(clip_len))) {
        alts <- c(alts, substr(seq, base_pos, base_pos))
        pos_vector <- c(pos_vector, start_pos + base_pos - 1)
      }
    }

    # get all bases in soft clip at the end of read
    if (cigar_df[nrow(cigar_df),"type"]=="S") {
      clip_len <- cigar_df[nrow(cigar_df),"size"]
      clip_start <- sum(as.numeric(as.character(cigar_df$size))) - as.numeric(as.character(clip_len)) + 1
      for (base_pos in clip_start:(clip_start + as.numeric(as.character(clip_len)))) {
        alts <- c(alts, substr(seq, base_pos, base_pos))
        pos_vector <- c(pos_vector, start_pos + base_pos - 1)
      }
    }
  }

  #return a data frame of the
  mutation_df <- as.data.frame(cbind(alt=alts, position=pos_vector))
  return(mutation_df)
}
