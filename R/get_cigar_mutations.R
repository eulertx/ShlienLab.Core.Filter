get_cigar_mutations <- function(bam_file=NULL, mut_file=NULL, chr=NULL, position=NULL, operation="S", offsets=10) {

  if(is.null(bam_file)) stop("Mandatory argument bam_file is missing")
  if(is.null(mut_file)) stop("Mandatory argument mut_file is missing")
  if(is.null(chr)) stop("Mandatory argument chr is missing")

  #load mutations file if file directory or data frame
  if (is.data.frame(mut_file)) {
    muts <- mut_file
  } else {
    snv <- load(mut_file)
    muts <- get(snv)
  }

  # update offsets to be non zero only if soft clips are called (clipped at the front)
  if (operation == "S") {
    offsets = offsets
  } else {
    offsets = 0
  }

  # filter only needed chromosome if chr included, change position to useable numbers
  chr_mut <- muts %>% filter(annovar_chr == as.character(chr))
  chr_mut$annovar_start <- as.numeric(as.character(chr_mut$annovar_start))

  # load bam file.
  s_reads <- load_bam_reads(bam_file=bam_file, chr_mut=chr_mut, chr=chr, position=position, operation=operation, offsets=offsets)

  # get unique soft clipping positions from reads
  s_clips_pos <- NULL
  if (nrow(s_reads) > 0) {
    for (read in 1:nrow(s_reads)) {
      s_frame <- cigar_to_data_frame(cigar=s_reads$cigar[read], position=s_reads$pos[read])

      # if the operation is hard clipping, mark the entire cigar read as hard clipped
      if (operation=="H") {
        H_read_pos <- data.frame(start=s_frame[1,"start"], end=s_frame[nrow(s_frame), "end"])
        s_clips_pos <- rbind(s_clips_pos, H_read_pos)
      } else {
        if (operation=="S" & s_frame[1,"type"]=="S") {
          s_reads[read,"pos"] <- s_reads[read,"pos"] - as.numeric(as.character(s_frame[1,"size"]))
          s_frame <- cigar_to_data_frame(cigar=s_reads$cigar[read], position=s_reads$pos[read])
        }
        s_frame <- s_frame %>% filter(type==operation)
        s_clips_pos <- rbind(s_clips_pos, s_frame[,c('start', 'end')])
      }
    }
    s_clips_pos$start <- as.numeric(as.character(s_clips_pos$start))
    s_clips_pos$end <- as.numeric(as.character(s_clips_pos$end))
  }

  # select only the mutations in the chromosome which occur in soft clippings
  f <- function(pos) any(ifelse(s_clips_pos$start<=pos & pos<=s_clips_pos$end, TRUE, FALSE))
  chr_mut$in_soft_clips <- sapply(chr_mut$annovar_start, f)

  chr_mut$index <- 1:nrow(chr_mut)

  # go through each positive, check where the mutation occurs (from md_snv_position)
  # and compare against the given position, to see if they are at the same position
  s_muts <- chr_mut %>% filter(in_soft_clips == TRUE)

  if (nrow(s_muts) > 0) {
    mut_in_hard_read <- 0

    # go through each read to see if mutation occurs in correct location
    for (read in 1:nrow(s_reads)) {
      mutation_pos <- md_snv_position(md=as.character(s_reads[read, "md"]),
                                      position=as.numeric(as.character(s_reads[read, "pos"])),
                                      cigar=as.character(s_reads[read, "cigar"]),
                                      seq=as.character(s_reads[read, "seq"]),
                                      operation=operation)

      if (is.data.frame(mutation_pos)) { #if position in right place, keep true
        if (nrow(merge(data.frame("alt"=as.character(s_muts[1,"annovar_alt"]), "position"=position), mutation_pos,
                       by=c("alt", "position"))) > 0) {
          mut_in_hard_read <- mut_in_hard_read + 1
        }
      } else if (mutation_pos == 'indel') { #indel, keep as true, just in case
        mut_in_hard_read <- mut_in_hard_read + 1
      }
    }
    chr_mut$in_soft_clips[s_muts$index[1]] <- mut_in_hard_read
  }

  # append column to initial data frame
  cigar_type <- paste("in", get_cigar_type(operation), sep="_")

  muts$in_S <- ifelse(muts[,"annovar_chr"]==chr, chr_mut$in_soft_clips, 0)

  # rename in_S column
  names(muts) <- c(head(names(muts), -1), eval(parse(text="cigar_type")))
  return(muts)
}
