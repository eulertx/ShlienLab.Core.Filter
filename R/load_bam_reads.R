load_bam_reads <- function(bam_file = NULL, chr_mut = NULL, chr = NULL, position = NULL, operation = "S", offsets = NULL) {

  if(is.null(bam_file)) stop("Mandatory argument bam_file is missing")
  if(is.null(chr)) stop("Mandatory argument chr is missing")
  if(is.null(offsets)) stop("Mandatory argument offsets is missing")
  if(is.null(chr_mut) & is.null(position)) stop("Function requires one of either position or chr_mut argument")

  # If a position isn't specified, look through all chromosome mutation positions
  if (is.null(position)) {
    positions <- chr_mut$annovar_start
    whole_bam <- NULL

    for (pos in 1:length(positions)) {
      what <- c('pos', 'cigar', 'mapq', 'seq')
      which <- get_which(chr=as.character(chr), position=positions[pos], offsets=offsets)
      param <- ScanBamParam(what=what, which=which, tag="MD")
      bam <- as.data.frame(scanBam(bam_file, param = param))
      names(bam) <- c('pos', 'mapq', 'cigar', 'seq', 'md')
      whole_bam <- rbind(whole_bam, bam)
    }
  }

  # If a position is specified, do not go through for loop
  else {
    what <- c('pos', 'cigar', 'mapq', 'seq')
    which <- get_which(chr=as.character(chr), position=position, offsets=offsets)
    param <- ScanBamParam(what=what, which=which, tag="MD")
    whole_bam <- as.data.frame(scanBam(bam_file, param = param))
    names(whole_bam) <- c('pos', 'mapq', 'cigar', 'seq', 'md')
  }

  # filter only reads with soft clippings (or other operation) in CIGAR string
  s_cigars <- whole_bam %>% filter(grepl(operation, cigar))
  s_cigars$cigar <- sapply(s_cigars$cigar, as.character)
  return(s_cigars)
}
