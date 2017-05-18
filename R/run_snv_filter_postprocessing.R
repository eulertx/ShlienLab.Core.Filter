run_snv_filter_postprocessing <- function(data=NULL, cfilter=NULL, min.dist.complex=0, max_number_of_flagged=2, pon=NULL, pon_max=2, hardfilter=FALSE) {
  if (is.null(data)) stop("Mandatory argument data is missing")
  if (is.null(cfilter)) stop("Mandatory argument cfilter is missing")
  if (is.null(pon)) stop("Mandatory argument pon is missing")

  data <- add.snvid(data=data)
  data.cfilter <- read.table(
    file=cfilter,
    header=TRUE,
    as.is=TRUE,
    sep='\t',
    quote="\""
    )
  data.cfilter <- data.cfilter %>% dplyr::mutate(
    snvid=paste(
      sep='_',
      Chromosome1,
      Position1,
      Position1,
      Ref.Allele,
      Alt.Allele
      )
    )
  data.cfilter <- data.cfilter %>% dplyr::select(
    snvid,
    in_centromere,
    normal_coverage_threshold,
    unique_mapping,
    high_depth,
    distance_to_low_complexity_1,
    distance_to_low_complexity_2,
    multi_mapping
    )

  # change the distance_to_low_complexity_1 values to "FLAGGED" or 0
  data.cfilter <- flag.complex.overlap.loci(
    data=data.cfilter,
    minimum=min.dist.complex
    )
  data.cfilter$flagged_count <- apply(
    X=data.cfilter[,c('in_centromere', 'normal_coverage_threshold', 'unique_mapping', 'high_depth', 'multi_mapping', 'dist_low_complexity')],
    MARGIN=1,
    FUN=add.flagged.count.column
    )

  data <- dplyr::left_join(x=data, y=data.cfilter, by=c('snvid'))

  # threshold for flagged items
  if (hardfilter == TRUE) {
    data <- data %>% filter(flagged_count < max_number_of_flagged)
    }
  # load the panel of normal count table to annotate the somatic SNV data
  load(pon)
  data <- left_join(x=data, y=data.frame(pon_count), by=c('snvid'))

  # NAs are introduced if a snvid has no entry in the pon_count vector
  data[is.na(data$pon_count),]$pon_count <- 0

  # threshold for recurrence count
  if (hardfilter == TRUE) {
    data <- data %>% filter(count < pon_max)
    }

  # cleanup the dataframe, remove extraneous columns
  data$distance_to_low_complexity_2 <- NULL

  return(data)
  }
