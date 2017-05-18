\name{run_snv_filter_postprocessing}
\alias{run_snv_filter_postprocessing}
\title{
Run the SNV filter postprocessing function.
}
\description{
A function that wraps the SNV filter postprocessing pipeline.
}
\usage{
run_snv_filter_postprocessing(data = NULL, cfilter = NULL, min.dist.complex = 0, max_number_of_flagged = 2, pon = NULL, pon_max = 2, hardfilter = FALSE)
}
\arguments{
  \item{data}{somatic mutation dataframe (required)}
  \item{cfilter}{output file from the cfilter pipeline (required)}
  \item{min.dist.complex}{minimum distance to complex region to be considered overlap (optional, default: 0)}
  \item{max_number_of_flagged}{maximum number of "FLAGGED" tag in all filters (optional, default=2)}
  \item{pon}{full path to the panel of normal count table (required)}
  \item{pon_max}{minimum number of recurrence within panel of normals (optional, default=2)}
  \item{hardfilter}{flag to issue when removing filtered mutations instead of annotating}
}
\value{
Returns a dataframe with hard-filtered mutations and annotations.
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\keyword{ snv }
\keyword{ filter }
