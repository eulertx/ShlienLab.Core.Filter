\name{get_position_reads}
\alias{get_position_reads}
\title{get_position_reads}
\description{
Produces a list of all the reads (CIGAR string and position) from a BAM file which overlap a certain chromosome and position, and whose CIGAR string contains a classification of the specified operation.
}
\usage{
get_position_reads(bam_file = NULL, chr = NULL, position = NULL, operation = "S", offsets = NULL)
}
\arguments{
  \item{bam_file}{
Path directory to where an aligned BAM file is, containing aligned reads in the correct format (required).
  }
  \item{chr}{
String or character of the chromosome to be focused on, (e.g. "1", 1, "X") (required).
  }
  \item{position}{
Number representing the location within the chromosome. Reads which overlap this position on the given chromosome will be returned (required).
  }
  \item{operation}{
Character of either one of: M, I, D, N, S, H, P, =, or X, representing the available operations within a CIGAR string. This parameter will determine which overlapping reads will be returned, if they contain this operation. If "H", then the entire read will be classified as hard clipped, not just the portion from the cigar string (default "S" for soft clippings).
  }
  \item{offsets}{
Positive integer dictating interval to the right of the position in which to access all reads from the BAM file (required).
  }
}
\value{
Returns a data frame containing a unique subset of the reads which overlap the given position and chromosome> The subset of the location's reads is specified by whether the read contains the CIGAR string operation, and whether that operation surrounds the given position.
}
\author{
L. Christine Schreiner
}
\keyword{ read }
\keyword{ reads }
