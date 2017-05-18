\name{load_bam_reads}
\alias{load_bam_reads}
\title{load_bam_reads}
\description{
Loads reads from a given BAM file which both overlap a given position on a chromosome and contain a part which is classified as the same type as the given CIGAR string operation. If the position is not given, then the function will loop through every mutation within the chr_mut data frame and extract the given position. The function will then load and return all reads which overlap one of those positions, and which also contain that given CIGAR string operation.
}
\usage{
load_bam_reads(bam_file = NULL, chr_mut = NULL, chr = NULL, position = NULL, operation = "S", offsets = NULL)
}
\arguments{
  \item{bam_file}{
Path directory to where an aligned BAM file is, containing aligned reads in the correct format (required).
  }
  \item{chr_mut}{
Mutation data frame containing mutations from only the chr chromosome in the correct format. Only required if the position argument is not included, in which case the mutations file will be used to find all positions at which a mutation occurs within the given chromosome, and loads all reads associated with those positions which contain the given CIGAR string.
  }
  \item{chr}{
String or character of the chromosome to be focused on, (e.g. "1", 1, "X") (required).
  }
  \item{position}{
If included, function will load the reads across this position with the given CIGAR classification. If not included, the function will look across all positions of the chromosome where mutations occur (which requires the passing of the chr_mut argument), and will load reads with the given CIGAR classification in all those positions.
  }
  \item{operation}{
Character of either one of: M, I, D, N, S, H, P, =, or X, representing the available operations within a CIGAR string. This parameter will determine what the circumstance is a given mutation to be found within. If operation is H for hard clipping, the entire read will be considered a hard clipping, not just the subset of the CIGAR string defined as "H", and will adjust the position to take the entire read into account, not just the hard clipped portion (default "S" for soft clippings).
  }
    \item{offsets}{
Positive integer dictating interval to the right of the position in which to access all reads from the BAM file (required).
  }
}
\value{
Returns a data frame containing all of the reads (unique position, mapping quality, CIGAR string, base sequence, and MD tag) from the given BAM file which overlap the chromosome and position stated. If the position is not stated, will return all BAM file reads in the given chromosome and in every position at which a mutation (from chr_mut) occurs.
}
\author{L. Christine Schreiner}

\keyword{ reads }
\keyword{ BAM }
\keyword{ CIGAR }
\keyword{ load }
