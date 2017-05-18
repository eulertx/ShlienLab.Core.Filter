\name{get_cigar_mutations}
\alias{get_cigar_mutations}
\title{get_cigar_mutation}
\description{
Calculates which mutations of a given RDA file and chromosome are found within a region classified as a specific CIGAR operation, given a BAM file of aligned reads and a character to determine the classification of alignment of CIGAR string.
}
\usage{
get_cigar_mutations(bam_file = NULL, mut_file = NULL, chr = NULL, position = NULL, operation = "S", offsets = 10)
}
\arguments{
  \item{bam_file}{
Path directory to where an aligned BAM file is, containing aligned reads in the correct format (required).
  }
  \item{mut_file}{
Path directory to where a RDA file is *or* previously loaded mutation data frame containing mutations in the correct format (required).
  }
  \item{chr}{
String or character of the chromosome to be focused on, (e.g. "1", 1, "X") (required).
  }
  \item{position}{
If included, function will only look at the reads across this position to determine whether there are any mutations in the given CIGAR classification. If not included, the function will look across all positions of the chromosome where mutations occur, to look for the given CIGAR classification of reads in all those positions.
  }
  \item{operation}{
Character of either one of: M, I, D, N, S, H, P, =, or X, representing the available operations within a CIGAR string. This parameter will determine what the circumstance is a given mutation to be found within (default "S" for soft clippings). If operation is H for hard clipping, the entire read will be considered a hard clipping, not just the subset of the CIGAR string defined as "H", and will adjust the position to take the entire read into account, not just the hard clipped portion.
  }
  \item{offsets}{
Positive integer dictating interval to the right of the position in which to access all reads from the BAM file. This will be used only if the operation is "S" for soft clipped. This allows the bam file to access any read which is soft clipped on the left side where the position in the BAM file is offset to ignore the soft clipped region (default 10).
  }
}
\value{
Returns a data frame similar to the contents of mut_file, but with an added column (value of FALSE or a natural number) stating whether each mutation is within the specified CIGAR operation as well as in the chromosome and position specified, and if so, how many reads of that operation the mutation occurs in. If the position is not included, then will calculate number of reads within the CIGAR operation for each mutation within the RDA file which occur in the correct chromosome.
}
\author{L. Christine Schreiner}

\keyword{ mutations }
\keyword{ cigar }
\keyword{ CIGAR }
\keyword{ bam }
\keyword{ soft }
\keyword{ hard }
\keyword{ clip }
