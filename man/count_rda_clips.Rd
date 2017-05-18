\name{count_rda_clips}
\alias{count_rda_clips}
\title{count_rda_clips}
\description{
For given BAM and RDA files, calls get_cigar_mutations on each row (mutation) of the file, checking if it occurs within whichever CIGAR string operation specified. Avoids get_cigar_mutations from being recursive, and focuses functionality on specific mutations and not entire RDA files.}
\usage{
count_rda_clips(bam = NULL, rda = NULL, operation = "S", offsets = 10)
}
\arguments{
  \item{bam}{
Path directory to BAM file for sample (required).
  }
  \item{rda}{
Path directory to where a RDA file is *or* previously loaded mutation data frame containing mutations in the correct format (required).
  }
  \item{operation}{
Character of either one of: M, I, D, N, S, H, P, =, or X, representing the available operations within a CIGAR string. This parameter will determine what the circumstance is a given mutation to be found within (default "S" for soft clippings). If operation is H for hard clipping, the entire read will be considered a hard clipping, not just the subset of the CIGAR string defined as "H", and will adjust the position to take the entire read into account, not just the hard clipped portion.
  }
  \item{offsets}{
Positive integer dictating interval to the right of the position in which to access all reads from the BAM file. This will be used only if the operation is "S" for soft clipped. This allows the bam file to access any read which is soft clipped on the left side where the position in the BAM file is offset to ignore the soft clipped region (default 10).
  }
}
\value{
Returns a data frame similar to the contents of mut_file, but with an added column (value of FALSE or a natural number) stating whether each mutation is within the specified CIGAR operation, and if so how many reads of that operation the mutation occurs in. The chromosome and position are extracted from each row, and used to calculate whether or not that specific row is within the operation provided.
}
\author{L. Christine Schreiner}

\keyword{ mutations }
\keyword{ cigar }
\keyword{ CIGAR }
\keyword{ bam }
\keyword{ soft }
\keyword{ hard }
\keyword{ clip }
