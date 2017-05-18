\name{cigar_to_data_frame}
\alias{cigar_to_data_frame}
\title{
  cigar_to_data_frame
}
\description{
  Creates a data frame from a CIGAR string which differentiates all the different CIGAR operations, and positions them all based on an inputted starting position.
}
\usage{
cigar_to_data_frame(cigar = NULL, position = NULL)
}
\arguments{
  \item{cigar}{
    String representing a CIGAR string, in the correct format (required).
  }
  \item{position}{
    Number representing the beginning position of the cigar string along a chromosome (required).
  }
}
\value{
  Returns a data frame representing the dissected CIGAR string. Has 4 columns: type, size, start, end. Type represents the type of CIGAR operation (one of M, I, D, N, S, H, P, =, or X). Size represents how long the type goes on for, from the start position to the end position (represented by start and end respectfully).
}
\author{
  L. Christine Schreiner
}
\keyword{ cigar }
\keyword{ CIGAR }
\keyword{ position }
