\name{get_which}
\alias{get_which}
\title{
get_which
}
\description{
Function to calculate the 'which' parameter when loading a BAM file through Rsamtools. Given a chromosome and position, will calculate the appropriate access of the subset of data from the file, without loading unnecessary data.  }
\usage{
get_which(chr = NULL, position = NULL, offsets = NULL)
}
\arguments{
  \item{chr}{
String representing a chromosome (e.g. "X", "12") (required).
  }
  \item{position}{
Integer representing the position in question, The function will return all the reads which overlap this position.
  }
  \item{offsets}{
Positive integer dictating interval to the right of the position in which to access all reads from the BAM file (required).
  }
}

\value{
Produces a range that can be used to create a param object to open a subset of a BAM file in a separate function.
}
\author{
L. Christine Schreiner
}
\keyword{ RSamtools }
\keyword{ which }
\keyword{ param }
