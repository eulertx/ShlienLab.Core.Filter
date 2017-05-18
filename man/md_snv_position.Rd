\name{md_snv_position}
\alias{md_snv_position}
\title{md_snv_position}
\description{
Calculates all the mutated bases and the positions thereof from the information given about a single read. Produces either a data frame of those mutations, or a string detailing why no mutations were found in that specific read.
}
\usage{
md_snv_position(md = NULL, position = NULL, cigar = NULL, seq = NULL, operation = "S")
}
\arguments{
  \item{md}{
String representing an MD tag of a specific read extracted from a BAM file.
  }
  \item{position}{
Number representing the start position of the read (this is adjusted automatically if the read is hard clipped)
  }
  \item{cigar}{
String representing the cigar string of the read, in the correct format.
  }
  \item{seq}{
String representing the sequence of bases in the read (TCAG).
  }
  \item{operation}{
Character representing a CIGAR string operation, either M, I, D, N, S, H, P, =, or X.
  }
}
\value{
Returns either a data frame or a string. If the read does not contain any mutations, returns the string "unmutated". If the read contains either an insertion or deletion in comparison to the reference genome, then this function returns the string "indel". If the read is mutated and does not contain any length differences, then this function returns a data frame with the columns *alt* and *position* which represent the alternative base, and at which position the SNV occurs.
}
\author{L. Christine Schreiner}

\keyword{ md }
\keyword{ snv }
\keyword{ hard }
\keyword{ clip }
