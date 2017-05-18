\name{flag.complex.overlap.loci}
\alias{flag.complex.overlap.loci}
\title{
Flag overlap of loci to complex regions of the genome.
}
\description{
A method that identifies loci and determines whether they are found within/near complex genomic regions.
}
\usage{
flag.complex.overlap.loci(data = NULL, minimum = 0, far.dist = 999999)
}
\arguments{
  \item{data}{A dataframe containing the cFilter output}
  \item{minimum}{maximum distance from a complex region, default: 0}
  \item{far.dist}{an integer used to represent a far distance to replace NA values, default: 999999}
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\keyword{ snv }
\keyword{ filter }
