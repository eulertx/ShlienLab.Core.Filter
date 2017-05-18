\name{add.flagged.count.column}
\alias{add.flagged.count.column}
\title{
Count number of "FLAGGED" in a vector.
}
\description{
Inputs a vector containing strings 0 or FLAGGED.  Returns the count of "FLAGGED" found in the vector.
}
\usage{
add.flagged.count.column(x)
}
\arguments{
  \item{x}{character vector containing 0 or FLAGGED}
}
\value{
Returns an integer representing the count of FLAGGED.
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\examples{
x <- c('FLAGGED', '0', 'FLAGGED', 'FLAGGED', '0', 'FLAGGED')
add.flagged.count.column(x=x)
}
\keyword{ snv }
\keyword{ snvid }
