\name{add.snvid}
\alias{add.snvid}
\title{
Add a SNVID column.
}
\description{A method for adding the SNVID to the mutation dataframe using chromosome, start, end, ref, alt.}
\usage{
add.snvid(data = NULL)
}
\arguments{
  \item{data}{A dataframe from the standard Shlien Lab SNV pipeline.}
}
\value{
Returns a dataframe containing appended SNV ID.
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\examples{
x = data.frame(annovar_chr=c(1, 2, 3), annovar_start=c(1000, 2000, 3000), annovar_end=c(1000, 2000, 3000), annovar_ref=c('A', 'G', 'C'), annovar_alt=c('C', 'C', 'T'))
add.snvid(data=x)
}
\keyword{ snv }
\keyword{ snvid }
