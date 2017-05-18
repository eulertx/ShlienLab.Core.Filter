\name{get.recurrent.snvids}
\alias{get.recurrent.snvids}
\title{
Get recurrent SNV IDs from dataframe.
}
\description{
Returns a dataframe summarizing recurrent SNV IDs.
}
\usage{
get.recurrent.snvids(data = NULL, id = "snvid")
}
\arguments{
  \item{data}{Dataframe containing standard Shlien Lab somatic mutation dataframe with SNV ID column appended (required)}
  \item{id}{column header representing the SNV ID (optional, default: 'snvid')}}
\value{
Returns a data table containing SNVID and count.
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\keyword{ snv }
\keyword{ recurrent }
\keyword{ snvid }
