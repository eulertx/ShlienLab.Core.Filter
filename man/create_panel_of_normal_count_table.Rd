\name{create_panel_of_normal_count_table}
\alias{create_panel_of_normal_count_table}
\title{
Create a panel of normal count table.
}
\description{
A function that imports a list of mutations in normal samples and aggregates into a single table of normal counts based on SNV ID.
}
\usage{
create_panel_of_normal_count_table(path = NULL, recurrence = 2, pattern = "annotated_filtered.rda$")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{Full path to the directory containing panel of normal R dataframes}
  \item{recurrence}{a recurrent threshold for reporting (optional, default: 2)}
  \item{pattern}{file pattern (optional, default: annotated_filtered.rda$)}
}
\value{Returns a dataframe containing panel of normal SNV IDs and counts. }
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\keyword{ snv }
\keyword{ pon }
\keyword{ filter }
