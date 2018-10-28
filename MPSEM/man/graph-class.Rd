\encoding{utf8}
\name{graph-class}
\alias{ graph-class }
\alias{ print.graph }
\docType{class}
\title{ Graph class and methods }
\description{
  Class and methods to handle graphs.
}
\usage{
  \method{print}{graph}(x, ...)
}
\arguments{
  \item{ x }{ An object of \code{\link{graph-class}} }
  \item{ ... }{ Additional parameters to be passed to the
    method. Currently ignored. }
}
\details{
  Prints user-relevant information about the graph: number of edges and
  vertices, edge and vertex labels, addition edge properties and vertex
  properties.
}
\value{
  A \code{\link{graph-class}} object contains:
  \item{ edge }{ a list whose first two unnamed members are the indices
    of the origin and destination vertices. Additional members must be
    named and are additional edge properties (e.g. length), and }
  \item{ vertex }{ a list that optionally contain vertex properties, if
    any (or an empty list if none). }
}
\references{
  Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic
  eigenvector maps (PEM): a framework to model and predict species
  traits. Meth. Ecol. Evol. In press.
}
\author{ Guillaume Guénard, Département de sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
\seealso{
  \code{\link{PEM.build}}
  \code{\link{PEM-class}}
}
