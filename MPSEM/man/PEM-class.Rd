\encoding{utf8}
\name{PEM-class}
\alias{ PEM-class }
\alias{ print.PEM }
\alias{ as.data.frame.PEM }
\alias{ predict.PEM }
\docType{class}
\title{ Class and methods for Phylogenetic Eigenvector Maps (PEM) }
\description{
  Class and methods to handle Phylogenetic Eigenvector Maps (PEM).
}
\usage{
  \method{print}{PEM}(x, ...)
  \method{as.data.frame}{PEM}(x, row.names = NULL, optional = FALSE, ...)
  \method{predict}{PEM}(object, targets, lmobject, newdata, 
  interval = c("none", "confidence", "prediction"), level = 0.95, ...)
}
\arguments{
  \item{ x }{ A \code{\link{PEM-class}} object containing a Phylogenetic
    Eigenvector Map. }
  \item{ row.names }{ Included for method consistency reason; ignored. }
  \item{ optional }{ Included for method consistency reason; ignored. }
  \item{ object }{ A \code{\link{PEM-class}} object. }
  \item{ targets }{ Output of \code{\link{getGraphLocations}}. }
  \item{ lmobject }{ An object of class \sQuote{lm} (see
    \code{\link{lm}} for details). }
  \item{ newdata }{ auxiliary trait values }
  \item{ interval }{ The kind of limits (confidence or prediction) to
    return with the predictions. \code{interval="none"}: do not return a
    confidence interval. }
  \item{ level }{ Probability of the confidence of prediction interval. }
  \item{ ... }{ Further parameters to be passed to other functions or
    methods (currently ignored). }
}
\details{
  The \code{\link{print}} method provides the number of eigenvectors,
  the number of observations these vectors are spanning, and their
  associated eigenvalues.

  The \code{\link{as.data.frame}} method extracts the eigenvectors from
  the object and allows one to use \code{\link{PEM-class}} objects as
  \code{data} parameter in function such as \code{\link{lm}} and
  \code{\link{glm}}.

  The \code{\link{predict}} object is a barebone interface meant to make
  predictions. It must be given species locations with respect to the
  phylogenetic graph (\code{target}), which are provided by function
  \code{\link{getGraphLocations}} and a linear model in the form of an
  object from \code{\link{lm}}. The user must provide auxiliary trait
  values if \code{lmobject} involves such trait.
}
\value{
  A \code{\link{PEM-class}} object contains:
  \item{ x }{ the \code{\link{graph-class}} object that was used to
    build the PEM (see \code{\link{PEM.build}}), }
  \item{ sp }{ a \code{\link{logical}} vector specifying which vertex is
    a tip, }
  \item{ B }{ the influence matrix for those vertices that are tips, }
  \item{ ne }{ the number of edges, }
  \item{ nsp }{ the number of tips, }
  \item{ Bc }{ the column-centred influence matrix, }
  \item{ means }{ the column means of \code{B} }
  \item{ dist }{ edge lengths, }
  \item{ a }{ the steepness parameter (see \code{\link{PEM.build}} for
    details), }
  \item{ psi }{ the relative evolution rate along the edges (see
    \code{\link{PEM.build}} for details),}
  \item{ w }{ edge weights, }
  \item{ BcW }{ the weighted and column-centred influence matrix, }
  \item{ d }{ the singular values of \code{BcW}, }
  \item{ u }{ the eigenvectors (left singular vectors) of \code{BcW},
    and }
  \item{ vt }{ the right singular vectors of \code{BcW}.}
  In addition to these standard component, function,
  \code{\link{PEM.fitSimple}} and \code{\link{PEM.forcedSimple}} add the
  following members, which are necessary to make predictions:
  \item{ S2 }{ the variance(s) of the response(s), }
  \item{ y }{ a copy of the response(s), and}
  \item{ opt }{ the list returned by \code{\link{optim}}, }
  as well as a copy of the estimated weighting parameters as edge
  properties.
}
\references{
  Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic
  eigenvector maps (PEM): a framework to model and predict species
  traits. Meth. Ecol. Evol. In press.
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
