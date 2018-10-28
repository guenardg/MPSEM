\encoding{utf8}
\name{lmutils}
\alias{ lmforwardsequentialsidak }
\alias{ lmforwardsequentialAICc }
\title{ Linear modelling utility functions }
\description{ Utility functions to build linear models using
  Phylogenetic Eigenvector Maps as their features.
}
\usage{
  lmforwardsequentialsidak(y, x, object, alpha=0.05)
  lmforwardsequentialAICc(y, x, object)
}
\arguments{
  \item{ y }{ a response variable }
  \item{ x }{ descriptors to be used as auxiliary traits }
  \item{ object }{ a \code{PEM} class object }
  \item{ alpha }{ the threshold above which to stop adding variables }
}
\value{
  Both functions return a \code{lm} class object.
}
\details{
  Function \code{\link{lmforwardsequentialsidak}}, performs a forward
  stepwise selection of the PEM eigenvectors until the familywise test
  of significance of the new variable to be included exceeds the
  threshold \code{alpha}. The familiwise type I error probability is
  obtained using the Holm-Sidak correction of the testwise
  probabilities, thereby correcting for type I error rate inflation due
  to multiple testing. \code{lmforwardsequentialAICc} carries out
  forward stepwise selection of the eigenvectors as long as the
  candidate model features a lower sample-size-corrected Akaike
  information criterion than the previous model. The final model should
  be regarded as overfit from the Neyman-Pearson (\emph{i.e.}
  frequentist) point of view, but it is the model that minimizes
  information loss from the standpoint of information theory.
}
\references{
  Burnham, K. P. & Anderson, D. R. 2002. Model selection and multimodel
  inference: a practical information-theoretic approach, 2nd
  ed. Springer-Verlag. xxvi + 488 pp.

  Holm, S. 1979. A simple sequentially rejective multiple test procedure.
  Scand. J. Statist. 6: 65-70.

  Sidak, Z. 1967. Rectangular confidence regions for means of multivariate
  normal distributions. J. Am. Stat. Ass. 62, 626-633.
}
\examples{
## No example has yet been produced.
}
\author{ Guillaume Guénard, Département de sciences biologiques
  Université de Montréal, Montréal, QC, Canada. }
