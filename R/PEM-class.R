## **************************************************************************
##
##    (c) 2010-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Class and method for Phylogenetic Eigenvector maps**
##
##    This file is part of MPSEM
##
##    MPSEM is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    MPSEM is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with MPSEM. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Class and Methods for Phylogenetic Eigenvector Maps (PEM)
#' 
#' @description Class and methods to handle Phylogenetic Eigenvector Maps (PEM).
#' 
#' @docType class
#' 
#' @name PEM-class
#'  
#' @param x A \code{\link{PEM-class}} object containing a Phylogenetic
#' Eigenvector Map.
#' @param row.names Included for method consistency reason; ignored.
#' @param optional Included for method consistency reason; ignored.
#' @param object A \code{\link{PEM-class}} object.
#' @param targets Output of \code{\link{getGraphLocations}}.
#' @param lmobject An object of class \sQuote{lm} (see \code{\link{lm}} for
#' details).
#' @param newdata Auxiliary trait values.
#' @param interval The kind of limits (confidence or prediction) to return with
#' the predictions; \code{interval="none"}: do not return a confidence interval.
#' @param level Probability of the confidence of prediction interval.
#' @param ... Additional parameters to be passed to the method. Currently
#' ignored.
#' 
#' @details The \code{\link{print}} method provides the number of eigenvectors,
#' the number of observations these vectors are spanning, and their associated
#' eigenvalues.
#' 
#' The \code{\link{as.data.frame}} method extracts the eigenvectors from the
#' object and allows one to use \code{\link{PEM-class}} objects as \code{data}
#' parameter in function such as \code{\link{lm}} and \code{\link{glm}}.
#' 
#' The \code{\link{predict}} object is a barebone interface meant to make
#' predictions. It must be given species locations with respect to the
#' phylogenetic graph (\code{target}), which are provided by function
#' \code{\link{getGraphLocations}} and a linear model in the form of an object
#' from \code{\link{lm}}. The user must provide auxiliary trait values if
#' \code{lmobject} involves such trait.
#' 
#' @format A \code{\link{PEM-class}} object contains:
#' \describe{
#'   \item{ x }{ The \code{\link{graph-class}} object that was used to
#'   build the PEM (see \code{\link{PEM.build}}). }
#'   \item{ sp }{ A \code{\link{logical}} vector specifying which vertex is
#'   a tip. }
#'   \item{ B }{ The influence matrix for those vertices that are tips. }
#'   \item{ ne }{ The number of edges. }
#'   \item{ nsp }{ The number of species (tips). }
#'   \item{ Bc }{ The column-centred influence matrix. }
#'   \item{ means }{ The column means of \code{B}. }
#'   \item{ dist }{ Edge lengths. }
#'   \item{ a }{ The steepness parameter (see \code{\link{PEM.build}} for
#'   details). }
#'   \item{ psi }{ The relative evolution rate along the edges (see
#'   \code{\link{PEM.build}} for details). }
#'   \item{ w }{ Edge weights. }
#'   \item{ BcW }{ The weighted and column-centred influence matrix. }
#'   \item{ d }{ The singular values of \code{BcW}. }
#'   \item{ u }{ The eigenvectors (left singular vectors) of \code{BcW}. }
#'   \item{ vt }{ The right singular vectors of \code{BcW}. }
#' }
#' In addition to these standard component, function,
#' \code{\link{PEM.fitSimple}} and \code{\link{PEM.forcedSimple}} add the
#' following members, which are necessary to make predictions:
#' \describe{
#'   \item{ S2 }{ The variances of responses (one value for each response).
#'   }
#'   \item{ y }{ A copy of the responses. }
#'   \item{ opt }{ The list returned by \code{\link{optim}}. }
#' }
#' The estimated weighting parameters are also given as an edge property.
#' 
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps (PEM): a framework to model and predict species traits. Meth. Ecol.
#' Evol. 4: 1120--1131
#' 
#' @seealso \code{\link{PEM.build}}, \code{\link{PEM-class}}
#' 
#' @importFrom stats qt
#' 
NULL
#' 
#' @describeIn PEM-class
#' 
#' Print method for PEM-class objects
#' 
#' @export
print.PEM <- function(x, ...) {
  cat("A phylogenetic eigenvector map (PEM) for ",x$nsp," species:\n")
  if(x$nsp >= 10L)
    cat(paste(rownames(x$u)[1L:8L],collapse=","),"...",
        rownames(x$u)[nrow(x$u)],"\n")
  else
    cat(paste(rownames(x$u),collapse=","),"\n")
  cat("obtained from the following phylogenetic graph:\n")
  print(x$x)
  return(invisible(NULL))
}
#' 
#' @describeIn PEM-class
#' 
#' Method \code{as.data.frame} for PEM-class objects
#' 
#' @export
as.data.frame.PEM <- function(x, row.names = NULL, optional = FALSE, ...) {
  return(as.data.frame(x$u))
}
#' 
#' @describeIn PEM-class
#' 
#' Predict method for PEM-class objects
#' 
#' @export
predict.PEM <- function (object, targets, lmobject, newdata,
                         interval = c("none", "confidence", "prediction"),
                         level = 0.95, ...) {
  if(missing(newdata)) newdata <- Locations2PEMscores(object, targets) else {
    if(nrow(targets$locations)!=nrow(newdata))
      stop("'newdata' has ",nrow(newdata),
           " rows but the number of target species is ",
           nrow(targets$locations),".")
    tmp <- Locations2PEMscores(object, targets)
    rownames(newdata) <- rownames(tmp$scores)
    tmp$scores <- cbind(newdata,tmp$scores)
    newdata <- tmp
    rm(tmp)
  }
  interval <- match.arg(interval)
  Residual.variance <- diag(t(lmobject$residuals) %*%
                              lmobject$residuals)/lmobject$df
  Xh <- cbind(1, as.matrix(newdata$scores[, attr(lmobject$terms, "term.labels"),
                                          drop = FALSE]))
  pred <- Xh %*% lmobject$coefficients
  if (interval == "none") return(pred)
  R <- qr.R(lmobject$qr)
  invXtX <- solve(t(R) %*% R)
  XhinvXtXtXh <- diag(Xh %*% invXtX %*% t(Xh))
  if (interval == "confidence")
    S <- sqrt(t((newdata$VarianceFactor/nrow(object$y)) + 
                  matrix(Residual.variance, length(Residual.variance), 
                         length(XhinvXtXtXh)) * XhinvXtXtXh))
  if (interval == "prediction") 
    S <- sqrt(t(newdata$VarianceFactor + matrix(Residual.variance, 
                                                length(Residual.variance),
                                                length(XhinvXtXtXh)) * 
                  (1 + XhinvXtXtXh)))
  return(list(values = pred,
              lower = pred + S * qt(0.5 * (1 - level), lmobject$df),
              upper = pred + S * qt(0.5 * (1 - level), lmobject$df,
                                    lower.tail = FALSE)))
}
##
