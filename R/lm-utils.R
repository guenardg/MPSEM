## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Linear modelling utility functions **
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
#' Linear Modelling Utility Functions
#' 
#' @description Utility functions to build linear models using Phylogenetic
#' Eigenvector Maps among their explanatory variables.
#' 
#' @name lm-utils
#' 
#' @aliases utils
#' 
#' @param formula an object of class "\code{\link{formula}}" (or one that can be
#' coerced to that class): a symbolic description of the model data to be
#' prepared. See ‘Details’ in \code{\link{lm}} for further information about
#' model specification.
#' @param data An optional data frame, list or environment (or object coercible
#' by \code{\link{as.data.frame}} to a data frame) containing the variables. If
#' not found in data, the variables are taken from environment(formula),
#' typically the environment from which the function is called.
#' @param ... Additional parameters to be passed to the method.
#' @param na.action A function (of the name of a function) for treating missing
#' values (\code{NA}s) (default: \code{na.pass}).
#' @param contrasts An optional list. See the \code{contrasts.arg} of
#' \code{\link{model.matrix.default}}. (default: \code{NULL}).
#' @param discard.intercept A \code{\link{logical}}; whether of not to discard
#' the intercept from the model matrix (default: \code{TRUE}).
#' @param obs A numeric vector of observations.
#' @param prd A numeric vector of model predictions.
#' 
#' @details Function \code{model.data} is useful to prepare data to be given as
#' response and auxiliary trait(s) to other functions such as
#' \code{\link{evolution.model.PEM}}. In general, the implicit constant term
#' (intercept) is not useful and can be explicitly discarded.
#' 
#' @return
#' \describe{
#'   \item{model.data}{A three-member list with member \code{$y} (a vector or
#'     matrix of response traits), member \code{$x} (a matrix auxiliary traits
#'     coded as numeric values), and member \code{$terms} (A model
#'     description).}
#'   \item{Psquare}{A numeric value.}
#' }
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution 4: 1120-1131
#' 
#' @importFrom stats model.response model.matrix model.frame na.pass
#' 
NULL
#' 
#' @describeIn lm-utils
#' 
#' Model Data Preparation
#' 
#' Transforms data from various types into a list containing the response
#' trait(s) and a strictly numeric auxiliary traits data matrix.
#' 
#' @export
model.data <- function (formula, data, ..., na.action = na.pass,
                        contrasts = NULL, discard.intercept = TRUE) {
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  x <- model.matrix(mt, mf, contrasts)
  y <- model.response(mf, "numeric")
  
  if(discard.intercept) {
    attr(mt, "intercept") <- 0L
    i <- which(attr(x,"assign") == 0L)
    if(length(i)) {
      asn <- attr(x,"assign")
      ctr <- attr(x,"contrasts")
      x <- x[,-i,drop=FALSE]
      asn <- asn[-i]
      attr(x,"assign") <- asn
      attr(x,"contrasts") <- ctr
    }
  }
  
  if(ncol(x)) {
    list(y=y, x=x, terms=mt)
  } else {
    list(y=y, x=NULL, terms=NULL)
  }
}
#' 
#' @describeIn lm-utils
#' 
#' Coefficient of Prediction
#' 
#' Calculates the prediction coefficient between observations and model
#' predictions.
#' 
#' @export
Psquare <- function (obs, prd)
  1 - sum((obs - prd)^2)/sum((obs - mean(obs))^2)
#' 
