## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Inverse Distance Weighting calculation **
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
#' @name invDistWeighting
#' 
#' @title Inverse Distance Weighting
#' 
#' @description A function to calculate inverse distance weights associated with
#' a vector of weights.
#' 
#' @param x A numeric vector containing a set of distances.
#' @param a A numeric value for the distance exponent (default \code{0}).
#' @param ... Any further arguments to be ignored by the function.
#' 
#' @returns A numeric vector containing a set of weights summing to 1.
#' 
#' @details Particular cases are when \code{a = 0} (the default), whereby all
#' the weights have equal values (1/n, where n is \code{length(x)}) and when
#' any \code{d[i] == 0}, whereby the weights are non-zero only for the
#' zero-valued distance(s). The weights are equal to \code{w[i] = d[i]^(-a)/S},
#' where S is the sum of \code{d[i]^(-a)} for all the \code{i}.
#' 
#' @return A numeric value.
#' 
#' @author \packageAuthor{MPSEM}
#' 
#' @examples
#' 
#' ## Example 1, equal weights (the function's default when a = 0):
#' d0 <- c(10,1,0.1,0.2,2,20,0.3)
#' w0 <- invDistWeighting(d0)
#' w0
#' sum(w0)
#' 
#' ## Example 2, inverse distance (a = 1):
#' w1 <- invDistWeighting(d0, a=1)
#' w1
#' sum(w1)
#' 
#' ## Example 3, inverse squared distance (a = 2):
#' w2 <- invDistWeighting(d0, a=2)
#' round(w2, 5)
#' sum(w2)
#' 
#' ## A case some of the distances are 0:
#' d1 <- c(10,0,0.1,0.2,0,20,0.3)
#' w3 <- invDistWeighting(d1, a=1)  ## or any a != 0
#' w3
#' sum(w3)
#' 
#' 
#' @export
invDistWeighting <- function(x, a = 0, ...) {
  
  n <- length(x)
  storage.mode(x) <- "double"
  storage.mode(a) <- "double"
  
  .C("invDistWeightingC", n, a, x, double(n), PACKAGE="MPSEM")[[4L]]
}
