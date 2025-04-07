## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Harmonic mean calculation **
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
#' @name harmonicMean
#' 
#' @title Harmonic Mean
#' 
#' @description Functions to calculate the harmonic mean and weighted harmonic
#' mean of a set of data.
#' 
#' @param x A numeric vector containing a set of observations.
#' @param w A numeric vector containing the observation weights.
#' 
#' @return A numeric value.
#' 
#' @details The (weighted) harmonic mean of a data set is the inverse of the
#' (weighted) mean of the inverse values. As such, any 0 value in the data (or
#' in the non-zero weighted data in the weighted case) will result in the
#' functions returning a value of 0. It is worth recall that whereas the value
#' of the arithmetic mean is strongly influenced by even a few extreme values,
#' the value of the harmonic mean is strongly influenced by even a few small
#' values. The function take any numeric values. However, the harmonic mean of
#' data containing both positive and negative values is unlikely to be stable.
#' 
#' @author \packageAuthor{MPSEM}
#' 
#' 
NULL
#' 
#' @describeIn harmonicMean
#' 
#' Harmonic Mean.
#' 
#' Calculates the (unweighted) harmonic mean.
#' 
#' @export
hmean <- function(x) {
  
  n <- length(x)
  storage.mode(x) <- "double"
  
  .C("hmeanC", n, x, double(1L), PACKAGE="MPSEM")[[3L]]
}
#' 
#' @describeIn harmonicMean
#' 
#' Weighted Harmonic Mean
#' 
#' Calculates the weighted harmonic mean.
#' 
#' @export
weighted.hmean <- function(x, w) {
  
  n <- length(x)
  storage.mode(x) <- "double"
  w <- rep(w, length.out = n)
  storage.mode(w) <- "double"
  
  .C("whmeanC", n, x, w, double(1L), PACKAGE="MPSEM")[[4L]]
}
#' 
