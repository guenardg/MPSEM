## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Phylogenetic Eigenvector Maps (PEM) Helper Functions **
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
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with MPSEM. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Phylogenetic Eigenvector Maps
#' 
#' @description Functions to calculate and manipulate Phylogenetic Eigenvector
#' Maps (PEM), which are sets of eigenfunctions describing the structure of a 
#' phylogenetic graph. Each computation function is briefly described in 
#' section \code{Functions} below.
#' 
#' @name PEM-functions
#' 
#' @param d A numeric vector of the evolutionary distances (\code{PEMweights})
#' or a character string specifying the row name of the edge table (located in
#' the \code{\link{graph-class}} object's attribute \code{edge}) where the
#' evolutionary distances (edge lengths) can be found.
#' @param a The steepness parameter describing whether changes occur, on
#' average: progressively long edges (a close to 0) or abruptly at vertices (a
#' close to 1; default: \code{0}).
#' @param psi Relative evolution rate along the edges. This parameter only
#' becomes relevant when multiple values are assigned to different portions of
#' the phylogeny (default: \code{1}).
#' 
#' @return The returned value depends on the function:
#' \describe{
#'   \item{PEMweights}{A set of numeric value to be used as weights during PEM
#'     calculation.}
#' }
#' 
#' @author \packageAuthor{MPSEM} --
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution. 4: 1120--1131
#' 
#' Makarenkov, V., Legendre, L. & Desdevise, Y. 2004. Modelling phylogenetic
#' relationships using reticulated networks. Zoologica Scripta 33: 89--96
#' 
#' Blanchet, F. G., Legendre, P. & Borcard, D. 2008. Modelling directional
#' spatial processes in ecological data. Ecological Modelling 215: 325--336
#' 
#' @seealso \code{\link{PEM-class}}
#' 
## @examples
#' 
NULL
#' 
#' @describeIn PEM-functions
#' 
#' PEM Weighting
#' 
#' A power function to obtain the edge weights used during PEM calculation.
#' 
#' @export
PEMweights <- function (d, a = 0, psi = 1) {
  
  nd <- length(d)
  a <- rep(a, length.out = nd)
  psi <- rep(psi, length.out = nd)
  
  .C(
    "PEMweightC",
    as.double(d),
    as.integer(nd),
    as.double(a),
    as.double(psi),
    res = double(nd),
    PACKAGE = "MPSEM"
  )$res
}
#' 
