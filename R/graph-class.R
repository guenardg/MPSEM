## **************************************************************************
##
##    (c) 2010-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Class and method for directed graphs**
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
#' Class and Method for Directed Graphs
#' 
#' @description Class and methods to handle graphs.
#' 
#' @docType class
#' 
#' @name graph-class
#' 
#' @param x An object of \code{\link{graph-class}}.
#' @param ... Additional parameters to be passed to the method. Currently
#' ignored.
#' 
#' @details Prints user-relevant information about the graph: number of edges
#' and vertices, edge and vertex labels, addition edge properties and vertex
#' properties.
#' 
#' @format A \code{graph-class} object contains:
#' \describe{
#'   \item{ edge }{ A list whose first two unnamed members are the indices of
#'   the origin and destination vertices. Additional members must be named and
#'   are additional edge properties (e.g. length). }
#'   \item{ vertex }{ A list that optionally contain vertex properties, if any
#'   (or an empty list if none). }
#' }
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
NULL
#' 
#' @describeIn graph-class
#' 
#' Print method for graph-class objects
#' 
#' @export
print.graph <- function(x, ...) {
  cat("\nA graph with",attr(x,"ev")[1],"edges and",attr(x,"ev")[2],"vertices.",
      "\n")
  if(!is.null(attr(x,"elabel")))
    cat("Edge labels:",paste(attr(x,"elabel")),"\n")
  if(!is.null(attr(x,"vlabel")))
    cat("Vertex labels:",paste(attr(x,"vlabel")),"\n")
  if(length(attr(x$edge,"names")>2)) {
    cat("Available edge information: ",
        paste(attr(x$edge,"names")[-(1:2)],collapse=", "),"\n")
  } else {
    cat("No available edge information\n")
  }
  if(length(attr(x$vertex,"names")>0)) {
    cat("Available vertex information: ",
        paste(attr(x$vertex,"names"),collapse=", "),"\n")
  } else {
    cat("No available vertex information\n")
  }
  cat("\n")
}
##
