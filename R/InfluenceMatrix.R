## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Influence matrix class Methods, and functions **
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
#' Influence Matrix
#' 
#' @description Functions and methods calculate and manipulate graph influence
#' matrix.
#' 
#' @docType class
#' 
#' @name InflMat-class
#' 
#' @aliases inflMat
#' 
#' @param x A \code{\link{graph-class}} or \code{InflMat-class} object.
#' @param value A vector or \code{\link{data.frame}} containing the values to be
#' given to the \code{InflMat-class} object.
#' @param ... Further arguments to be passed internally to other functions or
#' methods.
#' 
#' @return The returned value depends on the function:
#' \describe{
#' \item{InflMat}{A binary influence matrix of the graph with as many rows as
#' its number of vertices and as many columns as its number of edges.}
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
#' @seealso \code{\link{PEM-class}} \code{\link{PEM-functions}}
#' 
#' @examples  ## Exemplary graph:
#' data.frame(
#'   species = as.logical(c(0,0,1,0,0,0,0,0,0,0,1,1,1)),
#'   type = c(2,2,3,1,2,2,2,2,2,2,3,3,3),
#'   x = c(1,3,4,0,1.67,4,1,1.33,2.33,3.33,4.33,4,5),
#'   y = c(1,1,1,0,0.5,0,-1,0,0,-0.5,-1,-0.5,-0.5),
#'   row.names = sprintf("V%d",1:13)
#' ) %>%
#'   st_as_sf(
#'     coords=c("x","y"),
#'     crs = NA
#'   ) %>%
#'   graph %>%
#'   add.edge(
#'     from = c(1,2,1,5,4,4,5,9,4,8,9,4,7,7,6,6,9,10,10),
#'     to = c(2,3,5,2,1,5,9,2,8,9,6,7,8,9,3,13,10,12,11),
#'     data = data.frame(
#'       distance = c(4.2,4.7,3.9,3.0,3.6,2.7,4.4,3.4,3.6,3.3,
#'                    4.8,3.2,3.5,4.4,2.5,3.4,4.3,3.1,2.2),
#'       row.names = sprintf("E%d",1:19)
#'     )
#'   ) -> x
#' 
#' ## Calculation of the influence matrix:
#' y <- InflMat(x)
#' 
#' ## Obtain the number of edges:
#' nedge(y)
#' 
#' ## Obtain the edge names:
#' edgenames(y)
#' 
#' ## Obtain the edge data frame:
#' edge(y)
#' 
NULL
#' 
#' @describeIn InflMat-class
#' 
#' Influence Matrix
#' 
#' Calculates the influence matrix of a phylogenetic graph. The influence matrix
#' is a binary matrix whose rows and columns correspond to the vertices and
#' edges of the phylogenetic graph, respectively, and whose elements describe
#' whether a given edge had been taken by any ancestors of a vertex
#' (representing extinct of extant species) during evolution (value = 1) or not
#' (value = 0).
#' 
#' @export
InflMat <- function(x) {
  
  if(!inherits(x,"graph"))
    stop("Parameter 'x' must be of class 'graph'")
  
  nv <- nrow(x)
  ne <- nedge(x)
  edge <- edge(x)
  
  structure(
    matrix(
      .C(
        "InflMatC",
        as.integer(ne),
        as.integer(nv),
        as.integer(edge[[1L]]),
        as.integer(edge[[2L]]),
        B = integer(nv * ne),
        PACKAGE = "MPSEM"
      )$B,
      nrow = nv,
      ncol = ne,
      dimnames = list(
        rownames(x),
        rownames(edge)
      )
    ),
    edge = edge(x),
    class = "InflMat"
  )
}
#' 
#' @describeIn InflMat-class
#' 
#' Print Graph
#' 
#' A print method for \code{InflMat-class} objects.
#' 
#' @method print InflMat
#' 
#' @importFrom utils head tail
#' 
#' @export
print.InflMat <- function (x, ...) {
  
  nv <- nrow(x)
  ne <- ncol(x)
  edge <- edge(x)
  
  cat("\nAn Influence matrix involving", nv, "vertices and", ne, "edges\n")
  cat("--------------------------------------------------------------\n\n")
  
  if(ncol(x) > 12L) {
    cat(
      paste(
        c(
          "",
          head(colnames(x), 7L),
          paste("+", nrow(x) - 10L, " ...", sep=""),
          tail(colnames(x), 3L)
        ),
        collapse="\t"
      ),
      "\n"
    )
  } else cat(paste(colnames(x), collapse="\t"), "\n")
  
  i <- 1L
  
  while(!(i > nrow(x))) {
    
    if(ncol(x) > 12L) {
      cat(
        paste(
          c(
            rownames(x)[i],
            head(x[i,], 7L),
            paste("+", nrow(x) - 10L, " ...", sep=""),
            tail(x[i,], 3L)
          ),
          collapse="\t"
        ),
        "\n"
      )
    } else cat(paste(x[i,], collapse="\t"), "\n")
    
    if(nrow(x) > 12L) {
      if(i == 7L) {
        cat("+", nrow(x) - 10L, " ...\n", sep="")
        i <- nrow(x) - 2L
      } else i <- i + 1L
    } else i <- i + 1L
  }
  
  if(ncol(edge) > 2L) {
    cat("\nAvailable edge information: ",
        paste(colnames(edge)[-(1L:2L)], collapse = ", "), "\n")
  } else
    cat("No edge information available\n")
  
  cat("\n")
  
  invisible(NULL)
}
#' 
#' @describeIn InflMat-class
#' 
#' Number of Edges
#' 
#' Get the number of edges in an \code{InflMat-class} object.
#' 
#' @method nedge InflMat
#' 
#' @export
nedge.InflMat <- function(x) nrow(attr(x, "edge"))
#' 
#' @describeIn InflMat-class
#' 
#' Edge Extraction
#' 
#' Extracts the edges of an \code{InflMat-class} object.
#' 
#' @method edge InflMat
#' 
#' @export
edge.InflMat <- function(x) attr(x, "edge")
#' 
#' @describeIn InflMat-class
#' 
#' Edge Assignment
#' 
#' Assigns edges to an \code{InflMat-class} object.
#' 
#' @method edge<- InflMat
#' 
#' @export
`edge<-.InflMat` <- function(x, value) {
  
  if(is.null(value)) {
    
    attr(x, "edge") <- data.frame(from=integer(0L), to=integer(0L))
    
  } else if(is.data.frame(value) && !is.null(value[[1L]]) &&
            !is.null(value[[2L]])) {
    
    attr(x, "edge") <- value
    
  } else
    stop("The 'value' must be a data frame with at least two columns")
  
  x
}
#' 
#' @describeIn InflMat-class
#' 
#' Edge Names Extraction
#' 
#' Extracts the edge names of an \code{InflMat-class} object.
#' 
#' @method edgenames InflMat
#' 
#' @export
edgenames.InflMat <- function(x) rownames(attr(x, "edge"))
#' 
#' @describeIn InflMat-class
#' 
#' Edge Names Assignment
#' 
#' Assigns edge names to an \code{InflMat-class} object.
#' 
#' @method edgenames<- InflMat
#' 
#' @export
`edgenames<-.InflMat` <- function(x, value) {
  rownames(attr(x, "edge")) <- value
  x
}
#' 
