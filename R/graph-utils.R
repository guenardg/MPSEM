## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Directed Graph - Utility Functions **
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
#' Graph Utility Functions
#' 
#' @description A suite of graph utility functions.
#' 
#' @name graph-utils
#' 
#' @param x A \code{\link{graph-class}} object.
#' @param ...  Further argument(s) to be passed to the weighting function
#' described below.
#' @param mean One of character strings \code{"arithmetic"}, \code{"geometric"},
#' \code{"harmonic"}, or any unambiguous abbreviation thereof, specifying the
#' type of mean used for averaging the distances when a vertex has multiple
#' ascendent edges.
#' @param weighting A weighting function; it takes a set of distances as its
#' first argument and returns a set of weights summing to \code{1} (default:
#' \code{\link{invDistWeighting}}).
#' @param order An integer vector of the vertex indices.
#' @param shuffleEdge A Boolean. Whether to randomly shuffle the order that the
#' edges are stored in the \code{\link{graph-class}} object (\code{FALSE}).
#' 
#' @details A origin vertex is one having only outgoing edge(s) and no incoming
#' edge, whereas a terminal vertex is one having only incoming edge(s) and no
#' outgoing edge. A median vertex has one incoming edge and one outgoing edge. A
#' non-connected vertex has no edge, whereas a connected vertex may have
#' incoming edge(s), outgoing edge(s), or both.
#' 
#' Reordering a graph with a \code{processOrder} attribute will come with a
#' recalculation of the process order, whereas doing so on a graph with a
#' \code{dist} attribute cause the pairwise distance matrix to also be
#' reordered.
#' 
#' @return
#' \describe{
#' \item{getOrigin}{A vector of integer.}
#' \item{getConnected}{A vector of integer.}
#' \item{getNonConnected}{A vector of integer.}
#' \item{getTerminal}{A vector of integer.}
#' \item{getMedian}{A vector of integer.}
#' \item{isTree}{A \code{logical} stipulating whether the graph is a tree.}
#' \item{isDivergent}{A \code{logical} stipulating whether the graph has
#' divergence.}
#' \item{isLinear}{A \code{logical} stipulating whether the graph is a linear
#' sequence of vertices.}
#' \item{graphDist}{A pairwise distance matrix such as the one obtained from
#' function \code{\link[stats]{dist}}.}
##   \item{reorderGraph}{A \code{\link{graph-class}} object.}
#' }
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @seealso \code{\link{graph-class}}.
#' 
#' @examples ## Create and example graph with 10 vertices and 16 edges:
#' data.frame(
#'   species = rep(TRUE,10),
#'   x = c(2,3,2,4,3,4,2,1,1,0),
#'   y = c(-2,1,2,0,-0.5,-2,0,-1,1,0),
#'   row.names = sprintf("V%d",1:10)
#' ) %>%
#' graph %>%
#' add.edge(
#'   from = c(10,10,9,9,8,8,3,7,7,10,2,2,5,1,4,5),
#'   to = c(9,8,3,7,7,1,2,2,5,2,1,4,4,4,6,6),
#'   data = data.frame(
#'     distance = c(1,1,1,1,1,1,1,1,1,4,2,1,1,3,1,1),
#'     row.names = sprintf("E%d",1:16)
#'   )
#' ) -> x
#' 
#' getOrigin(x)         ## The graph has a single origin vertex.
#' getConnected(x)      ## All the vertices
#' getNonConnected(x)   ## are connected.
#' getTerminal(x)       ## The graph has a single terminal vertex.
#' getMedian(x)         ## The graph has a single median vertex.
#' isTree(x)            ## The graph is not a tree.
#' isDivergent(x)       ## The graoh has divergences.
#' isLinear(x)          ## The graph is not a linear vertex sequence.
#' 
#' ## The average pairwise distances between the vertices:
#' graphDist(x)
#' 
#' ## Reordering of the vertices:
#' xr <- reorderGraph(x, order=c(5:1,8,6,7,10,9))
#' xr
#' 
#' getOrigin(xr)     ## Same origin vertex, but at a different index.
#' getTerminal(xr)   ## Same terminal vertex, but at a different index.
#' graphDist(xr)     ## Same distances, but in a different order.
#' 
#' ## Comparison between distances obtained using various means and weighting
#' ## parameter:
#' cmpr <- function(x, y) 200*(x - y)/(x + y)
#' 
#' cmpr(graphDist(x), graphDist(x, mean="geo"))
#' cmpr(graphDist(x), graphDist(x, mean="har"))
#' cmpr(graphDist(x), graphDist(x, mean="ari", a=1))
#' cmpr(graphDist(x), graphDist(x, mean="geo", a=1))
#' cmpr(graphDist(x), graphDist(x, mean="har", a=1))
#' cmpr(graphDist(x), graphDist(x, mean="ari", a=2))
#' cmpr(graphDist(x), graphDist(x, mean="geo", a=2))
#' cmpr(graphDist(x), graphDist(x, mean="har", a=2))
#' 
#' 
NULL
#' 
#' 
#' @describeIn graph-utils
#' 
#' Get Origin Vertex
#' 
#' Obtain the origin vert(ex/ices) of a directed graph; an origin vertex is one
#' with no incoming edge.
#' 
#' @export
getOrigin <- function(x) {
  
  tmp <- rep(FALSE, nrow(x))
  tmp[attr(x,"edge")[[1L]]] <- TRUE
  tmp[attr(x,"edge")[[2L]]] <- FALSE
  
  tmp <- which(tmp)
  names(tmp) <- rownames(x)[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Get Connected Vertex
#' 
#' Obtain the connected vert(ex/ices) of a graph.
#' 
#' @export
getConnected <- function(x) {
  
  tmp <- rep(FALSE, nrow(x))
  tmp[attr(x,"edge")[[1L]]] <- TRUE
  tmp[attr(x,"edge")[[2L]]] <- TRUE
  
  tmp <- which(tmp)
  names(tmp) <- rownames(x)[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Get Non-connected Vertex
#' 
#' Obtain the non-connected connected vert(ex/ices) of a graph.
#' 
#' @export
getNonConnected <- function(x) {
  
  tmp <- rep(TRUE, nrow(x))
  tmp[attr(x,"edge")[[1L]]] <- FALSE
  tmp[attr(x,"edge")[[2L]]] <- FALSE
  
  tmp <- which(tmp)
  names(tmp) <- rownames(x)[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Get Terminal Vertex
#' 
#' Obtain the terminal vert(ex/ices) of a directed graph; a terminal vertex is
#' one with no outgoing edge.
#' 
#' @export
getTerminal <- function(x) {
  
  tmp <- rep(TRUE, nrow(x))
  tmp[attr(x,"edge")[[1L]]] <- FALSE
  
  tmp <- which(tmp)
  names(tmp) <- rownames(x)[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Get Terminal Vertex
#' 
#' Obtain the median vert(ex/ices) of a directed graph; a median vertex is one
#' with one incoming edge and one outgoing edge.
#' 
#' @export
getMedian <- function(x) {
  
  tmp <- rep(FALSE, nrow(x))
  
  for(i in 1L:nrow(x))
    if((sum(edge(x)[[1L]] == i) == 1L) && (sum(edge(x)[[2L]] == i) == 1L))
      tmp[i] <- TRUE
  
  tmp <- which(tmp)
  names(tmp) <- rownames(x)[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Tree Test
#' 
#' Testing whether the graph is a tree.
#' 
#' @export
isTree <- function(x) {
  
  for(i in 1L:nrow(x))
    if(sum(attr(x,"edge")[[2L]] == i) > 1L)
      return(FALSE)
  
  TRUE
}
#' 
#' @describeIn graph-utils
#' 
#' Divergence Test
#' 
#' Testing whether the graph has divergence.
#' 
#' @export
isDivergent <- function(x) {
  
  for(i in 1L:nrow(x))
    if(sum(attr(x,"edge")[[1L]] == i) > 1L)
      return(TRUE)
  
  FALSE
}
#' 
#' @describeIn graph-utils
#' 
#' Linearity Test
#' 
#' Testing whether the graph is a linear sequence.
#' 
#' @export
isLinear <- function(x) isTree(x) && !isDivergent(x)
#' 
#' 
#' @describeIn graph-utils
#' 
#' Graph Distance Matrix
#' 
#' Obtain a matrix of the (average) graph distance among the vertices.
#' 
#' @export
graphDist <- function (x, ..., mean = c("arithmetic","geometric","harmonic"),
                       weighting = invDistWeighting) {
  
  mean <- match.arg(mean)
  
  edge <- edge(x)
  
  if(is.null(edge$distance)) 
    stop("Edges must have a 'distance' property.")
  if(length(getOrigin(x)) > 1L) 
    stop("This procedure is only suitable for single-origin graphs.")
  if(length(getNonConnected(x))) 
    stop("All vertices must be connected.")
  
  nv <- nrow(x)
  
  ord <- attr(x, "processOrder")
  if(is.null(ord)) 
    ord <- getProcessOrder(x)
  
  diss <- numeric(nv * (nv - 1L)/2L)
  
  for (k in 2L:nv) {
    
    asc <- match(edge[[1L]][edge[[2L]] == ord[k]], ord)
    tstp <- edge$distance[edge[[2L]] == ord[k]]
    diss[dst_idx(nv, ord[asc], ord[k])] <- tstp
    
    others <- (1L:k)[-c(asc, k)]
    
    if(length(others)) {
      
      dd <- diss[dst_idx(nv, ord[asc[1L]], ord[others])] + tstp[1L]
      
      if(length(asc) > 1L) {
        
        w <- weighting(tstp, ...)
        
        switch(
          mean,
          arithmetic = {
            dd <- w[1L]*dd
            for(i in 2L:length(asc))
              dd <- dd + w[i]*(diss[dst_idx(nv, ord[asc[i]], ord[others])] + tstp[i])
          },
          geometric = {
            dd <- w[1L]*log(dd)
            for(i in 2L:length(asc))
              dd <- dd + w[i]*log(diss[dst_idx(nv, ord[asc[i]], ord[others])] + tstp[i])
            dd <- exp(dd)
          },
          harmonic = {
            dd <- w[1L]/dd
            for(i in 2L:length(asc))
              dd <- dd + w[i]/(diss[dst_idx(nv, ord[asc[i]], ord[others])] + tstp[i])
            dd <- 1/dd
          }
        )
      }
      
      diss[dst_idx(nv, ord[k], ord[others])] <- dd
    }
    
  }
  
  structure(
    diss,
    Size = nv,
    Labels = rownames(x),
    Diag = FALSE, 
    Upper = FALSE,
    method = "patristic",
    class = "dist", 
    call = match.call()
  )
}
#' 
#' @describeIn graph-utils
#' 
#' Reorder Vertices
#' 
#' Reorder the vertices of a directed graph.
#' 
#' @export
reorderGraph <- function(x, order, shuffleEdge = FALSE) {
  
  nv <- nrow(x)
  
  if(length(order) != nv)
    stop("Argument 'order' is not of the correct size.")
  
  if(is.character(order))
    order <- match(order, rownames(x))
  
  if(!all(order %in% 1L:nv))
    stop("Incorrect vertex reference(s) in 'order'")
  
  if(!all(1L:nv %in% order))
    stop("Missing / duplicated vertex reference(s) in 'order'")
  
  edge <- attr(x,"edge")
  ne <- nrow(edge)
  
  graph(
    data.frame(
      if(ncol(x))
        lapply(x, function(x, o) x[o], o=order),
      row.names = rownames(x)[order]
    )
  ) -> out
  
  if(shuffleEdge) {
    
    so <- sample(ne, replace=FALSE)
    
    add.edge(
      x = out,
      from = match(edge[[1L]], order)[so],
      to = match(edge[[2L]], order)[so],
      data = data.frame(
        lapply(edge[-(1L:2L)], function(x, o) x[o], o=so),
        row.names = rownames(edge)[so]
      )
    ) -> out
    
  } else
    add.edge(
      x = out,
      from = match(edge[[1L]], order),
      to = match(edge[[2L]], order),
      data = data.frame(
        edge[-(1L:2L)],
        row.names = rownames(edge)
      )
    ) -> out
  
  if(!is.null(attr(x,"processOrder")))
    attr(out,"processOrder") <- getProcessOrder(out)
  
  if(!is.null(attr(x,"dist")))
    attr(out,"dist") <- graphDist(out)
  
  out
}
#' 
