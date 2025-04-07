## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Graph purging functions **
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
#' MPSEM Graph Purging Functions
#' 
#' @description A set of functions for purging possibly uninformative edges from
#' MPSEM graphs.
#' 
#' @name graph-purge
#' 
#' @param x A \code{graph-class} object.
#' @param combine A function for combining the distances of redundant edges,
#' should they occur. The default function calculates the inverse of the sum of
#' the inverse distances.
#' @param ... Further argument to be internally passed to the function given as
#' argument \code{combine} or function \code{\link{graphDist}}.
#' 
#' @details Function \code{purge.terminal} or \code{purge.median} will only
#' purge the vertices that are not marked as species. Removal follows certain
#' constraints. Function \code{purge.terminal} removes the terminal vertices,
#' irrespective of the number of incoming edges. Function \code{purge.median}
#' removes the median vertices (vertices with a single incoming edge and a
#' single outgoing edge) by joining the incoming and outgoing edges into a
#' single edge (having the name of the incoming edge). When joining is done,
#' care is taken that any edge having the same origin vertex and destination
#' vertex be consolidated. The default function for consolidating the distances
#' is the inverse of the sum of the inverse distances. For all the other edge
#' characteristics, values of the incoming edge involved in the last removal is
#' taken.
#' 
#' 
#' 

#' The latter are vertices that have only a single incoming edge and a single
#' outgoing edge, whereas the former are vertices that have no outgoing edges.
#' These vertices are generally uninformative for phylogenetic modelling and do
#' not carry known trait values; thus making them safe for removal.
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @return The purged \code{graph-class} object, possibly with attributes
#' \code{removedVertex} (whenever vertices has to be removed) and/or
#' \code{removedEdge} (whenever edges had to be removed).
#' 
#' @examples
#' ## A 16-vertex graph with 24 edges:
#' data.frame(
#'   species = as.logical(c(1,0,0,0,1,0,0,0,1,1,1,1,1,0,0,0)),
#'   x = c(1,3,4,0,1.67,4,1,1.33,2.33,3.33,4.33,4,5,5,5,2.33),
#'   y = c(1,1,1,0,0.5,0,-1,0,0,-0.5,-1,-0.5,-0.5,1,0.5,-1),
#'   row.names = sprintf("V%d",1:16)
#' ) %>%
#'   st_as_sf(
#'     coords=c("x","y"),
#'     crs = NA
#'   ) %>%
#'   graph %>%
#'   add.edge(
#'     from = c(1,2,1,5,4,4,5,9,4,8,9,4,7,7,3,6 ,9 ,10,10,3 ,3 ,7 ,9, 10),
#'     to =   c(2,3,5,2,1,5,9,2,8,9,6,7,8,9,6,13,10,12,11,14,15,16,16,16),
#'     data = data.frame(
#'       distance = c(4.2,4.7,3.9,3.0,3.6,2.7,4.4,3.4,3.6,3.3,4.8,3.2,3.5,
#'                    4.4,2.5,3.4,4.3,3.1,2.2,2.1,0.9,1.0,2.1,0.9),
#'       row.names = sprintf("E%d",1:24)
#'     )
#'   ) -> gr1
#' 
#' ## Plotting the exemplary graph:
#' plot(gr1)
#' 
#' ## Purging the terminal vertices:
#' tmp <- purge.terminal(gr1)
#' plot(tmp)
#' attr(tmp,"removedVertex")
#' attr(tmp,"removedEdge")
#' 
#' ## Purging the median vertices:
#' tmp2 <- purge.median(tmp)
#' plot(tmp2)
#' attr(tmp2,"removedVertex")
#' attr(tmp2,"removedEdge")
#' 
NULL
#' 
#' @describeIn graph-purge
#' 
#' Purge Terminal Vertices
#' 
#' Attempts to purge the terminal vertices of a graph that are not marked as
#' species.
#' 
#' @export
purge.terminal <- function (x, ...) {
  
  if(!inherits(x, "graph")) 
    stop("Argument 'x' is not a graph-class object")
  
  if(is.null(x$species)) 
    stop("'x' has no vertex property called 'species'")
  
  if(all(x$species)) {
    
    attr(x,"removedVertex") <- integer(0L)
    attr(x,"removedEdge") <- integer(0L)
    
    return(x)
  }
  
  removedVertex <- logical(nrow(x))
  removedEdge <- logical(nedge(x))
  
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=11L
    for(i in which(!removedVertex))
      if(!x$species[i])
        if(!any(!removedEdge & (edge(x)[[1L]] == i))) {
          
          removedEdge[!removedEdge & (edge(x)[[2L]] == i)] <- TRUE
          removedVertex[i] <- TRUE
          
          end <- FALSE
          break
        }
  }
  
  removedVertex <- which(removedVertex)
  removedEdge <- which(removedEdge)
  
  if(length(removedEdge))
    x <- rm.edge(x, removedEdge)
  
  if(length(removedVertex))
    x <- rm.vertex(x, removedVertex)
  
  if(length(removedVertex) || length(removedVertex)) {
    
    if(!is.null(attr(x,"processOrder")))
      attr(x,"processOrder") <- getProcessOrder(x)
    
    if(!is.null(attr(x,"dist")))
      attr(x,"dist") <- graphDist(x, ...)
  }
  
  attr(x,"removedVertex") <- removedVertex
  attr(x,"removedEdge") <- removedEdge
  
  x
}
#' 
#' @describeIn graph-purge
#' 
#' Purge Median Vertices
#' 
#' Attempts to purge the median vertices of a graph that are not marked as
#' species, connecting the two end vertices.
#' 
#' @export
purge.median <- function(x, combine = function(x, ...) 1/sum(1/x), ...) {
  
  if(!inherits(x, "graph")) 
    stop("Argument 'x' is not a graph-class object")
  
  if(is.null(x$species)) 
    stop("'x' has no vertex property called 'species'")
  
  if(all(x$species)) {
    
    attr(x,"removedVertex") <- integer(0L)
    attr(x,"removedEdge") <- integer(0L)
    
    return(x)
  }
  
  if(is.null(edge(x)$distance)) 
    stop("'x' has no edge property called 'distance'")
  
  removedVertex <- logical(nrow(x))
  removedEdge <- logical(nedge(x))
  
  end <- FALSE
  while (!end) {
    end <- TRUE
    
    for(i in which(!removedVertex))
      if(!x$species[i]) {
        
        down <- which(!removedEdge & (edge(x)[[1L]] == i))
        
        if(length(down) == 1L) {
          
          up <- which(!removedEdge & (edge(x)[[2L]] == i))
          
          if(length(up) == 1L) {
            
            removedVertex[i] <- TRUE
            edge(x)[up,2L] <- edge(x)[down,2L]
            edge(x)$distance[up] <- sum(edge(x)$distance[c(up, down)])
            removedEdge[down] <- TRUE
            
            ## Begin check for edge duplicates
            for(j in 1L:nedge(x))
              if(!removedEdge[j]) {
                
                which(
                  !removedEdge & (edge(x)[[1L]][j] == edge(x)[[1L]]) &
                    (edge(x)[[2L]][j] == edge(x)[[2L]])
                ) -> w
                
                if(length(w) > 1L) {
                  
                  edge(x)$distance[w[1L]] <- combine(edge(x)$distance[w], ...)
                  removedEdge[w[-1L]] <- TRUE
                }
              }
            ## End check for edge duplicates
            
            end <- FALSE
            break
          }
        }
      }
  }
  
  removedVertex <- which(removedVertex)
  removedEdge <- which(removedEdge)
  
  if(length(removedEdge))
    x <- rm.edge(x, removedEdge)
  
  if(length(removedVertex))
    x <- rm.vertex(x, removedVertex)
  
  if(length(removedVertex) || length(removedVertex)) {
    
    if(!is.null(attr(x,"processOrder")))
      attr(x,"processOrder") <- getProcessOrder(x)
    
    if(!is.null(attr(x,"dist")))
      attr(x,"dist") <- graphDist(x, ...)
  }
  
  attr(x,"removedVertex") <- removedVertex
  attr(x,"removedEdge") <- removedEdge
  
  x
}
#' 
