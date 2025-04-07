## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Directed graph functions **
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
#' MPSEM Graph Manipulation Functions
#' 
#' @description A set of primitive functions for creating and munipulating
#' MPSEM graphs.
#' 
#' @name graph-functions
#' 
#' @param x A \code{graph-class} object.
#' @param data A \code{\link{data.frame}} providing data (and labels) to be
#' assigned to the vertices or edges. It can be an empty data frame with row
#' names providing the vertex labels.
#' @param from An integer or character vector. References to the origin of the
#' edges to be added.
#' @param to An integer or character vector. The destinations of the edges to be
#' added (vertex references).
#' @param id An integer or character vector. The identity of vertex or edge to
#' be removed (references).
#' 
#' @details A new graph is populated with vertices using function
#' \code{graph()}. The function must be provided with a data frame. This data
#' frame may be a 0 column data frame as long as its \code{row.names} attribute
#' contains the number of elements necessary to reference the vertices.
#' Additional vertices can be added later with function \code{add.vertex()}. The
#' graph thus created contains no edge; the latter are added using function
#' \code{add.edge()}. Edges and vertices are removed using functions
#' \code{rm.edge()} and \code{rm.vertex()}, respectively.
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution 4: 1120-1131
#' 
#' Makarenkov, V., Legendre, L. & Desdevise, Y. 2004. Modelling phylogenetic
#' relationships using reticulated networks. Zoologica Scripta 33: 89-96
#' 
#' Blanchet, F. G., Legendre, P. & Borcard, D. 2008. Modelling directional
#' spatial processes in ecological data. Ecological Modelling 215: 325-336
#' 
#' @examples
#' 
#' ## Populate a graph with 7 vertices labeled A-G and edge properties x and y:
#' data.frame(
#'   quantitative = c(1.1,7.2,7.2,4.1,5.5,6.9,3.3),
#'   factor = factor(c("A","A","A","B","B","B","B")),
#'   row.names = c("A","B","C","D","E","F","G")
#' ) %>%
#'   graph -> x
#' 
#' ## Note from package magrittr:
#' ## x %>% f(y, z)    is equivalent to f(x, y, z)
#' ## x %<>% f(y, z)   is equivalent to x <- f(x, y, z)
#' 
#' ## Add three vertices without descriptors:
#' x %>% add.vertex(
#'   data = data.frame(row.names = c("H","I","J"))
#' )
#' 
#' ## This is another way to add vertices:
#' x %<>% add.vertex(
#'   data = data.frame(
#'     factor = factor(c("C","C","C")),
#'     ordered = ordered(c(1,2,2)),
#'     row.names = c("H","I","J")
#'   )
#' )
#' 
#' ## Adding 10 edges, labeled E1-E10 and with properties d and r, to the graph:
#' x %<>% add.edge(
#'   from = c("A","B","B","C","C","D","D","E","E","F"),
#'   to = c("A","C","D","E","F","F","G","H","I","J"),
#'   data = data.frame(
#'     distance = c(1,5,2,3,3,2,5,1,1,1),
#'     reversible = c(TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE),
#'     row.names = paste("E",1:10,sep="")
#'   )
#' )
#' 
#' ## Adding three more edges, this time without variable 'reversible', but
#' ## adding a variable called 'factor':
#' x %<>% add.edge(
#'   from = c("E","F","G"),
#'   to = c("A","B","C"),
#'   data.frame(
#'     distance = c(2,3,1),
#'     factor = factor(c("S","S","G")),
#'     row.names = c("E1","E11","E23")
#'   )
#' )
#' 
#' ## Removing two edges (E3 and E5):
#' x %<>% rm.edge(id=c("E3","E5"))
#' 
#' ## Removing vertices B, F, and G with their associated edges:
#' x %<>% rm.vertex(id=c("B","F","G"))
#' 
#' 
NULL
#' 
#' 
#' @describeIn graph-functions
#' 
#' Create Graph
#' 
#' Create a graph and populates it with vertices.
#' 
#' @export
graph <- function(data = data.frame()) {
  
  if(!is.data.frame(data))
    stop("'Argument 'data' must be a data.frame")
  
  cl <- class(data)
  
  if(!("graph" %in% cl))
     cl <- c("graph", cl)
  
  structure(
    data,
    edge = data.frame(from=integer(0L), to=integer(0L)),
    class = cl
  )
}
#' 
#' @describeIn graph-functions
#' 
#' Add Vertices
#' 
#' Add vertices to an existing graph.
#' 
#' @export
add.vertex <- function(x, data = data.frame()) {
  
  if(!inherits(x, "graph"))
    stop("Argument 'x' must be a graph-class object")
  
  if(!is.data.frame(data))
    stop("Values for vertices must be provided as a data frame")
  
  if(!ncol(data)) {

    for(i in colnames(x))
      data[[i]] <- rep(NA, nrow(data))
    
  } else {
    
    for(i in colnames(x)[!(colnames(x) %in% colnames(data))])
      data[[i]] <- rep(NA, nrow(data))
    
    for(i in colnames(data)[!(colnames(data) %in% colnames(x))])
      x[[i]] <- rep(NA, nrow(x))
    
    data <- data[,colnames(x)]
    
  }

  rbind(x, data)
}
#' 
#' @describeIn graph-functions
#' 
#' Add Edges
#' 
#' Add edges to a graph.
#' 
#' @export
add.edge <- function(x, from, to, data = data.frame()) {
  
  if(!inherits(x, "graph"))
    stop("Argument 'x' must be a graph-class object")
  
  if(length(from) != length(to))
    stop("The length of 'from' (", length(from), ") does not equal the length ",
         "of 'to' (", length(to), ")")
  
  if(!is.data.frame(data))
    stop("Values for edges must be provided as a data frame")
  
  if(nrow(data) != length(from))
    stop("Argument 'data' has ", nrow(data), " rows, but ", length(from),
         " edges are being added to the graph")
  
  if(is.character(from)) {
    
    safe <- from
    from <- match(from, rownames(x))
    
    if(any(is.na(from)))
      stop("Unknown origin vertices (",
           paste(safe[which(is.na(from))],collapse=","),").")
    
  } else
    if(any(from > nrow(x)))
      stop("Unknown origin vertices (",
           paste(from[from > nrow(x)],collapse=","),").")
  
  if(is.character(to)) {
    
    safe <- to
    to <- match(to, rownames(x))
    
    if(any(is.na(to)))
      stop("Unknown destination vertices (",
           paste(safe[which(is.na(to))],collapse=","),").")
    
  } else
    if(any(to > nrow(x)))
      stop("Unknown destination vertices (",
           paste(to[to > nrow(x)],collapse=","),").")
  
  data <- cbind(from, to, data)
  
  edge <- attr(x,"edge")
  
  if(!nrow(edge)) {
    
    edge <- rbind(edge, data)
    
  } else {
    
    for(i in colnames(edge)[!(colnames(edge) %in% colnames(data))])
      data[[i]] <- rep(NA, nrow(data))
    
    for(i in colnames(data)[!(colnames(data) %in% colnames(edge))])
      edge[[i]] <- rep(NA, nrow(edge))
    
    edge <- rbind(edge, data[,colnames(edge)])
    
  }
  
  attr(x,"edge") <- edge
  
  x
}
#' 
#' @describeIn graph-functions
#' 
#' Remove Edges
#' 
#' Remove edges from a graph.
#' 
#' @export
rm.edge <- function(x, id) {
  
  if(!length(id))
    return(x)
  
  if(!inherits(x, "graph"))
    stop("Argument 'x' must be a graph-class object")
  
  edge <- attr(x,"edge")
  
  if(is.character(id)) {
    
    safe <- id
    id <- match(id, rownames(edge))
    
    if(any(is.na(id)))
      stop("Unknown edge(s) (", paste(safe[which(is.na(id))], collapse=","),
           ")")
    
  } else
    if(any(id > nrow(edge)))
      stop("Unknown edge(s) (", paste(id[id > nrow(edge)], collapse=","), ")")
  
  attr(x,"edge") <- edge[-id,]
  
  x
}
#' 
#' @describeIn graph-functions
#' 
#' Remove Vertices
#' 
#' Remove vertices from a graph.
#' 
#' @export
rm.vertex <- function(x, id) {
  
  if(!length(id))
    return(x)
  
  if(!inherits(x, "graph"))
    stop("Argument 'x' must be a graph-class object.")
  
  if(is.character(id)) {
    
    safe <- id
    id <- match(id, rownames(x))
    
    if(any(is.na(id)))
      stop("Unknown vertices (", paste(safe[which(is.na(id))],collapse=","),
           ").")
    
  } else {
    
    if(any(id > nrow(x)))
      stop("Unknown vertices (",paste(id[id > nrow(x)],collapse=","), ").")
    
  }
  
  ## Removing the edges involved with the vertices:
  which(!is.na(match(attr(x,"edge")[[1L]], id)) |
          !is.na(match(attr(x,"edge")[[2L]], id))) -> sel
  x <- rm.edge(x, id=sel)
  
  ## Re-indexing the vertices:
  mask <- rep(NA, nrow(x))
  mask[-id] <- 1L:(nrow(x) - length(id))
  attr(x,"edge")[[1L]] <- mask[attr(x,"edge")[[1L]]]
  attr(x,"edge")[[2L]] <- mask[attr(x,"edge")[[2L]]]
  
  ## Returing with the edges removed:
  `[.data.frame`(x,-id,,drop=FALSE)
}
#' 
