## **************************************************************************
##
##    (c) 2010-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Directed graph functions**
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
#' Graph Manipulation Functions
#' 
#' @description A set of primitive functions for creating and munipulating
#' graphs.
#' 
#' @name graph-functions
#' 
#' @param x A \code{graph-class} object.
#' @param n The number of vertex to populate a new graph (\code{pop.graph}) or
#' to add to an existing graph (\code{add.vertex}).
#' @param vertex A list of vertex properties.
#' @param edge A list of edge properties.
#' @param label Labels to be given to edges or vertices.
#' @param from The origins of the edge to be added (vertex labels or indices).
#' @param to The destinations of the edge to be added (vertex labels or
#' indices).
#' @param id Indentity (label or index) of vertex or edge to be removed.
#' @param tp Phylogenetic tree object of class \sQuote{phylo}, as defined in
#' \code{\link{ape-package}}.
#' 
#' @details A new graph can be populated with \code{n} vertices using function
#' \code{pop.graph} and vertices can be added later with function
#' \code{add.vertex}. The graphs so created contain no edges; the latter are
#' added using function \code{add.edge}. Vertices and edges are removed using
#' functions \code{rm.vertex} and \code{rm.edge}, respectively.
#' 
#' Function \code{collapse.vertex} allows one to remove a vertex while
#' reestablishing the connections between the vertices located above and below
#' that vertex using a new set of edges.
#' 
#' Function \code{Phylo2DirectedGraph} uses the graph functions to convert a
#' rooted phylogenetic tree of class \sQuote{phylo} (see
#' \code{\link{ape-package}}) to a \code{\link{graph-class}} object. It recycles
#' tip labels and creates default node labels, if they were absent from the
#' \sQuote{phylo} object, and uses them as vertex labels. The resulting acyclic
#' graph can then be edited to represent cases that do not have a tree topology.
#' 
#' @returns A \code{\link{graph-class}} object. Objects returned by
#' \code{\link{Phylo2DirectedGraph}} have a \code{\link{numeric}} edge property
#' called \sQuote{distance} featuring branch lengths and a \code{link{logical}}
#' vertex property called \sQuote{species} specifying whether a vertex is a tree
#' tip or an internal node.
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps (PEM): a framework to model and predict species traits. Meth. Ecol.
#' Evol. 4: 1120--1131
#' 
#' Makarenkov, V., Legendre, L. & Desdevise, Y. 2004. Modelling phylogenetic
#' relationships using reticulated networks. Zool. Scr. 33: 89--96
#' 
#' Blanchet, F. G., Legendre, P. & Borcard, D. 2008. Modelling directional
#' spatial processes in ecological data. Ecol. Model. 215: 325--336
#' 
#' @seealso \code{\link{graph-class}}.
#' 
#' @importFrom ape is.rooted
#' 
#' @examples
#' ## Populate a graph with 7 vertices labeled A-G having properties x and y:
#' gr <- pop.graph(n=7,
#'                 vertex=list(x=rnorm(7,0,1),y=rnorm(7,0,1)),
#'                 label=c("A","B","C","D","E","F","G"))
#' gr
#' 
#' ## Adding 3 vertices H, I, and J with property x (y is absent) and a new
#' ## property z (type character), which is unknown for A-G:
#' gr <- add.vertex(x=gr,
#'                  n=3,
#'                  label=c("H","I","J"),
#'                  vertex=list(x=rnorm(3,0,1),z=c("A","B","C")))
#' gr
#' gr$vertex
#' 
#' ## Adding 10 edges, labeled E1-E10 and with properties a and b, to the graph:
#' gr <- add.edge(x=gr,
#'                from=c("A","B","B","C","C","D","D","E","E","F"),
#'                to=c("A","C","D","E","F","F","G","H","I","J"),
#'                edge=list(a=rnorm(10,0,1),b=rnorm(10,0,1)),
#'                label=paste("E",1:10,sep=""))
#' gr
#' gr$edge
#' 
#' ## Removing edges 2, 4, and 7 from the graph:
#' print(rm.edge(gr,id=c(2,4,7)))
#' 
#' ## Removing vertices 1, 3, 7, and 10 from the graph:
#' print(rm.vertex(gr,id=c(1,3,7,10)))
#' # Notice that the edges that had one of the removed vertex as their
#' # origin or destination are also removed:
#' print.default(rm.vertex(gr,id=c(1,3,7,10)))
#' 
#' ## Vertex collapsing.
#' x <- pop.graph(n=9,label=c("A","B","C","D","E","F","G","H","I"))
#' x <- add.edge(x,from=c("A","A","B","B","C","C","D","D","E","E"),
#'               to=c("B","C","D","E","E","I","F","G","G","H"),
#'               label=paste("E",1:10,sep=""),
#'               edge=list(length=c(1,2,3,2,1,3,2,2,1,3)))
#' print.default(x)
#' for(i in c("A","B","C","D","E","F","G","H","I"))
#'   print(collapse.vertex(x,id=i))
#' 
#' if(require(ape)) {
#'   tree1 <- read.tree(
#'     text=paste(
#'       "(((A:0.15,B:0.2)N4:0.15,C:0.35)N2:0.25,((D:0.25,E:0.1)N5:0.3,",
#'       "(F:0.15,G:0.2)N6:0.3)N3:0.1)N1;",sep=""))
#'   x <- Phylo2DirectedGraph(tree1)
#'   print(x)
#' }
#' 
NULL
#' 
#' @describeIn graph-functions
#' 
#' Creates a graph and populates it with vertices.
#' 
#' @export
pop.graph <- function(n, vertex=list(), label=NULL) {
  if(!is.list(vertex))
    stop("Parameter vertex must be a list.")
  if(length(vertex))
    for(i in 1L:length(vertex))
      if(length(vertex[[i]]) != n)
        stop("Vertex property '",names(vertex)[i],"' has length ",
             length(vertex[[i]])," but the graph has ",n," vertices.")
  if(!is.null(label)) {
    if(is.character(label)) {
      if(length(label) != n)
        stop(length(label),"labels are provided, but",n,"are required.")
    } else {
      stop("Labels should be of type character.")
    }
  } else {
    label = as.character(1L:n)
  }
  return(structure(list(edge=list(numeric(0L),numeric(0L)),
                        vertex=vertex),
                   ev=c(0L,n),
                   class="graph",
                   elabel=character(0L),
                   vlabel=label))
}
#' 
#' @describeIn graph-functions
#' 
#' Adds vertices to an existing graph.
#' 
#' @export
add.vertex <- function(x,n,vertex=list(),label=NULL) {
  if(class(x) != "graph")
    stop("Parameter x must be of class graph.")
  if(!is.list(vertex))
    stop("Values for vertices must be provided as a list.")
  if(length(vertex))
    for(i in 1L:length(vertex))
      if(length(vertex[[i]]) != n)
        stop("Vertex property '",names(vertex)[i],"' has length ",
             length(vertex[[i]])," but the ",n," vertices are to be added.")
  if(!is.null(label)) {
    if(is.character(label)) {
      if(length(label) != n)
        stop(length(label)," labels are provided, but ",n," are required.")
    } else {
      stop("Labels should be of type character.")
    }
  } else {
    label = as.character(attr(x,"ev")[2L]+(1L:n))
  }
  for (i in names(x$vertex))
    x$vertex[[i]] <-
      if(is.null(vertex[[i]])) c(x$vertex[[i]],rep(NA,n)) else c(x$vertex[[i]],vertex[[i]])
  for (i in names(vertex))
    x$vertex[[i]] <-
      if(is.null(x$vertex[[i]])) c(rep(NA,attr(x,"ev")[2L]),vertex[[i]]) else x$vertex[[i]]
  attr(x,"ev")[2L] <- attr(x,"ev")[2L]+n
  attr(x,"vlabel") <- c(attr(x,"vlabel"),label)
  return(x)
}
#' 
#' @describeIn graph-functions
#' 
#' Adds edges to a graph.
#' 
#' @export
add.edge <- function(x,from,to,edge=list(),label=NULL) {
  if(class(x) != "graph")
    stop("Parameter x must be of class graph.")
  if(length(from) != length(to))
    stop("Number of origins(",length(from),") mismatch that of destinations (",
         length(to),").")
  if(!is.list(edge))
    stop("Values for edges must be provided as a list.")
  if(length(edge))
    for(i in 1L:length(edge))
      if(length(edge[[i]]) != length(from))
        stop("Edge property '",names(edge)[i],"' has length ",length(edge[[i]]),
             " but the ",length(from)," edges are to be created.")
  if(!is.null(label)) {
    if(is.character(label)) {
      if(length(label) != length(from))
        stop(length(label)," labels are provided, but ",
             length(from)," are required.")
    } else {
      stop("Labels should be of type character.")
    }
  } else {
    label = as.character(attr(x,"ev")[1L]+(1L:length(from)))
  }
  if(is.character(from)) {
    safe <- from
    from <- match(from,attr(x,"vlabel"))
    if(any(is.na(from)))
      stop("Unknown origin vertices (",
           paste(safe[which(is.na(from))],collapse=","),").")
  } else {
    if(any(from > attr(x,"ev")[2L]))
      stop("Unknown origin vertices (",
           paste(from[from > attr(x,"ev")[2L]],collapse=","),").")
  }
  if(is.character(to)) {
    safe <- to
    to <- match(to,attr(x,"vlabel"))
    if(any(is.na(to)))
      stop("Unknown destination vertices (",
           paste(safe[which(is.na(to))],collapse=","),").")
  } else {
    if(any(to > attr(x,"ev")[2L]))
      stop("Unknown destination vertices (",
           paste(to[to > attr(x,"ev")[2L]],collapse=","),").")
  }
  x$edge[[1L]] <- c(x$edge[[1L]],from)
  x$edge[[2L]] <- c(x$edge[[2L]],to)
  for (i in names(x$edge)) {
    if(i != "") {
      x$edge[[i]] <-
        if(is.null(edge[[i]])) c(x$edge[[i]],rep(NA,length(from))) else c(x$edge[[i]],edge[[i]])
    }
  }
  for (i in names(edge)) {
    x$edge[[i]] <-
      if(is.null(x$edge[[i]])) c(rep(NA,attr(x,"ev")[1]),edge[[i]]) else x$edge[[i]]
  }
  attr(x,"ev")[1L] <- attr(x,"ev")[1L]+length(from)
  attr(x,"elabel") <- c(attr(x,"elabel"),label)
  return(x)
}
#' 
#' @describeIn graph-functions
#' 
#' Removes edges from a graph.
#' 
#' @export
rm.edge <- function(x,id) {
  if(class(x) != "graph")
    stop("Parameter x must be of class graph.")
  if(is.character(id)) {
    safe <- id
    id <- match(id,attr(x,"elabel"))
    if(any(is.na(id)))
      stop("Unknown edge(s) (",paste(safe[which(is.na(id))],collapse=","),").")
  } else {
    if(any(id > attr(x,"ev")[1L]))
      stop("Unknown edge(s) (",paste(id[id > attr(x,"ev")[1L]],collapse=","),
           ").")
  }
  for (i in 1L:length(x$edge))
    x$edge[[i]] <- x$edge[[i]][-id]
  attr(x,"ev")[1L] <- attr(x,"ev")[1L]-length(id)
  attr(x,"elabel") <- attr(x,"elabel")[-id]
  return(x)
}
#' 
#' @describeIn graph-functions
#' 
#' Removes vertices from a graph.
#' 
#' @export
rm.vertex <- function(x,id) {
  if(class(x) != "graph")
    stop("Parameter x must be of class graph.")
  if(is.character(id)) {
    safe <- id
    id <- match(id,attr(x,"vlabel"))
    if(any(is.na(id)))
      stop("Unknown vertex(es) (",paste(safe[which(is.na(id))],collapse=","),
           ").")
  } else {
    if(any(id > attr(x,"ev")[2L]))
      stop("Unknown vertex(es) (",paste(id[id > attr(x,"ev")[2L]],collapse=","),
           ").")
  }
  x <- rm.edge(x,id=which(!is.na(match(x$edge[[1L]],id)) |
                            !is.na(match(x$edge[[2L]],id))))
  mask <- rep(NA,attr(x,"ev")[2L])
  mask[-id] <- 1L:(attr(x,"ev")[2L]-length(id))
  x$edge[[1L]] <- mask[x$edge[[1L]]] ; x$edge[[2L]] <- mask[x$edge[[2L]]]
  for (i in names(x$vertex))
    x$vertex[[i]] <- x$vertex[[i]][-id]
  attr(x,"ev")[2L] <- attr(x,"ev")[2L]-length(id)
  attr(x,"vlabel") <- attr(x,"vlabel")[-id]
  return(x)
}
#' 
#' @describeIn graph-functions
#' 
#' Removes vertices from a graph while also removing their associated edges.
#' 
#' @export
collapse.vertex <- function(x,id) {
  if(class(x) != "graph")
    stop("Parameter x must be of class graph.")
  if(is.character(id)) {
    safe <- id
    id <- match(id,attr(x,"vlabel"))
    if(any(is.na(id)))
      stop("Unknown vertex(es) (",paste(safe[which(is.na(id))],collapse=","),
           ").")
  } else {
    if(any(id > attr(x,"ev")[2L]))
      stop("Unknown vertex(es) (",paste(id[id > attr(x,"ev")[2L]],collapse=","),
           ").")
  }
  for(i in id) {
    up <- which(!is.na(match(x$edge[[2L]],i))) ; lup <- length(up)
    down <- which(!is.na(match(x$edge[[1L]],i))) ; ldown <- length(down)
    # If the vertex is not an intermediary between other vertex, simply remove
    # it with its edges.
    if(!(lup&ldown))
      x <- rm.vertex(x,i)
    else {
      from <- x$edge[[1L]][up] ; to <- x$edge[[2L]][down]
      from <- rep(from,each=ldown) ; to <- rep(to,lup)
      # Prevents the edge that already exist to be recreated.
      lstrip <- lup*ldown ; strip <- rep(FALSE,lstrip)
      for (j in 1L:lstrip)
        strip[j] <- any((from[j] == x$edge[[1L]]) & (to[j] == x$edge[[2L]]))
      # If all the intermediary connections already exist, simply remove the
      # vertex with its edges.
      if(all(strip))    
        x <- rm.vertex(x,i)
      else {
        if(!is.null(attr(x,"elabel"))) {
          if(!is.null(attr(x,"vlabel")))
            newlab <- paste(attr(x,"vlabel")[from[!strip]],
                            attr(x,"vlabel")[to[!strip]],sep="->")
          else
            newlab <- paste("V#",from[!strip],"->V#",to[!strip],sep="")
        } else
          newlab <- NULL
        if(!is.null(x$edge$length)) {
          ll <- list(length=rep(x$edge$length[up],each=ldown)[!strip] +
                       rep(x$edge$length[down],each=lup)[!strip])
        } else
          ll <- list()
        x <- add.edge(x,from[!strip],to[!strip],ll,newlab)
        x <- rm.vertex(x,i)
      }
    }
  }
  return(x)
}
#' 
#' @describeIn graph-functions
#' 
#' Transforms a phylogenetic tree into a directed graph.
#' 
#' @export
Phylo2DirectedGraph <- function(tp) {
  if(!is.rooted(tp))
    warning("The tree is not rooted. Direction taken from the first edge.")
  if(is.null(tp$node.label))
    tp$node.label <- paste("n",1:tp$Nnode,sep="")
  x <- pop.graph(n=tp$Nnode+length(tp$tip.label),
                 label=c(tp$tip.label,tp$node.label),
                 vertex=list(species=c(rep(TRUE,length(tp$tip.label)),
                                       rep(FALSE,tp$Nnode))))
  x <- add.edge(x,from=tp$edge[,1L],to=tp$edge[,2L],
                label=c(paste("E",1L:nrow(tp$edge),sep="")),
                edge=list(distance=tp$edge.length))
  if(!is.null(tp$root.edge))
    warning("The root edge has been omitted from the phylogenetic graph.")
  return(x)
}
##
