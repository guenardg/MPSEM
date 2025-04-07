## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Class and method for directed graphs **
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
#' @description Class and methods to handle MPSEM graphs.
#' 
#' @docType class
#' 
#' @name graph-class
#' 
#' @param x A \code{graph-class} object.
#' @param value A vector or \code{\link{data.frame}} containing the values to be
#' given to the \code{graph-class} object.
#' @param y An optional set of numeric value available at the vertices to be
#' shown using marker of different sizes or background colors.
#' @param ylim A length 2 numeric vector giving minimum and maximum bounds for
#' \code{y}. It is only pertinent when \code{y} is given. When omitted, it is
#' estimated from the range of \code{y} values.
#' @param use.geometry A logical value specifying whether to plot the graph
#' using its specified geometry, if available (default: \code{TRUE}). When no
#' geometry is available or if \code{use.geometry = FALSE}, the graph is plotted
#' using a back bone tree (obtained using the \code{as.phylo} method) with added
#' supplementary edges.
#' @param pch An integer. Graphical parameter \code{pch} used internally by
#' function \code{\link{points}} (default: 21, a round marker).
#' @param bg Either a single character string specifying a single background
#' color or a vector of colors used as a color scale to display the values of
#' \code{y} (default: \code{"white"}).
#' @param cex.min A numeric value giving the marker size used for the smallest
#' value of \code{y} or for all markers when \code{y} is \code{NULL} (default:
#' \code{2}).
#' @param cex.max A numeric value. When \code{y} is not \code{NULL}, the marker
#' size used for the largest value of \code{y} (default: \code{cex.min}).
#' @param cex.lab A numeric value specifying the size of the vertex labels
#' (default: \code{cex.min/3}).
#' @param axes A logical; whether to show the axis graduation of the plot
#' (default: \code{FALSE}).
#' @param xlab A character string; the title of the abscissa of the plot
#' (default: \code{""}).
#' @param ylab A character string; the title of the ordinates of the plot
#' (default: \code{""}).
#' @param edge.color A list of two colors or two color vectors (mode: character)
#' specifying the color (or colors) of the edges. The first elements is used for
#' the backbone tree when not using the geometry or all edges when using the
#' geometry. The second element is only used to specify the color (or colors) of
#' the supplementary edges when not using the geometry. The default is
#' \code{list("black","red")}.
#' @param length A numeric value giving the size of the arrow heads (in inches,
#' Default: \code{0.05}).
#' @param code An integer determining kind of arrows to be drawn (see
#' \code{\link{arrows}} for the details, default: \code{2}).
#' @param show.vertex.labels A logical specifying whether to plot the vertex
#' labels (default: \code{TRUE}).
#' @param direction One of character strings \code{"rightwards"} (the default),
#' \code{"leftwards"}, \code{"upwards"}, or \code{"downwards"}, or any
#' unambiguous thereof, specifying the direction of the plot, when not using the
#' geometry.
#' @param name A literal character string or a \code{\link{name}}.
#' @param target An integer or character vector containing the indices or names
#' to be specified as target vertices.
#' @param ... Additional parameters to be internally passed to other functions
#' or methods.
#' 
#' @details Prints user-relevant information about the graph: number of edges
#' and vertices, edge and vertex labels, addition edge properties and vertex
#' properties.
#' 
#' Function \code{as.graph()} uses the MPSEM graph functions to convert a
#' phylogenetic tree of class \sQuote{phylo} (see
#' \code{\link[ape]{ape-package}}) to a \code{\link{graph-class}} object. It
#' recycles tip labels. It also creates default node labels if they were absent
#' from the \sQuote{phylo} object, and uses them as vertex labels. The resulting
#' graph can then be edited to represent cases that do not strictly conform to
#' a tree topology.
#' 
#' The method \code{as.phylo()} convert a \code{graph-class} object into a
#' class \sQuote{phylo} objects as defined in package ape. The function ensures
#' that all vertices have only a single incoming edge. When multiple incoming
#' edges are found, the one originating from the vertex that is the closer from
#' the origin (or root) of the graph is kept and any other is discarded. The
#' discarded edges are gathered into a data frame which is appended to the
#' returned class \sQuote{phylo} object as an attribute called \code{discarded}.
#' Another attribute called \code{subst} is also appended to the returned
#' object; it contains the indices of the graph vertices as they appear in the
#' returned class \sQuote{phylo} object.
#' 
#' The plot method follows two possible approaches. The graph may be assigned
#' a geometry (see \code{\link[sf]{sf-package}}), in which case the latter is
#' used to plot the graph provided that argument \code{use.geometry = TRUE},
#' which is the default. If the graph has no geometry (or if it has one, but
#' argument \code{use.geometry = FALSE}), the graph is plotted as an augmented
#' tree (or as a regular tree if it is a non-reticulated graph). This tree is
#' obtained using the \code{as.phylo} method for graph-class objects herein
#' described. #' Package ape's internal functions are used to determine the
#' coordinates of the plot. For a reticulated graph, the discarded edges are
#' added to the backbone tree obtained from the \code{as.phylo} method, thereby
#' resulting in an augmented tree. The \code{plot} method also has the
#' possibility to display trait values at the vertices using maker sizes.
#' 
#' Evolutionary distances (ie., edge lengths) are abstract measurements that
#' may represent time (eg., ka, Ma), the number of generations, or any relevant
#' metric of the evolutionary processes occurring along the edges of the
#' evolutionary graph. This information can be included in the graph as an edge
#' property (see `Format`).
#' 
#' @format Minimally, a \code{graph-class} object contains a
#' \code{\link{data.frame}} object that contains vertex properties (or a
#' zero-column data.frame when there is no vertex property), to which it adds a
#' \code{edge} attribute. The \code{edge} attribute is an other
#' \code{\link{data.frame}} object with its first two integer columns called
#' \code{from} and \code{to} that contain the indices of the origin and
#' destination vertices, respectively, for each of the edges in addition to any
#' number of edge properties.
#' 
#' A vertex property commonly found is a \code{graph-class} object is a
#' \code{\link{logical}} property called `species` that identifies which vertex
#' has trait values. Also, an edge property commonly found in a
#' \code{graph-class} object is a numeric property called `distance` that gives
#' the evolutionary distance (or length) associated with the edges. These two
#' additional properties are mandatory for calculating a \code{\link{PEM}}, but
#' they may bear different names than the ones previously mentioned.
#' 
#' @return
#' \describe{
#'   \item{print}{\code{NULL}.}
#'   \item{nedge}{An integer; the number of edges in the graph.}
#'   \item{edge}{A data frame of the edge information, namely the origin and
#'   destination of each edge (integers), and, possibly, additional edge
#'   features.}
#'   \item{edgenames}{A vector of type \sQuote{character} giving the name of
#'   each of the edges.}
#'   \item{geometry}{A geometry object (e.g., points, lines, polygons).}
#'   \item{as.phylo}{An object of class \sQuote{phylo}, possibly with an
#'   attribute called \sQuote{discarded} containing the edges that had been
#'   discarded in order to coerce the graph into a tree and (in every cases) an
#'   attribute called \sQuote{subst} containing an integer vector giving the
#'   indices of the vertices on the tree plot.}
#'   \item{plot}{When the geometry is used, the graph object given through
#'   argument \code{x}, when the geometry is not used, the object of class
#'   \sQuote{phylo} produced internally be the \code{as.phylo} method, with an
#'   attribute called \sQuote{xy} containing the display coordinates of the tree
#'   plot.}
#'   \item{locate}{A list with two elements: a \code{graph-class} object
#'   (\code{$x}) and a data frame (\code{$location}). The \code{graph-class}
#'   object is the residual graph after the targets have been removed and has
#'   two \code{\link{logical}} attributes; one called "removedVertex" used to
#'   identify any vertex of the input graph (argument \code{x}) that is no
#'   longer found in the residual graph, and one called "removedEdge" used to
#'   identify any edge of the input graph that is no longer found in the
#'   residual graph. The data frame has three columns; a first column called
#'   `ref` giving the index of the vertex or edge where each target can be
#'   found, a second column called `dist` giving the (evolutionary) distance
#'   along the edge where the target can be found (\code{NA} for a target
#'   located at a vertex) and `lca` giving the distance between the latest
#'   common ancestor of the target in the graph and the target itself.}
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
#' @seealso \code{\link{graph-functions}}, \code{\link[ape]{plot.phylo}}, and
#' \code{\link[sf]{plot.sf}}.
#' 
#' @importFrom sf st_as_sf
#' 
#' @examples ## Create an exemplary graph:
#' data.frame(
#'   species = rep(TRUE,13),
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
#' ## The object appears as follows:
#' x
#' 
#' ## The number of vertices is the number of rows in the data frame:
#' nrow(x)
#' 
#' ## The number of vertex properties is the number of columns in the data
#' ## frame:
#' ncol(x)
#' 
#' ## Methods defined for data frames are applicables to graph-class objects:
#' dim(x)
#' rownames(x)
#' colnames(x)
#' dimnames(x)
#' labels(x)
#' names(x)
#' as.data.frame(x)
#' 
#' ## MPSEM defines new generics and implements methods for them:
#' nedge(x)       ## To obtain the number of edges
#' edge(x)        ## To obtain the edge data frame
#' edgenames(x)   ## To obtain the names of the edges
#' geometry(x)    ## To obtain any geometry associated with the edges
#' 
## ## The last three methods are two-way. For instances:
## edge(x)$weight <- PEMweights(edge(x)$distance, a=0.5, psi=2.3)
## edge(x)
#' 
#' edgenames(x) <- NULL
#' edgenames(x)
#' 
#' edgenames(x) <- sprintf("E_%d",1:nedge(x))
#' edgenames(x)
#' edge(x)
#' 
#' ## Plotting the graph:
#' plot(x)
#' plot(x, use.geometry=FALSE)
#' 
#' geom <- geometry(x)
#' geometry(x) <- NULL    ## Removing the geometry.
#' plot(x)
#' 
#' geometry(x) <- geom    ## Putting the geometry back into place.
#' plot(x)
#' 
#' ## The graph is transformed or coerced into a phylogenetic tree as follows:
#' phy <- as.phylo(x)
#' plot(phy, show.node.label=TRUE)
#' 
#' ## A phylogenetic tree with 7 tips and 6 internal nodes:
#' tree1 <- read.tree(
#'   text=paste(
#'   "(((A:0.15,B:0.2)N4:0.15,C:0.35)N2:0.25,((D:0.25,E:0.1)N5:0.3,",
#'   "(F:0.15,G:0.2)N6:0.3)N3:0.1)N1:0.1;", sep=""))
#' 
#' ## Default: excluding the root vertex:
#' as.graph(tree1)
#' 
#' ## Including the root vertex:
#' as.graph(tree1, includeRoot = TRUE)
#' 
#' 
NULL
#' 
#' 
#' @describeIn graph-class
#' 
#' Print Graph
#' 
#' A print method for graph-class objects.
#' 
#' @method print graph
#' 
#' @importFrom utils head tail
#' 
#' @export
print.graph <- function (x, ...) {
  
  nv <- nrow(x)
  ne <- nedge(x)
  edge <- edge(x)
  
  cat("\nA graph with", nv, "vertices and", ne, "edges\n")
  cat("----------------------------------\n")
  
  vl <- rownames(x)
  
  if(length(vl)) {
    
    cat("Vertex labels: ")
    
    if(length(vl) > 10L) {
      cat(paste(c(head(vl,5L), paste("... +", length(vl) - 10L, "more ..."),
                  tail(vl,3L)), collapse=", "))
    } else cat(paste(vl, collapse=", "))
    
    cat("\n")
  }
  
  el <- rownames(edge)
  
  if(length(el)) {
    
    cat("Edge labels: ")
    
    if(length(el) > 10L) {
      cat(paste(c(head(el, 5L), paste("... +", length(el) - 10L, "more ..."),
                  tail(el, 3L)), collapse=", "))
    } else cat(paste(el, collapse=", "))
    
    cat("\n")
  }
  
  if(ncol(x) > 0L) {
    cat("Vertex information: ", paste(colnames(x), collapse = ", "), "\n")
  } else cat("No available vertex information\n")
  
  if(ncol(edge) > 2L) {
    cat("Available edge information: ",
        paste(colnames(edge)[-(1L:2L)], collapse = ", "), "\n")
  } else
    cat("No edge information available\n")
  
  cat("\n")
  
  invisible(NULL)
}
#' 
#' @describeIn graph-class
#' 
#' Number of Edges
#' 
#' Get the number of edges in a graph.
#' 
#' @method nedge graph
#' 
#' @export
nedge.graph <- function(x) nrow(attr(x, "edge"))
#' 
#' @describeIn MPSEM-generics
#' 
#' Edge Extraction
#' 
#' Extracts the edges of a graph-class object.
#' 
#' @method edge graph
#' 
#' @export
edge.graph <- function(x) attr(x, "edge")
#' 
#' @describeIn MPSEM-generics
#' 
#' Edge Assignment
#' 
#' Assigns edges to a graph-class object.
#' 
#' @method edge<- graph
#' 
#' @export
`edge<-.graph` <- function(x, value) {
  
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
#' @describeIn MPSEM-generics
#' 
#' Edge Names Extraction
#' 
#' Extracts the edge names of a graph-class object.
#' 
#' @method edgenames graph
#' 
#' @export
edgenames.graph <- function(x) rownames(attr(x, "edge"))
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Edge Names Assignment
#' 
#' Assigns edge names to a graph-class object.
#' 
#' @method edgenames<- graph
#' 
#' @export
`edgenames<-.graph` <- function(x, value) {
  rownames(attr(x, "edge")) <- value
  x
}
#' 
#' @describeIn graph-class
#' 
#' Extract Geometry
#' 
#' Extracts a geometry from a \code{graph-class} object.
#' 
#' @method geometry graph
#' 
#' @export
geometry.graph <- function(x)
  if(is.null(attr(x,"sf_column"))) NULL else x[[attr(x,"sf_column")]]
#' 
#' @describeIn graph-class
#' 
#' Assign Geometry
#' 
#' Assigns a geometry to a \code{graph-class} object.
#' 
#' @method geometry<- graph
#' 
#' @export
`geometry<-.graph` <- function(x, value) {
  
  if(is.null(value)) {
    
    if(inherits(x, "sf")) {
      x[[attr(x, "sf_column")]] <- NULL
      attr(x, "sf_column") <- NULL
      class(x) <- class(x)[class(x) != "sf"]
    }
    
  } else {
    
    if(nrow(x) != length(value))
      stop("The number of coordinates provided (", length(value), ") does not ",
           "equal the number of vertices (", nrow(x),")")
    
    if(inherits(x, "sf")) {
    
      x[[attr(x, "sf_column")]] <- value
    
    } else {
    
      attr(x, "sf_column") <- "geometry"
      x[[attr(x, "sf_column")]] <- value
      cl <- class(x)
      class(x) <- c(cl[1L:which(cl == "graph")], "sf",
                    cl[(which(cl == "graph") + 1L):length(cl)])
      
    }
    
  }
  
  x
}
#' 
#' 
#' @describeIn graph-class
#' 
#' Insert or Replace a Vertex Property
#' 
#' Insert or replace a vertex property into a \code{graph-class} object.
#' 
#' @export
`$<-.graph` <- function(x, name, value)
  `[[<-.data.frame`(x, name, value=value)
#' 
#' 
#' @describeIn graph-class
#' 
#' Transformation to a Tree
#' 
#' An \code{\link[ape]{as.phylo}} method to transforms a graph-class object into
#' a "phylo" class object.
#' 
#' @method as.phylo graph
#' 
#' @importFrom ape as.phylo
#' 
#' @export
as.phylo.graph <- function(x, ...) {
  
  nv <- nrow(x)
  
  po <- attr(x, "processOrder")
  if(is.null(po))
    po <- getProcessOrder(x)
  
  x <- reorderGraph(x, po)
  
  o <- getOrigin(x)
  
  if(length(o) > 1L)
    stop("The graph has to have a single origin to be transformable ",
         "into a tree.")
  
  ## The graph is coerced into as tree, if necessary.
  if(!isTree(x)) {
    
    discarded <- NULL
    
    d <- attr(x,"dist")
    if(is.null(d))
      d <- graphDist(x)
    
    do <- numeric(nv)
    do[-o] <- d[dst_idx(nv, o)]
    
    for(i in 1L:nv) {
      wh <- which(edge(x)[[2L]] == i)
      if(length(wh) > 1L) {
        
        wh <- wh[-which.min(do[edge(x)[[1L]][wh]])]
        
        as.data.frame(
          lapply(edge(x), function(x) x[wh]),
          row.names = edgenames(x)[wh],
        ) -> df
        
        discarded <- rbind(discarded, df)
        
        x <- rm.edge(x, wh)
      }
    }
    
  } else
    discarded <- NULL
  
  sw <- integer(nv)
  tip <- getTerminal(x)
  
  sw[tip] <- 1L:length(tip)
  sw[-tip] <- (length(tip) + 1L):nv
  
  if(!is.null(discarded)) {
    discarded[[1L]] <- sw[discarded[[1L]]]
    discarded[[2L]] <- sw[discarded[[2L]]]
  }
  
  subst <- integer(nv)
  subst[sw] <- po
  
  structure(
    list(
      edge = cbind(sw[edge(x)[[1L]]], sw[edge(x)[[2L]]]),
      tip.label = rownames(x)[tip],
      node.label = rownames(x)[-tip],
      Nnode = nv - length(tip),
      edge.length = edge(x)$distance
    ),
    discarded = discarded,
    subst = subst,
    class = "phylo"
  )
}
#' 
#' @describeIn graph-class
#' 
#' Transformation From a Tree Into a Graph
#' 
#' An \code{\link{as.graph}} method to transforms a "phylo" class object into a
#' \code{graph-class} object.
#' 
#' @method as.graph phylo
#' 
#' @export
as.graph.phylo <- function(x, ...) {
  
  args <- list(...)
  
  includeRoot <- if(is.null(args$includeRoot)) FALSE else args$includeRoot
  rootKnown <- if(is.null(args$rootKnown)) TRUE else args$rootKnown
  
  root <- includeRoot && !is.null(x$root.edge)
  
  if(is.null(x$node.label))
    x$node.label <- paste("N", 1L:x$Nnode, sep="")
  
  graph(
    data.frame(
      species = c(rep(TRUE, length(x$tip.label)),
                  rep(FALSE, x$Nnode), if(root) rootKnown),
      row.names = c(x$tip.label, x$node.label, if(root) "ROOT")
    )
  ) -> out
  
  add.edge(
    out,
    from = c(x$edge[,1L], if(root) nrow(out)),
    to = c(x$edge[,2L], if(root) x$edge[1L,1L]),
    data = data.frame(
      distance = c(x$edge.length, if(root) x$root.edge),
      row.names = sprintf("E%d", 1L:(nrow(x$edge) + as.integer(root)))
    )
  )
}
#' 
#' @describeIn graph-class
#' 
#' Plot Graph
#' 
#' A plotting method for graph-class objects.
#' 
#' @method plot graph
#' 
#' @importFrom graphics arrows segments points text
#' 
#' @export
plot.graph <- function(x, y, ylim, use.geometry = TRUE, pch = 21L, bg = "white",
                       cex.min = 2, cex.max = cex.min, cex.lab = cex.min/3,
                       axes = FALSE, xlab = "", ylab = "",
                       edge.color = list("black","red"), length = 0.05,
                       code = 2L, show.vertex.labels = TRUE,
                       direction = c("rightwards", "leftwards", "upwards",
                                     "downwards"),
                       ...) {
  
  .midArrows <- function(x0, y0, x1, y1, length, code, ...) {
    arrows(x0=x0, y0=y0, x1=0.5*(x0 + x1), y1=0.5*(y0 + y1), length=length,
           code=code, ...)
    segments(x0=0.5*(x0 + x1), y0=0.5*(y0 + y1), x1=x1, y1=y1, ...)
  }
  
  if(use.geometry && !is.null(geometry(x))) {
    
    xy <- st_coordinates(geometry(x))
    
    plot(xy, asp=1, type="n", xlab=xlab, ylab=ylab, axes=FALSE, ...)
    
    .midArrows(
      x0 = xy[edge(x)[[1L]],1L],
      y0 = xy[edge(x)[[1L]],2L],
      x1 = xy[edge(x)[[2L]],1L],
      y1 = xy[edge(x)[[2L]],2L],
      col = edge.color[[1L]],
      length = length,
      code = code,
      ...
    )
    
    if(missing(y)) {
      
      points(x=xy[,1L], y=xy[,2L], pch=pch, cex=cex.min, bg=bg[1L], xpd=TRUE,
             ...)
      
    } else {
      
      if(missing(ylim))
        ylim <- range(y)
      
      cc <- (y - ylim[1L])/(ylim[2L] - ylim[1L])
      cc[cc < 0] <- 0
      cc[cc > 1] <- 1
      
      points(
        x = xy[,1L],
        y = xy[,2L],
        pch = pch,
        cex = cex.min + (cex.max - cex.min)*cc,
        bg = bg[1L + floor(cc*(length(bg) - 1L))],
        xpd = TRUE,
        ...
      )
      
    }
    
    if(show.vertex.labels)
      text(
        x = xy[,1L],
        y = xy[,2L],
        labels = rownames(x),
        cex = cex.lab,
        ...
      )
    
    invisible(x)
    
  } else {
    
    direction <- match.arg(direction)
    
    tre <- as.phylo(x)
    
    xy <- getPhyloXY(tre, direction, FALSE)
    attr(tre,"xy") <- xy
    
    plot(xy, type="n", xlab=xlab, ylab=ylab, axes=FALSE, ...)
    
    .midArrows(
      x0 = xy[tre$edge[,1L],1L],
      y0 = xy[tre$edge[,1L],2L],
      x1 = xy[tre$edge[,2L],1L],
      y1 = xy[tre$edge[,2L],2L],
      col = edge.color[[1L]],
      length = length,
      code = code,
      ...
    )
    
    if(!is.null(attr(tre,"discarded")))
      .midArrows(
        x0 = xy[attr(tre,"discarded")[[1L]],1L],
        y0 = xy[attr(tre,"discarded")[[1L]],2L],
        x1 = xy[attr(tre,"discarded")[[2L]],1L],
        y1 = xy[attr(tre,"discarded")[[2L]],2L],
        col = edge.color[[2L]],
        length = length,
        code = code,
        ...
      )
    
    if(missing(y)) {
      
      points(x=xy[,1L], y=xy[,2L], pch=pch, cex=cex.min, bg=bg[1L], xpd=TRUE,
             ...)
      
    } else {
      
      if(missing(ylim))
        ylim <- range(y)
      
      cc <- (y[attr(tre,"subst")] - ylim[1L])/(ylim[2L] - ylim[1L])
      cc[cc < 0] <- 0
      cc[cc > 1] <- 1
      
      points(
        x = xy[,1L],
        y = xy[,2L],
        pch = pch,
        cex = cex.min + (cex.max - cex.min)*cc,
        bg = bg[1L + floor(cc*(length(bg) - 1L))],
        xpd = TRUE,
        ...
      )
      
    }
    
    if(show.vertex.labels)
      text(
        x = xy[,1L],
        y = xy[,2L],
        labels = c(tre$tip.label,tre$node.label),
        cex = cex.lab,
        ...
      )
    
    invisible(tre)
    
  }
}
#'
#' 
#' @describeIn graph-class
#' 
#' Locate Graph Targets
#' 
#' Calculate target locations and residual graph from an initial graph and a set
#' of targets.
#' 
#' @method locate graph
#' 
#' @export
locate.graph <- function(x, target, ...) {
  
  if(is.character(target)) {
    
    idx <- match(target, rownames(x))
    if(any(is.na(idx)))
      stop("Unknown target(s): ", paste(target[is.na(idx)], collapse=", "))
    
    data.frame(
      row.names = target,
      ref = idx,
      dist = as.double(rep(NA, length(target))),
      ladist = double(length(target))
    ) -> ttab
    
    target <- idx
    
  } else {
    
    outrange <- which((target < 0) | (target > nrow(x)))
    if(length(outrange))
      stop("Unknown target(s): ", paste(target[outrange], collapse=", "))
    
    data.frame(
      row.names = names(target),
      ref = target,
      dist = as.double(rep(NA, length(target))),
      ladist = double(length(target))
    ) -> ttab
  }
  
  if(is.null(x$species))
    stop("'x' has no vertex property called 'species'")
  
  if(!all(x$species[target]))
    stop("Non-species Target(s): ",
         paste(target[!x$species[target]], collapse=", "))
  
  ## Sanity check: is there any terminal vertex not marked as species?
  if(!all(x$species[getTerminal(x)]))
    stop("Sanity check failed: the graph has one or more terminal vertices ",
         "not marked as species.\nFunction purge.terminal() can be used to ",
         "discard them.")
  
  ## Sanity check: is there any median vertex not marked as a species?
  if(!all(x$species[getMedian(x)]))
    stop("Sanity check failed: the graph has one or more median vertices not ",
         "marked as species.\nFunction purge.median() can be used to discard ",
         "them automatically.")
  
  edge <- edge(x)
  
  if(is.null(edge$distance))
    stop("'x' has no edge property called 'distance'")
  
  vrm <- logical(nrow(x))
  erm <- logical(nrow(edge))
  trm <- logical(length(target))
  
  ## Step 1: removing any target that is terminal vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=2L
    for(i in which(!trm))
      if(!any(!erm & (edge[[1L]] == target[i]))) {
        
        up <- which(!erm & (edge[[2L]] == target[i]))
        
        if(length(up) == 1L) {
          
          erm[up] <- TRUE
          vrm[target[i]] <- TRUE
          trm[i] <- TRUE
          
          s <- which(is.na(ttab$dist) & ttab$ref == edge[up,2L])
          ttab$ref[s] <- edge[up,1L]
          ttab$ladist[s] <- ttab$ladist[s] + edge$distance[up]
          
          end <- FALSE
          break
        }
      }
  }
  
  ## Step 2: removing any target that is a median vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=2L
    for(i in which(!trm)) {
      
      down <- which(!erm & (edge[[1L]] == target[i]))
      
      if(length(down) == 1L) {
        
        up <- which(!erm & (edge[[2L]] == target[i]))
        
        if(length(up) == 1L) {
          
          if(!any(!erm & (edge[,1L] == edge[up,1L]) &
                  (edge[,2L] == edge[down,2L]))) {
            
            vrm[target[i]] <- TRUE
            
            edge[up,2L] <- edge[down,2L]
            dup <- edge$distance[up]
            edge$distance[up] <- sum(edge$distance[c(up,down)])
            
            erm[down] <- TRUE
            trm[i] <- TRUE
            
            ## Problem here...
            ## s <- which(ttab$ref == edge[down,1L])
            ## ttab$ref[s] <- up
            ## ttab$dist[s] <- ifelse(is.na(ttab$dist[s]), dup, ttab$dist[s] + dup)
            ## This is the old code.
            
            ## Treat vertices differently from edges
            ## Vertices:
            s <- which(is.na(ttab$dist) & (ttab$ref == edge[down,1L]))
            ttab$ref[s] <- up
            ttab$dist[s] <- dup
            
            ## Edges (if any):
            s <- which(!is.na(ttab$dist) & (ttab$ref == down))
            ttab$ref[s] <- up
            ttab$ladist[s] <- ttab$ladist[s] + ttab$dist[s]
            ttab$dist[s] <- dup
            
            end <- FALSE
            break
          }
        }
      }
    }
  }
  
  ## Step 3: Purging any new terminal vertex (ie., a previously non-terminal
  ## vertex that have become a terminal vertex following the removal of terminal
  ## target species) that is not marked as a species.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=10L
    for(i in which(!(vrm | x$species)))
      if(!any(!erm & (edge[[1L]] == i))) {
        
        up <- which(!erm & (edge[[2L]] == i))
        
        if(length(up) == 1L) {
          
          s <- which(is.na(ttab$dist) & ttab$ref == edge[up,2L])
          ttab$ref[s] <- edge[up,1L]
          ttab$ladist[s] <- ttab$ladist[s] + edge$distance[up]
          
          erm[up] <- TRUE
          vrm[edge[up,2L]] <- TRUE
          
          end <- FALSE
          break
        }
      }
  }
  
  ## Step 4: purging any new median vertex (ie., a previously non-median vertex
  ## that have become a median vertex following the previous removal of terminal
  ## vertice; either target species at step 1 or vertices not marked as a
  ## species at step 2-3) that is not-marked as a species vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=3L
    for(i in which(!(vrm | x$species))) {
      
      down <- which(!erm & (edge[[1L]] == i))
      
      if(length(down) == 1L) {
        
        up <- which(!erm & (edge[[2L]] == i))
        
        if(length(up) == 1L) {
          
          if(!any(!erm & (edge[,1L] == edge[up,1L]) &
                  (edge[,2L] == edge[down,2L]))) {
            
            vrm[i] <- TRUE
            edge[up,2L] <- edge[down,2L]
            dup <- edge$distance[up]
            edge$distance[up] <- sum(edge$distance[c(up,down)])
            erm[down] <- TRUE
            
            ## Treat vertices differently from edges
            ## Vertices:
            s <- which(is.na(ttab$dist) & (ttab$ref == edge[down,1L]))
            ttab$ref[s] <- up
            ttab$dist[s] <- dup
            
            ## Edges (if any):
            ## Problem here: !is.na(dist), but ttab$ref is a vertex!
            ## s <- which(!is.na(ttab$dist) & (ttab$ref == edge[down,1L]))
            ## Old code above...
            s <- which(!is.na(ttab$dist) & (ttab$ref == down))
            ttab$ref[s] <- up
            ttab$ladist[s] <- ttab$ladist[s] + ttab$dist[s]
            ttab$dist[s] <- dup
            
            end <- FALSE
            break
          }
        }
      }
    }
  }
  
  ## Vertices that cannot be removed simply lose their species status.
  x$species[target[!trm]] <- FALSE
  
  ## Recalculating the new edge indices:
  mask <- rep(NA, nrow(edge))
  mask[!erm] <- 1L:sum(!erm)
  
  ## Changing the edge indices from the target table:
  ttab$ref[!is.na(ttab$dist)] <- mask[ttab$ref[!is.na(ttab$dist)]]
  
  ## Removing the edges marked for removal:
  edge <- edge[!erm,]
  
  ## Recalculating the new vertex indices:
  mask <- rep(NA, nrow(x))
  mask[!vrm] <- 1L:sum(!vrm)
  
  ## Changing the vertex indices from the target table:
  ttab$ref[is.na(ttab$dist)] <- mask[ttab$ref[is.na(ttab$dist)]]
  
  ## Reindexing the vertices:
  edge[[1L]] <- mask[edge[[1L]]]
  edge[[2L]] <- mask[edge[[2L]]]
  
  ## Removing the vertices marked for removal and reassign the edges:
  out <- x[!vrm,,drop=FALSE]
  edge(out) <- edge
  class(out) <- c("graph",class(out)[-which(class(out) == "graph")])
  
  if(!is.null(attr(out,"processOrder")))
    attr(out,"processOrder") <- getProcessOrder(out)
  
  if(!is.null(attr(out,"dist")))
    attr(out,"dist") <- graphDist(out)
  
  attr(out,"removedVertex") <- vrm
  attr(out,"removedEdge") <- erm
  
  list(
    x = out,
    location = ttab
  )
}
#' 
