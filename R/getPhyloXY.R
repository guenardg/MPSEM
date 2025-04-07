## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Internal: Get the Vertex Coordinates from a Phylo (ape) Object **
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
## Get the vertex (nodes and tips) display coordinates of a phylo object
## (package ape). The function was obtained by taking the source code from
## package ape's function plot.phylo and stripping all the unnecessary code,
## keeping only the coordinates.
#'
#' @importFrom stats reorder
#' @importFrom ape reorder.phylo node_depth_edgelength node_height
#'
getPhyloXY <- function(
    x, direction = c("rightwards", "leftwards", "upwards", "downwards"),
    root.edge = FALSE) {
  
  if(!inherits(x, "phylo"))
    stop("Argument 'x' must be a phylo object (from R package ape)")
  
  Ntip <- length(x$tip.label)
  if (Ntip < 2L) {
    warning("Found fewer than 2 terminal vertices")
    return(NULL)
  }
  
  direction <- match.arg(direction)
  
  if(root.edge && is.null(x$root.edge)) {
    warning("Argument 'root.edge' = TRUE, but there is not root edge")
    root.edge <- FALSE
  }
  
  .nodeHeight <- function(edge, Nedge, yy)
    .C("node_height",
       as.integer(edge[,1L]),
       as.integer(edge[,2L]),
       as.integer(Nedge),
       as.double(yy),
       PACKAGE = "MPSEM")[[4L]]
  
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, edge.length)
    .C("node_depth_edgelength",
       as.integer(edge[,1L]),
       as.integer(edge[,2L]),
       as.integer(Nedge),
       as.double(edge.length),
       double(Ntip + Nnode),
       PACKAGE = "MPSEM")[[5L]]
  
  Nedge <- nrow(x$edge)
  Nnode <- x$Nnode
  
  if(any(x$edge < 1L) || any(x$edge > Ntip + Nnode))
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  
  ROOT <- Ntip + 1L
  
  horizontal <- direction %in% c("rightwards", "leftwards")
  
  x <- reorder.phylo(x)
  ## attr(x, "order")
  
  yy <- numeric(Ntip + Nnode)
  TIPS <- x$edge[x$edge[,2L] <= Ntip,2L]
  yy[TIPS] <- 1L:Ntip
  
  z <- reorder(x, order = "postorder")
  yy <- .nodeHeight(z$edge, Nedge, yy)
  xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, z$edge.length)
  
  if(!horizontal) {
    tmp <- yy
    yy <- xx
    xx <- tmp - min(tmp) + 1
  }
  
  if(root.edge) {
    if(direction == "rightwards") 
      xx <- xx + x$root.edge
    if(direction == "upwards") 
      yy <- yy + x$root.edge
  }
  
  x.lim <- c(ifelse(horizontal, 0, 1), max(xx[1L:Ntip]))
  
  if(direction == "leftwards") 
    xx <- x.lim[2L] - xx
  
  y.lim <- c(ifelse(horizontal, 1, 0), max(yy[1L:Ntip]))
  
  if(direction == "downwards") 
    yy <- y.lim[2L] - yy
  
  if(root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  
  data.frame(x = xx, y = yy)
}
