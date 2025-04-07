## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Distance-based directed graph **
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
#' Distance-Based Directed Graph
#' 
#' Calculates a distance-based directed graph from a dissimilarity matrix, a
#' threshold value, and an origin (or root) vertex.
#' 
#' @name dstGraph
#' 
#' @param d A dissimilarity matrix such as the one obtained from
#' \code{\link[stats]{dist}} or \code{\link[ape]{dist.dna}}.
#' @param th Numeric. A threshold value for dissimilarity. Vertices are
#' considered as connected whenever their pairwise dissimilarity value is
#' smaller or equal to that value.
#' @param origin Integer. Index of the origin vertex from which the edges are
#' directed.
#' @param stretch Numeric (optional). When a vertex is unreachable, stretch the
#' threshold value for the shortest edge connecting it to the rest of the graph
#' up to that value.
#' 
#' @return A \code{\link{graph-class}} object.
#' 
#' @details
#' The algorithm
#' 
#' Beginning on a user-defined origin vertex, the algorithm proceeds by
#' connecting all vertices within a given dissimilarity value from the ones that
#' have already been connected, until all the vertices that can be reached has
#' been reached. Optionally, the dissimilarity threshold value can be stretched
#' for the vertices that are unreachable. Vertices that cannot be reached in any
#' way are reported by the function.
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
#' @import magrittr
#' 
#' @examples
#' 
#' ## We set the seed to obtain a consistent example, but you are welcome to
#' ## experiment with other graphs.
#' set.seed(7653401)
#' 
#' ## Here, the dissimilarity matrix is generated from the Euclidean distance of
#' ## a two-dimensional plot for the sake of simplicity. In practice, the matrix
#' ## will come from DNA data using a dissimilarity method such as those
#' ## implemented by  ape packages's function dist.dna().
#' 
#' N <- 100
#' coords <- cbind(x=runif(N,-1,1), y=runif(N,-1,1))
#' rownames(coords) <- sprintf("N%d",1:N)
#' dst <- dist(coords)
#' 
#' ## Calculate the distance-based graph:
#' gr <- dstGraph(d=dst, th=0.25, origin=15)
#' 
#' ## This graph have unconnected vertices.
#' 
#' ## Plotting the graph with colors indicating the order of the edges:
#' plot(coords, type="n", asp=1)
#' col <- head(rainbow(max(gr$order) + 1), max(gr$order))
#' arrows(
#'   x0 = coords[edge(gr)[[1]],1],
#'   x1 = coords[edge(gr)[[2]],1],
#'   y0 = coords[edge(gr)[[1]],2],
#'   y1 = coords[edge(gr)[[2]],2],
#'   length = 0.05,
#'   col = col[gr$order[edge(gr)[[2]]]]
#' )
#' points(coords, pch=21, bg="black", cex=0.25)
#' 
#' ## Try again raising the threshold to help in connecting all the vertices:
#' gr <- dstGraph(d=dst, th=0.28, origin=15)
#' 
#' ## It helped, but does not entirely solve the matter.
#' 
#' plot(coords, type="n", asp=1)
#' col <- head(rainbow(max(gr$order) + 1), max(gr$order))
#' arrows(
#'   x0 = coords[edge(gr)[[1]],1],
#'   x1 = coords[edge(gr)[[2]],1],
#'   y0 = coords[edge(gr)[[1]],2],
#'   y1 = coords[edge(gr)[[2]],2],
#'   length = 0.05,
#'   col = col[gr$order[edge(gr)[[2]]]]
#' )
#' points(coords, pch=21, bg="black", cex=0.25)
#' 
#' ## Try again while stretching the threshold for the unconnected vertices:
#' gr <- dstGraph(d=dst, th=0.28, origin=15, stretch=0.5)
#' 
#' ## All the vertices are now connected.
#' 
#' plot(coords, type="n", asp=1)
#' col <- head(rainbow(max(gr$order) + 1), max(gr$order))
#' arrows(
#'   x0 = coords[edge(gr)[[1]],1],
#'   x1 = coords[edge(gr)[[2]],1],
#'   y0 = coords[edge(gr)[[1]],2],
#'   y1 = coords[edge(gr)[[2]],2],
#'   length = 0.05,
#'   col = col[gr$order[edge(gr)[[2]]]]
#' )
#' points(coords, pch=21, bg="black", cex=0.25)
#' 
#' ## This is the influence matrix of that directed graph:
#' tmp <- InflMat(gr)
#' 
#' ## An image plot of this influence matrix:
#' image(t(tmp[nrow(tmp):1L,]), col=gray(c(1,0)), asp=1)
#' 
## ## Generate a PEM for that graph:
## pem_gr <- PEM.build(gr)
## 
## ## Plotting the different eigenvectors one-by-one on the graph plot:
## for(i in 1:ncol(pem_gr$u)) {
##   v <- pem_gr$u[,i]
##   plot(coords, type="n", asp=1, main=colnames(pem_gr$u)[i])
##   col <- head(rainbow(max(gr$order) + 1), max(gr$order))
##   arrows(
##     x0 = coords[edge(gr)[[1]],1],
##     x1 = coords[edge(gr)[[2]],1],
##     y0 = coords[edge(gr)[[1]],2],
##     y1 = coords[edge(gr)[[2]],2],
##     length = 0.05,
##     col = col[gr$order[edge(gr)[[2]]]]
##   )
##   
##   points(coords, pch=21, bg=gray(c(0.9,0.1))[1 + (sign(v) + 1)/2],
##          cex=10*abs(v))
##   
##   if(is.null(locator(1))) break
## }
#' 
#' @export
dstGraph <- function(d, th, origin, stretch) {
  
  n <- attr(d, "Size")
  
  ## All edges associated with distances > the threshold are written off from
  ## the start by assigning NA to their distance value.
  if(!missing(stretch)) {
    ## If edge stretching is allowed:
    dc <- d
    d[d > th] <- NA
    
    for(i in 1L:n) {
      idx <- dst_idx(n,i)
      if(all(is.na(d[idx]))) {
        wh <- which(dc[idx] <= stretch)
        d[idx[wh]] <- dc[idx[wh]]
      }
    }
    
  } else {
    ## If edge stretching is not allowed:
    d[d > th] <- NA
  }
  
  graph(
    data.frame(
      order = integer(n),
      row.names = attr(d,"Label")
    )
  ) -> gr

  ord <- 1L
  oo <- origin
  
  while(length(oo)) {
    ## i=oo[1L]
    for(i in oo) {
      idx <- dst_idx(n,i)
      dd <- d[idx]
      wh <- which(!is.na(dd))
      dd <- dd[wh]
      vv <- (1L:n)[-i][wh]
      gr$order[vv] <- ord
      nvv <- length(vv)
      if(nvv)
        add.edge(
          gr,
          from = rep(i,nvv),
          to = vv,
          data = data.frame(
            distance = dd,
            row.names = sprintf("E%d", nedge(gr) + (1L:nvv))
          )
        ) -> gr
      d[idx[wh]] <- NA          ## These edges should be written off.
    }
    oo <- which(gr$order == ord)
    ord <- ord + 1L
  }
  
  unconnected <- gr$order == 0
  unconnected[origin] <- FALSE  ## The origin vertex will appear as unconnected
  
  attr(gr,"processOrder") <- getProcessOrder(gr)
  
  if(any(unconnected)) {
    message(
      sum(unconnected),
      " vertices are not reachable from the specified origin:\n",
      sprintf("%s\t",rownames(gr)[unconnected])
    )
  }
  
  ## By default, all the vertices have observations.
  gr$species <- TRUE
  
  gr
}
