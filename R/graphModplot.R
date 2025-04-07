## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Edit the Plot of a Graph Object **
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
#' Graph Plot Editor
#' 
#' @name graphModplot
#' 
#' @description Interactively modify the plot of a graph object by manually
#' changing the location where vertices are being placed.
#' 
#' @param x A \code{\link{graph-class}} object.
#' @param col Color of the arrows representing the edges. Default is
#' \code{"grey80"}.
#' @param bg Background colors for the origin, intermediate, and terminal
#' vertices, respectively. Default is \code{c("red","black","blue")}.
#' @param pch Plotting 'character' (see \code{\link[graphics]{points}}) for the
#' details. Default is \code{21} (a solid dot with a colored background).
#' @param length Length of the edge of the arrow head (see
#' \code{\link[graphics]{arrows}}). The default value is \code{0.05}.
#' @param pt.cex The relative point size used for the vertex markers. The
#' default value is \code{0.75}.
#' @param ... Additional parameters to be passed to other functions of methods.
#' 
#' @details After calling the function, a mouse-click on the display will
#' select the nearest vertex. Which one has been selected is indicated by the
#' background of its marker turning white, with its incoming edges appearing as
#' red dashed segments and its outgoing edges appearing as blue dashed segments.
#' A second mouse-click will cause the vertex to be relocated at the location
#' of the mouse pointer.
#' 
#' @returns A point geometry object.
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @seealso \code{\link{graph-class}}.
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
#' ## Original coordinates:
#' plot(x)
#' 
#' ## Edit the vertex locations manually:
#' x <- graphModplot(x)
#' 
#' ## Plot with the new coordinates
#' plot(x)
#' 
#' 
#' @importFrom graphics locator segments
#' @importFrom sf st_coordinates
#' 
#' @export
graphModplot <- function(x, col = "grey80", bg = c("red","black","blue"),
                         pch = 21L, length = 0.05, pt.cex = 0.75, ...) {
  
  if(is.null(geometry(x)))
    stop("The graph must have a geometry")
  
  nv <- nrow(x)
  
  if(is.null(x$type)) {
    x$type <- rep(2L, nv)
    x$type[getOrigin(x)] <- 1L
    x$type[getTerminal(x)] <- 3L
  }
  
  par <- par(no.readonly=TRUE)
  on.exit(par(par))
  par(mar=c(1,1,1,1))
  
  coord <- st_coordinates(geometry(x))
  
  repeat {
    
    plot(NA, xlim=range(coord[,1L]), ylim=range(coord[,2L]), type="n", asp=1,
         axes=FALSE, ...)
    
    arrows(
      x0 = coord[edge(x)[[1L]],1L],
      x1 = coord[edge(x)[[2L]],1L],
      y0 = coord[edge(x)[[1L]],2L],
      y1 = coord[edge(x)[[2L]],2L],
      col = col,
      length = length,
      ...
    )
    
    points(
      x = coord[,1L],
      y = coord[,2L],
      pch = pch,
      bg = bg[x$type],
      cex = pt.cex,
      ...
    )
    
    if(!is.null(xy <- locator(1L))) {
      
      i <- which.min(sqrt((coord[,1L] - xy$x)^2 + (coord[,2L] - xy$y)^2))
      
      inc <- which(edge(x)[[2L]] == i)
      
      if(length(inc))
        segments(
          x0 = coord[edge(x)[[1L]][inc],1L],
          x1 = coord[edge(x)[[2L]][inc],1L],
          y0 = coord[edge(x)[[1L]][inc],2L],
          y1 = coord[edge(x)[[2L]][inc],2L],
          col = "red",
          lwd = 3,
          lty = 2L
        )
      
      dec <- which(edge(x)[[1L]] == i)
      
      if(length(dec))
        segments(
          x0 = coord[edge(x)[[1L]][dec],1L],
          x1 = coord[edge(x)[[2L]][dec],1L],
          y0 = coord[edge(x)[[1L]][dec],2L],
          y1 = coord[edge(x)[[2L]][dec],2L],
          col = "blue",
          lwd = 3,
          lty = 2L
        )
      
      points(x=coord[i,1L], y=coord[i,2L], pch=pch, bg="white", cex=2*pt.cex)
      
      if(!is.null(xy <- locator(1L)))
        coord[i,] <- c(xy$x,xy$y)
    } else
      break
  }
  
  geometry(x) <- st_as_sf(as.data.frame(coord), coords=c("X","Y"))$geometry
  
  cl <- class(x)
  wh <- which(cl == "graph")
  if(wh != 1L) {
    cl <- c(cl[wh],cl[-wh])
    class(x) <- cl
  }
  
  x
}
#' 
