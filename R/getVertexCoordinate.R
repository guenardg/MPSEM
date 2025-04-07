## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Internal: Get Vertex Coordinates **
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
## A function to calculate vertex coordinates suitable for plotting.
## This is an internal function; it as no formal R documentation.
getVertexCoordinate <- function(x, ...) {
  
  nv <- nrow(x)
  o <- getOrigin(x)
  
  coord <- data.frame(x=numeric(nv), y=numeric(nv))
  
  type <- rep(2L, nv)
  type[o] <- 1L
  type[getTerminal(x)] <- 3L
  
  imat <- InflMat(x)
  
  if(!is.null(edge(x)$distance))
    imat <- t(t(imat)*sqrt(edge(x)$distance))
  
  rawc <- svd(scale(imat, scale=FALSE), nu=2L, nv=0L)$u
  
  ord <- getVertexMinOrder(x)
  
  ## i=1L
  for(i in 1L:max(ord)) {
    
    wh <- which((ord == i) & (type != 3L))
    n <- length(wh)
    coord$x[wh] <- seq(i - 0.25, i + 0.75, length.out=n)[order(rawc[wh,1L])]
    coord$y[wh] <- seq(0, 1, length.out=n)[order(rawc[wh,2L])]
  }
  
  i <- max(ord) + 1L
  wh <- which(type == 3L)
  n <- length(wh)
  
  coord$x[wh] <- seq(i - 0.25, i + 0.75, length.out=n)[order(rawc[wh,2L])]
  coord$y[wh] <- seq(0, 1, length.out=n)[order(rawc[wh,1L])]
  
  st_as_sf(coord, coords=c("x","y"), ...)$geometry
}
