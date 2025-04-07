## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Internal: Get the Minimum Order of the Vertices of a Graph **
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
## Get the mininum order of a vertex with respect to the origin(s). Though
## multiple may be present, such a situation may involve problems and the
## default behavior is to warn whenever multiple origins are present in order
## to help in solving programming issues. Setting warnMultipleOrigins to FALSE
## will turn off multiple origin warning.
## This is an internal function; it as no formal R documentation.
getVertexMinOrder <- function(x, warnMultipleOrigins=TRUE) {
  
  nv <- nrow(x)
  o <- getOrigin(x)
  
  if(warnMultipleOrigins && (length(o) != 1L))
    warning("The graph has multiple origins!")
  
  isProcessed <- logical(nv)
  ord <- integer(nv)
  isProcessed[o] <- TRUE
  
  while(!all(isProcessed))
    for(i in which(!isProcessed))
      if(all(isProcessed[attr(x,"edge")[[1L]][attr(x,"edge")[[2L]] == i]])) {
        ord[i] <- min(ord[attr(x,"edge")[[1L]][attr(x,"edge")[[2L]] == i]]) + 1L
        isProcessed[i] <- TRUE
      }
  
  ord
}
