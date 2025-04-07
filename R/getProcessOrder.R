## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Internal: Get the Vertex Processing Order **
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
## Get a vertex processing order. This order is suitable to perform calculations
## that require that the calculations had first been performed in all the
## preceding vertices.
## This is an internal function; it as no formal R documentation.
getProcessOrder <- function(x) {
  
  ord <- integer(nrow(x))
  is_proc <- rep(TRUE, nrow(x))
  edge <- edge(x)
  
  is_proc[edge[[2L]]] <- FALSE
  wh <- which(is_proc)
  ord[1L:length(wh)] <- wh
  i <- length(wh) + 1L
  
  while(!all(is_proc)) {
    for(j in which(!is_proc)) {
      if(all(is_proc[edge[[1L]][edge[[2L]] == j]])) {
        ord[i] <- j
        is_proc[j] <- TRUE
        i <- i + 1L
      }
    }
  }
  
  ord
}
