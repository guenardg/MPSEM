## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Internal function for pairwise distance vector indexing **
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
dst_idx <- function(n, i, j) {
  
  if(missing(j)) {
    if(length(i) != 1L)
      stop("Only a single 'i' can be selected when argument 'j' is omitted.")
    j <- (1L:n)[-i]
  }
  
  nn <- ifelse(length(i) > length(j), length(i), length(j))
  
  .C(
    "dstIdxC",
    as.integer(n),
    length(i),
    length(j),
    nn,
    as.integer(i),
    as.integer(j),
    integer(nn),
    PACKAGE = "MPSEM"
  )[[7L]]
}
