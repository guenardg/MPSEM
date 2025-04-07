## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Global variable definitions and residual code **
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
#' 
#' @useDynLib MPSEM, .registration = TRUE
#' 
## Global variable definition:
globalVariables(c(
  ".",                                  ## Used by package magrittr.
  "candidate", "included", "pemModel"   ## Used by pemlm.
))
## 
PEMvar <- function(d, a = 0, psi = 1) {
  
  nd <- length(d)
  a <- rep(a, length.out = nd)
  psi <- rep(psi, length.out = nd)
  
  .C(
    "PEMvarC",
    as.double(d),
    as.integer(nd),
    as.double(a),
    as.double(psi),
    res = double(nd),
    PACKAGE = "MPSEM"
  )$res
}
## 
NULL
##
