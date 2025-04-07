## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Sequence concatenation **
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
#' @name concatenate
#' 
#' @title Sequence Concatenation
#' 
#' @description A function to transform a \code{\link{raw}} matrix, used to
#' store sets of single nucleotides, into a vector of strings, optionally
#' removing specific value(s).
#' 
#' @param x A \code{\link{raw}} vector storing the molecular traits (single
#' nucleotides, gaps, and so on).
#' @param discard A vector of characters to exclude from the results (default:
#' \code{NULL}, meaning nothing is excluded).
#' 
#' @returns A vector of character strings.
#' 
#' @details That function concatenates the rows of a \code{\link{raw}} matrix,
#' which is used for storing the molecular trait values, into a vector of
#' strings. If argument \code{discard} is given a set of values, they will be
#' excluded from the resulting strings.
#' 
#' @author \packageAuthor{MPSEM}
#' 
#' @import magrittr
#' 
#' @examples ## Define a raw vector for storing nuceotide values:
#' c(Sequence_1 = "ATCG-TTTCG--CCCCA--TTA--TTAA-GTAA-GTAATCTTTCA",
#'   Sequence_2 = "TTGGCTTCC--TC--CAT-TTC--TTCA-GT-ACG-ATTCTTTTA",
#'   Sequence_3 = "TTCG-TACC-T-T---A-ATAA--T-AA-GTCTTGTAATCGTTCA") %>%
#'   sapply(charToRaw) %>%
#'   t -> sqn
#' 
#' ## Display the raw sequence:
#' sqn
#' 
#' ## Transforming the sequence to character strings
#' concatenate(sqn)
#' 
#' ## Transforming the sequence to character strings without the gaps:
#' concatenate(sqn, discard="-")
#' 
#' ## Clean-up:
#' rm(sqn)
#' 
#' @export
concatenate <- function(x, discard = NULL) {
  apply(
    x,
    1L,
    function(x, discard, linebreak, sep) {
      if(!is.null(discard)) x <- x[!(x %in% charToRaw(discard))]
      rawToChar(x)
    },
    discard = discard
  )
}
