## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Show DNA Sequence **
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
#' @name show.sequence
#' 
#' @title Show DNA Sequences
#' 
#' @description A function displaying a set of DNA sequences in the form of a
#' mosaic of colored tiles.
#' 
#' @param x A vector containing the DNA sequences (ATCG-; character).
#' @param xlim The range of nucleotides to show (a length 2 numeric vector).
#' @param ylim The range of sequences to show (a length 2 numeric vector).
#' @param xlab An optional label for the x-axis (character).
#' @param text Whether to show the letters associated with the DNA bases
#' (default: \code{FALSE}; logical).
#' @param col A named vector of color values with which to show the nucleotides.
#' @param ... Further arguments (graphical parameters) to be passed internally
#' to function \code{\link{plot}} and \code{\link{text}}.
#' 
#' @author \packageAuthor{MPSEM}
#' 
#' @returns \code{NULL} (invisibly).
#' 
#' @examples ## A set of exemplary DNA sequences:
#' sqn <- c(Sequence_0 = "ATCG-TTT-G--C--CA--TTA--TTAA-GTGC-GTGGTCT-TCA",
#'          Sequence_1 = "ATCG-TTTCG--CCCCA--TTA--TTAA-GTAA-GTAATCTTTCA",
#'          Sequence_2 = "TTGGCTTCC--TC--CAT-TTC--TTCA-GT-ACG-ATTCTTTTA",
#'          Sequence_3 = "TTCG-TACC-T-T---A-ATAA--T-AA-GTCTTGTAATCGTTCA",
#'          Sequence_4 = "TTCG-TACC-T-T-T-C-ATAA--T-AA-GTCTGGTGGTCGTTCA",
#'          Sequence_5 = "TTCG-TACG-T-T-T-C-ATAT--T-AA-GTCTGGTGGTCGTTCA",
#'          Sequence_6 = "TTCG-TACG-T-T-T-GATTTT--T-AA-GTCTGGT---CGTTCA",
#'          Sequence_7 = "TTAG-TACG-T-T-T-GATTTT--T-TT-G---GGT---CGTTTT",
#'          Sequence_8 = "TTAG-TACG-T-T-T-GATTTT--T-TT-G---GA----CGTTTT",
#'          Sequence_9 = "TTAG-TACG-TATCT-GAT--TAAT-TT-G---GA----CG--TA")
#' 
#' ## With all the defaults:
#' show.sequence(sqn)
#' 
#' ## Displays the nucleotides and the gaps:
#' show.sequence(sqn, text=TRUE, cex=0.75)
#' 
#' ## Clean up:
#' rm(sqn)
#' 
#' @importFrom grDevices dev.flush dev.hold
#' @importFrom graphics axis par rect
#' 
#' @export
show.sequence <- function(
    x, xlim, ylim, xlab = "Position", text = FALSE,
    col = c(A="red",G="yellow",C="green",`T`="blue", `-`="grey"), ...) {
  
  par(mar=c(4,7,1,1))
  
  if(missing(xlim))
    xlim <- c(1, max(nchar(x)))
  
  if(missing(ylim))
    ylim <- c(length(x), 0)
  
  dev.hold()
  
  plot(NA, xlim=xlim, ylim=ylim, axes=FALSE, xlab=xlab, ylab="", ...)
  
  axis(1L)
  
  axis(2L, at=1:length(x), labels = names(x), las=1)
  
  for(i in 1L:length(x)) {
    cc <- unlist(strsplit(x[i],""))
    for(j in 1L:length(cc)) {
      rect(j - 1, i - 0.5, j, i + 0.5, col = col[cc[j]], border = col[cc[j]])
      if(text) text(j - 0.5, i, cc[j], ...)
    }
  }
  
  dev.flush()
  
  invisible(NULL) 
  
}
