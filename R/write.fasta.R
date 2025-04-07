## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Write sequence to a text file in FASTA format **
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
#' @name write.fasta
#' 
#' @title Write Sequences to a FASTA Text File
#' 
#' @description A function to write a vector of strings representing nucleotide
#' sequences into a text file, with optional line breaks
#' 
#' @param x A \code{\link{character}} vector storing the nucleotide sequences.
#' @param file A character string describing the connection (see the description
#' of argument \code{description} in function \code{\link{file}} help page).
#' @param linebreak Interval for line breaking to improve file readability
#' (default: \code{NULL}, meaning no line breaking).
#' @param sep Character used for line breaking (default: \code{"\n"}).
#' 
#' @returns \code{NULL} (invisibly).
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
#' tmp <- concatenate(sqn)
#' tmp
#' 
#' ## Transforming the sequence to character strings without the gaps:
#' concatenate(sqn, discard="-")
#' 
#' ## Write the sequences into a text file:
#' write.fasta(tmp, file="Sequences.fst", linebreak=15)
#' 
#' ## Clean-up:
#' file.remove("Sequences.fst")
#' rm(sqn, tmp)
#' 
#' @export
write.fasta <- function(x, file, linebreak = NULL, sep = "\n") {
  
  con <- file(file, open="w", blocking=TRUE)
  
  if(is.null(names(x)))
    names(x) <- sprintf("SEQ%d",1L:length(x))
  
  if(!is.null(linebreak)) {
    for(i in 1L:length(x)) {
      cat(sprintf(">%s\n", names(x)[i]), file=con, append=TRUE)
      ss <- seq(1L, nchar(x[i]), linebreak)
      cat(
        sprintf(
          "%s\n",
          paste(
            paste(substring(x[i],ss[-length(ss)],ss[-1L] - 1L), collapse=sep),
            substring(x[i],ss[length(ss)],nchar(x[i])), sep=sep
          )
        ),
        file=con,
        append=TRUE
      )
    }
  } else {
    for(i in 1L:length(x)) {
      cat(sprintf(">%s\n", names(x)[i]), file=con, append=TRUE)
      cat(sprintf("%s\n", x[i]), file=con, append=TRUE)
    }
  }
  
  close(con)
  
  invisible(NULL)
}
