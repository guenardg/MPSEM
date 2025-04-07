## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Random (Directed) Graph Generator **
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
#' Random Graph Generator
#' 
#' @description A function for generating a random directed graph.
#' 
#' @name randomGraph
#' 
#' @param NV An integer. The number of vertices to generate.
#' @param NC A function with a \code{...} argument returning the maximum number
#' outgoing edges from each of the vertices; effectively determining the maximum
#' number of children vertices for any given vertex.
#' @param NP A function with a \code{...} argument returning the maximum number
#' of incoming edges to each of the vertices; effectively determining the
#' maximum number of parent vertices for any given vertex.
#' @param timestep A function with a \code{...} argument returning the amount of
#' time associated with the edges.
#' @param maxDist A function with a \code{...} argument returning the maximum
#' distances, in terms of time, allowed between any two parents vertices.
#' @param ... Any arguments to be passed internally to the functions given as
#' arguments \code{NC}, \code{NP}, \code{timestep}, \code{maxDist}, or
#' \code{weighting}.
#' @param mean One of character strings \code{"arithmetic"}, \code{"geometric"},
#' \code{"harmonic"}, or any unambiguous abbreviation thereof, specifying the
#' type of mean used for averaging the distances when a vertex has multiple
#' ascendant edges.
#' @param weighting A weighting function; it takes a set of distances as its
#' first argument and returns a set of weights summing to \code{1} (default:
#' \code{\link{invDistWeighting}}).
#' @param verbose A Boolean. Whether or not to print messages associated with
#' the graph simulation process (default: \code{TRUE}).
#' @param saveDist A Boolean. Whether or not to save the graph distance matrix
#' as an attribute to the returned \code{\link{graph-class}} object (default:
#' \code{TRUE}).
#' 
#' @details Details contents...
#' 
#' @returns A \code{\link{graph-class}} object.
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @seealso \code{\link{graph-class}}.
#' 
#' @examples
#' ## Setting the RNG seed to obtain consistent examples:
#' set.seed(2182955)
#' 
#' ## A linear evolutionary sequence with random edge lengths between 2 and 5:
#' randomGraph(
#'   NV = 100,
#'   NC = function(...) 1,
#'   NP = function(...) 1,
#'   timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
#'   maxDist = function(...) NULL,
#'   ts_min = 2,
#'   ts_max = 5
#' )
#' 
#' ## As above, but allowing for dichotomic splitting.
#' randomGraph(
#'   NV = 100,
#'   NC = function(...) 2,
#'   NP = function(...) 1,
#'   timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
#'   maxDist = function(...) NULL,
#'   ts_min = 2,
#'   ts_max = 5
#' )
#' 
#' ## A random evolutionary graph with random numbers of children and parents per
#' ## node, random time steps, and a random maximum distance between the parents:
#' randomGraph(
#'   NV = 250,
#'   NC = function(lambda_child, ...) 1 + rpois(1, lambda_child),
#'   NP = function(lambda_parent, ...) 1 + rpois(1, lambda_parent),
#'   timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
#'   maxDist = function(max_anc, ...) runif(1, 0, max_anc),
#'   lambda_child = 2.5,
#'   lambda_parent = 4,
#'   ts_min = 2,
#'   ts_max = 5,
#'   max_anc = 4
#' )
#' 
#' @export
randomGraph <- function (NV, NC, NP, timestep, maxDist, ...,
                         mean = c("arithmetic","geometric","harmonic"),
                         weighting = invDistWeighting, verbose = TRUE,
                         saveDist = TRUE) {
  
  mean <- match.arg(mean)
  
  x <- graph(
    data.frame(
      species = rep(TRUE, NV),
      ncld = c(NC(...), rep(NA, NV - 1L)),
      row.names = sprintf("V%d", 1L:NV)
    )
  )
  
  diss <- numeric(NV * (NV - 1L)/2L)
  
  ecnt <- 1L
  
  for(k in 2L:NV) {
    
    asc <- which(!is.na(x$ncld) & x$ncld > 0L)
    
    if(length(asc) > 1L) {
      np <- NP(...)
      md <- maxDist(...)
      s <- sample(asc, 1L)
      r <- asc[!(asc %in% s)]
      r <- r[diss[dst_idx(NV, s, r)] <= md]
      asc <- c(s, if (length(r) > (np - 1L)) sample(r, np - 1L) else r)
    }
    
    np <- length(asc)
    
    tstp <- numeric(np)
    
    for(i in 1L:np)
      tstp[i] <- timestep(...)
    
    for(i in 1L:np) {
      x <- add.edge(
        x,
        from = asc[i],
        to = k,
        data = data.frame(
          distance = tstp[i],
          row.names = sprintf("E%d", ecnt)
        )
      )
      
      ecnt <- ecnt + 1L
      x$ncld[asc[i]] <- x$ncld[asc[i]] - 1L
      
    }
    
    x$ncld[k] <- NC(...)
    
    others <- (1L:k)[-c(asc, k)]
    
    if(length(others)) {
      
      dd <- diss[dst_idx(NV, asc[1L], others)] + tstp[1L]
      
      if(length(asc) > 1L) {
        
        w <- weighting(tstp, ...)
        
        switch(
          mean,
          arithmetic = {
            dd <- w[1L]*dd
            for(i in 2L:length(asc))
              dd <- dd + w[i]*(diss[dst_idx(NV, asc[i], others)] + tstp[i])
          },
          geometric = {
            dd <- w[1L]*log(dd)
            for(i in 2L:length(asc))
              dd <- dd + w[i]*log(diss[dst_idx(NV, asc[i], others)] + tstp[i])
            dd <- exp(dd)
          },
          harmonic = {
            dd <- w[1L]/dd
            for(i in 2L:length(asc))
              dd <- dd + w[i]/(diss[dst_idx(NV, asc[i], others)] + tstp[i])
            dd <- 1/dd
          }
        )
      }
      
      diss[dst_idx(NV, k, others)] <- dd
    }
    
    diss[dst_idx(NV, asc, k)] <- tstp
    
    if(verbose) {
      message(sprintf("%d edges:\n", length(asc)))
      message(sprintf("\t%d -> %d (%f)\n", asc, k, tstp))
      message("-----------------------\n")
    }
  }
  
  attr(x, "processOrder") <- 1L:NV
  
  if(saveDist)
    structure(
      diss,
      Size = NV,
      Labels = rownames(x),
      Diag = FALSE,
      Upper = FALSE,
      method = "patristic",
      class = "dist",
      call = match.call()
    ) -> attr(x, "dist")
  
  x$species <- TRUE
  
  if(verbose)
    message(
      sum(x$species), " terminal vertices: ",
      paste(rownames(x)[x$species], collapse = ", ")
    )
  
  x
}
#' 
