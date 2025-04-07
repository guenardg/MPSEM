## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Molecular Evolution Simulator (Markov Process) **
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
#' @name molEvolSim
#' @aliases DNArate drawDNASequence drawEvolRate simulateSequence
#' 
#' @title Molecular Evolution Simulator
#' 
#' @description Functions to simulate the evolution of DNA sequences following a
#' Markov substitution process.
#' 
#' @param model A character string identifying one of eight DNA evolution
#' models, namely \code{"JC69"}: Jukes and Cantor (1969); \code{"K80"}: Kimura
#' (1980; aka K2P); \code{"K81"}: Kimura (1981, aka K3P); \code{"F81"}:
#' Felsenstein (1981); \code{"HKY85"}: Hasegawa, Kishino and Yano (1985);
#' \code{"T92"}: Tamura (1992); \code{"TN93"}: Tamura and Nei (1993);
#' \code{"GTR"}: Tavaré (1986); or any unambiguous abbreviation thereof
#' (default: \code{"JC69"}).
#' @param piGap Numeric; the equilibrium frequency of the gaps (a value
#' between 0 and 1, default: 0).
#' @param insertionRate Numeric; the rate of nucleotide insertion (default: 0).
#' @param deletionRate Numeric; the rate of nucleotide deletion (default: 0).
#' @param pi Optional. A length 4 numeric between 0 and 1, summing to 1 and
#' giving the equilibrium base frequencies or the probabilities for drawing one
#' of the four DNA bases (A, G, C, or T, in that order). When missing, the
#' function assumes equal equilibrium base frequencies
#' (\code{pi = rep(0.25, 4)}).
#' @param piGC Optional. A numeric between 0 and 1; the GC relative frequency
#' with respect to AT. When missing, the function assumes equal equilibrium base
#' frequencies (\code{piGC = 0.5}).
#' @param par Optional. A set with 1 to 6 evolution rate parameters (numeric;
#' see the details below).
#' @param Q A 5 x 5 shift intensity matrix (transition and transversion) such as
#' the one obtained from function \code{DNAevol} (see details).
#' @param step The simulation time step (arbitrary units of time).
#' @param rho An evolution rate factor to be used on top of the shift intensity
#' matrix.
#' @param NN An integer; the number of locations to be generated. A locations is
#' either one of the four nucleotides or a gap.
#' @param gamma.shape A numeric; shape parameter of the beta distribution used
#' to draw the nucleotide (and gaps) evolution rates.
#' @param gamma.scale A numeric; scale parameter of the beta distribution used
#' to draw the nucleotide (and gaps) evolution rates.
#' @param x A \code{\link{graph-class}} object.
#' @param sqn A raw vector or matrix containing one or more initial DNA
#' sequence.
#' @param rate A numeric vector of mean evolution rates.
#' @param contrib A function that determine the contribution of the parent
#' vertices as a function of the evolutionary distances (passed as its first
#' argument. It may have any number of argument, each with a default value, and 
#' must accept arbitrary arguments (...). Default:
#' \code{function(x, a = 0, ...) (x^-a)/sum(x^-a)}.
#' @param ... Any named arguments to be internally passed to other functions or
#' methods.
#' 
#' @return
#' \describe{
#'   \item{ \code{DNArate} }{A 5 x 5 shift intensity matrix (transition and
#'   transversion)}
#'   \item{ \code{molEvolSim} }{...}
#'   \item{ \code{drawDNASequence} }{A raw vector of length \code{NN} containing
#'   the ASCII values for characters '-' (0x2d), 'A' (0x41), 'C' (0x43), 'G'
#'   (0x47), and 'T' (0x54) representing random nucleotides to be used as the
#'   seed sequence for a DNA fragment.}
#'   \item{ \code{drawEvolRate} }{A numeric vectors of length \code{NN}
#'   containing the evolution rate for each of the nucleotides in the sequence.}
#'   \item{ \code{simulateSequence} }{A raw matrix with as many rows as the
#'   number of vertices in the graph and as many columns as the number of
#'   nucleotides in the simulated sequence. The elements of the matrix are
#'   ASCII values for characters '-' (0x2d), 'A' (0x41), 'C' (0x43), 'G' (0x47),
#'   and 'T' (0x54) representing the simulated nucleotides.}
#' }
#' 
#' @details The molecular evolution model is based on a set of locations that
#' can take one of five states, namely gap ('-'), adenine ('A'), guanine ('G')
#' cytosine ('C'), or Thymine ('T'). A seed sequence is evolved one location at
#' a time. Changes from one nucleotides to another appear as nucleotide
#' transitions or transversions (shifts), whereas changes from a gap to one of
#' the nucleotides appear as an insertion and changes from one the nucleotides
#' to a gap as a deletion.
#' 
#' The changes are simulated as a simple Markov process, using a shift
#' probability matrix, which is calculated using a shift intensity matrix
#' \code{Q}. The off-diagonal elements of this matrix can be calculated from a
#' DNA evolution model by \code{DNArate()}, The diagonal elements are defined by
#' construct. Multiple shift intensity matrices can be employed for various
#' simulated sequences, it necessary.
#' 
#' Gap-only locations are discarded at the outset of the process, yielding sets
#' of sequences that are more or less shorter than the prescribed number of
#' nucleotides depending on the gap frequency used in the initial sequence. The
#' initial (root) sequence can be drawn from a uniform distribution with
#' user-defined frequencies (i.e., using function
#' \code{link{drawDNASequence}}), whereas each location's evolution rate can be
#' drawn from a gamma distribution (i.e., using function
#' \code{link{drawEvolRate}}).
#' 
#' The molecular evolution simulator can be instantiated multiple times, by
#' calling \code{molEvolSim} multiple times and storing the results into a list.
#' For instance, if every single location is assigned its own evolution rate, a
#' simulator is implemented for every single location. Then, member function
#' \code{$evolve(N)} is called to evolve a location 'N', with the returned value
#' being the new location value. When the time step of evolution rate change,
#' member function \code{$recalculate(step, rho)} is called to update the
#' shift probability matrix with time step 'step' and evolution rate 'rho'.
#' The shift matrix is obtained using member function \code{$getMt()}.
#' 
#' The gap opening and closure rates (arguments \code{gepOpen} and
#' \code{gspClose}, respectively) correspond to the total of the four changes
#' (from A, G, C, or T to a gap in the case of the gap opening rate and from a
#' gap to A, G, C, or T, in the case of a gap closure).
#' 
#' Models \code{"JC69"}, \code{"K80"}, and \code{"K81"} do not require
#' equilibrium base frequencies (argument \code{pi}), since they assume equal
#' equilibrium base frequencies, whereas model \code{"T92"} uses parameter
#' \code{piGC} to calculate these frequencies. Therefore, these models will
#' disregard any values provided to argument \code{pi}. For the other model,
#' unequal equilibrium base frequencies may be provided through argument
#' \code{pi}. Otherwise, equal equilibrium base frequencies are assumed in all
#' cases. Argument \code{par} is disregarded by models \code{"JC69"} and
#' \code{"F81"}. The number of values provided to that argument varies depending
#' on the model. \code{"K80"}, \code{"HKY85"}, and \code{"T92"} require a single
#' value, \code{"K81"} and \code{"TN93"}, require two; and \code{"GTR"} requires
#' six. In all cases, omitting \code{par} will yield a model equivalent to
#' \code{"F81"} (or \code{"JC69"} with equal equilibrium base frequencies).
#' Permissible values for argument \code{par} are positive values and are
#' subjected to further constraints that depend on the model.
#' 
#' The \code{"JC69"} and \code{"F81"} models take no \code{par} value. The
#' value of the single parameter (\code{par[1]}) of models \code{"K80"},
#' \code{"HKY85"}, and \code{"T92"}, as well as the two parameters of model
#' \code{"TN93"}, have to be between 0 and 1. For the \code{"K81"} model, the
#' two parameters also must have a sum smaller than, or equal to 1. Finally, the
#' six parameters of the \code{"GTR"} model must sum to 6, but can be
#' individually larger than 1.
#' 
#' Function \code{simulateSequence()} is assigned a random sequence at the
#' origin vertex (or vertices) of the evolutionary graph through its argument
#' \code{sqn}, and also takes a set of evolution rate (one for each location)
#' through its argument \code{rate}. The DNA sequence evolution is simulated as
#' a Markov process described by the transition intensity matrix given as
#' argument \code{Q}.
#' 
#' @references
#' Jukes, T.H. & Cantor, C.R. (1969). Evolution of Protein Molecules.
#' New York: Academic Press pp. 21-132. doi:10.3389/fgene.2015.00319
#' 
#' Kimura, M. 1980. A simple method for estimating evolutionary rates of base
#' substitutions through comparative studies of nucleotide sequences. Journal of
#' Molecular Evolution 16(2): 111-120. doi:10.1007/BF01731581
#' 
#' Kimura, M. 1981. Estimation of evolutionary distances between homologous
#' nucleotide sequences. Proceedings of the National Academy of Sciences of the
#' United States of America 78(1): 454-458. doi:10.1073/pnas.78.1.454
#' 
#' Felsenstein, J. 1981. Evolutionary trees from DNA sequences: a maximum
#' likelihood approach. Journal of Molecular Evolution 17(6): 368-376.
#' doi:10.1007/BF01734359
#' 
#' Hasegawa, M.; Kishino, H. & Yano, T. 1985. Dating of the human-ape splitting
#' by a molecular clock of mitochondrial DNA. Journal of Molecular Evolution
#' 22(2): 160-174. doi:10.1007/BF02101694
#' 
#' Tamura, K. 1992. Estimation of the number of nucleotide substitutions when
#' there are strong transition-transversion and G+C-content biases. Molecular
#' Biology and Evolution. 9(4): 678-687.
#' doi:10.1093/oxfordjournals.molbev.a040752
#' 
#' Tamura, K. & Nei, M. 1993. Estimation of the number of nucleotide
#' substitutions in the control region of mitochondrial DNA in humans and
#' chimpanzees. Molecular Biology and Evolution 10(3): 512-526.
#' doi:10.1093/oxfordjournals.molbev.a040023
#' 
#' Tavaré S. 1986. Some Probabilistic and Statistical Problems in the Analysis
#' of DNA Sequences. Lectures on Mathematics in the Life Sciences. 17: 57-86.
#' 
#' @author \packageAuthor{MPSEM}
#' 
#' @importFrom expm expm
#' @importFrom stats rgamma
#' 
#' @examples ## Examples of various molecular evolution models
#' 
#' ## Jukes and Cantor (1969):
#' DNArate("JC69", piGap=0.30, insertionRate=0.02, deletionRate=0.02)
#' 
#' ## Kimura 1980:
#' DNArate("K80", piGap=0.30, insertionRate=0.02, deletionRate=0.02)
#' DNArate("K80", piGap=0.30, insertionRate=0.02, deletionRate=0.02, par=0)
#' DNArate("K80", piGap=0.30, insertionRate=0.02, deletionRate=0.02, par=2/3)
#' DNArate("K80", piGap=0.30, insertionRate=0.02, deletionRate=0.02, par=1)
#' 
#' ## Kimura 1981:
#' DNArate("K81", piGap=0.30, insertionRate=0.02, deletionRate=0.02)
#' DNArate("K81", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         par=c(0,1/3))
#' DNArate("K81", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         par=c(1/3,0))
#' DNArate("K81", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         par=c(1/3,0))
#' DNArate("K81", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         par=c(0.1,0.4))
#' 
#' ## Felsenstein 1981:
#' DNArate("F81", piGap=0.30, insertionRate=0.02, deletionRate=0.02)
#' DNArate("F81", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         pi=c(0.5,0.15,0.10,0.25))
#' DNArate("F81", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         pi=c(0.15,0.5,0.12,0.23))
#' DNArate("F81", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         pi=c(0.15,0.5,0.05,0.3))
#' 
#' ## Hasegawa, Kishino, and Yano (1985)
#' DNArate("HKY85", piGap=0.30, insertionRate=0.02, deletionRate=0.02)
#' DNArate("HKY85", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         pi=c(0.5,0.15,0.10,0.25))
#' DNArate("HKY85", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         pi=c(0.5,0.15,0.10,0.25), par=0.1)
#' DNArate("HKY85", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         pi=c(0.5,0.15,0.10,0.25), par=0.5)
#' DNArate("HKY85", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         pi=c(0.5,0.15,0.10,0.25), par=0.9)
#' DNArate("HKY85", piGap=0.30, insertionRate=0.02, deletionRate=0.02, par=0.9)
#' 
#' ## Hasegawa, Kishino, and Yano (1985)
#' DNArate("HKY85")
#' DNArate("HKY85", pi = c(0.5, 0.15, 0.10, 0.25))
#' DNArate("HKY85", par = 0.1, pi = c(0.5, 0.15, 0.10, 0.25))
#' DNArate("HKY85", par = 0.5, pi = c(0.5, 0.15, 0.10, 0.25))
#' DNArate("HKY85", par = 0.9, pi = c(0.5, 0.15, 0.10, 0.25))
#' DNArate("HKY85", par = 0.5)
#' 
#' ## Tamura (1992)
#' DNArate("T92", piGap=0.30, insertionRate=0.02, deletionRate=0.02)
#' DNArate("T92", piGap=0.30, insertionRate=0.02, deletionRate=0.02, par=0.1)
#' DNArate("T92", piGap=0.30, insertionRate=0.02, deletionRate=0.02, piGC=0.25)
#' DNArate("T92", piGap=0.30, insertionRate=0.02, deletionRate=0.02, piGC=1/3,
#'         par=0.2)
#' 
#' ## Tamura & Nei (1993)
#' DNArate("TN93", piGap=0.30, insertionRate=0.02, deletionRate=0.02)
#' DNArate("TN93", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         par=c(0.1,0.2))
#' DNArate("TN93", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         par = c(1/3,1/2))
#' 
#' ## Generalized time-reversible (GTR; Tavaré 1986)
#' DNArate("GTR", piGap=0.30, insertionRate=0.02, deletionRate=0.02)
#' DNArate("GTR", piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         par=c(1.5,1,0.5,0,2,1))
#' DNArate("GTR", , piGap=0.30, insertionRate=0.02, deletionRate=0.02,
#'         pi=c(0.5,0.25,0.15,0.10), par=c(1.5,1,0.5,0,2,1))
#' 
#' ## The transition intensity matrix from a Kimura (1980) model:
#' Q <- DNArate(model="K80", piGap = 0.3, deletionRate=0.1, insertionRate=0.1)
#' Q
#' 
#' ## Implement the molecular evolution simulator for a single nucleotide:
#' molEvolSim(Q = Q, step = 1, rho = 5) -> em1
#' 
#' ## Get the transition probability matrix as follows:
#' em1$getMt()
#' 
#' ## A vector of raw as examples of initial traits:
#' tr <- charToRaw("-AGCT")
#' 
#' ## Set the RNG seed.
#' set.seed(28746549L)
#' 
#' ## Simulate molecular evolution from:
#' rawToChar(em1$evolve(tr[1L]))    ## a gap.
#' rawToChar(em1$evolve(tr[2L]))    ## an adenine base.
#' rawToChar(em1$evolve(tr[3L]))    ## a guanine base.
#' rawToChar(em1$evolve(tr[4L]))    ## a cytosine base.
#' rawToChar(em1$evolve(tr[5L]))    ## a thymine base.
#' 
#' ## Recalculate the probabilities for a lower mean evolution rate (one tenth
#' ## the previous one):
#' em1$recalculate(1, 0.1)
#' 
#' em1$getMt()        ## The recalculated transition probability matrix.
#' 
#' ## Simulate molecular evolution from:
#' rawToChar(em1$evolve(tr[1L]))    ## a gap.
#' rawToChar(em1$evolve(tr[2L]))    ## an adenine base.
#' rawToChar(em1$evolve(tr[3L]))    ## a guanine base.
#' rawToChar(em1$evolve(tr[4L]))    ## a cytosine base.
#' rawToChar(em1$evolve(tr[5L]))    ## a thymine base.
#' 
#' ## Base changes are now less probable.
#' 
#' ## Simulate the evolution of a sequence with 100 base pairs for 250 generations.
#' Ngeneration <- 250
#' Nnucleotide <- 100
#' 
#' ## This is the matrix holding the sequences:
#' seq <- matrix(raw(), Ngeneration, Nnucleotide)
#' 
#' ## Drawing the initial sequence:
#' seq[1,] <- drawDNASequence(Nnucleotide, piGap = 0.25, pi = c(0.4,0.4,0.1,0.1))
#' 
#' ## Each site has its own mean evolution rate, which are drawn as follows:
#' erate <- drawEvolRate(Nnucleotide, gamma.shape = 5, gamma.scale = 5e-03)
#' 
#' ## Using the Hasegawa, Kishino, and Yano (1985) model:
#' DNArate(
#'   model = "HKY85",
#'   piGap = 0.25,
#'   deletionRate = 0.1,
#'   insertionRate = 0.1,
#'   pi = c(0.4, 0.4, 0.1, 0.1),
#'   par = 0.25
#' ) -> Q
#' 
#' ## Instantiating a molecular evolution models for each site using the single
#' ## change rate matrix, a constant time step of 1 (Ma), and individual mean
#' ## evolution rates drawn previously:
#' em <- list()
#' for(j in 1:Nnucleotide)
#'   em[[j]] <- molEvolSim(Q = Q, step = 1, rho = erate[j])
#' 
#' ## These for loops call the $evolve() function to evolve each site for each
#' ## generation as follows:
#' for(i in 2:Ngeneration)
#'   for(j in 1:Nnucleotide)
#'     seq[i, j] <- em[[j]]$evolve(seq[i - 1, j])
#' 
#' ## Sequences with the gaps (perfect alignment):
#' concatenate(seq) %>%
#'   show.sequence
#' 
#' ## Sequences with the gaps removed (prior to multiple sequence alignment):
#' concatenate(seq, discard="-") %>%
#'   show.sequence
#' 
#' ### Examples for molEvolSim() here...
#' 
#' ## Clean up:
#' rm(Q, em1, tr, Ngeneration, Nnucleotide, seq, erate, em, i, j)
#' 
NULL
#' 
#' @describeIn molEvolSim
#' 
#' @title DNA Evolution Shift Rate
#' 
#' @description Calculates the shift rate matrix associated with one of eight
#' DNA evolution models.
#' 
#' @export
DNArate <- function(
    model = c("JC69","K80","K81","F81","HKY85","T92","TN93","GTR"), piGap = 0,
    deletionRate = 0, insertionRate = 0, pi, piGC, par) {
  
  model <- match.arg(model)
  
  matrix(
    0,
    5L, 5L,
    dimnames = list(
      c("-","A","G","C","T"),
      c("-","A","G","C","T")
    )
  ) -> Q
  
  rate <- c(alpha=NA, beta=NA, delta=NA, gamma=NA, epsilon=NA, eta=NA)
  
  if((piGap < 0) || (piGap > 1))
    stop("Argument 'piGap' must be between 0 and 1")
  
  if((deletionRate < 0) || (deletionRate > 1))
    stop("Argument 'deletionRate' must be between 0 and 1")
  
  if((insertionRate < 0) || (insertionRate > 1))
    stop("Argument 'insertionRate' must be between 0 and 1")
  
  if(missing(pi)) {
    
    pi <- rep(0.25, 4L)
    
  } else {
    
    if(any(pi < 0) && any(pi > 1))
      stop("The value of the proportion parameters (argument 'prop') must be ",
           "greater then 0 and no greater than 1 (see documentation).")
    
    if(length(pi) != 4L) {
      warning("Argument 'prop' had to be recycled/padded to lenght 4!")
      pi <- rep(pi, length.out=4L)
    }
    
    if(sum(pi) != 1) {
      warning("The sum of the proportion parameters given through argument ",
              "'pi', sum to ", sum(pi), " rather than 1, and were made to sum ",
              "to 1 (see documentation).")
      pi <- pi/sum(pi)
    }
    
  }
  
  if(missing(piGC)) piGC <- 0.5 else if((piGC < 0) || (piGC > 1))
    stop("The value of the GC content bias parameter (argument 'piGC') must ",
          "be between 0 and 1 (see documentation).")
  
  if(model == "JC69") {
    
    rate[] <- 1                 ## Equal evolution rates for all changes.
    pi[] <- 0.25                ## Force equal equilibrium base frequencies.
    
  }
  
  if(model == "K80") {
    
    if(missing(par)) par <- 1/3 else {
      if((par[1L] < 0) || (par[1L] > 1))
        stop("The value of transition-transversion rate disparity parameter ",
             "(argument 'par[1]') must be between 0 and 1 (see documentation).")
    }
    
    rate[c(1L,6L)] <- 3*par[1L]        ## Any transversions
    rate[2L:5L] <- 1.5*(1 - par[1L])   ## Any transitions
    
    pi[] <- 0.25            ## Force equal equilibrium base frequencies.
    
  }
  
  if(model == "K81") {
    
    if(missing(par)) par <- rep(1/3, 2L) else {
      if(length(par) != 2L)
        stop("When present, argument 'par' must be of length 2.")
      if((par[1L] < 0) || (par[2L] < 0) || ((par[1L] + par[2L]) > 1))
        stop("The values provided to both parameters 'alpha' (argument ",
              "par[1]) and 'beta' (argument par[2] must be between 0 and 1 ",
             " and their sums must be <= 1 (see documentation).")
    }
    
    rate[c(1L,6L)] <- 3*par[1L]                 ## transversions
    rate[c(2L,5L)] <- 3*par[2L]                 ## SW transitions, conserves AK
    rate[c(3L,4L)] <- 3*(1 - par[1L] - par[2L]) ## AK transitions, conserves SW
    
    pi[] <- 0.25            ## Force equal equilibrium base frequencies.
    
  }
  
  if(model == "F81") {
    
    rate[] <- 1          ## Equal evolution rates for all transitions.
    ## Uses the provided equilibrium base frequencies.
    
  }
  
  if(model == "HKY85") {
    
    if(missing(par)) par <- 1/3 else {
      if((par[1L] < 0) || (par[1L] > 1))
        stop("The value of transition-transversion rate disparity parameter ",
             "(argument 'par[1]') must be between 0 and 1 (see documentation).")
    }
    
    rate[c(1L,6L)] <- 3*par[1L]        ## Any transversions
    rate[2L:5L] <- 1.5*(1 - par[1L])   ## Any transitions
    ## Uses the provided equilibrium base frequencies.
    
  }
  
  if(model == "T92") {
    
    if(missing(par)) par <- 1/3 else {
      if((par[1L] < 0) || (par[1L] > 1))
        stop("The value of transition-transversion rate disparity parameter ",
             "(argument 'par[1]') must be between 0 and 1 (see documentation).")
    }
    
    rate[c(1L,6L)] <- 3*par[1L]       ## Any transversions
    rate[2L:5L] <- 1.5*(1 - par[1L])  ## Any transitions
    pi[c(1L,4L)] <- 0.5*(1 - piGC)    ## A and T frequencies
    pi[c(2L,3L)] <- 0.5*piGC          ## G and C frequencies
    
  }
  
  if(model == "TN93") {
    
    if(missing(par)) par <- rep(1/3, 2L) else {
      if(length(par) != 2L)
        stop("When present, argument 'par' must be of length 2.")
      if(any(par < 0) && any(par > 1))
        stop("The value of any of the two transition-transversion rate ",
             "disparity parameters (argument 'K[1]' and 'K[2]') must be ",
             "between 0 and 1 (see documentation).")
    }
    
    rate[c(1L,6L)] <- 3*par[1L:2L]
    rate[2L:5L] <- 1 - 0.5*par[1L] + 0.5*par[2L]
    ## Uses the provided equilibrium base frequencies.
    
  }
  
  if(model == "GTR") {
    
    if(missing(par)) par <- rep(1, 6L) else {
      
      if(length(par) < 6L)
        stop("When present, argument 'par' must be of length 6.")
      
      if(sum(par) != 6) {
        warning("The parameters of the generalized time-reversible model given ",
                "through argument 'par', sum to ", sum(par), " rather than ",
                "6, and were made to sum to 6 (see documentation).")
        par <- 6*par/sum(par)
      }
      
    }
    
    rate[] <- par
    ## Uses the provided equilibrium base frequencies.
    
  }
  
  ## Gap shift rate used to model insertions and deletions:
  piN <- 1 - piGap
  Q[1L,-1L] <- insertionRate*pi*piN
  Q[1L, 1L] <- -sum(Q[1L,-1L])
  Q[-1L,1L] <- piGap*deletionRate
  
  ##  Calculation of the transition intensities:
  
  c(-(sum(rate[c(1L,2L,4L)]*piN*pi[c(2L,3L,4L)]) + Q[2L,1L]),
    rate[c(1L,2L,4L)]*piN*pi[c(2L,3L,4L)]) -> Q[2L,-1L]       ## From A
  
  c(rate[c(1L)]*piN*pi[c(1L)],
    -(sum(rate[c(1L,3L,5L)]*piN*pi[c(1L,3L,4L)]) + Q[3L,1L]),
    rate[c(3L,5L)]*piN*pi[c(3L,4L)]) -> Q[3L,-1L]             ## From G
  
  c(rate[c(2L,3L)]*piN*pi[c(1L,2L)],
    -(sum(rate[c(2L,3L,6L)]*piN*pi[c(1L,2L,4L)]) + Q[4L,1L]),
    rate[c(6L)]*piN*pi[c(4L)]) -> Q[4L,-1L]                   ## From C
  
  c(rate[c(4L,5L,6L)]*piN*pi[c(1L,2L,3L)],
    -(sum(rate[c(4L,5L,6L)]*piN*pi[c(1L,2L,3L)]) + Q[5L,1L])
  ) -> Q[5L,-1L]                                              ## From T
  
  Q
}
#' 
#' @describeIn molEvolSim
#' 
#' @title DNA Molecular Evolution Simulator
#' 
#' @description Instantiate a DNA (molecular) evolution simulator from a shift
#' intensity matrix over a given evolutionary time.
#' 
#' @export
molEvolSim <- function(Q, step, rho) {
  Mt <- expm(step*rho*Q)
  nbytes <- charToRaw(paste(rownames(Q),collapse=""))
  list(
    recalculate = function(step, rho)
      Mt <<- expm(step*rho*Q),
    evolve = function(N)
      sample(x=nbytes, size=1L, prob=Mt[which(N == nbytes),]),
    getMt = function() Mt
  )
}
#' 
#' @describeIn molEvolSim
#' 
#' @title Random Sequence Generator
#' 
#' @description Generates a random sequence of gaps or DNA bases.
#' 
#' @export
drawDNASequence <- function(NN, piGap = 0, pi = rep(0.25,4L)) {
  
  if((piGap < 0) || (piGap > 1))
    stop("Argument 'piGap' must be between 0 and 1")
  
  if(any(pi < 0) && any(pi > 1))
    stop("The value of the proportion parameters (argument 'prop') must be ",
         "greater then 0 and no greater than 1 (see documentation).")
  
  if(length(pi) != 4L) {
    warning("Argument 'prop' had to be recycled/padded to lenght 4!")
    pi <- rep(pi, length.out=4L)
  }
  
  if(sum(pi) != 1) {
    warning("The sum of the proportion parameters given through argument ",
            "'pi', sum to ", sum(pi), " rather than 1, and were made to sum ",
            "to 1 (see documentation).")
    pi <- pi/sum(pi)
  }
  
  sample(
    x = charToRaw("-AGCT"),
    size = NN,
    prob = c(piGap,(1 - piGap)*pi),
    replace = TRUE
  )
}
#' 
#' @describeIn molEvolSim
#' 
#' @title Random Evolution Rate Generator
#' 
#' @description Generates a set of random nucleotide evolution rate from a gamma
#' distribution.
#' 
#' @export
drawEvolRate <- function(NN, gamma.shape = 5, gamma.scale = 5e-04)
  rgamma(
    NN,
    shape = gamma.shape,
    scale = gamma.scale
  )
#' 
#' @describeIn molEvolSim
#' 
#' @title Sequence Evolution Simulator
#' 
#' @description Generates a set of DNA sequences along the edge of a
#' phylogenetic network by evolving filial sequences from parental ones
#' following a random Markov process.
#' 
#' @export
simulateSequence <- function(
    x, Q, sqn, rate,
    contrib = function(x, a = 0, ...) (x^-a)/sum(x^-a), ...
  ) {
  
  nv <- nrow(x)
  
  origin <- getOrigin(x)
  norig <- length(origin)
  
  if(norig != NCOL(sqn))
    stop("The graph has ",norig," origins, but 'seq' as ",NCOL(sqn),
         " sequences.")
  
  NN <- NROW(sqn)
  
  if(NN != length(rate)) {
    
    warning(
      "The number of 'rate' values (",length(rate),
      ") did not match the length ",
      "of the DNA sequences (",NN,") and had to be adjusted."
    )
    
    rate <- rep(rate, length.out=NN)
  }
  
  timestep <- 1
  
  em <- list()
  
  for(i in 1L:NN)
    em[[i]] <- molEvolSim(Q, timestep, rate[i])
  
  out <- matrix(raw(), nv, NN, dimnames = list(rownames(x), names(sqn)))
  
  out[origin,] <- sqn
  
  ord <- attr(x,"processOrder")
  if(is.null(ord))
    ord <- getProcessOrder(x)
  
  edge <- attr(x, "edge")
  
  for(k in ord[(norig + 1L):nv]) {
    
    av <- edge[[1L]][edge[[2L]] == k]
    da <- edge$distance[edge[[2L]] == k]
    ctb <- contrib(da, ...)
    s <- sample(length(ctb), 1L, FALSE, ctb)
    i <- av[s]
    
    if(da[s] != timestep) {
      
      timestep <- da[s]
      
      for(j in 1L:NN)
        em[[j]]$recalculate(timestep, rate[j])
    }
    
    for(j in 1L:NN)
      out[k,j] <- em[[j]]$evolve(out[i,j])
  }
  
  out
}
##
