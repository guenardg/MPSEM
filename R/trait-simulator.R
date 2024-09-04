## **************************************************************************
##
##    (c) 2010-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Trait value simulator**
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
#' Simulate the Evolution of a Quantitative Trait
#' 
#' @description Functions to simulate the evolution of a quantitative trait
#' along a phylogenetic tree inputted as an object of class \sQuote{phylo}
#' (package \link[ape]{ape}) or a \code{\link{graph-class}} object.
#' 
#' @name trait-simulator
#' 
#' @aliases trait simulation
#' 
#' @param tp A rooted phylogenetic tree of class \sQuote{phylo} (see package
#' \link[ape]{ape}).
#' @param tw Transition matrix giving the probability that the optimum trait
#' value changes from one state (row) to another (column) at vertices. 
#' All rows must sum to 1.
#' @param anc Ancestral state of a trait (at the root).
#' @param p Number of variates to generate.
#' @param root Root node of the tree.
#' @param d Phylogenetic distances (edge lengths).
#' @param a Selection rate in function (\code{\link{OUvar}}) or steepness in
#' (\code{\link{PEMvar}}).
#' @param theta Adaptive evolution rate, i.e. mean trait shift by natural
#' selection.
#' @param sigma Neutral evolution rate, i.e. mean trait shift by drift.
#' @param psi Mean evolution rate.
#' @param opt An index vector of optima at the nodes.
#' @param x A \code{\link{graph-class}} object.
#' @param variance Variance function: \code{\link{OUvar}}, \code{\link{PEMvar}},
#' or any other suitable user-defined function.
#' @param distance The name of the member of \sQuote{x$edge} where edge lengths
#' can be found.
#' @param ... Additional parameters for the specified variance function.
#' 
#' @details Function \code{EvolveOptimMarkovTree} allows one to simulate the
#' changes of optimum trait values as a Markov process. The index whereby the
#' process starts, at the tree root, is set by parameter \code{anc}; this is the
#' ancestral character state. From the root onwards to the tips, the optimum is
#' given the opportunity to change following a multinomial random draw with
#' transition probabilities given by the rows of matrix \code{tw}. The integers
#' thus obtained can be used as indices of a vector featuring the actual optimum
#' trait values corresponding to the simulated selection regimes. 
#'
#' The resulting
#' optimum trait values at the nodes are used by \code{\link{TraitOUsimTree}} as
#' its argument \code{opt} to simulate trait values at nodes and tips.
#' 
#' Function \code{\link{TraitVarGraphSim}} uses a graph variance function
#' (either \code{OUvar} or \code{PEMvar}) to reconstruct a covariance matrix,
#' used to generate covariates drawn from a multi-normal distribution.
#' 
#' @return Functions \code{\link{EvolveOptimMarkovTree}} and
#' \code{\link{TraitOUsimTree}} return a matrix whose rows represent the
#' vertices (nodes and tips) of the phylogenetic tree and whose columns stand
#' for the \code{n} different trials the function was asked to perform. 
#' 
#' For \code{EvolveQTraitTree}, the elements of the matrix are integers,
#' representing the selection regimes prevailing at the nodes and tips, whereas
#' for \code{\link{TraitOUsimTree}}, the elements are simulated quantitative
#' trait values at the nodes and tips. These functions are implemented in C
#' language and therefore run swiftly even for large (10000+ species) trees.
#' 
#' Function \code{\link{TraitVarGraphSim}} returns \code{p} phylogenetic signals.
#' It is implemented using a rotation of a matrix of standard normal random
#' (mean=0, variance=1) deviates. The rotation matrix is itself obtained by
#' Choleski factorization of the trait covariance matrix expected for a given
#' set of trees, variance function, and variance function parameters.
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Butler, M. A. & King, A. A. 2004. Phylogenetic comparative analysis: a
#' modeling approach for adaptive evolution. American Naturalist 164: 683-695
#' 
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits.  Methods in Ecology 
#' and Evolution 4: 1120-1131
#' 
#' @importFrom stats rnorm
#' 
#' @examples
#' opt <- c(-2,0,2) # Three trait optima: -2, 0, and 2
#' ## Transition probabilities:
#' transit <- matrix(c(0.7,0.2,0.2,0.2,0.7,0.1,0.1,0.1,0.7),
#'                   length(opt),length(opt),dimnames=list(from=opt,to=opt))
#' 
#' ## In this example, the trait has a probability of 0.7 to stay at a given
#' ## optimum, a probability of 0.2 for the optimum to change from -2 to 0,
#' ## from 0 to -2, and from 2 to -2, and a probability of 0.1 for the
#' ## optimum to change from -2 to 2, from 0 to 2, and from 2 to 0.
#' nsp <- 25  # A random tree for 25 species.
#' tree2 <- rtree(nsp,tip.label=paste("Species",1:nsp,sep=""))
#' tree2$node.label=paste("N",1:tree2$Nnode,sep="")  # Node labels.
#' 
#' ## Simulate 10 trials of optimum change.
#' reg <- EvolveOptimMarkovTree(tp=tree2,tw=transit,p=10,anc=2)
#' y1 <- TraitOUsimTree(tp=tree2,a=0,sigma=1,
#'                      opt=opt[reg[,1]],p=10)    ## Neutral
#' y2 <- TraitOUsimTree(tp=tree2,a=1,sigma=1,
#'                      opt=opt[reg[,1]],p=10)    ## Few selection.
#' y3 <- TraitOUsimTree(tp=tree2,a=10,sigma=1,
#'                      opt=opt[reg[,1]],p=10)    ## Strong selection.
#' 
#' ## Display optimum change with colours.
#' displayOUprocess <- function(tp,trait,regime,mvalue) {
#'   layout(matrix(1:2,1,2))
#'   n <- length(tp$tip.label)
#'   ape::plot.phylo(tp,show.tip.label=TRUE,show.node.label=TRUE,root.edge=FALSE,
#'                   direction="rightwards",adj=0,
#'                   edge.color=rainbow(length(trait))[regime[tp$edge[,2]]])
#'   plot(y=1:n,x=mvalue[1:n],type="b",xlim=c(-5,5),ylab="",xlab="Trait value",yaxt="n",
#'        bg=rainbow(length(trait))[regime[1:n]],pch=21) 
#'   text(trait[regime[1:n]],y=1:n,x=5,col=rainbow(length(trait))[regime[1:n]])
#'   abline(v=0)
#' }
#' 
#' displayOUprocess(tree2,opt,reg[,1],y1[,1])  # Trait evolve neutrally,
#' displayOUprocess(tree2,opt,reg[,1],y2[,1])  # under weak selection,
#' displayOUprocess(tree2,opt,reg[,1],y3[,1])  # under strong selection.
#' 
#' x <- Phylo2DirectedGraph(tree2)
#' y4 <- TraitVarGraphSim(x, variance = OUvar, p=10, a=5)
#' 
#' DisplayTreeEvol <- function(tp,mvalue) {
#'   layout(matrix(1:2,1,2))
#'   n <- length(tp$tip.label)
#'   ape::plot.phylo(tp,show.tip.label = TRUE, show.node.label = TRUE,
#'                   root.edge = FALSE, direction = "rightwards", adj = 0)
#'   plot(y=1:n, x=mvalue[1:n], type="b", xlim=c(-5,5), ylab="",
#'        xlab="Trait value", yaxt="n", pch=21)
#'   abline(v=0)
#' }
#' 
#' ## Iteratively displays the simulated traits.
#' ## Left-click on the display area to go to the next plot.
#' ## To terminate: right-click (WIndows, X11), esc key (Mac), or hit the
#' ## "finish" button (RStudio).
#' 
#' for(i in 1:10) {
#'   DisplayTreeEvol(tree2,y4[i,])
#'   if(is.null(locator(1)))
#'     break  ## Terminate: 
#' }
#' 
#' @useDynLib MPSEM, .registration = TRUE
#' 
NULL
#' 
#' @describeIn trait-simulator
#' 
#' Trait Optima Simulator
#' 
#' Simulates the evolution of trait optima along a phylogeny as a Markov
#' process.
#' 
#' @export
EvolveOptimMarkovTree <- function(tp, tw, anc, p=1, root=tp$edge[1,1]) {
  nn <- length(tp$tip.label)+tp$Nnode
  if(nrow(tw) != ncol(tw))
    stop("Transition probability matrix (tw) must be a square matrix")
  if(anc > nrow(tw))
    stop("Ancestral state (anc) not defined in the transition probability matrix (tw).")
  if(any((abs(rowSums(tw)-1)) > sqrt(.Machine$double.eps)))
    warning("The sum of transition probabilities is not systematically 1.")
  if(root > nn)
    stop("Invalid parameter root.")
  res <- t(matrix(.C("EvolveQC",
                     as.integer(tp$edge[,1]),
                     as.integer(tp$edge[,2]),
                     as.integer(nrow(tp$edge)),
                     as.integer(nn),
                     nv = double(p*nn),
                     as.double(t(tw)),
                     as.integer(nrow(tw)),
                     as.integer(anc),
                     as.integer(p),
                     as.integer(root))$nv,p,nn))
  if(!is.null(tp$node.label)) {
    rownames(res) <- c(tp$tip.label,tp$node.label)
  } else {
    rownames(res) <- c(tp$tip.label,rep("",tp$Nnode))
  }
  colnames(res) <- paste("Trial",1:p,sep="_")
  return(res)
}
#' 
#' @describeIn trait-simulator
#' 
#' Trait Value Simulator
#' 
#' Simulates the evolution of trait values along a phylogeny as a
#' Ornstein–Uhlenbeck process.
#' 
#' @export
TraitOUsimTree <- function(tp, a, sigma, opt, p=1, root=tp$edge[1,1]) {
  nn <- length(tp$tip.label)+tp$Nnode
  if(root > nn)
    stop("Invalid parameter root.")
  if(length(opt) != nn)
    stop("Optima don't match the number of nodes.")
  res <- matrix(.C("OUsim",
                   as.integer(tp$edge[,1]),
                   as.integer(tp$edge[,2]),
                   as.integer(nrow(tp$edge)),
                   as.integer(nn),
                   as.double(tp$edge.length),
                   as.double(a[1]),
                   as.double(sigma[1]),
                   as.double(opt),
                   as.integer(p),
                   as.integer(root),
                   out=double(p*nn))$out,nn,p)
  if(!is.null(tp$node.label)) {
    rownames(res) <- c(tp$tip.label,tp$node.label)
  } else {
    rownames(res) <- c(tp$tip.label,rep("",tp$Nnode))
  }
  colnames(res) <- paste("Trial",1:p,sep="_")
  return(res)
}
#' 
#' @describeIn trait-simulator
#' 
#' Ornstein–Uhlenbeck Variance Calculator
#' 
#' Calculates the expected covariance matrix for a trait evolving following an
#' Ornstein–Uhlenbeck process. This function is meant to be used with function
#' \code{TraitVarGraphSim}.
#' 
#' @export
OUvar <- function(d, a=0, theta=1, sigma=1) {
  nd <- length(d)
  w <- numeric(nd)
  a <- rep(a, length.out = nd)
  theta <- rep(theta, length.out = nd)
  sigma <- rep(sigma, length.out = nd)
  nz <- a != 0
  w[nz] <- (theta[nz]*(1-exp(-a[nz]*d[nz])))**2 +
    (sigma[nz]**2)*(1-exp(-2*a[nz]*d[nz]))/(2*a[nz])
  w[!nz] <- (sigma[!nz]**2)*d[!nz]
  return(w)
}
#' 
#' @describeIn trait-simulator
#' 
#' Phylogenetic Eigenvector Maps Variance Calculator
#' 
#' Calculates the covariance on the basis of the covariance model (power
#' function) associated used in calculating Phylogenetic Eigenvector Maps. This
#' function is meant to be used with function \code{TraitVarGraphSim}.
#' 
#' @export
PEMvar <- function(d, a=0, psi=1) {
  nd <- length(d)
  a <- rep(a, length.out = nd)
  psi <- rep(psi, length.out = nd)
  return(.C("PEMvarC",
            as.double(d),
            as.integer(nd),
            as.double(a),
            as.double(psi),
            res=double(nd))$res)
}
#' 
#' @describeIn trait-simulator
#' 
#' Covariance-based Trait Evolution Simulator.
#' 
#' Simulates trait evolution as covariates drawn from a multi-normal
#' distribution whose covariance is estimated using an external function
#' (functions \code{OUvar}, \code{PEMvar} provided with the package or any
#' user-provided function).
#' 
#' @export
TraitVarGraphSim <- function(x, variance, distance="distance", p=1, ...) {
  if(attr(x,"class") != "graph")
    stop("Parameter 'x' must be of class 'graph'")
  if(is.null(x$edge[[distance]]))
    stop("There is no property '",distance,"' for the edges of the graph.")
  if(is.null(x$vertex$species)) {
    B <- PEMInfluence(x,mroot = FALSE)
  } else {
    B <- PEMInfluence(x,mroot = FALSE)[x$vertex$species,]
  }
  n <- nrow(B)
  if(!missing(variance)) {
    fargs <- as.list(match.call())
    wfform <- formals(variance)
    wfform[[1L]] <- x$edge[[distance]]
    for(i in names(wfform)[-1L]) {
      if (!is.null(fargs[[i]])) {
        wfform[[i]] <- fargs[[i]]
      } else if(!is.null(x$edge[[i]])) {
        wfform[[i]] <- x$edge[[i]]
      }
    }
    return(matrix(rnorm(p*n,0,1),p,n)%*%chol(B%*%diag(do.call(what=variance,args=as.list(wfform)))%*%t(B)))
  } else {
    return(matrix(rnorm(p*n,0,1),p,n)%*%chol(B%*%t(B)))
  }
}
##
