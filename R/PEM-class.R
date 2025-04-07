## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Class and method for Phylogenetic Eigenvector maps **
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
#' Class and Methods for Phylogenetic Eigenvector Maps (PEM)
#' 
#' @description Class and methods to calculate and manipulate Phylogenetic
#' Eigenvector Maps (PEM).
#' 
#' @docType class
#' 
#' @name PEM-class
#'  
#' @param x A \code{\link{graph-class}} object, a \code{\link{PEM-class}}
#' object, or a matrix containing auxiliary traits.
#' @param a A numeric vector of linear coefficients to estimate the steepness in
#' the trait evolution model; it must have as many values as the number of
#' columns in argument \code{mm_a} (default: \code{-Inf}).
#' @param mm_a A numeric matrix with as many rows as the number of
#' edges and as many columns as the number of values of argument \code{a}
#' (default: an all-ones single-column matrix).
#' @param psi A numeric vector of linear coefficients to estimate the evolution
#' rate in the trait evolution model; it must have as many values as the number
#' of columns in argument \code{mm_psi} (default: \code{NULL}, meaning an
#' homogeneous evolution rate throughout the graph).
#' @param mm_psi A numeric matrix with as many rows as the number of
#' edges and as many columns as the number of values of argument \code{psi}
#' (default: a zero-column matrix).
#' @param d A character string; the name of the (numeric) edge property
#' containing the phylogenetic distances (default: \code{"distance"}).
#' @param sp A character string; the name of the (logical) vertex property
#' containing vertex property (default: \code{"species"}).
#' @param tol A numeric; the lowest singular value above which to retain a
#' singular vector as part of a PEM (default: \code{.Machine$double.eps^0.5}).
#' @param row.names Included for method consistency reason; ignored.
#' @param optional Included for method consistency reason; ignored.
#' @param object A \code{\link{PEM-class}} object.
#' @param y A vector or matrix of (numeric) trait values to estimate the
#' parameters of the trait evolution model.
#' @param newdata A three-column data frame giving the locations of a set of
#' targets in the graph such as the one provided by function
#' \code{\link{locate.graph}} (element \code{$location}; see Details below).
#' @param ... Additional parameters to be passed to the method. Currently
#' ignored.
#' 
#' @details The class is implemented using function \code{PEM}. It provides
#' methods \code{\link{as.matrix}} and \code{\link{as.data.frame}} to extract
#' the phylogenetic vectors in the forms of a matrix or a data frame,
#' respectively, an \code{update} method to recalculate an existing PEM with
#' new parameters (ie., arguments \code{a} and \code{psi}), an
#' \code{\link{evolution.model}} method to empirically estimate parameter values
#' from a data set of traits (ie., argument \code{y} and \code{x}), and a
#' \code{\link{predict}} method to estimate PEM scores for arbitrary graph
#' locations (using argument \code{newdata}).
#' 
#' Graph locations are given as a three-column data frame:
#' \describe{
#'   \item{ref}{the index of the vertex or edge where each target can be found,}
#'   \item{dist}{the phylogenetic distance along the edge where each target can
#'     be found, or \code{NA} for a target located at a vertex,}
#'   \item{lca}{the phylogenetic distance between the latest common ancestor of
#'     the target in the graph and the target itself.}
#' }
#' One such data frame is one of the second element of the list produced by
#' method \code{\link{locate.graph}} (called "location").
#' 
#' @format A \code{\link{PEM-class}} object contains:
#' \describe{
#'   \item{d}{the name of the edge property containing the phylogenetic
#'     distances (see argument \code{d} above),}
#'   \item{sp}{the name of the vertex property identifying which which of the
#'     vertices are bearing trait values (see argument \code{sp} above),}
#'   \item{tol}{the tolerance value (see argument \code{tol} above),}
#'   \item{mm}{a list with the model matrices (see arguments \code{mm_a} and
#'     \code{mm_psi} above) that are used to estimate the steepness and relative
#'     evolution rate throughout the graph manifold,}
#'   \item{nsp}{the number of vertices that bear trait values,}
#'   \item{B}{the influence matrix (see \code{\link{inflMat}}),}
#'   \item{print}{function; prints the PEM object,}
#'   \item{graph}{function; returns the graph object used to calculate the PEM,}
#'   \item{pem}{function; returns a list of components (ie., vectors and
#'     matrices) involved in PEM calculation,}
#'   \item{par}{function; returns the actual model parameter values as a list,}
#'   \item{nodal}{function; returns the steepness and evolution rates at any
#'     given vertex of the graph,}
#'   \item{update}{function; updates the PEM with a new set of parameters,}
#'   \item{var}{function; returns the PVCV matrix: phylogenetic variance and
#'     covariances matrix among the vertices with trait values,}
#'   \item{inv_var}{function; returns the inverse PVCV matrix,}
#'   \item{logdet}{function; returns the natural logarithm of the PVCV matrix
#'     determinant,}
#'   \item{dev}{function; calculates the deviance of a target trait, or many
#'     such traits, with respect to the PVCV matrix, optionally conditional to
#'     one or more auxiliary trait(s),}
#'   \item{S2}{function; calculates the variance(s) of one or more trait(s)
#'     with respect to the PVCV matrix or the residual variance(s) of said
#'     trait(s) conditioned to one or more auxiliary trait(s),}
#'   \item{predictVertex}{function; calculates the PEM prediction scores at the
#'     graph's vertices, and}
#'   \item{predictEdge}{function; calculates the PEM prediction scores at the
#'     graph's edge.}
#' }
#' 
#' @examples
#' ## Synthetic example
#' 
#' ## This example describes the phyogeny of 7 species (A to G) in a tree with 6
#' ## nodes, presented in Newick format, read by function
#' ## read.tree of package ape.
#' 
#' t1 <- read.tree(text=paste(
#'   "(((A:0.15,B:0.2)N4:0.15,C:0.35)N2:0.25,((D:0.25,E:0.1)N5:0.3,",
#'   "(F:0.15,G:0.2)N6:0.3)N3:0.1)N1:0.1;",sep=""))
#' t1
#' summary(t1)
#' 
#' ## Example of a made up set of trait values for the 7 species
#' y <- c(A=-1.1436265,B=-0.3186166,C=1.9364105,D=1.7164079,E=1.0013993,
#'        F=-1.8586351,G=-2.0236371)
#' 
#' ## The phylogenetic tree is turned into an evolutionary graph as follows:
#' x <- as.graph(t1)
#' 
#' ## The graph:
#' x
#' 
#' ## This is the edge information of that simple graph:
#' edge(x)
#' 
#' ## Calculate the (binary) influence matrix; E1 to E12 are the tree edges:
#' IM1 <- InflMat(x)
#' IM1
#' 
#' ## This is the edges associated with the influence matrix:
#' edge(IM1)
#' 
#' ## This is the influence matrix for the tree leaves (ie., the graph's 
#' ## terminal vertices):
#' IM1[x$species,]
#' 
#' ## Here is a suite of weighting function profiles:
#' seq(0,1.5,0.001) %>%
#'   plot(y=PEMweights(., a=0), ylim=c(0,1.7), type="l", xlab="distance",
#'        ylab="weight")
#' 
#' seq(0,1.5,0.001) %>%
#'   lines(y=PEMweights(., a=0.5), col="red")
#' 
#' seq(0,1.5,0.001) %>%
#'   lines(y=PEMweights(., a=0.5, psi=1.5), col="green")
#' 
#' seq(0,1.5,0.001) %>%
#'   lines(y=PEMweights(., a=0.9), col="blue")
#' 
#' ## A Phylogenetic eigenvector maps (PEM) is build with the default parameters
#' ## as follows:
#' PEM1 <- PEM(x)
#' 
#' ## The PEM object:
#' PEM1
#' 
#' ## The as.matrix methods returns the phylogenetic eigenvectors as a matrix:
#' as.matrix(PEM1)
#' 
#' ## The as.data.frame methods returns the phylogenetic eigenvectors as a data
#' ## frame:
#' as.data.frame(PEM1)
#' 
#' ## The update method recalculates the PEM with new parameters:
#' update(PEM1, a=log(0.25/(1 - 0.25)), psi=NULL)
#' 
#' ## Note: The way 'a' is estimated is based on an inverse-logit-transformed
#' ## linear sub-model. This approach is used to facilitate the bounding of the
#' ## steepness parameter in the [0,1] interval. Therefore, to set a given
#' ## steepness value, one has to provide its logit-transformed value. As a
#' ## reminder, the logit transformation is f(x) = log(x/(1 - x)).
#' 
#' ## To show the edges of the graph within the PEM object:
#' edge(PEM1$graph())
#' 
#' ## A second PEM is initialized with a = 0.2 as follows:
#' PEM2 <- PEM(x, a=log(0.2/(1 - 0.2)))
#' PEM2
#' 
#' ## The edges of PEM2 graph:
#' edge(PEM2$graph())
#' 
#' ## A third PEM is initialized with a = 1 as follows:
#' PEM3 <- PEM(x, a=Inf)
#' PEM3
#' 
#' ## The edges of PEM3 graph:
#' edge(PEM3$graph())
#' 
#' ## Estimation on the steepness parameter empirically using data set y
#' ## Initial value must be finite (not -Inf nor Inf, meaning a > 0 and a < 1):
#' PEM4 <- PEM(x, a=0)
#' 
#' ## Method evolution.model carry out the estimation:
#' opt <- evolution.model(PEM4, y=y)
#' opt
#' 
#' ## The edges of PEM4 graph after the estimation:
#' edge(PEM4$graph())
#' 
#' ## Graph locations for target species X, Y, and Z not found in the original
#' ## data set:
#' read.tree(
#'   text = paste(
#'     "((X:0.45,((A:0.15,B:0.2)N4:0.15,(C:0.25,Z:0.2)NZ:0.1)N2:0.05)NX:0.2,",
#'     "(((D:0.25,E:0.1)N5:0.05,Y:0.25)NY:0.25,(F:0.15,G:0.2)N6:0.3)N3:0.1)N1;",
#'     sep=""
#'   )
#' ) -> tree
#' tree
#' 
#' ## The tree is converted into a graph as follows:
#' x2 <- as.graph(tree)
#' x2
#' 
#' ## The target X, Y, and Z are extracted from graph x2 as follows:
#' loc <- locate(x2, target = c("X","Y","Z"))
#' loc
#' 
#' ## The list contains the residual graph and a data frame with the locations
#' ## of the target in the residual graph.
#' 
#' ## Building a PEM on the residual graph (assuming a = 0):
#' PEM5 <- PEM(loc$x)
#' 
#' ## For making predictions, we need to do any of the following:
#' ## 1. run evolution.model(); doing this operation will change the parameter
#' ##    values,
#' ## 2  estimate the deviance using $dev(y, ...)
#' ## 3. estimate the variance using $S2(y, ...)
#' 
#' ## Before any of these operations is done, no variance estimate is available:
#' PEM5$S2()
#' 
#' ## Presenting a variable to the PEM makes its variance estimate available
#' ## later on. For instance, presenting 'y' to the PEM5 while estimating its
#' ## deviance as follows:
#' PEM5$dev(y)
#' 
#' ## This operation makes the variance estimate for 'y' available for later
#' ## access:
#' PEM5$S2()
#' 
#' ## The predict method is used to generate the predictors as follows:
#' scr <- predict(PEM5, loc$location)
#' scr
#' 
#' ## The 'vf' attribute, which can be accessed as follows:
#' attr(scr,"vf")
#' ## is the amount of variance associated with the evolutionary distance
#' ## between each target species and its latest common ancestor in the
#' ## evolutionary graph.
#' 
#' ## Given a model built using the PEM as follows:
#' lm1 <- lm(y ~ U_2 + U_3 + U_5, data=as.data.frame(PEM5))
#' 
#' ## Predictions are obtained as follows:
#' ypred <- predict(lm1, as.data.frame(scr))
#' ypred
#' 
#' ## Estimate the ancestral trait values
#' 
#' ## The ancestors of the species of the tree are the vertices without extant
#' ## species:
#' data.frame(
#'   ref = which(!x$species),
#'   dist = rep(NA_real_, sum(!x$species)),
#'   lca = rep(0, sum(!x$species)),
#'   row.names = rownames(x)[!x$species]
#' ) -> anc
#' 
#' ## The scores of the ancestors are obtained as follows:
#' src_anc <- predict(PEM1, anc)
#' ## Note the absence of the 'vf' attributes since no variable has been
#' ## presented to PEM1.
#' 
#' ## Predictions are obtained from the previous linear model as follows:
#' predict(lm1, as.data.frame(src_anc))
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution 4: 1120-1131
#' 
#' @seealso \code{\link{PEM-functions}}
#' 
#' @importFrom stats qt
#' 
NULL
#' 
#' @describeIn PEM-class
#' 
#' Create PEM-class
#' 
#' A function creating a PEM-class object.
#' 
#' @export
PEM <- function(x, ..., a = -Inf, mm_a = matrix(1, nedge(x), 1L),
                psi = NULL, mm_psi = matrix(NA, nedge(x), 0L),
                d = "distance", sp = "species", tol = .Machine$double.eps^0.5) {
  
  if(!inherits(x, "graph"))
    stop("Parameter 'x' must be a graph-class object")
  
  if(is.null(edge(x)[[d]]))
    stop("Value of argument 'd' (", d, ") corresponds to no edge property")
  
  if(!is.numeric(edge(x)[[d]]))
    stop("Edge property '", d, "' is not numeric")
  
  if(is.null(x[[sp]]))
    stop("Value of argument 'sp' (", sp, ") corresponds to no vertex property.")
  
  if(!is.logical(x[[sp]]))
    stop("Vertex property '", sp, "' is not logical")
  
  if(any(!is.matrix(mm_a), !is.matrix(mm_psi)))
    stop("Argument 'mm_a' and 'mm_psi', must be matrices")
  
  ne <- nedge(x)
  nv <- nrow(x)
  
  if(any(ne != c(nrow(mm_a),nrow(mm_psi))))
    stop("Model matrices given as arguments 'mm_a' and 'mm_psi' must have ",
         ne, " rows")
  
  na <- length(a)
  npsi <- length(psi)
  
  if(na != ncol(mm_a))
    stop("The model matrix 'mm_a' must have ", na, " column(s)")
  
  if(npsi != ncol(mm_psi))
    stop("The model matrix 'mm_psi' must have ", npsi, " column(s)")
  
  ## Calculate the nodal values for the model matrices
  aggregateOnVertex(
    x,
    mm_a,
    function(x) ifelse(any(as.logical(x)), 1, 0),
    default = c(1, rep(0, na - 1L))
  ) -> mm_a_aov
  
  aggregateOnVertex(
    x,
    mm_psi,
    function(x) ifelse(any(as.logical(x)), 1, 0),
    default = rep(0, npsi)
  ) -> mm_psi_aov
  
  ## Calculate the values:
  edge(x)$a <- as.double(1/(1 + exp(-(mm_a %*% a))))
  edge(x)$psi <- if(npsi) as.double(2/(1 + exp(-(mm_psi %*% psi)))) else 1
  
  nsp <- sum(x[[sp]])
  
  matrix(
    .C("InflMatC", ne, nv, as.integer(edge(x)[[1L]]), as.integer(edge(x)[[2L]]),
       B = integer(nv * ne), PACKAGE = "MPSEM")$B,
    nrow = nv,
    ncol = ne
  ) -> B
  
  .C(
    "PEMbuildC",
    ne,
    nsp,
    Bc = as.double(B[x[[sp]],]),      ## [[3]] Bc
    means = double(ne),               ## [[4]] Column means of Bc
    as.double(edge(x)[[d]]),
    as.double(edge(x)$a),
    as.double(edge(x)$psi),
    w = double(ne),                   ## [[8]] The weights
    BcW = double(nsp*ne),             ## [[9]] BcW for the SVD
    PACKAGE = "MPSEM"
  )[c(3L,4L,8L,9L)] -> pem
  
  attr(pem$Bc, "dim") <- c(nsp, ne)
  attr(pem$BcW, "dim") <- c(nsp, ne)
  
  list(rownames(x)[x[[sp]]], edgenames(x)) ->
    dimnames(pem$Bc) -> dimnames(pem$BcW)
  
  pem <- c(pem, La.svd(pem$BcW, nsp, nsp))
  
  sel <- pem$d >= tol
  pem$d <- pem$d[sel]
  pem$u <- pem$u[,sel,drop = FALSE]
  pem$vt <- pem$vt[sel,,drop = FALSE]
  
  S2val <- NULL
  
  ## Embedded functions:
  
  nodal <- function(v) {
    
    matrix(
      data = NA,
      nrow = length(v),
      ncol = 2L,
      dimnames = list(rownames(x)[v], c("a","psi"))
    ) -> out
    
    for(i in 1L:length(v)) {
      out[i,1L] <- as.double(1/(1 + exp(-sum(mm_a_aov[v[i],] * a))))
      out[i,2L] <- as.double(2/(1 + exp(-sum(mm_psi_aov[v[i],] * psi))))
    }
    
    out
  }
  
  update <- function(a, psi, ...) {
    
    if(length(a) != na)
      stop("That PEM requires ", na, " selection parameter(s) (argument 'a')")
    
    if(length(psi) != npsi)
      stop("That PEM requires ", npsi, " evolution rate parameter(s) ",
           "(argument 'psi')")
    
    a <<- a
    psi <<- psi
    
    edge(x)$a <<- as.double(1/(1 + exp(-(mm_a %*% a))))
    edge(x)$psi <<- if(npsi) as.double(2/(1 + exp(-(mm_psi %*% psi)))) else 1
    
    .C(
      "PEMupdateC",
      ne,
      nsp,
      as.double(pem$Bc),
      as.double(edge(x)[[d]]),
      as.double(edge(x)$a),
      as.double(edge(x)$psi),
      as.double(pem$w),           ## [[7]]  The new weights.
      as.double(pem$BcW),         ## [[8]]  The new weighted influence matrix.
      PACKAGE = "MPSEM"
    )[c(7L,8L)] -> tmp
    
    pem$w <<- tmp[[1L]]
    pem$BcW[] <<- tmp[[2L]]
    tmp <- La.svd(pem$BcW, nsp, nsp)
    sel <- tmp$d > tol
    pem$d <<- tmp$d[sel]
    pem$u <<- tmp$u[, sel, drop = FALSE]
    pem$vt <<- tmp$vt[sel, , drop = FALSE]
    
    invisible(NULL)
  }
  
  var <- function(...)
    pem$u %*% diag(pem$d^2) %*% t(pem$u)
  
  inv_var <- function(...)
    pem$u %*% diag(pem$d^(-2)) %*% t(pem$u)
  
  logdet <- function(...)
    sum(log(pem$d^2))
  
  dev <- function(y, x = NULL, ...) {
    
    invS <- inv_var()
    
    logdetS <- logdet()
    
    res <- scale(y, scale=FALSE)
    
    if(!is.null(x)) {
      
      BX <- MASS::ginv(t(x) %*% invS %*% x) %*% t(x) %*% invS %*% res
      
      res <- res - x %*% BX
    }
    
    cte <- nsp + (nsp - 1)*log(2*pi)
    
    S2 <- double(ncol(res))
    
    dev <- 0
    
    for(i in 1L:ncol(res)) {
      
      S2[i] <- as.double(t(res[,i,drop = FALSE]) %*% invS %*%
                            res[,i,drop = FALSE]/nsp)
      
      dev <- dev + cte + nsp*log(S2[i]) + logdetS
    }
    
    names(S2) <- colnames(res)
    
    S2val <<- S2
    
    dev
  }
  
  S2 <- function(y, x = NULL, ...) {
    
    if(missing(y))
      return(S2val)
    
    invS <- inv_var()
    
    res <- scale(y, scale=FALSE)
    
    if(!is.null(x)) {
      
      BX <- MASS::ginv(t(x) %*% invS %*% x) %*% t(x) %*% invS %*% res
      
      res <- res - x %*% BX
    }
    
    S2 <- double(ncol(res))
    
    for(i in 1L:ncol(res))
      S2[i] <- drop(t(res[,i,drop = FALSE]) %*% invS %*%
                      res[,i,drop = FALSE])/nsp
    
    names(S2) <- colnames(res)
    
    S2val <<- S2
    
    S2
  }
  
  predictVertex <- function(v, lca) {
    
    if(is.character(v)) {
      
      i <- match(v, rownames(x))
      
      if(any(is.na(i)))
        warning("Unknown vertex: ", paste(v[is.na(i)], collapse=", "))
      
      v <- i[!is.na(i)]
    }
    
    if(missing(lca)) {
      lca <- rep(0, length(v))
    } else
      lca <- rep(lca, length.out=length(v))
    
    ntgt <- length(v)
    ex <- edge(x)
    ed <- ex[[d]]
    nsv <- length(pem$d)
    
    matrix(
      .C(
        "PEMLoc2Scores",
        ne,
        as.double(pem$means * pem$w),
        ntgt,
        as.double(t(apply(B[v,,drop=FALSE], 1L, function(x, y) x*y, y=ed))),
        as.double(ex$a),
        as.double(ex$psi),
        nsv,
        as.double(pem$d),
        as.double(pem$vt),
        double(nsv * ntgt),
        PACKAGE = "MPSEM"
      )[[10L]],
      nrow = ntgt,
      ncol = nsv,
      dimnames = list(
        rownames(x)[v],
        colnames(pem$u)
      )
    ) -> scr
    
    a_psi <- nodal(v)
    
    list(score=scr, vf=PEMvar(lca, a_psi[,1L], a_psi[,2L]))
  }
  
  predictEdge <- function(e, dst, lca) {
    
    if(is.character(e)) {
      
      nms <- names(e)
      i <- match(e, edgenames(x))
      
      if(any(is.na(i)))
        stop("Unknown vertex: ", paste(e[is.na(i)], collapse=", "))
      
      e <- i
      names(e) <- nms
      
    } else
      if(any(e > nedge(x)))
        stop("Edge index (", paste(e[e > nedge(x)], collapse=", "),
             ") higher than the number of edges")
    
    if(missing(dst)) {
      dst <- rep(0, length(e))
    } else
      dst <- rep(dst, length.out=length(e))
    
    if(missing(lca)) {
      lca <- rep(0, length(e))
    } else
      lca <- rep(lca, length.out=length(e))
    
    ntgt <- length(e)
    ex <- edge(x)
    ed <- ex[[d]]
    nsv <- length(pem$d)
    loc <- matrix(NA, ntgt, nedge(x))
    
    for(i in 1L:ntgt) {
      
      j <- e[i]
      loc[i,] <- B[ex[[1L]][j],] * ed
      
      if(dst[i] > ed[j])
        warning("Edge ", rownames(ex)[j], ", d = ",dst[i],
                " is higher than the maximum distance of ",
                round(ed[j], 6L))
      
      loc[i,j] <- dst[i]
    }
    
    nsv <- length(pem$d)
    
    matrix(
      .C(
        "PEMLoc2Scores",
        ne,
        as.double(pem$means * pem$w),
        ntgt,
        as.double(loc),
        as.double(edge(x)$a),
        as.double(edge(x)$psi),
        nsv,
        as.double(pem$d),
        as.double(pem$vt),
        double(nsv * ntgt),
        PACKAGE = "MPSEM"
      )[[10L]],
      nrow = ntgt,
      ncol = nsv,
      dimnames = list(
        names(e),
        colnames(pem$u)
      )
    ) -> scr
    
    list(score=scr, vf=PEMvar(lca, ex$a[e], ex$psi[e]))
  }
  
  structure(
    list(
      d = d,
      sp = sp,
      tol = tol,
      mm = list(a=mm_a, psi=mm_psi),
      nsp = nsp,
      B = B,
      print = function(...) {
        cat("A phylogenetic eigenvector map (PEM) for", nsp, "species:\n")
        if(nsp >= 11L)
          cat(paste(head(rownames(pem$u), 9L), collapse = ", "), "..., ", 
              tail(rownames(pem$u), 1L), "\n")
        else cat(paste(rownames(pem$u), collapse = ","), "\n")
        cat("obtained from the following phylogenetic graph:\n")
        print(x)
        invisible(NULL)
      },
      graph = function(...) x,
      pem = function(...) pem,
      par = function(...) list(a=a, psi=psi),
      nodal = nodal,
      update = update,
      var = var,
      inv_var = inv_var,
      logdet = logdet,
      dev = dev,
      S2 = S2,
      predictVertex = predictVertex,
      predictEdge = predictEdge
    ),
    class = "PEM"
  )
}
#' 
#' @describeIn PEM-class
#' 
#' Print PEM-class
#' 
#' A print method for PEM-class objects.
#' 
#' @method print PEM
#' 
#' @export
print.PEM <- function(x, ...) x$print(...)
#' 
#' @describeIn PEM-class
#' 
#' Method \code{as.matrix} for PEM-class Objects
#' 
#' A method to extract the phylogenetic eigenvectors from a PEM-class object.
#' 
#' @method as.matrix PEM
#' 
#' @export
as.matrix.PEM <- function(x, ...) {
  
  out <- x$pem(...)$u
  gr <- x$graph(...)
  
  dimnames(out) <- list(
    rownames(gr)[gr[[x$sp]]],
    sprintf("U_%d", 1L:ncol(out))
  )
  
  out
}
#' 
#' @describeIn PEM-class
#' 
#' Method \code{as.data.frame} for PEM-class Objects
#' 
#' A method to extract the phylogenetic eigenvectors from a PEM-class object.
#' 
#' @method as.data.frame PEM
#' 
#' @export
as.data.frame.PEM <- function(x, row.names = NULL, optional = FALSE, ...) {
  
  out <- as.data.frame(x$pem(...)$u)
  gr <- x$graph(...)
  if(is.null(row.names))
    row.names <- rownames(gr)[gr[[x$sp]]]
  
  dimnames(out) <- list(row.names, sprintf("U_%d", 1L:ncol(out)))
  
  out
}
#' 
#' @describeIn PEM-class
#' 
#' Update Method for PEM-class Objects
#' 
#' An updating method to recalculate a PEM with new smoothness and evolution
#' rate parameters.
#' 
#' @method update PEM
#' 
#' @export
update.PEM <- function(object, a, psi, ...)
  object$update(a=a, psi=psi, ...)
#' 
#' @describeIn PEM-class
#' 
#' Evolution Model Method for PEM-class Objects
#' 
#' A method to estimate the smoothness and evolution rate parameters of a
#' PEM-class object empirically on the basis of an observed quantitative trait
#' and, optionally, optional trait(s).
#' 
#' @importFrom stats optim
#' 
#' @method evolution.model PEM
#' 
#' @export
evolution.model.PEM <- function(object, y, ..., x = NULL) {
  
  if(!is.matrix(y))
    y <- cbind(y = y)
  
  if(object$nsp != nrow(y))
    stop("'y' has ", nrow(y), " observations, but 'object' has ", object$nsp,
         " species")
  
  if(!is.null(x)) {
    
    if(!is.matrix(x))
      x <- cbind(aux = x)
    
    if(object$nsp != nrow(x))
      stop("'x' has ", nrow(x), " observations, but 'object' has ", object$nsp,
           " species")
  }
  
  objf <- function(par, y, x, w, na, npsi, ..., verbose = FALSE) {
    
    w$update(a=par[1L:na], psi=if(npsi) par[na + (1L:npsi)] else c(), ...)
    
    out <- w$dev(y, x, ...)
    
    if(verbose) {
      pp <- w$par()
      
      cat("a:", paste(pp$a, collapse=", "), ";",
          "psi: ", paste(pp$psi, collapse=", "), ";",
          "value: ", out, "\n")
    }
    
    out
  }
  
  p <- object$par()
  
  if(any(!is.finite(p$a))) {
    
    warning("One or more non-finite initial value for parameter 'a'; ",
            "assumed to be 0")
    
    p$a[!is.finite(p$a)] <- 0
  }
  
  optim(
    par = unlist(p),
    fn = objf,
    method = "BFGS",
    y = y,
    x = x,
    w = object,
    na = length(p$a),
    npsi = length(p$psi)
  ) -> opt
  
  list(
    a = opt$par[1L:length(p$a)],
    psi = if(length(p$psi)) opt$par[length(p$a) + (1L:length(p$psi))] else c()
  ) -> opt$par
  
  opt
}
#' 
#' @describeIn PEM-class
#' 
#' Predict Method for PEM-class Objects
#' 
#' A predict method to obtain PEM prediction scores for arbitrary graph
#' locations.
#' 
#' @method predict PEM
#' 
#' @export
predict.PEM <- function(object, newdata, ...) {
  
  if(!is.data.frame(newdata))
    stop("Argument 'newdata' must be a data frame")
  
  if(is.null(newdata[[1L]]) || !is.numeric(newdata[[1L]]))
    stop("Missing or non-numeric 'ref' column in 'newdata'")
  
  if(is.null(newdata[[2L]]) || !is.numeric(newdata[[2L]]))
    stop("Missing or non-numeric 'dist' column in 'newdata'")
  
  if(is.null(newdata[[3L]]) || !is.numeric(newdata[[3L]]))
    stop("Missing or non-numeric 'lca' column in 'newdata'")
  
  matrix(
    data = NA,
    nrow = nrow(newdata),
    ncol = ncol(object$pem()$u),
    dimnames = list(
      rownames(newdata),
      sprintf("U_%d", 1L:ncol(object$pem()$u))
    )
  ) -> out
  
  S2val <- object$S2()
  if(!is.null(S2val)) {
    vf <- double(NROW(newdata))
    names(vf) <- rownames(newdata)
  }
    
  
  wh <- is.na(newdata[,2L])
  
  if(any(wh)) {
    
    object$predictVertex(
      v = newdata[wh,1L],
      lca = newdata[wh,3L]
    ) -> tmp
    
    out[wh,] <- tmp$score
    
    if(!is.null(S2val))
      vf[wh] <- tmp$vf
  }
  
  if(any(!wh)) {
    
    object$predictEdge(
      e = newdata[!wh,1L],
      dst = newdata[!wh,2L],
      lca = newdata[!wh,3L]
    ) -> tmp
    
    out[!wh,] <- tmp$score
    
    if(!is.null(S2val))
      vf[!wh] <- tmp$vf
  }
  
  if(!is.null(S2val))
    if(length(S2val) > 1L) {
      attr(out,"vf") <- matrix(vf, ncol=1L) %*% matrix(S2val, nrow=1L)
      dimnames(attr(out,"vf")) <- list(names(vf),names(S2val))
    } else {
      attr(out,"vf") <- S2val * vf
      names(attr(out,"vf")) <- names(vf)
    }
  
  out
}
#' 
