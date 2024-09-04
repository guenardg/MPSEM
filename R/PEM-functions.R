## **************************************************************************
##
##    (c) 2010-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Phylogenetic Eigenvector Maps (PEM) functions**
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
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with MPSEM. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Phylogenetic Eigenvector Maps
#' 
#' @description Functions to calculate and manipulate Phylogenetic Eigenvector
#' Maps (PEM), which are sets of eigenfunctions describing the structure of a 
#' phylogenetic graph. Each computation function is briefly described in 
#' section \code{Functions} below.
#' 
#' @name PEM-functions
#' 
#' @param x A \code{\link{graph-class}} object containing a phylogenetic graph.
#' @param w A \code{\link{graph-class}} object containing a phylogenetic graph.
#' @param object A \code{\link{PEM-class}} object containing a Phylogenetic
#' Eigenvector Map.
#' @param y One or many response variable(s) in the form of a single numeric 
#' vector or a \code{\link{matrix}}, respectively.
#' @param mroot Boolean (TRUE or FALSE) specifying whether multiple roots are
#' allowed.
#' @param d The name of the member of \code{x$edge} where the phylogenetic
#' distances (edge lengths) can be found.
#' @param a The steepness parameter describing whether changes occur, on
#' average: progressively long edges (a close to 0) or abruptly at vertices (a
#' close to 1).
#' @param psi Relative evolution rate along the edges (default: 1). This
#' parameter is only relevant when multiple values are assigned to different
#' portions of the phylogeny.
#' @param sp Name of the member of \code{x$vertex} where a \code{\link{logical}}
#' vertex property can be found, specifying which vertices are species (see
#' \code{\link{graph-class}}).
#' @param tol Eigenvalue threshold indicating that eigenvectors as usable.
#' @param lower Lower limit for the L-BFGS-B optimization algorithm
#' implemented in \code{\link{optim}}.
#' @param upper Upper limit for the L-BFGS-B optimization algorithm
#' implemented in \code{\link{optim}}.
#' @param tpall First parameter of function \code{getGraphLocations}:
#' Phylogenetic tree object with class \sQuote{phylo} (package \link[ape]{ape})
#' containing all species (model and target) used in the study.
#' @param targets Name of the target species to extract using the tree
#' \code{tpall}.
#' @param gsc The output of \code{getGraphLocations}.
#' 
#' @details Functions \code{\link{PEMInfluence}} and \code{\link{PEMweights}}
#' are used internally by \code{\link{PEM.build}} to create a binary matrix
#' referred to as an \sQuote{influence matrix} and weight its columns. That
#' matrix has a row for each vertex (or node) of graph \sQuote{x} and a column
#' for each of its edges. The elements of the influence matrix are 1 whenever
#' the vertex associated with a row is located in the tree, either directly or
#' indirectly downward the edge associated with a column. That function is
#' implemented in \code{C} language using recursive function calls. Although
#' \code{\link{PEMInfluence}} allows one to use multiple roots as its default
#' argument, it is called within \code{PEM.build} with \code{mroot = FALSE}.
#' User must therefore make sure that the graph provided to \code{PEMap} is
#' single-rooted.
#' 
#' Function \code{\link{PEM.build}} is used to produce a phylogenetic
#' eigenvector map, while function \code{\link{PEM.updater}} allows one to
#' re-calculate a \code{\link{PEM-class}} object with new weighting function
#' parameters. Function \code{\link{PEM.fitSimple}} performs a maximum
#' likelihood estimation of parameters \code{a} and \code{psi} assuming single
#' values for the whole tree, whereas function \code{\link{PEM.forcedSimple}}
#' allows one to impose values to arguments \code{a} and \code{psi} of a
#' \code{\link{PEM-class}} object, while making the function produce the same
#' details as \code{\link{PEM.fitSimple}} would have produced; these details are
#' necessary to make predictions.
#' 
#' Functions \code{\link{getGraphLocations}} returns the coordinates of a
#' species in terms of its position with respect to the influence matrix while
#' function \code{\link{Locations2PEMscores}} transforms these coordinates into
#' sets of scores that can be used to make predictions. Function
#' \code{\link{getAncGraphLocations}} produces the same output as
#' \code{\link{getGraphLocations}}, but for the ancestral species (i.e. the
#' nodes of the phylogeny) in order to estimate ancestral trait values.
#' 
#' @returns Function \code{\link{PEMInfluence}} returns the influence matrix of
#' graph \code{x} and function \code{\link{PEMweights}} returns weights
#' corresponding to the distances. Functions \code{\link{PEM.build}},
#' \code{\link{PEM.fitSimple}} and \code{\link{PEM.forcedSimple}} return a
#' \code{\link{PEM-class}} object. Function \code{\link{getGraphLocations}}
#' returns a list whose first member is an influence coordinate matrix whose
#' rows refer to the target species and columns refer to the edges. The second
#' member contains the lengths of the terminal edges connecting each target
#' species to the rest of the phylogeny.
#' 
#' Function \code{\link{Locations2PEMscores}} returns a list whose first member
#' is a PEM score matrix whose rows refer to the target species and columns
#' refer to the eigenvectors. The second member contains the variance associated
#' with the terminal edges connecting the target species to the phylogeny.
#' 
#' @author \packageAuthor{MPSEM} --
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution. 4: 1120--1131
#' 
#' Makarenkov, V., Legendre, L. & Desdevise, Y. 2004. Modelling phylogenetic
#' relationships using reticulated networks. Zoologica Scripta 33: 89--96
#' 
#' Blanchet, F. G., Legendre, P. & Borcard, D. 2008. Modelling directional
#' spatial processes in ecological data. Ecological Modelling 215: 325--336
#' 
#' @seealso \code{\link{PEM-class}}
#' 
#' @importFrom ape is.rooted drop.tip
#' @importFrom stats optim na.omit
#' @importFrom MASS ginv
#' 
#' @examples
#' ### Synthetic example
#' 
#' ## This example describes the phyogeny of 7 species (A to G) in a tree with 6
#' ## nodes, presented in Newick format, read by function
#' ## read.tree of package ape.
#' 
#' t1 <- read.tree(text=paste(
#'             "(((A:0.15,B:0.2)N4:0.15,C:0.35)N2:0.25,((D:0.25,E:0.1)N5:0.3,",
#'             "(F:0.15,G:0.2)N6:0.3)N3:0.1)N1;",sep=""))
#' t1                 # Summary of the structure of the tree
#' summary(t1)
#' 
#' x <- Phylo2DirectedGraph(t1)
#' 
#' ## Calculate the (binary) influence matrix; E1 to E12 are the tree edges
#' ## Edge E12 comes from the tree origin
#' PEMInfluence(x)
#' PEMInfluence(x)[x$vertex$species,]
#' 
#' ## Building phylogenetic eigenvector maps
#' PEM1 <- PEM.build(x)
#' PEM2 <- PEM.build(x, a = 0.2)
#' PEM3 <- PEM.build(x, a = 1)
#' PEM4 <- PEM.updater(PEM3,a=0.5)
#' 
#' ## Print summary statistics about PEM1
#' print(PEM1)
#' 
#' ## Extract the eigenvectors (species A--G, 6 eigenvectors)
#' as.data.frame(PEM4)
#' 
#' ## Example of a made up set of trait values for the 7 species
#' y <- c(A=-1.1436265,B=-0.3186166,C=1.9364105,D=1.7164079,E=1.0013993,
#'        F=-1.8586351,G=-2.0236371)
#' 
#' ## Estimate a single steepness parameter for the whole tree
#' PEMfs1 <- PEM.fitSimple(y=y, x=NULL, w=x, d="distance", sp="species",
#'                         lower=0, upper=1)
#' PEMfs1$optim       # Optimisation results
#' 
#' ## Force neutral evolution over the whole tree
#' PEMfrc1 <- PEM.forcedSimple(y=y,x=NULL,w=x,d="distance",sp="species",a=0)
#' PEMfrc1$x$edge$a   # Steepness parameter forced on each individual edge
#' 
#' ## Graph locations for target species X, Y, and Z not found in the original
#' ## data set
#' tpAll <- read.tree(text=paste("((X:0.45,((A:0.15,B:0.2)N4:0.15,",
#'                               "(C:0.25,Z:0.2)NZ:0.1)N2:0.05)NX:0.2,",
#'                               "(((D:0.25,E:0.1)N5:0.05,Y:0.25)NY:0.25,",
#'                               "(F:0.15,G:0.2)N6:0.3)N3:0.1)N1;",sep=""))
#' tpAll
#' summary(tpAll)     # Summary of the structure of the tree
#' 
#' grloc <- getGraphLocations(tpAll, c("X","Y","Z"))
#' grloc
#' 
#' PEMfs2 <- PEM.fitSimple(y=y, x=NULL, w=grloc$x, d="distance", sp="species",
#'                         lower=0, upper=1)
#' PEMfs2
#'
#' ## Same as for PEMfs1$optim
#' PEMfs2$optim
#' 
#' ## Get the PEM scores from the species graph locations:
#' PEMsc1 <- Locations2PEMscores(PEMfs2, grloc)
#' lm1 <- lm(y ~ V_2 + V_3 + V_5, data=PEMfs2)
#' 
#' ## Making prdictions for the species in locations `grloc`
#' ## using linear model `lm1`:
#' ypred <- predict(object=PEMfs2, targets=grloc, lmobject=lm1, interval="none")
#' 
#' ## Removing species X, Y, and Z from the tree in `tpAll`:
#' tpModel <- drop.tip(tpAll, c("X","Y","Z"))
#' 
#' ## Plot the results
#' layout(t(c(1,1,2)))
#' par(mar=c(6,2,2,0.5)+0.1)
#' plot(tpModel, show.tip.label=TRUE, show.node.label=TRUE, root.edge = TRUE,
#'      srt = 0, adj=0.5, label.offset=0.08, font=1, cex=1.5, xpd=TRUE)
#' edgelabels(paste("E", 1:nrow(tpModel$edge), sep=""),
#'            edge=1:nrow(tpModel$edge), bg="white", font=1, cex=1)
#' points(x=0.20,y=2.25,pch=21,bg="black")
#' lines(x=c(0.20,0.20,0.65), y=c(2.25,0.55,0.55), xpd=TRUE, lty=2)
#' text("X",x=0.69, y=0.55, xpd=TRUE, font=1, cex=1.5)
#' points(x=0.35, y=4.5,pch=21,bg="black")
#' lines(x=c(0.35,0.35,0.6), y=c(4.5,5.47,5.47), xpd=TRUE, lty=2)
#' text("Y", x=0.64, y=5.47, xpd=TRUE, font=1, cex=1.5)
#' points(x=0.35, y=3, pch=21, bg="black")
#' lines(x=c(0.35,0.35,0.55), y=c(3,3.5,3.5), xpd=TRUE, lty=2)
#' text("Z", x=0.59, y=3.5, xpd=TRUE, font=1, cex=1.5)
#' text(c("NX","NY","NZ"), x=c(0.20,0.35,0.35), y=c(2.25,4.5,3)+0.3*c(1,-1,-1),
#'      font=1, cex=1)
#' add.scale.bar(length=0.1, cex=1.25)
#' par(mar=c(3.75,0,2,2)+0.1)
#' plot(x=y, y=1:7, ylim=c(0.45,7), xlim=c(-4,4), axes=FALSE, type="n", xlab="")
#' axis(1, label=c("-4","-2","0","2","4"), at=c(-4,-2,0,2,4))
#' abline(v=0)
#' 
#' ## Plot the observed values
#' points(x=y, y=1:7, xlim=c(-2,2), pch=21, bg="black")
#' text("B)", x=-3.5, y=7, cex=1.5, xpd=TRUE)
#' text("Trait value", x=0, y=-0.5, cex=1.25, xpd=TRUE)
#' 
#' ## Plot the predicted values
#' points(x=ypred, y=c(0.5,5.5,3.5), pch=23, bg="white", cex=1.25)
#' 
#' ## Estimate the ancestral trait values
#' ANCloc <- getAncGraphLocations(x)
#' PEMfsAnc <- PEM.fitSimple(y=y, x=NULL, w=ANCloc$x, d="distance",
#'                           sp="species", lower=0, upper=1)
#' PEMfsAnc$optim
#' 
#' ## Get the PEM scores from the species graph locations:
#' PEManc1 <- Locations2PEMscores(PEMfsAnc, ANCloc)
#' 
#' ## Making predictions for the ancestral species whose locations are found in
#' ## `ANCloc` using the linear model `lm1`:
#' y_anc <- predict(object=PEMfsAnc, targets=ANCloc, lmobject=lm1,
#'                  interval="confidence")
#' 
#' @useDynLib MPSEM, .registration = TRUE
#' 
NULL
#' 
#' @describeIn PEM-functions
#' 
#' Influence Matrix
#' 
#' Calculates the influence matrix of a phylogenetic graph. The influence matrix
#' is a binary matrix whose rows and columns correspond to the vertices and
#' edges of the phylogenetic graph, respectively, and whose elements describe
#' whether a given edge had been taken by any ancestors of a vertex
#' (representing extinct of extant species) during evolution (value = 1) or not
#' (value = 0).
#' 
#' @export
PEMInfluence <- function(x,mroot=TRUE) {
  if(attr(x,"class") != "graph")
    stop("Parameter 'x' must be of class 'graph'")
  from <- x$edge[[1L]] ; to <- x$edge[[2L]] ; ev <- attr(x,"ev")
  roots <- (1L:ev[2L])[is.na(match(1L:ev[2L],to))] ; nroots <- length(roots)
  ## If several roots, add a guide vertex and nroots edges connecting each root
  ## from the guide.
  if(nroots != 1L) {
    if(!mroot)
      stop("Multiple (",nroots,
           ") roots are detected but multiple rooting is not allowed ",
           "(i.e. mroot == FALSE).")
    warning("Multiple (",nroots,") roots detected.")
    from <- c(rep(1L,nroots),from+1L)
    to <- c(roots,to)+1L
    ev <- ev+c(nroots,1L)
  }
  B <- matrix(.C("PEMInfMat",
                 as.integer(from),
                 as.integer(to),
                 as.integer(ev[1L]),
                 as.integer(ev[2L]),
                 B = integer(ev[2L]*ev[1L]))$B,ev[2L],ev[1L])
  if(nroots != 1L)
    B <- B[-1L,]    ## Strip the guidance vertex if many roots.
  if(!is.null(attr(x,"vlabel")))
    rownames(B) <- attr(x,"vlabel")
  if(!is.null(attr(x,"elabel"))) {
    if(nroots == 1L) {  ## Edge labels should be handled with 
      colnames(B) <- attr(x,"elabel")
    } else {
      colnames(B) <- c(paste("V->",attr(x,"vlabel")[roots],sep=""),
                       attr(x,"elabel"))
    }
  }
  return(B)
}
#' 
#' @describeIn PEM-functions
#' 
#' PEM Weighting
#' 
#' A power function to obtain the edge weights used during PEM calculation.
#' 
#' @export
PEMweights <- function (d, a=0, psi=1) {
  nd <- length(d)
  a <- rep(a, length.out = nd)
  psi <- rep(psi, length.out = nd)
  return(.C("PEMweightC",
            as.double(d),
            as.integer(nd),
            as.double(a),
            as.double(psi),
            res = double(nd))$res)
}
#' 
#' @describeIn PEM-functions
#' 
#' PEM Building
#' 
#' Calculates a PEM with parameters given by arguments a and psi.
#' 
#' @export
PEM.build <- function(x, d="distance", sp="species", a=0, psi=1,
                      tol=.Machine$double.eps^0.5) {
  if(attr(x,"class") != "graph")
    stop("Parameter 'x' must be of class 'graph'")
  if(is.null(x$edge[[d]]))
    stop("There is no property '",d,"' to be used as edges lengths.")
  if(is.null(x$vertex[[sp]]))
    stop("There is no property '",sp,"' to indicate species vertices.")
  nsp <- sum(x$vertex[[sp]])
  ### All graphs sould be single-rooted
  ev <- as.integer(attr(x, "ev"))  # Just to be sure.
  a <- rep(a,length.out=ev[1L])
  psi <- rep(psi,length.out=ev[1L])
  out <- list(x=x,sp=x$vertex[[sp]])
  out[["B"]] <-
    matrix(.C("PEMInfMat",
              as.integer(x$edge[[1L]]),
              as.integer(x$edge[[2L]]),
              ev[1L],ev[2L],
              B = integer(ev[2L]*ev[1L]))$B,ev[2L],ev[1L])[x$vertex[[sp]],]
  out <- c(out,.C("PEMbuildC",
                  ne=ev[1L],nsp=nsp,
                  Bc=as.double(out$B),
                  means=double(ev[1L]),
                  dist=as.double(x$edge[[d]]),
                  a=as.double(a),
                  psi=as.double(psi),
                  w=double(ev[1L]),
                  BcW=double(nsp*ev[1L])))
  attr(out$Bc,"dim") <- c(nsp,ev[1L])
  attr(out$BcW,"dim") <- c(nsp,ev[1L])
  dimnames(out$Bc) <- dimnames(out$BcW) <-
    list(attr(x,"vlabel")[x$vertex[[sp]]],attr(x,"elabel"))
  out <- c(out, La.svd(out$BcW,nsp,nsp))
  sel <- out$d >= tol
  out$d <- out$d[sel]
  out$u <- out$u[,sel,drop=FALSE]
  out$vt <- out$vt[sel,,drop=FALSE]
  rownames(out$vt) <- colnames(out$u) <- paste("V",1L:sum(sel),sep="_")
  rownames(out$u) <- attr(x,"vlabel")[x$vertex[[sp]]]
  colnames(out$vt) <- attr(x,"elabel")
  attr(out,"class") <- "PEM"
  return(out)
}
#' 
#' @describeIn PEM-functions
#' 
#' PEM Update
#' 
#' Update a PEM with new parameters given by arguments a and psi.
#' 
#' @export
PEM.updater <- function(object,a,psi=1,tol=.Machine$double.eps^0.5) {
  nsp <- object$nsp
  ne <- object$ne
  object$a <- rep(a,length.out=ne)
  object$psi <- rep(psi,length.out=ne)
  out <- .C("PEMupdateC",
            as.integer(ne),
            as.integer(nsp),
            as.double(object$Bc),
            as.double(object$dist),
            as.double(object$a),
            as.double(object$psi),
            w=as.double(object$w),
            BcW=as.double(object$BcW))
  object$w <- out$w
  object$BcW[] <- out$BcW
  out <- La.svd(object$BcW,nsp,nsp)
  sel <- out$d > tol
  object$d <- out$d[sel]
  object$u[] <- out$u[, sel, drop = FALSE]
  object$vt[] <- out$vt[sel, , drop = FALSE]
  return(object)
}
#' 
#' @describeIn PEM-functions
#' 
#' Fitting a PEM to Data while Estimating Global Steepness
#' 
#' Fits a PEM to a data set estimating the selection (steepness) parameter using
#' gradient descent. The selection and evolution rate (psi = 1) are assumed to
#' be homogeneous for the whole phylogenetic network.
#' 
#' @export
PEM.fitSimple <- function(y, x, w, d="distance", sp="species", lower=0, upper=1,
                          tol=.Machine$double.eps^0.5) {
  object <- PEM.build(w,d=d,sp=sp,a=lower)
  nsp <- sum(w$vertex[[sp]])
  if(is.matrix(y)) {
    if(nrow(y)!=nsp)
      stop("The number of rows in 'y' must match the number of species in ",
           "phylogenetic graph 'w'")
    p <- ncol(y)
  } else {
    if(length(y)!=nsp)
      stop("The number of elements in 'y' must match the number of species ",
           "in phylogenetic graph 'w'")
    y <- cbind(y)
    p <- 1L
  }
  f <- function(par,y,x,object,nsp,p) {
    object <- PEM.updater(object,a=par[1L],psi=1,tol=tol)
    Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
    logdetS <- sum(log(object$d^2))
    if(is.null(x)) {
      BX <- mean(y)
      res <- y-BX
    } else {
      BX <- ginv(t(x)%*%Sinv%*%x) %*% t(x) %*% Sinv %*% y
      res <- y-x%*%BX
      BX <- c(mean(res),BX)
      res <- res-BX[1]
    }
    dev <- 0
    for(i in 1L:p) {
      S2 <- as.numeric(t(res[,i,drop=FALSE])%*%Sinv%*%res[,i,drop=FALSE]/nsp)
      dev <- dev + nsp + (nsp-1L)*log(2*pi) + nsp*log(S2) + logdetS
    }
    return(dev)
  }
  opt <- optim(par=0,f,method="L-BFGS-B",lower=lower,upper=upper,y=y,x=x,
               object=object,nsp=nsp,p=p)
  if(opt$convergence)
    warning("No optimum found... Message from optim() - ",opt$message,
            ". Status = ",opt$convergence)
  object <- PEM.updater(object,a=opt$par[1L],psi=1,tol=tol)
  Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
  if(is.null(x)) {
    BX <- mean(y)
    res <- y-BX
  } else {
    BX <- MASS::ginv(t(x)%*%Sinv%*%x) %*% t(x) %*% Sinv %*% y
    res <- y-x%*%BX
    BX <- c(mean(res),BX)
    res <- res-BX[1]
  }
  S2 <- numeric(p)
  for(i in 1L:p)
    S2[i] <- as.numeric(t(res[,i,drop=FALSE])%*%Sinv%*%res[,i,drop=FALSE]/nsp)
  object$S2 <- S2
  object$y <- y
  object$optim <- opt
  object$x$edge[["a"]] <- rep(opt$par[1L],attr(object$x,"ev")[1L])
  return(object)
}
#' 
#' @describeIn PEM-functions
#' 
#' Fitting a PEM to Data while Forcing Global Steepness
#' 
#' Fits a PEM to a data set forcing a user-provided selection (steepness)
#' parameter. The selection and evolution rate (psi = 1) are assumed to be
#' homogeneous for the whole phylogenetic network.
#' 
#' @export
PEM.forcedSimple <- function(y, x, w, d="distance", sp="species", a=0, psi=1,
                             tol=.Machine$double.eps^0.5) {
  object <- PEM.build(w,d=d,sp=sp,a=a,psi=psi,tol=tol)
  nsp <- sum(w$vertex[[sp]])
  if(is.matrix(y)) {
    if(nrow(y)!=nsp)
      stop("The number of rows in 'y' must match the number of species in ",
           "phylogenetic graph 'w'")
    p <- ncol(y)
  } else {
    if(length(y)!=nsp)
      stop("The number of elements in 'y' must match the number of species ",
           "in phylogenetic graph 'w'")
    y <- cbind(y) ; p <- 1L
  }
  Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
  if(is.null(x)) {
    BX <- mean(y) ; res <- y-BX
  } else {
    BX <- MASS::ginv(t(x)%*%Sinv%*%x) %*% t(x) %*% Sinv %*% y
    res <- y-x%*%BX
    BX <- c(mean(res),BX)
    res <- res-BX[1]
  }
  S2 <- numeric(p)
  for(i in 1L:p)
    S2[i] <- as.numeric(t(res[,i,drop=FALSE])%*%Sinv%*%res[,i,drop=FALSE]/nsp)
  object$S2 <- S2
  object$y <- y
  object$optim <- list()
  object$optim$par <- c(a,psi)
  object$x$edge[["a"]] <- rep(a,attr(object$x,"ev")[1L])
  object$x$edge[["psi"]] <- rep(psi,attr(object$x,"ev")[1L])
  return(object)
}
#' 
#' @describeIn PEM-functions
#' 
#' Get Phylogenetic Graph Locations
#' 
#' Takes a phylogenetic tree and a list of species to be removed, and produce a
#' phylogenic graph without these species together with the locations of the
#' removed species on that graph (i.e., the location where the removed species
#' would be found should they be inserted again in the phylogenetic graph).
#' 
#' @export
getGraphLocations <- function(tpall, targets) {
  oldnlab <- tpall$node.label
  tpall$node.label <- newnlab <- paste("N",1:tpall$Nnode,sep="")
  tpmodel <- drop.tip(tpall, targets)
  tpmodel$root.edge <- tpall$root.edge
  xmodel <- Phylo2DirectedGraph(tpmodel)
  Bmodel <- PEMInfluence(xmodel)
  loc <- matrix(NA, length(targets), ncol(Bmodel),
                dimnames = list(targets,colnames(Bmodel)))
  dtt <- rep(NA, length(targets))
  if(length(targets))
    for (i in 1L:length(targets)) {
      ## Target species need to be striped one-by-one.
      tptargetonly <-
        if (length(targets) == 1L) tpall else drop.tip(tpall, targets[-i])
      tptargetonly$root.edge <- tpall$root.edge
      Vnames <- c(tptargetonly$tip.label,tptargetonly$node.label)
      ## The vertices names.
      VBidx <- which(tptargetonly$edge[,2L] == which(Vnames == targets[i]))
      ## The index of the vertex immediatly below.
      VB <- tptargetonly$edge[VBidx,1L]
      ## VB: the vertex immediatly below
      VBName <- Vnames[VB]
      ## The name of VB
      dBX <- tptargetonly$edge.length[VBidx]
      ## Distance between VB and VX
      VA <- tptargetonly$edge[tptargetonly$edge[,2L] == VB,1L]
      ## VA: the vertex below the one immediatly below and
      VAName <- Vnames[VA]
      ## its name. That vertex may not exist.
      dAB <- tptargetonly$edge.length[tptargetonly$edge[,2L] == VB]
      ## Distance between VB and VA (if exists)
      VC <- tptargetonly$edge[tptargetonly$edge[,1L] == VB,2L]
      ## The vertices above the vertex immediatly below and
      VCName <- Vnames[VC]
      ## their names. 2 or more (in case of tri+ chotomy).
      dBC <- tptargetonly$edge.length[tptargetonly$edge[,1L] == VB]
      ## Distances above VB.
      VC <- VC[VCName!=targets[i]]
      dBC <- dBC[VCName!=targets[i]]
      VCName <- Vnames[VC]
      ## Strip the target itself from the VC list.
      if(length(VA)==0) {
        ## If the target species is beyond the root:
        loc[i,] <- 0
        ## In all case: location = the root.
        dtt[i] <- dBX+min(dBC)
        ## Distance between the off-root species and the remaining ones.
      } else {
        ## If the target species is not beyond the root:
        dtt[i] <- dBX
        ## In all cases: dtt == the distance between the target and the vertex
        ## immediatly below.
        if(length(VC)>1) {
          ## When VB is more than dichotomic (several Vertices C):
          loc[i,] <- Bmodel[VBName,] * xmodel$edge$distance
          ## Coordinates are those of the still-existing vertex B.
        } else {
          ## When VB collapses (it was dichotomic):
          loc[i,] <- Bmodel[VAName,] * xmodel$edge$distance
          ## Copy the coordinates of vertex VA.
          loc[i,Bmodel[VAName,] != Bmodel[VCName,]] <- dAB
          ## Change the length of that edge connecting VA to the single VC
          ## (0 for VA) to dAB.
        }
      }
    }
  if(!is.null(oldnlab)) {
    oldnlab <- oldnlab[match(attr(xmodel,"vlabel"),newnlab)]
    attr(xmodel,"vlabel")[!is.na(oldnlab)] <- na.omit(oldnlab)
  }
  return(list(x = xmodel, locations = loc, LCA2target = dtt))
}
#' 
#' @describeIn PEM-functions
#' 
#' Get Ancestral Species Location
#' 
#' Get the location on the phylogenetic graph of the immediate ancestors for a
#' list of species. The species of the list remain in the resulting phylogenetic
#' graph. This function is useful for estimating the ancestral state of a trait.
#' 
#' @export
getAncGraphLocations <- function(x, tpall) {
  if(missing(x) && !missing(tpall))
    x <- Phylo2DirectedGraph(tpall)
  loc <- PEMInfluence(x)[!x$vertex$species,] * x$edge$distance
  return(list(x = x, locations = loc, LCA2target = numeric(nrow(loc))))
}
#' 
#' @describeIn PEM-functions
#' 
#' PEM Score Calculation
#' 
#' Calculates the scores of an extant or ancestral species on a phylogenetic
#' eigenvector map (i.e., its value on the eigenvectors of the map) from its
#' location on the phylogenetic graph used to build that map.
#' 
#' @export
Locations2PEMscores <- function(object, gsc) {
  if(object$ne != ncol(gsc$locations))
    stop("Numbers of edge coordinates mismatch: gsc$locations has ",
         ncol(gsc$locations)," columns while the phylogenetic graph has ",
         object$ne," edges.")
  if(is.null(object$y)||is.null(object$optim))
    stop("PEM has no attached data/optimization parameters.")
  ntgt <- nrow(gsc$locations) ; nsv <- length(object$d)
  out <- list()
  out[["scores"]] <- matrix(.C("PEMLoc2Scores",
                               as.integer(object$ne),
                               as.double(object$means*object$w),
                               as.integer(ntgt),
                               as.double(gsc$locations),
                               as.double(object$a),
                               as.double(object$psi),
                               as.integer(nsv),
                               object$d,
                               object$vt,
                               scores=double(nsv*ntgt))$scores,ntgt,nsv)
  dimnames(out[["scores"]]) <- list(rownames(gsc$locations),colnames(object$u))
  VarianceFactor <- PEMvar(gsc$LCA2target,a=object$optim$par[1L],psi=1)
  VarianceFactor <-
    matrix(object$S2,length(object$S2),1L) %*%
    matrix(VarianceFactor,1L,length(VarianceFactor))
  dimnames(VarianceFactor) <- list(colnames(object$y),rownames(gsc$locations))
  out[["VarianceFactor"]] <- VarianceFactor
  return(out)
}
##
