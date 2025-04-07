## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Class and method for PEM linear models **
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
#' PEM Linear Model Interface
#' 
#' @description An interface to estimate a linear model 
#' 
#' @docType class
#' 
#' @name pemlm-class
#' 
#' @param formula An object of class "\code{\link{formula}}" (or one that can be
#' coerced to that class): a symbolic description of the model to be fitted.
#' (see `details` in \code{\link{lm}} for further details on model
#' specification.
#' @param data An optional data frame, list or environment (or object coercible
#' by \code{\link{as.data.frame}} to a data frame) containing the variables in
#' the model. If not found in data, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{pemlm} is called.
#' @param pem A \code{\link{PEM-class}} object.
#' @param contrasts An optional list. See the \code{contrasts.arg} of
#' \code{\link{model.matrix.default}}.
#' @param x,object A \code{pemlm-class} object.
#' @param newdata An optional data frame containing auxiliary traits used for
#' making predictions.
#' @param newloc A data frame containing the graph locations of species targeted
#' for making predictions.
#' @param se.fit A \code{\link{logical}} specifying whether to return standard
#' errors (default: \code{FALSE}).
#' @param interval One of \code{"none"}, \code{"confidence"}, or
#' \code{"prediction")} specifying the type of interval associated with the
#' predictions.
#' @param level Tolerance/confidence level (default: \code{0.95}).
#' @param ... Additional arguments to be passed to the low level regression
#' fitting functions.
#' 
#' @details ...
#' 
#' @format A \code{pemlm-class} object contains:
#' \describe{
#'   \item{auxModel}{A \code{\link{lm}}-class object.}
#'   \item{resetPEM}{A function ...}
#'   \item{getPEM}{A function ...}
#'   \item{aic}{A function ...}
#'   \item{getIncluded}{A function ...}
#'   \item{getCandidate}{A function ...}
#'   \item{scanCandidate}{A function ...}
#'   \item{forward}{A function ...}
#'   \item{promote}{A function ...}
#'   \item{pemModel}{A function ...}
#' }
#' 
#' @return
#' \describe{
#'   \item{pemlm}{A \code{pemlm-class} object.}
#'   \item{print.pemlm}{\code{NULL} (invisibly).}
#'   \item{summary.pemlm}{A \code{summary.lm-class} object.}
#'   \item{anova.pemlm}{An \code{anova-class} data frame.}
#'   \item{predict.pemlm}{...}
#' }
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution 4: 1120-1131
#' 
#' @examples ## Example here...
#' 
#' @importFrom stats anova delete.response predict is.empty.model lm.fit .getXlevels as.formula terms .checkMFClasses na.fail
#' 
NULL
#' 
#' @describeIn pemlm-class
#' 
#' PEM Linear Model
#' 
#' Calculate a PEM-based linear model for estimating trait values.
#' 
#' @export
pemlm <- function(formula, data, pem, ..., contrasts = NULL) {
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  
  ## if(is.matrix(y))
  ##   stop("This function does not support multivariate regression natively.")
  
  ny <- NROW(y)
  
  if(ny != pem$nsp)
    stop("The number of species in 'pem' (", pem$nsp,") is not equal to the ",
         "number of observations (", ny,")")
  
  if(is.empty.model(mt))
    stop("The model is empty.")
  
  mm <- model.matrix(mt, mf, contrasts)
  u <- as.matrix(pem)
  
  resetPEM <- function() {
    
    assign("included", integer(0L), inherits=TRUE)
    assign("candidate", 1L:ncol(u), inherits=TRUE)
    assign("pemModel", NULL, inherits=TRUE)
    
    invisible(NULL)
  }
  
  resetPEM()
  
  getAux <- function(...) {
    
    z <- lm.fit(mm, y, offset = NULL, singular.ok = TRUE, ...)
    
    class(z) <- c(if(is.matrix(y)) "mlm","lm")
    z$contrasts <- attr(mm, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    z$model <- mf
    
    z
  }
  
  auxModel <- getAux(...)
  
  aic <- function(s, ..., k = 2, corrected = TRUE) {
    
    z <-  lm.fit(cbind(mm, if(!missing(s)) u[,s,drop=FALSE]), y, offset = NULL,
                 singular.ok = TRUE, ...)
    
    res <- z$residuals
    N <- NROW(res)
    M <- NCOL(res)
    logl <- -0.5 * (N*M*(log(2*pi) + M - log(N*M) + log(sum(res^2))))
    aic <- -2*logl + k*M*(z$rank + 1L)
    
    if(corrected) {
      k1 <- length(z$coefficients)
      aic <- aic + 2*k1*(k1 + 1)/(N*M - k1 - 1)
    }
    
    aic
  }
  
  scanCandidate <- function(..., k = 2, corrected = TRUE) {
    
    out <- numeric(length(candidate))
    
    for(i in 1L:length(candidate))
      out[i] <- aic(c(included,candidate[i]), ..., k=k, corrected=corrected)
    
    out
  }
  
  getModel <- function(...) {
    
    x <- cbind(mm, u[,included,drop=FALSE])
    asn <- attr(mm, "assign")
    attr(x, "assign") <- c(asn, rep(max(asn) + 1L, length(included)))
    attr(x, "contrasts")$U <- "pem"
    
    z <- lm.fit(x, y, offset = NULL, singular.ok = TRUE, ...)
    class(z) <- c(if(is.matrix(y)) "mlm", "lm")
    z$contrasts <- attr(x, "contrasts")
    
    z$xlevels <- .getXlevels(mt, mf)
    z$xlevels$U <- rownames(pem)
    
    z$call <- cl
    
    mtat <- attributes(mt)
    str2lang(
      paste(
        "list(",
        paste(c(as.character(mtat$variables)[-1L],"U"), collapse=","),
        ")",
        sep=""
      )
    ) -> mtat$variables
    
    mtat$factors <- rbind(cbind(mtat$factors, U=0L), U=0L)
    mtat$factors["U","U"] <- 1L
    mtat$term.labels <- c(mtat$term.labels,"U")
    mtat$order <- c(mtat$order,1L)
    
    str2lang(
      paste(
        "list(",
        paste(c(as.character(mtat$predvars)[-1L],"U"), collapse=","),
        ")",
        sep=""
      )
    ) -> mtat$predvars
    
    mtat$dataClasses <- c(mtat$dataClasses, U="pem")
    z$terms <- as.formula(paste(paste(deparse(mt), collapse=""), "U", sep=" + "))
    attributes(z$terms) <- mtat
    z$model <- mf
    
    z
  }
  
  forward <- function(..., k = 2, corrected = TRUE) {
    
    a <- aic(included, ..., k=k, corrected=corrected)
    
    repeat{
      
      b <- scanCandidate(..., k=k, corrected=corrected)
      
      if(min(b) < a) {
        wh <- which.min(b)
        assign("included", c(included, candidate[wh]), inherits=TRUE)
        assign("candidate", candidate[-wh], inherits=TRUE)
        a <- b[wh]
      } else
        break
    }
    
    assign("pemModel", getModel(...), inherits=TRUE)
    
    included
  }
  
  promote <- function(wh, ...) {
    
    idx <- match(wh, candidate)
    
    if(any(is.na(idx))) {
      warning("Unknown or unavailable eigenfunction(s): ",
              paste(wh[is.na(idx)],collapse=", "))
      idx <- idx[!is.na(idx)]
    }
    
    if(length(idx)) {
      assign("included", c(included, candidate[idx]), inherits=TRUE)
      assign("candidate", candidate[-idx], inherits=TRUE)
    }
    
    assign("pemModel", getModel(...), inherits=TRUE)
    
    included
  }
  
  structure(
    list(
      auxModel = auxModel,
      resetPEM = resetPEM,
      getPEM = function() pem,
      aic = aic,
      getIncluded = function() included,
      getCandidate = function() candidate,
      scanCandidate = scanCandidate,
      forward = forward,
      promote = promote,
      pemModel = function() pemModel
    ),
    class = "pemlm"
  )
}
#' 
#' @describeIn pemlm-class
#' 
#' Print PEM Linear Model
#' 
#' A print method for pemlm-class objects.
#' 
#' @method print pemlm
#' 
#' @export
print.pemlm <- function(x, ...) {
  
  cat("A PEM-based linear model\n-------------------\n\n")
  cat("Auxiliary linear model:\n")
  print(x$auxModel)
  incl <- x$getIncluded()
  cat("Includes",if(!length(incl)) "no" else length(incl),"eigenfunction(s)\n")
  if(length(incl)) {
    cat("PEM-based linear model:\n")
    print(x$pemModel())
  }
  
  invisible(NULL)
}
#' 
#' @describeIn pemlm-class
#' 
#' Summarize pemlm Model
#' 
#' A summary method for \code{pemlm-class} objects.
#' 
#' @method summary pemlm
#' 
#' @export
summary.pemlm <- function(object, ...) {
  out <- if(is.null(object$pemModel())) object$auxModel else object$pemModel()
  summary(out, ...)
}
#' 
#' @describeIn pemlm-class
#' 
#' Anova of a pemlm Model
#' 
#' An anova method for \code{pemlm-class} objects.
#' 
#' @method anova pemlm
#' 
#' @export
anova.pemlm <- function(object, ...) {
  out <- if(is.null(object$pemModel())) object$auxModel else object$pemModel()
  anova(out, ...)
}
#' 
#' @describeIn pemlm-class
#' 
#' Predict Trait Values
#' 
#' A predict method for \code{pemlm-class} objects.
#' 
#' @method predict pemlm
#' 
#' @export
predict.pemlm <- function(object,
                          newdata = data.frame(row.names = rownames(newloc)),
                          newloc, se.fit = FALSE,
                          interval = c("none","confidence","prediction"),
                          level = 0.95, ...) {
  
  tt <- delete.response(terms(object$auxModel))
  m <- model.frame(tt, newdata, na.action=na.fail,
                   xlev=object$auxModel$xlevels)
  
  if(!is.null(cl <- attr(tt, "dataClasses")))
    .checkMFClasses(cl, m)
  
  X <- model.matrix(tt, m, contrasts.arg=object$auxModel$contrasts)
  
  pm <- object$pemModel()
  incl <- object$getIncluded()
  n <- NROW(pm$residuals)
  m <- NCOL(pm$residuals)
  p <- pm$rank
  
  if(p) {
    p1 <- seq_len(p)
    piv <- pm$qr$pivot[p1]
  } else
    stop("The model has no rank.")
  
  U <- predict.PEM(object$getPEM(), newloc)
  XU <- cbind(X, U[, incl, drop = FALSE])[, piv, drop = FALSE]
  
  if(p < NCOL(XU))
    stop("The fit is rank-deficient; consider removing model terms.")
  
  interval <- match.arg(interval)
  
  if(se.fit || interval != "none") {
    
    ## Verify if the vf is suitable, otherwise, vf <- NULL
    if(!is.null(vf <- attr(U,"vf")))
      if(NROW(vf) != NROW(X) || NCOL(vf) != m) {
        warning("Unsuitable variance factor vector or matrix.")
        vf <- NULL
      }
    
    df <- pm$df.residual
    XRinv <- XU[,piv,drop=FALSE] %*% qr.solve(qr.R(pm$qr)[p1, p1])
    
    if(m > 1L) {
      res.var <- colSums(pm$residuals^2)/df
      ip <- drop(XRinv^2 %*% matrix(res.var, p, m, TRUE, list(NULL, names(res.var))))
    } else {
      res.var <- sum(pm$residuals^2)/df
      ip <- drop(XRinv^2 %*% rep(res.var, p))
    }
  }
  
  if(m > 1L) {
    
    predictor <- drop(XU %*% pm$coefficients[piv,])
    
    if(interval != "none") {
      
      tfrac <- qt((1 - level)/2, df)
      
      if(interval == "confidence") {
        hwid <- tfrac * sqrt(ip + if(!is.null(vf)) vf/n)
      } else if(interval == "prediction") {
        matrix(res.var, nrow(X), m, TRUE,
               list(rownames(X), names(res.var))) -> tmp
        hwid <- tfrac * sqrt(ip + tmp + if(!is.null(vf)) vf)
      }
      
      array(
        c(predictor, predictor + hwid, predictor - hwid),
        dim = c(NROW(predictor),NCOL(predictor),3L),
        dimnames = list(
          rownames(predictor),
          colnames(predictor),
          c("fit","lwr","upr")
        )
      ) -> predictor
    }
    
    if(se.fit) {
      se <- sqrt(ip + if(!is.null(vf)) vf/n)
      list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
    } else
      predictor
    
  } else {
    
    predictor <- drop(XU %*% pm$coefficients[piv])
    
    if(interval != "none") {
      
      tfrac <- qt((1 - level)/2, df)
      
      switch(
        interval,
        confidence = sqrt(ip + if(!is.null(vf)) vf/n),
        prediction = sqrt(ip + res.var + if(!is.null(vf)) vf)
      ) * tfrac -> hwid
      
      predictor <- cbind(predictor, predictor + hwid %o% c(1, -1))
      colnames(predictor) <- c("fit", "lwr", "upr")
    }
    
    if(se.fit) {
      se <- sqrt(ip + if(!is.null(vf)) vf/n)
      list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
    } else
      predictor
  }
}
#' 
