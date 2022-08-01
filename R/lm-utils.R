## **************************************************************************
##
##    (c) 2010-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Linear modelling utilisy functions**
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
#' Linear Modelling Utility Functions
#' 
#' @description Utility functions to build linear models using Phylogenetic
#' Eigenvector Maps among their explanatory variables.
#' 
#' @name lm-utils
#' 
#' @aliases utils
#' 
#' @param y A response variable.
#' @param x Descriptors (numeric of \code{\link{factor}}) to be used as
#' auxiliary traits.
#' @param object A \code{\link{PEM-class}} object.
#' @param alpha The p-value threshold above which the function will stop adding
#' variables.
#' 
#' @details Function \code{\link{lmforwardsequentialsidak}}, performs a forward
#' stepwise selection of the PEM eigenvectors until the familywise test of
#' significance of the new variable to be included exceeds the p-value 
#' threshold \code{alpha}. The familiwise type I error probability is obtained
#' using the Holm-Sidak correction of the testwise probabilities, thereby
#' correcting for type I error rate inflation due to multiple testing.
#' 
#' Function \code{lmforwardsequentialAICc} carries out forward stepwise
#' selection of the eigenvectors as long as the candidate model features a 
#' sample-size-corrected Akaike information criterion lower than the previous
#' model. The final model should be regarded as overfitted from the
#' Neyman-Pearson (\emph{i.e.} frequentist) point of view, but this is the model
#' that minimizes information loss from the standpoint of information theory.
#' 
#' @return An \code{\link{lm}-class} object.
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Burnham, K. P. & Anderson, D. R. 2002. Model selection and multimodel
#' inference: a practical information-theoretic approach, 2nd ed.
#' Springer-Verlag. xxvi + 488 pp.
#' 
#' Holm, S. 1979. A simple sequentially rejective multiple test procedure.
#' Scand. J. Statist. 6: 65-70
#' 
#' Sidak, Z. 1967. Rectangular confidence regions for means of multivariate
#' normal distributions. J. Am. Stat. Ass. 62, 626-633
#' 
#' @importFrom stats lm AIC anova as.formula
#' 
NULL
#' 
#' @describeIn lm-utils
#' 
#' Forward Stepwise Regression AICc
#' 
#' Forward stepwise variable addition using the sample-size-corrected Akaike
#' Information Criterion.
#' 
#' @export
lmforwardsequentialAICc <- function (y, x, object) {
  AuxTrait <- !missing(x)
  included <- numeric(0)
  candidates <- 1L:ncol(object$u)
  df1 <- if(AuxTrait) {x <- cbind(x) ; data.frame(x,object)} else object
  p1 <- if(AuxTrait) paste("y~",paste(colnames(x),collapse="+"),sep="") else "y~"
  while (TRUE) {
    p2 <- paste(if(length(included)) paste(if(AuxTrait)"+",paste(paste("V_", included, sep = ""), collapse = "+")) else if(AuxTrait) "" else "1" ,sep="")
    lm1 <- lm(as.formula(paste(p1,p2,sep="")), data = df1)
    k1 <- length(lm1$coef)
    AICc1 <- AIC(lm1) + (2 * k1 * (k1 + 1)/(length(y) - k1 - 1))
    AICc2 <- rep(NA, ncol(object$u))
    for (i in candidates) {
      lm2 <- lm(as.formula(paste(p1,p2,"+ V_",i,sep = "")), data = df1)
      k2 <- length(lm2$coef)
      AICc2[i] <- AIC(lm2) + (2 * k2 * (k2 + 1)/(length(y) - k2 - 1))
    }
    if (min(AICc2, na.rm = TRUE) < AICc1) {
      included <- c(included,candidates[candidates==which.min(AICc2)])
      candidates <- candidates[candidates!=which.min(AICc2)]
    } else {
      lm1$AICc <- AICc1
      return(lm1)
    }
  }
}
#' 
#' @describeIn lm-utils
#' 
#' Forward Stepwise Regression Sidak
#' 
#' Forward stepwise variable addition using a Sidak multiple testing
#' corrected alpha error threshold as the stopping criterion.
#' 
#' @export
lmforwardsequentialsidak <- function (y, x, object, alpha = 0.05) {
  AuxTrait <- !missing(x)## AuxTrait <- TRUE ## AuxTrait <- FALSE
  included <- numeric(0)
  candidates <- 1L:ncol(object$u)
  df1 <- if(AuxTrait) {x <- cbind(x) ; data.frame(x,object)} else object
  p1 <- if(AuxTrait) paste("y~",paste(colnames(x),collapse="+"),sep="") else "y~"
  while (TRUE) {
    p2 <- paste(if(length(included)) paste(if(AuxTrait)"+",paste(paste("V_", included, sep = ""), collapse = "+")) else if(AuxTrait) "" else "1" ,sep="")
    pval <- rep(NA, ncol(object$u))
    lm1 <- lm(as.formula(paste(p1,p2,sep="")), data = df1)
    for (i in candidates) {
      # i <- candidates[1L]
      lm2 <- lm(as.formula(paste(p1,p2,"+ V_",i,sep = "")), data = df1)
      aovcomp <- anova(lm1, lm2)
      pval[i] <- 1 - (1 - aovcomp[["Pr(>F)"]][2L])^(length(candidates) - length(included))
    }
    if (min(pval, na.rm = TRUE) < alpha) {
      included <- c(included, candidates[candidates == which.min(pval)])
      candidates <- candidates[candidates != which.min(pval)]
    } else {
      aovlm1 <- anova(lm1)
      lm1[["Familiwise"]] <- aovlm1[["Pr(>F)"]]
      idx <- match(colnames(object$u)[included],rownames(aovlm1))
      lm1[["Familiwise"]][idx] <- 1-(1-aovlm1[["Pr(>F)"]][idx])^(ncol(object$u):(ncol(object$u)-length(included)+1L))
      return(lm1)
    }
  }
}
##
