## **************************************************************************
##
##    (c) 2010-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Package MPSEM description**
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
#' \packageTitle{MPSEM}
#' 
#' @name MPSEM-package
#' 
#' @description \packageDescription{MPSEM}
#' 
#' @details Phylogenetic eignevector maps (PEM) is a method for using phylogeny
#' to model features of organism, most notably quantitative traits. It consists
#' in calculating sets of explanatory variables (eigenvectors) that are meant to
#' represent different patterns in trait values that are likely to have been
#' inducted by evolution. These patterns are used to model the data, using a
#' linear model for instance.
#' 
#' If one in interested in a \sQuote{target} species (i.e. a species for which
#' the trait value is unknown), and provided that we know the phylogenetic
#' relationships between that species and those of the model, the method allows
#' us to obtain the scores of that new species on the phylogenetic
#' eigenfunctions underlying a PEM. These scores are used to make empirical
#' predictions of trait values for the target species on the basis of those
#' observed for the species used in the model.
#' 
#' Functions \code{\link{PEM.build}}, \code{\link{PEM.updater}},
#' \code{\link{PEM.fitSimple}}, and \code{\link{PEM.forcedSimple}} allow one to
#' build, update (i.e. recalculate with alternative weighting parameters) as well
#' as to estimate or force arbitrary values for the weighting function
#' parameters.
#' 
#' Functions \code{\link{getGraphLocations}} and
#' \code{\link{Locations2PEMscores}} allow one to make predictions using method
#' \code{\link{predict.PEM}} and a linear model. To obtain this linear model,
#' one can use either function \code{\link{lm}} or auxiliary functions
#' \code{\link{lmforwardsequentialsidak}} or
#' \code{\link{lmforwardsequentialAICc}}, which perform forward-stepwise
#' variable addition on the basis of either familiwise type I error rate or the
#' Akaike Information Criterion (AIC), respectively.
#' 
#' The package provides low-level utility functions for performing operations on
#' graphs (see \link{graph-functions}), calculate influence matrix
#' (\code{\link{PEMInfluence}}), and simulate trait values (see
#' \link{trait-simulator}).
#' 
#' A phylogenetic modelling tutorial using \code{MPSEM} is available as a
#' package vignette. See example below.
#' 
#' The DESCRIPTION file:
#' \packageDESCRIPTION{MPSEM}
#' \packageIndices{MPSEM}
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution 4: 1120-1131
#' 
#' @seealso
#' Makarenkov, V., Legendre, P. & Desdevise, Y. 2004. Modelling phylogenetic
#' relationships using reticulated networks. Zoologica Scripta 33: 89-96
#' 
#' Blanchet, F. G., Legendre, P. & Borcard, D. 2008. Modelling directional
#' spatial processes in ecological data. Ecological Modelling 215: 325-336
#' 
#' @examples
#' ## To view MPSEM tutorial
#' vignette("MPSEM", package="MPSEM")
#'
NULL
##
