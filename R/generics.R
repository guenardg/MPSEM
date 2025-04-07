## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Method Definition **
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
#' Generic Function (Methods)
#' 
#' @description A set of non-standard generic functions used by package MPSEM.
#' 
#' @name MPSEM-generics
#' 
#' @param object A \code{\link{graph-class}} object or an object of any class
#' that can be converted into a \code{\link{graph-class}} object.
#' @param x A \code{\link{graph-class}} object or an object of any class
#' that can be converted into a \code{\link{graph-class}} object, or a numeric
#' vector or numeric matrix of auxiliary trait(s) to help in estimating the
#' parameters of the trait evolution model. In the second case, the default is
#' \code{NULL}, which entails that no auxiliary trait is used.
#' @param value A two-column \code{\link{data.frame}} giving the coordinates
#' (e.g., x, y, z) to be applied to the \code{graph-class} object.
#' @param target A numeric or character vector specifying one or more target
#' species in graph \code{x}.
#' @param ... Further arguments of the particular conversion method that is
#' being implemented.
#' @param y A numeric vector or numeric matrix of trait(s) against which to
#' estimate the parameters of the trait evolution model.
#' 
#' @return The return value depends on the implementation of the methods for a
#' particular S3 class.
#' 
NULL
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Number of Edges
#' 
#' A method to access the number of edges
#' 
#' @export
nedge <- function(x) UseMethod("nedge")
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Edge Extraction
#' 
#' A method to extract edges from an object.
#' 
#' @export
edge <- function(x) UseMethod("edge")
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Edge Assignment
#' 
#' A method to assign edges to an object.
#' 
#' @export
`edge<-` <- function(x, value)  UseMethod("edge<-")
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Edge Names Extraction
#' 
#' A method to extract edge names from an object.
#' 
#' @export
edgenames <- function(x) UseMethod("edgenames")
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Edge Assignment
#' 
#' A method to assign edges to an object.
#' 
#' @export
`edgenames<-` <- function(x, value)  UseMethod("edgenames<-")
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Geometry Extraction
#' 
#' A method to obtain an object's geometry.
#' 
#' @export
geometry <- function(x)  UseMethod("geometry")
#' 
#' @describeIn MPSEM-generics
#' 
#' Geometry Assignment
#' 
#' A method to assign a geometry to an object.
#' 
#' @export
`geometry<-` <- function(x, value)  UseMethod("geometry<-")
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Labels Assignment
#' 
#' A method to assign labels to an object.
#' 
#' @export
`labels<-` <- function(object, value)  UseMethod("labels<-")
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Graph-class Conversion
#' 
#' A method to convert an object into a graph-class object.
#' 
#' @export
as.graph <- function(x, ...)  UseMethod("as.graph")
#' 
#' 
#' @describeIn MPSEM-generics
#' 
#' Locate Targets
#' 
#' A method to locate a set of targets species from a graph-class object.
#' 
#' @export
locate <- function(x, target, ...) UseMethod("locate")
#'
#' 
#' @describeIn MPSEM-generics
#' 
#' Estimate Evolution Model
#' 
#' Estimate the parameters of an evolution model on a PEM2-class object.
#' 
#' @export
evolution.model <- function(object, y, ..., x = NULL)
  UseMethod("evolution.model")
#'
