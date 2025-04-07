## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Trait Evolution Simulator **
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
#' @name traitEvolSim
#' 
#' @title Trait Evolution Simulator
#' 
#' @description Functions to simulate the evolution of traits as an
#' Ornstein-Uhlenbeck process, optionally with trait optimal value evolving
#' following a Markov process.
#' 
#' @param name Name of the trait whose evolution is to be created.
#' @param sigma Variance parameter of the trait (default: \code{sigma = 1}).
#' @param step Simulation step size (default: \code{step = 1}).
#' @param optima Trait optimal value(s) (default: \code{optima = NULL}).
#' @param alpha Trait selection rate (default: \code{alpha = 0}, meaning purely
#' neutral trait evolution)
#' @param transition A square, n-by-n, transition intensity matrix (default:
#' \code{transition = NULL}, which is suitable for single trait optima).
#' @param x A \code{QTEM-class} or \code{\link{graph-class}} object.
#' @param tem A \code{QTEM-class} object from function \code{traitEvolSim} or a
#' list of such objects (see details).
#' @param state One or more initial trait state(s) (default: 1).
#' @param value One or more initial trait value(s).
#' @param weighting A weighting function that determine the contribution of the
#' parent vertices as a function of the evolutionary distances (passed as its
#' first argument. It may have any number of argument, each with a default
#' value, and must accept arbitrary arguments (...). Default:
#' \code{\link{invDistWeighting}}.
#' @param ... further arguments to be passed to other functions and methods
#' 
#' @return Function \code{traitEvolSim} returns a QTEM-class object, which
#' consists of 8 member functions:
#' \describe{
#'   \item{ \code{$getName} }{It takes no argument and returns the name of the
#'   trait.}
#'   \item{ \code{$getStep} }{It takes no argument and returns the size of the
#'   simulation time step.}
#'   \item{ \code{$setStep} }{It takes the time step as an argument and sets,
#'   the simulation time step, including the recalculation of the transition
#'   intensity matrix in the case of a trait with multiple optima.}
#'   \item{ \code{$getOptima} }{It has no argument and returns the trait
#'   optimum values.}
#'   \item{ \code{$getTransition} }{It takes no argument and returns the
#'   transition density matrix.}
#'   \item{ \code{$getProb} }{It takes no argument and returns the transition
#'   probability matrix.}
#'   \item{ \code{$updateState} }{It takes the index of the actual trait
#'   optimum and returns the updated index (allows shifts in trait optimum).}
#'   \item{ \code{updateValue} }{It takes the value of the trait and,
#'   optionally, the index of the optimum state, at returns the updated trait
#'   value.}
#'   \item{ \code{$dumpConfig} }{It takes no argument and returns the trait
#'   evolution simulator parameters as a list.}
#'   \item{ \code{$print} }{It takes no argument and prints the parameters of
#'   the trait simulator on the screen (returns \code{NULL} invisibly).}
#' }
#' 
#' Function \code{simulateTrait} returns a two-column \code{\link{data.frame}}:
#' \describe{
#'   \item{ \code{state} }{The simulated index of the trait optimum at the
#'   vertices of the graph.}
#'   \item{ \code{value} }{The simulated trait value at the vertices of the
#'   graph.}
#' }
#' 
#' @details The trait evolution simulator is based on a random walk process of
#' one of categories: Weiner process and the Ornstein-Uhlenbeck process. The
#' Weiner process, also known as 'Brownian motion' is a random neutral process
#' with unbounded traits values. It is characterized using a single variance
#' parameter called \code{sigma} that dictates to what extent the trait value
#' fluctuates randomly in time.
#' 
#' The Ornstein-Uhlenbeck process is a random process with attractors (optima),
#' with values varying about these optima (or single optimum). It is notable
#' that whenever multiple optima are possible, only a single optimum is
#' effective at any given time. The Ornstein-Uhlenbeck process involve two more
#' values: the trait optimum effective at a particular time and the selection
#' rate (\code{alpha}). The optima represents the value toward which the
#' trait is attracted, whereas the value of \code{alpha} dictates to what
#' extent such an attraction occurs. When \code{alpha = 0}, the optimum stops
#' being an attractor and the Ornstein-Uhlenbeck process becomes purely a Weiner
#' process, whereas when \code{alpha} is very large, the trait value detracts
#' only slightly from its optimal value.
#' 
#' The trait simulation interface enabled to handle traits evolving neutrally
#' (pure Weiner process), and according to an Ornstein-Uhlenbeck process with
#' one or many optima. When many optima are used, transitions among the values
#' are simulated as a Markov process. For that purpose, the simulator needs to
#' be provided with a transition density matrix.
#' 
#' When multiple selection regimes are to be simulated, argument \code{tem} is
#' given a list of quantitative trait evolution models instead of a single
#' model. In that case, function \code{simulateTrait} will expect the graph
#' (argument \code{x}) to have a edge property called \sQuote{tem}. This
#' property contains the indices of the trait evolution models for each edge of
#' the graph. In practice, to have reasonable chances to resolve multiple
#' evolution regimes, the number of selection regimes has to be small with
#' respect to the number of vertices.
#' 
#' 
#' @author \packageAuthor{MPSEM}
#' 
#' @importFrom stats rnorm
#' @importFrom expm expm
#' 
NULL
#' 
#' @rdname traitEvolSim
#' 
#' @examples
#' ## Setting the RNG (for example consistency):
#' set.seed(2182955)
#' 
#' ## A list to contain the trait simulator objects:
#' tem <- list()
#' 
#' ## The first simulator is non-neutral with three optima:
#' traitEvolSim(
#'   name = "Trait 1",
#'   sigma = 1.5,
#'   alpha = 0.15,
#'   optima = c(30,50,80),
#'   transition = matrix(c(NA,0.1,0.0,0.1,NA,0.1,0.0,0.1,NA), 3L, 3L)
#' ) -> tem[[1]]
#' tem[[1]]
#' 
#' ## The second simulator is neutral:
#' traitEvolSim(
#'   name = "Trait 2",
#'   sigma = 2.5,
#'   step = 1
#' ) -> tem[[2]]
#' tem[[2]]
#' 
#' ## The third simulator is non-neutral with a single optimum:
#' traitEvolSim(
#'   name = "Trait 3",
#'   alpha = 0.05,
#'   optima = 15
#' ) -> tem[[3]]
#' tem[[3]]
#' 
#' ## The fourth simulator is also non-neutral with a single optimum (but having
#' ## a different value than the latter):
#' traitEvolSim(
#'   name = "Trait 4",
#'   sigma = 2.0,
#'   alpha = 0.25,
#'   optima = -25
#' ) -> tem[[4]]
#' tem[[4]]
#' 
#' ## Each simulator has a set of embedded member functions that can be called
#' ## directly.
#' 
#' ## Member $getName() returns the name of the trait as follows:
#' unlist(lapply(tem, function(x) x$getName()))
#' 
#' ## Member $getStep() returns the step sizes for the traits as follows:
#' unlist(lapply(tem, function(x) x$getStep()))
#' 
#' ## Member $setStep() sets the step sizes as follows:
#' lapply(tem, function(x) x$setStep(0.1))
#' 
#' ## and returns the recalculated transition probability matrix of simulators
#' ## having a transition intensity matrix (i.e., for non-neutral simulators
#' ## with multiple trait optima).
#' 
#' ## This is the modified step sizes:
#' unlist(lapply(tem, function(x) x$getStep()))
#' 
#' ## Member $getOptima() returns the simulator's optimum (if available) or a
#' ## vector of optima (if multiple optima are available) as follows:
#' lapply(tem, function(x) x$getOptima())
#' 
#' ## Member $getTransition() returns the transition intensity matrix, when
#' ## available (NULL otherwise) as follows:
#' lapply(tem, function(x) x$getTransition())
#' 
#' ## Member $getProb() returns the transition intensity matrix, whenever
#' ## available (NULL otherwise) as follows:
#' lapply(tem, function(x) x$getProb())
#' 
#' ## When multiple optima are available, member $updateState() enables to
#' ## simulate the transition from one optimal trait state (the one given as the
#' ## argument) to another (the one which is returned by the function):
#' state <- 1
#' newstate <- tem[[1]]$updateState(state)
#' newstate
#' 
#' ## Member $updateValue() simulates the evolution of the trait from one value
#' ## to another. For a non-neutral with multiple optima, The trait state is
#' ## provided using argument state (default: 1, which is the only applicable
#' ## value for a single optimum).
#' oldvalue <- 31.5
#' newvalue <- tem[[1]]$updateValue(oldvalue, state=1)
#' newvalue
#' 
#' ## Member $dumpConfig() returns the configuration list, which can be used
#' 
#' cfg <- lapply(tem, function(x) x$dumpConfig())
#' 
#' ## The trait evolution simulators can be re-instantiated as follows:
#' lapply(
#'   cfg,
#'   function(x)
#'     traitEvolSim(
#'       name = x$name,
#'       sigma = x$sigma,
#'       step = x$step,
#'       alpha = x$alpha,
#'       optima = x$optima,
#'       transition = x$transition
#'     )
#' ) -> tem_clone
#' 
#' ## Clean up:
#' rm(cfg, tem_clone)
#' 
#' ## Simulate trait evolution using the four simulators described previously:
#' 
#' ## Set step size to 0.05
#' lapply(tem, function(x, a) x$setStep(a), a = 0.05)
#' 
#' ## Results list:
#' res <- NULL
#' trNms <- lapply(tem, function(x) x$getName())
#' res$state <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
#' res$optim <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
#' res$trait <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
#' rm(trNms)
#' 
#' ## res$state contains the trait state at the simulation time
#' ## res$optim contains the trait's optimum value at the simulation time
#' ## res$trait contains the trait value at the simulation time
#' 
#' ## Setting the optimal state of the four traits at the beginning of the
#' ## simulation period:
#' res$state[1L,] <- c(2,NA,1,1)  ## NBL trait #2 evolves neutrally.
#' 
#' ## Getting the trait optima at the beginning of the simulation period:
#' unlist(
#'   lapply(
#'     tem,
#'     function(x)
#'       x$getOptima()[res$state[1L,x$getName()]]
#'   )
#' ) -> res$optim[1L,]
#' 
#' ## Setting the initial trait values:
#' res$trait[1L,] <- c(50,0,15,-25)
#' 
#' ## The state of the simulation at the beginning of the simulation:
#' head(res$state)
#' head(res$optim)
#' head(res$trait)
#' 
#' ## Setting RNG state to obtain 
#' set.seed(1234567)
#' 
#' ## This loop simulates time steps #2 through #1001
#' for(i in 2L:1001L) {
#'   
#'   ## Simulate the evolution of the trait states (if relevant, which it is
#'   ## only for the first trait):
#'   unlist(
#'     lapply(
#'       tem,
#'       function(x)
#'         x$updateState(res$state[i - 1L,x$getName()])
#'     ) 
#'   )-> res$state[i,]
#'   
#'   ## Obtain the optimal trait value (relevant for all traits but the second,
#'   ## trait, which evolves neutrally):
#'   unlist(
#'     lapply(
#'       tem,
#'       function(x)
#'         x$getOptima()[res$state[i,x$getName()]]
#'     )
#'   ) -> res$optim[i,]
#'   
#'   ## Simulate the evolution of the trait value:
#'   unlist(
#'     lapply(
#'       tem,
#'       function(x)
#'         x$updateValue(
#'           state = res$state[i,x$getName()],
#'           value = res$trait[i - 1L,x$getName()]
#'         )
#'     )
#'   ) -> res$trait[i,]
#' }
#' 
#' ## Plot the results:
#' par(mar=c(4,4,1,1))
#' plot(NA, xlim=c(0,0.05*1000), ylim=range(res$trait), xlab="Time",
#'      ylab="Trait value", las=1L)
#' lines(x=0.05*(0:1000), y=res$trait[,1L], col="black")
#' if(!is.na(tem[[1L]]$getOptima())[1L])
#'   lines(x=0.05*(0:1000), y=res$optim[,1L], col="black", lty=3L)
#' lines(x=0.05*(0:1000), y=res$trait[,2L], col="red")
#' if(!is.na(tem[[2L]]$getOptima())[1L])
#'   lines(x=0.05*(0:1000), y=res$optim[,2L], col="red", lty=3L)
#' lines(x=0.05*(0:1000), y=res$trait[,3L], col="blue")
#' if(!is.na(tem[[3L]]$getOptima())[1L])
#'   lines(x=0.05*(0:1000), y=res$optim[,3L], col="blue", lty=3L)
#' lines(x=0.05*(0:1000), y=res$trait[,4L], col="green")
#' if(!is.na(tem[[4L]]$getOptima())[1L])
#'   lines(x=0.05*(0:1000), y=res$optim[,4L], col="green", lty=3L)
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
#' ) -> gr_lin
#' 
#' ## Simulate a trait (Ornstein Uhlenbeck):
#' simulateTrait(
#'   x = gr_lin,
#'   tem = tem[[1]],
#'   state = 2,
#'   value = 50,
#'   a = 1
#' ) -> simTrait
#' 
#' ## Showing the trait values with the optima (OU process):
#' x <- cumsum(c(0,attr(gr_lin,"edge")$distance))
#' plot(x=x, y=simTrait$value, type="l", las=1)
#' lines(x=x, y=c(30,50,80)[simTrait$state], lty=3)
#' 
#' ## Simulate an other trait (Brownian motion):
#' simulateTrait(
#'   x = gr_lin,
#'   tem = tem[[2]],
#'   value = 10,
#'   a = 1
#' ) -> simTrait
#' 
#' ## Showing the trait values:
#' plot(x=x, y=simTrait$value, type="l", las=1)
#' 
#' ## A distance-based network:
#' N <- 100
#' coords <- cbind(x=runif(N,-1,1), y=runif(N,-1,1))
#' rownames(coords) <- sprintf("N%d",1:N)
#' dst <- dist(coords)
#' gr_dst <- dstGraph(d=dst, th=0.35, origin=15)
#' 
#' ## Simulate a trait (Brownian motion) on the distance-based network:
#' simulateTrait(
#'   x = gr_dst,
#'   tem = tem[[2]],
#'   value = 0,
#'   a = 1
#' ) -> simTrait
#' 
#' ## Showing the results of the simulation:
#' gp <- par(no.readonly = TRUE)
#' par(mar=c(5,5,1,1))
#' 
#' plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1, las=1, xlab="X", ylab="Y")
#' points(x=coords[,"x"], y=coords[,"y"], pch=21, cex=abs(simTrait$value)/3,
#'        bg=gray(0.5 - 0.5*sign(simTrait$value)))
#' 
#' par(gp)
#' 
#' @export
traitEvolSim <- function(name, sigma = 1, step = 1, alpha = 0, optima = NULL,
                         transition = NULL) {
  
  if(sigma <= 0)
    stop("Argument 'sigma' value must be non-zero and positive!")
  
  if(step < 0)
    stop("Argument 'step' value must be positive!")
  
  if(alpha < 0)
    stop("Argument 'alpha' value must be positive!")
  
  if(is.null(optima) && (alpha != 0))
    stop("No optimal value provided for a non-neutral process!")
  
  nstate <- length(optima)
  
  if(nstate > 1L) {
    
    if(is.null(transition) || !is.matrix(transition))
      stop("Many optima provided without a transition intensity matrix!")
    
    if((NROW(transition) != length(optima)) ||
       (NROW(transition) != NCOL(transition)))
      stop("Inadequate transition matrix: ", NROW(transition), " by ",
           NCOL(transition), " for ", length(optima), " optima.")
    
    diag(transition) <- 0
    diag(transition) <- -rowSums(transition)
    prob <- expm(transition)
  } else
    prob <- NULL
  
  updateState <- function(state) {
    if(missing(state))
      stop("An initial state must be provided!")
    if(nstate > 1L)
      state <- sample(x=nstate, size=1L, prob=prob[state,])
    state
  }
  
  updateValue <- function(value, state = 1L) {
    if(missing(value))
      stop("No initial value!")
    if(alpha) {
      w <- exp(-alpha*step)
      s <- sigma*sqrt((1 - exp(-2*alpha*step))/(2*alpha))
      value <- rnorm(1L, w*value + (1 - w)*optima[state], s)
    } else {
      s <- sigma*sqrt(step)
      value <- rnorm(1L, value, s)
    }
    value
  }
  
  structure(
    list(
      getName = function() name,
      getStep = function() step,
      setStep = function(step) {
        step <<- step
        if(!is.null(transition))
          prob <<- expm(transition)
      },
      getOptima = function() if(is.null(optima)) NA else optima,
      getTransition = function() transition,
      getProb = function() prob,
      updateState = updateState,
      updateValue = updateValue,
      dumpConfig = function()
        list(
          name = name,
          sigma = sigma,
          step = step,
          alpha = alpha,
          optima = optima,
          transition = transition
        ),
      print = function() {
        cat("\nQuantitative trait evolution simulator (Ornstein-Uhlenbeck)",
            "\n----------------------------------------------\n")
        cat("Name:", name, "\n")
        if(alpha) {
          cat("Trait evolving non-neutrally (alpha = ",alpha,")\n",sep="")
        } else
          cat("Trait evolving neutrally (alpha = 0)\n")
        cat("Sigma = ",sigma,"\n",sep="")
        if(!is.null(optima))
          if(nstate > 1L) {
            cat("Trait optima:",paste(optima,collapse=", "),"\n")
          } else
            cat("Trait optimum:",optima,"\n")
        if(!is.null(transition)) {
          cat("Transition intensity matrix:\n")
          print(transition)
        }
        cat("\n")
        invisible(NULL)
      }
    ),
    class = "QTEM"
  )
}
#' 
#' @rdname traitEvolSim
#' 
#' @method print QTEM
#' 
#' @export
print.QTEM <- function(x, ...)
  x$print()
#' 
#' @rdname traitEvolSim
#' 
#' @export
simulateTrait <- function(x, tem, state, value, weighting = invDistWeighting,
                          ...) {
  
  if(inherits(tem, "QTEM")) {
    mtem <- FALSE
  } else {
    if(!all(unlist(lapply(tem, inherits, what="QTEM"))))
      stop("Argument 'tem' is neither a quantitative trait evolution model ",
           "nor a list thereof")
    if(is.null(edge(x)$tem)) {
      warning("Multiple trait evolution models are available, but the graph ",
              "has no edge property called 'tem'. Only the first element of ",
              "the list of trait evolution models will be used.")
      tem <- tem[[1L]]
      mtem <- FALSE
    } else {
      if(min(edge(x)$tem) < 1)
        stop("Indices of the trait evolution model were smaller than 1")
      if(max(edge(x)$tem) < length(tem))
        stop("Indices of the trait evolution model were larger than the ",
             "number of models in 'tem' (", length(tem), " models in argument ",
             "'tem', but ", max(edge(x)$tem)," indices found in edge ",
             "property 'tem')")
      mtem <- TRUE
    }
  }
  
  origin <- getOrigin(x)
  norig <- length(origin)
  
  if(missing(state)) {
    state <- rep(1L, norig)
  } else 
    if(norig != length(state)) {
      warning("The graph has ", norig, " origins, but 'state' is length ",
              length(state), ".")
      state <- rep(state, length.out = norig)
    }
  
  if(norig != length(value))
    stop("The graph has ", norig, " origins, but 'value' is length ",
         length(value), ".")
  
  out <- list(state = integer(nrow(x)), value = numeric(nrow(x)))
  out$state[origin] <- state
  out$value[origin] <- value
  
  ord <- attr(x,"processOrder")
  if(is.null(ord))
    ord <- getProcessOrder(x)
  
  edge <- attr(x,"edge")
  
  for(k in ord[(norig + 1L):nrow(x)]) {
    
    av <- edge[[1L]][edge[[2L]] == k]
    da <- edge$distance[edge[[2L]] == k]
    ctb <- weighting(da, ...)
    s <- sample(length(ctb), 1L, FALSE, ctb)
    i <- av[s]
    
    if(mtem) {
      
      if(da[s] != tem[[edge$tem[s]]]$getStep())
        tem[[edge$tem[s]]]$setStep(da[s])
      
      mergedValue <- sum(ctb*out$value[av])
      out$state[k] <- tem[[edge$tem[s]]]$updateState(out$state[i])
      out$value[k] <- tem[[edge$tem[s]]]$updateValue(mergedValue, out$state[k])
    } else {
      
      if(da[s] != tem$getStep())
        tem$setStep(da[s])
      
      mergedValue <- sum(ctb*out$value[av])
      out$state[k] <- tem$updateState(out$state[i])
      out$value[k] <- tem$updateValue(mergedValue, out$state[k])
    }
  }
  
  data.frame(out, row.names=rownames(x))
}
##
