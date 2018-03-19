gamma_dist_prob <- function(array,mean,sd) {
    theta <- sd * sd / mean
    k <- mean / theta
    return(1.0 - (1.0-pgamma(array+1,k,scale=theta))/(1.0-pgamma(array,k,scale=theta)))
}

#' @name spop
#' @export spop
#' @importFrom methods new
#' @importFrom stats pgamma rbinom
#' @title An S4 class to represent an age-structured population
#'
#' @description
#' \itemize{
#' \item \code{spop} implements the deterministic and stochastic age-structured population dynamics models described in Erguler et al. 2016 and 2017
#' \item \code{add} introduces a batch of individuals with a given age (default: 0)
#' \item \code{iterate} iterates the population for one day and calculates (overwrites) the number of dead individuals and the number of individuals designated to complete their development
#' \item \code{developed} reads the number of individuals designated to complete their development
#' \item \code{dead} reads the number of dead individuals after each iteration
#' \item \code{size} reads the total number of individuals
#' }
#'
#' @details
#' This is an R implementation of the age-structured population dynamics models described in Erguler et al. 2016 and 2017. The \code{spop} class records the number and age of individuals and implements two processes to exit from the population: development and death. The two processes act upon the population sequentially; survival is imposed prior to development. If the population survives for one day, then, it is allowed to grow and complete its development. Survival and development are defined either with an age-independent daily probability, or an age-dependent gamma-distributed probability.
#' \itemize{
#' \item \code{death} defines an age-independent daily probability of death
#' \item \code{death_mean} and \code{death_sd} define an age-dependent daily probability of death. \code{death_sd=0} indicates fixed life time of \code{death_mean}
#' \item \code{dev} defines an age-independent daily probability of development
#' \item \code{dev_mean} and \code{dev_sd} define an age-dependent daily probability of development. \code{dev_sd=0} indicates fixed development time of \code{dev_mean}
#' }
#' 
#' @examples
#' # Generate a population with stochastic dynamics
#' s <- spop(stochastic=TRUE)
#' # Add 1000 20-day-old individuals
#' add(s) <- data.frame(number=1000,age=20)
#' 
#' # Iterate one day without death and assume development in 20 (+-5) days
#' iterate(s) <- data.frame(dev_mean=20,dev_sd=5,death=0)
#' print(developed(s))
#'
#' # Iterate another day assuming no development but age-dependent survival
#' # Let each individual survive for 20 days (+-5)
#' iterate(s) <- data.frame(death_mean=20,death_sd=5,dev=0)
#' print(dead(s))
#' # Note that the previous values of developed and dead will be overwritten by this command
#' 
#' # Generate a deterministic population and observe the difference
#' s <- spop(stochastic=FALSE)
#' add(s) <- data.frame(number=1000,age=20)
#'
#' iterate(s) <- data.frame(dev_mean=20,dev_sd=5,death=0)
#' print(developed(s))
#'
#' iterate(s) <- data.frame(death_mean=20,death_sd=5,dev=0)
#' print(dead(s))
#'
spop <- setClass("spop", slots=list(stochastic="logical",
                                    pop="data.frame",
                                    developed="numeric",
                                    dead="numeric"))
setMethod("initialize",
          "spop",
          function(.Object, stochastic=TRUE) {
              .Object@stochastic <- stochastic
              .Object@pop <- data.frame(age=numeric(),number=numeric())
              .Object@developed <- 0
              .Object@dead <- 0
              return(.Object)
          })

#' @name add<-
#' @rdname add-spop
#' @exportMethod add<-
#' @title Add batch
#' @description
#' Introduce a batch of individuals with a given age
#' @param x spop class instant
#' @param value \code{data.frame} with \code{name} and \code{age} fields
setGeneric(name="add<-",
           def=function(x,value) standardGeneric("add<-"))
#' @rdname add-spop
#' @aliases add<-,spop-method
setMethod("add<-",
          c("spop","data.frame"),
          function(x, value) {
              if (!("age" %in% colnames(value))) value$age <- 0
              if (!("number" %in% colnames(value))) {
                  return(x)
              }
              i <- x@pop$age == value$age
              if (any(i))
                  x@pop$number[i] <- x@pop$number[i] + value$number
              else
                  x@pop[nrow(x@pop)+1,] <- list(value$age,value$number)
              return(x)
          })

#' @name iterate<-
#' @rdname iterate-spop
#' @exportMethod iterate<-
#' @title Iterate population
#' @description
#' Iterate the population for one day
#' @param x spop class instant
#' @param value \code{data.frame} with the following fields
#' \itemize{
#' \item \code{dev}
#' \item \code{dev_mean}, \code{dev_sd}
#' \item \code{death}
#' \item \code{death_mean}, \code{death_sd}
#' }
setGeneric(name="iterate<-",
           def=function(x,value) standardGeneric("iterate<-"))
#' @rdname iterate-spop
#' @aliases iterate<-,spop-method
setMethod("iterate<-",
          c("spop","data.frame"),
          function(x, value) {
              if (nrow(x@pop)==0) {
                  x@developed <- 0
                  x@dead <- 0
                  return(x)
              }
              if (!("dev" %in% colnames(value))) {
                  if (("dev_mean" %in% colnames(value)) && ("dev_sd" %in% colnames(value)))
                      if (value$dev_sd == 0)
                          dev <- as.numeric(x@pop$age >= value$dev_mean - 1.0)
                      else
                          dev <- gamma_dist_prob(x@pop$age,value$dev_mean,value$dev_sd)
                  else {
                      warning(sprintf("Error in development probability"))
                      return(x)
                  }
              } else
                  dev <- value$dev
              if (!("death" %in% colnames(value))) {
                  if (("death_mean" %in% colnames(value)) && ("death_sd" %in% colnames(value)))
                      if (value$death_sd == 0)
                          death <- as.numeric(x@pop$age >= value$death_mean - 1.0)
                      else
                          death <- gamma_dist_prob(x@pop$age,value$death_mean,value$death_sd)
                  else {
                      warning(sprintf("Error in probability of death"))
                      return(x)
                  }
              } else
                  death <- value$death
              if (x@stochastic) {
                  k <- rbinom(length(x@pop$number),x@pop$number,death)
                  x@pop$number <- x@pop$number - k
                  d <- rbinom(length(x@pop$number),x@pop$number,dev)
                  x@pop$number <- x@pop$number - d
              } else {
                  k <- x@pop$number * death
                  x@pop$number <- x@pop$number * (1.0 - death)
                  d <- x@pop$number * dev
                  x@pop$number <- x@pop$number * (1.0 - dev)
              }
              x@pop$age <- x@pop$age + 1
              x@pop <- x@pop[x@pop$number>0,]
              x@developed <- sum(d)
              x@dead <- sum(k)
              return(x)
          })

#' @name developed
#' @rdname developed-spop
#' @exportMethod developed
#' @title Read developed
#' @description
#' Read the number of individuals designated to complete their development
#' @param x spop class instant
setGeneric(name="developed",
           def=function(x) standardGeneric("developed"))
#' @rdname developed-spop
#' @aliases developed,spop-method
setMethod("developed",
          "spop",
          function(x) {
              return(x@developed)
          })

#' @name dead
#' @rdname dead-spop
#' @exportMethod dead
#' @title Read dead
#' @description
#' Read the number of dead individuals after each iteration
#' @param x spop class instant
setGeneric(name="dead",
           def=function(x) standardGeneric("dead"))
#' @rdname dead-spop
#' @aliases dead,spop-method
setMethod("dead",
          "spop",
          function(x) {
              return(x@dead)
          })

#' @name size
#' @rdname size-spop
#' @exportMethod size
#' @title Read size
#' @description
#' Read the total number of individuals
#' @param x spop class instant
setGeneric(name="size",
           def=function(x) standardGeneric("size"))
#' @rdname size-spop
#' @aliases size,spop-method
setMethod("size",
          "spop",
          function(x) {
              return(sum(x@pop$number))
          })
