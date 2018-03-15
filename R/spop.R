gamma_dist_prob <- function(array,mean,sd) {
    theta <- sd * sd / mean
    k <- mean / theta
    return(1.0 - (1.0-pgamma(array+1,k,scale=theta))/(1.0-pgamma(array,k,scale=theta)))
}

#' An S4 class to represent an age-structured population
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
#' \dontrun{
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
#' }
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

setGeneric(name="add<-",
           def=function(object,value) standardGeneric("add<-"))
setMethod("add<-",
          c("spop","data.frame"),
          function(object, value) {
              if (!("age" %in% colnames(value))) value$age <- 0
              if (!("number" %in% colnames(value))) {
                  return(object)
              }
              i <- object@pop$age == value$age
              if (any(i))
                  object@pop$number[i] <- object@pop$number[i] + value$number
              else
                  object@pop[nrow(object@pop)+1,] <- list(value$age,value$number)
              return(object)
          })
setGeneric(name="iterate<-",
           def=function(object,value) standardGeneric("iterate<-"))
setMethod("iterate<-",
          c("spop","data.frame"),
          function(object, value) {
              if (nrow(object@pop)==0) {
                  object@developed <- 0
                  object@dead <- 0
                  return(object)
              }
              if (!("dev" %in% colnames(value))) {
                  if (("dev_mean" %in% colnames(value)) && ("dev_sd" %in% colnames(value)))
                      if (value$dev_sd == 0)
                          dev <- as.numeric(object@pop$age >= value$dev_mean - 1.0)
                      else
                          dev <- gamma_dist_prob(object@pop$age,value$dev_mean,value$dev_sd)
                  else {
                      warning(sprintf("Error in development probability"))
                      return(object)
                  }
              } else
                  dev <- value$dev
              if (!("death" %in% colnames(value))) {
                  if (("death_mean" %in% colnames(value)) && ("death_sd" %in% colnames(value)))
                      if (value$death_sd == 0)
                          death <- as.numeric(object@pop$age >= value$death_mean - 1.0)
                      else
                          death <- gamma_dist_prob(object@pop$age,value$death_mean,value$death_sd)
                  else {
                      warning(sprintf("Error in probability of death"))
                      return(object)
                  }
              } else
                  death <- value$death
              if (object@stochastic) {
                  k <- rbinom(length(object@pop$number),object@pop$number,death)
                  object@pop$number <- object@pop$number - k
                  d <- rbinom(length(object@pop$number),object@pop$number,dev)
                  object@pop$number <- object@pop$number - d
              } else {
                  k <- object@pop$number * death
                  object@pop$number <- object@pop$number * (1.0 - death)
                  d <- object@pop$number * dev
                  object@pop$number <- object@pop$number * (1.0 - dev)
              }
              object@pop$age <- object@pop$age + 1
              object@pop <- object@pop[object@pop$number>0,]
              object@developed <- sum(d)
              object@dead <- sum(k)
              return(object)
          })
setGeneric(name="developed",
           def=function(object) standardGeneric("developed"))
setMethod("developed",
          "spop",
          function(object) {
              return(object@developed)
          })
setGeneric(name="dead",
           def=function(object) standardGeneric("dead"))
setMethod("dead",
          "spop",
          function(object) {
              return(object@dead)
          })
setGeneric(name="size",
           def=function(object) standardGeneric("size"))
setMethod("size",
          "spop",
          function(object) {
              return(sum(object@pop$number))
          })
