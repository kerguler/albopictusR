prob_list <- c('gamma','nbinom')

#' @name nbinom_dist_prob
#' @export nbinom_dist_prob
#' @importFrom stats pnbinom
#' @title Negative binomial probability of death or development
#' @description
#' Negative binomial-distributed probability of death or development happening in the next iteration
#' @param xr age of individuals
#' @param mn mean age of death or development
#' @param std standard deviation of the age of death or development
nbinom_dist_prob <- function(xr,mn,std) {
    p <- mn / (std * std)
    r <- mn * p / (1.0 - p)
    return(1.0 - (1.0-pnbinom(xr+1,size=r,prob=p))/(1.0-pnbinom(xr,size=r,prob=p)))
}

#' @name gamma_dist_prob
#' @export gamma_dist_prob
#' @importFrom stats pgamma
#' @title Gamma probability of death or development
#' @description
#' Gamma-distributed probability of death or development happening in the next iteration
#' @param xr age of individuals
#' @param mn mean age of death or development
#' @param std standard deviation of the age of death or development
gamma_dist_prob <- function(xr,mn,std) {
    theta <- std * std / mn
    k <- mn / theta
    return(1.0 - (1.0-pgamma(xr+1,k,scale=theta))/(1.0-pgamma(xr,k,scale=theta)))
}

#' @name spop
#' @export spop
#' @importFrom methods new
#' @importFrom stats pgamma rbinom pnbinom
#' @title An S4 class to represent an age-structured population
#'
#' @description
#' \itemize{
#' \item \code{spop} implements the deterministic and stochastic age-structured population dynamics models described in Erguler et al. 2016 and 2017
#' \item \code{add} introduces a batch of individuals with a given age, completed development cycles, and degree of development (default: 0)
#' \item \code{iterate} iterates the population for one day and calculates (overwrites) the number of dead individuals and the number of individuals designated to complete their development
#' \item \code{devtable} reads the number, age, and development cycle of individuals designated to complete their development
#' \item \code{developed} reads the total number of individuals designated to complete their development
#' \item \code{dead} reads the number of dead individuals after each iteration
#' \item \code{size} reads the total number of individuals
#' }
#'
#' @details
#' This is an R implementation of the age-structured population dynamics models described in Erguler et al. 2016 and 2017. The \code{spop} class records the number and age of individuals and implements two processes to exit from the population: development and death. The two processes act upon the population sequentially; survival is imposed prior to development. If the population survives for one day, then, it is allowed to grow and complete its development. Survival and development are defined either with an age-independent daily probability, or an age-dependent gamma- or negative binomial-distributed probability.
#' \itemize{
#' \item \code{stochastic}: a logical value indicating a deterministic or a stochastic population dynamics
#' \item \code{prob}: a character string indicating the basis of age-dependent survival or development (gamma: gamma-distributed, nbinom: negative binomial-distributed)
#' }
#' 
#' @examples
#' # Generate a population with stochastic dynamics
#' s <- spop(stochastic=TRUE)
#' # Add 1000 20-day-old individuals
#' add(s) <- data.frame(number=1000,age=20)
#' 
#' # Iterate one day without death and assume development in 20 (+-5) days (gamma-distributed)
#' iterate(s) <- data.frame(dev_mean=20,dev_sd=5,death=0)
#' print(developed(s))
#'
#' # Iterate another day assuming no development but age-dependent survival
#' # Let each individual survive for 20 days (+-5) (gamma-distributed)
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
                                    prob="character",
                                    pop="data.frame",
                                    developed="numeric",
                                    dead="numeric",
                                    devtable="data.frame",
                                    prob_fun="function"))
setMethod("initialize",
          "spop",
          function(.Object, stochastic=TRUE, prob='gamma') {
              .Object@prob <- prob
              if (.Object@prob == "nbinom")
                  .Object@prob_fun <- nbinom_dist_prob
              else if (.Object@prob == "gamma")
                  .Object@prob_fun <- gamma_dist_prob
              else {
                  cat(sprintf("ERROR: prob can be one of the following:\n"))
                  cat(sprintf("-> %s\n",paste(prob_list,collapse="\n-> ")))
                  return(NULL)
              }                  
              .Object@stochastic <- stochastic
              .Object@pop <- data.frame(age=numeric(),devcycle=numeric(),development=numeric(),number=numeric())
              .Object@developed <- 0
              .Object@dead <- 0
              .Object@devtable <- data.frame()
              return(.Object)
          })

#' @name add<-
#' @rdname add-spop
#' @exportMethod add<-
#' @title Add batch
#' @description
#' Introduce a batch of individuals with a given age
#' @param x spop class instant
#' @param value \code{data.frame} with \code{age}, \code{devcycle}, \code{development}, and \code{number} fields
setGeneric(name="add<-",
           def=function(x,value) standardGeneric("add<-"))
#' @rdname add-spop
#' @aliases add<-,spop-method
setMethod("add<-",
          c("spop","data.frame"),
          function(x, value) {
              if (nrow(value)==0) return(x);
              if (!("age" %in% colnames(value))) value$age <- 0
              if (!("devcycle" %in% colnames(value))) value$devcycle <- 0
              if (!("development" %in% colnames(value))) value$development <- value$age
              if (!("number" %in% colnames(value)) || value$number <= 0) {
                  return(x)
              }
              for (val in 1:nrow(value)) {
                  i <- (x@pop$age == value[val,]$age) & (x@pop$devcycle == value[val,]$devcycle) & (x@pop$development == value[val,]$development)
                  if (any(i))
                      x@pop$number[i] <- x@pop$number[i] + value[val,]$number
                  else
                      x@pop[nrow(x@pop)+1,] <- list(value[val,]$age,
                                                    value[val,]$devcycle,
                                                    value[val,]$development,
                                                    value[val,]$number)
              }
              return(x)
          })

iterfun <- function(x, value, pause) {
    if (nrow(x@pop)==0) {
        x@developed <- 0
        x@dead <- 0
        x@devtable <- data.frame()
        return(x)
    }
    if (!("dev" %in% colnames(value))) {
        if (("dev_mean" %in% colnames(value)) && (!("dev_sd" %in% colnames(value)) || (("dev_sd" %in% colnames(value)) && value$dev_sd == 0)))
            dev <- as.numeric(x@pop$development >= value$dev_mean - 1.0)
        else if (("dev_mean" %in% colnames(value)) && (("dev_sd" %in% colnames(value)) && (value$dev_sd > 0)))
            dev <- x@prob_fun(x@pop$development,value$dev_mean,value$dev_sd)
        else {
            warning(sprintf("Error in development probability"))
            return(x)
        }
    } else
        dev <- value$dev
    dev[!is.finite(dev)] <- 1.0
    if (!("death" %in% colnames(value))) {
        if (("death_mean" %in% colnames(value)) && (!("death_sd" %in% colnames(value)) || (("death_sd" %in% colnames(value)) && value$death_sd == 0)))
            death <- as.numeric(x@pop$age >= value$death_mean - 1.0)
        else if (("death_mean" %in% colnames(value)) && (("death_sd" %in% colnames(value)) && (value$death_sd > 0)))
            death <- x@prob_fun(x@pop$age,value$death_mean,value$death_sd)
        else {
            warning(sprintf("Error in probability of death"))
            return(x)
        }
    } else
        death <- value$death
    death[!is.finite(death)] <- 1.0
    if (x@stochastic) {
        k <- rbinom(length(x@pop$number),x@pop$number,death)
        x@pop$number <- x@pop$number - k
        d <- rbinom(length(x@pop$number),x@pop$number,dev)
        x@pop$number <- x@pop$number - d
    } else {
        k <- x@pop$number * death
        x@pop$number <- x@pop$number - k
        d <- x@pop$number * dev
        x@pop$number <- x@pop$number - d
    }
    if (!pause) {
        x@pop$age <- x@pop$age + 1
        x@pop$development <- x@pop$development + 1
    }
    x@devtable <- x@pop
    x@devtable$number <- d 
    x@devtable$devcycle <- x@devtable$devcycle + 1
    x@devtable$development <- 0
    x@devtable <- x@devtable[x@devtable$number>0,]
    
    x@pop <- x@pop[x@pop$number>0,]
    
    x@developed <- sum(d)
    x@dead <- sum(k)
    return(x)
}

#' @name iterate<-
#' @rdname iterate-spop
#' @exportMethod iterate<-
#' @title Iterate population
#' @description
#' Iterate the population for one day
#' @param x spop class instant
#' @param value \code{data.frame} with the following fields
#' \itemize{
#' \item \code{death}: age-independent daily probability of death
#' \item \code{death_mean} and \code{death_sd}: age-dependent daily probability of death (\code{death_sd=0} indicates fixed life time (defined by \code{death_mean}))
#' \item \code{dev}: age-independent daily probability of development
#' \item \code{dev_mean} and \code{dev_sd}: age-dependent daily probability of development (\code{dev_sd=0} indicates fixed development time (defined by \code{dev_mean}))
#' }
setGeneric(name="iterate<-",
           def=function(x,value) standardGeneric("iterate<-"))
#' @rdname iterate-spop
#' @aliases iterate<-,spop-method
setMethod("iterate<-",
          c("spop","data.frame"),
          function(x, value) {
              return(iterfun(x,value,FALSE))
          })

#' @name perturbate<-
#' @rdname perturbate-spop
#' @exportMethod perturbate<-
#' @title Iterate population without incrementing age or development
#' @description
#' Iterate the population for one day keeping age and development fixed
#' @param x spop class instant
#' @param value \code{data.frame} with the following fields
#' \itemize{
#' \item \code{death}: age-independent daily probability of death
#' \item \code{death_mean} and \code{death_sd}: age-dependent daily probability of death (\code{death_sd=0} indicates fixed life time (defined by \code{death_mean}))
#' \item \code{dev}: age-independent daily probability of development
#' \item \code{dev_mean} and \code{dev_sd}: age-dependent daily probability of development (\code{dev_sd=0} indicates fixed development time (defined by \code{dev_mean}))
#' }
setGeneric(name="perturbate<-",
           def=function(x,value) standardGeneric("perturbate<-"))
#' @rdname perturbate-spop
#' @aliases perturbate<-,spop-method
setMethod("perturbate<-",
          c("spop","data.frame"),
          function(x, value) {
              return(iterfun(x,value,TRUE))
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

#' @name devtable
#' @rdname devtable-spop
#' @exportMethod devtable
#' @title Read devtable
#' @description
#' Read the number, age, and development cycles of individuals completing their development after each iteration
#' @param x spop class instant
setGeneric(name="devtable",
           def=function(x) standardGeneric("devtable"))
#' @rdname devtable-spop
#' @aliases devtable,spop-method
setMethod("devtable",
          "spop",
          function(x) {
              return(x@devtable)
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
