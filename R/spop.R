gamma_dist_prob <- function(array,mean,sd) {
    theta <- sd * sd / mean
    k <- mean / theta
    return(1.0 - (1.0-pgamma(array+1,k,scale=theta))/(1.0-pgamma(array,k,scale=theta)))
}

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
           def=function(object,value,...) standardGeneric("add<-"))
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
           def=function(object,value,...) standardGeneric("iterate<-"))
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
                      dev <- gamma_dist_prob(object@pop$age,value$dev_mean,value$dev_sd)
                  else {
                      warning(sprintf("Error in development probability"))
                      return(object)
                  }
              } else
                  dev <- value$dev
              if (!("death" %in% colnames(value))) {
                  if (("death_mean" %in% colnames(value)) && ("death_sd" %in% colnames(value)))
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
           def=function(object,...) standardGeneric("developed"))
setMethod("developed",
          "spop",
          function(object) {
              return(object@developed)
          })
setGeneric(name="size",
           def=function(object,...) standardGeneric("size"))
setMethod("size",
          "spop",
          function(object) {
              return(sum(object@pop$number))
          })
setGeneric(name="dead",
           def=function(object,...) standardGeneric("dead"))
setMethod("dead",
          "spop",
          function(object) {
              return(object@dead)
          })
