\name{spop}
\alias{spop}
\alias{add}
\alias{iterate}
\alias{developed}
\alias{dead}
\title{ Age-Structured Population Dynamics Model }
\description{
    \code{spop} implements the deterministic and stochastic age-structured population dynamics models described in Erguler et al. 2016 and 2017
    
    \code{add} introduces a batch of individuals with a given age (default: 0)
    
    \code{iterate} iterates the population for one day and calculates (overwrites) the number of dead individuals and the number of individuals designated to complete their development
    
    \code{developed} reads the number of individuals designated to complete their development
    
    \code{dead} reads the number of dead individuals after each iteration
}
\author{
\itemize{
    \item Kamil Erguler, Associate Research Scientist, EEWRC, The Cyprus Institute \email{k.erguler@cyi.ac.cy}
    }
}

\details{
    This is an R implementation of the age-structured population dynamics models described in Erguler et al. 2016 and 2017. The spop class records the number and age of individuals and implements two processes to exit from the population. One of the processes is development and the other one is death. The two processes act upon the population sequentially, and survival is imposed prior to development. If the population survives for one day, then, it is allowed to grow and complete its development. Survival and development are defined either with an age-independent daily probability, or an age-dependent gamma-distributed probability.
    \itemize{
    \item death defines an age-independent daily probability of death
    \item death_mean and death_sd define an age-dependent daily probability of death. death_sd=0 indicates fixed life time of death_mean
    \item dev defines an age-independent daily probability of development
    \item dev_mean and dev_sd define an age-dependent daily probability of development. dev_sd=0 indicates fixed development time of dev_mean
    }
}
\references{
    Kamil Erguler, Stephanie E Smith-Unna, Joanna Waldock, Yiannis Proestos, George K Christophides, Jos Lelieveld, Paul E Parham. Large-scale modelling of the environmentally-driven population dynamics of temperate Aedes albopictus (Skuse). PLOS ONE, 2016
    
    Kamil Erguler, Nastassya L Chandra, Yiannis Proestos, Jos Lelieveld, George K Christophides, Paul E Parham. A large-scale stochastic spatiotemporal model for Aedes albopictus-borne chikungunya epidemiology. PLOS ONE, 2017
}
\examples{
    # Generate an integer-valued population with stochastic dynamics
    s <- spop(stochastic=TRUE)
    # Add 1000 20-day-old individuals to the population
    add(s) <- data.frame(number=1000,age=20)

    # Iterate the population for one day with no death process
    # and assume individuals develop in 20 (+-5) days
    iterate(s) <- data.frame(dev_mean=20,dev_sd=5,death=0)
    # Observe the ones ready to pass to the next stage
    print(developed(s))

    # Iterate the population for another day, this time, assuming no development,
    # but the life span of each individual to be 20 days (+-5)
    # Note that the previous developed and dead values will be overwritten by this command
    iterate(s) <- data.frame(death_mean=20,death_sd=5,dev=0)
    print(dead(s))

    # This time, we iterate a deterministic population and observe the difference
    s <- spop(stochastic=FALSE)
    add(s) <- data.frame(number=1000,age=20)
    
    iterate(s) <- data.frame(dev_mean=20,dev_sd=5,death=0)
    print(developed(s))
    
    iterate(s) <- data.frame(death_mean=20,death_sd=5,dev=0)
    print(dead(s))
}
