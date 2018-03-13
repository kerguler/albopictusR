# albopictus
Age-structured population dynamics model

## Installation
The R package can be installed from the command line,

    R CMD install albopictus_x.x.tar.gz

to be loaded easily at the R command prompt.

    library(albopictus)

## Usage

Generate an integer-valued population with stochastic dynamics

    s <- spop(stochastic=TRUE)

Add 1000 20-day-old individuals to the population

    add(s) <- data.frame(number=1000,age=20)

Iterate the population for one day with no death process and assume individuals develop in 20 (+-5) days

    iterate(s) <- data.frame(dev_mean=20,dev_sd=5,death=0)

Observe the ones ready to pass to the next stage

    print(developed(s))

Iterate the population for another day, this time, assuming no development, but the life span of each individual to be 20 days (+-5). Note that the previous developed and dead values will be overwritten by this command

    iterate(s) <- data.frame(death_mean=20,death_sd=5,dev=0)
    print(dead(s))

This time, we iterate a deterministic population and observe the difference

    s <- spop(stochastic=FALSE)
    add(s) <- data.frame(number=1000,age=20)
    
    iterate(s) <- data.frame(dev_mean=20,dev_sd=5,death=0)
    print(developed(s))
    
    iterate(s) <- data.frame(death_mean=20,death_sd=5,dev=0)
    print(dead(s))

For more information, please contact me at [k.erguler@cyi.ac.cy](mailto:k.erguler@cyi.ac.cy).
