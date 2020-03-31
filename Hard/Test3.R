# test 3  

# Linear Regression(lm)
rm(list=ls())
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d <- d[ , c("Divorce","MedianAgeMarriage","Marriage") ] 

summary(lm(Divorce ~ MedianAgeMarriage + Marriage, data=d))

# Bayesian approach
model <- map2stan( 
  alist(
    Divorce ~ dnorm( mu , sigma ) ,
    mu <- a + bA*MedianAgeMarriage + bR*Marriage ,
    a ~ dnorm(0,100),
    bR ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    sigma ~ dcauchy(0,2)
  ) ,
  data=d )
precis(model)

# Trace plot of our Markov Chain
plot(model)

# MCMC HM Implementation in R
rm(list=ls())
library(raster)

mcmcR <- function(samples, success, trials){
  n <- samples
  chain <- sapply(rep(NA, n), c)
  chain[1] <- 0.5
  for(x in 2:n){
    nextState <- rnorm(1, chain[x-1], sd = 0.16) # Monte Carlo
    
    # estimate coin toss p
    # dbeta(1,1) Almost uniform prior
    # dbinom(4, 10) 4 heads our of 10 tosses
    # Metropolis-Hasting Step
    likelihood_ratio <- 
      (dbeta(nextState,1,1) * dbinom(x=success,size=trials,prob=clamp(nextState,0,1)))/
      (dbeta(chain[x-1],1,1) * dbinom(x=success,size=trials,prob=clamp(chain[x-1],0,1)))
    
    acceptance_prob <- min(likelihood_ratio, 1)
    if(acceptance_prob == 1){
      chain[x] = nextState
    } else if(acceptance_prob > runif(1)){
      chain[x] = nextState
    } else {
      chain[x] = chain[x-1]
    }
  }
  return(chain)
}

# Posterior density from mcmcR
par(mfrow=c(1,1))
chain <- mcmcR(10000, 4, 10)
plot(chain, type="l")
graphics::hist(chain, freq=F, ylim=c(0, 3), ylab="", main="Samples from Posterior")
lines(density(chain), type="l",col = "red")
lines(density(rbeta(1000000, 5, 7)), type="l", col="blue")

# Posterior density from mcmcCpp
chainCpp <- mcmcCpp(10000, 4, 10)
plot(chainCpp, type="l")
plot(density(chainCpp), col="red", ylim=c(0, 3), main="Samples from Posterior(C++)")
lines(density(rbeta(1000000, 5, 7)), type="l", col="blue")

# Benchmark
library(microbenchmark)
library(ggplot2)
benchmark_result <- microbenchmark(mcmcR(10000, 4, 10), mcmcCpp(10000, 4, 10))
ggplot2::autoplot(benchmark_result)
