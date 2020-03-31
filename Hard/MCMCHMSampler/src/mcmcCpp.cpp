#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mcmcCpp(int length, int success, int trials) {
  NumericVector chain(length);
  
  chain[0] = 0.5;

  for(int x = 1; x<length; x++){
    double nextState = rnorm(1, chain[x-1], 0.16)[0];

    double numvalue = nextState;
    numvalue = numvalue<0.001?0.001:numvalue;
    numvalue = numvalue>0.999?0.999:numvalue;

    // double denumvalue = chain[x-1];
    // denumvalue = denumvalue<0.001?0.001:denumvalue;
    // denumvalue = denumvalue>0.999?0.999:denumvalue;

    double numerator = R::dbeta(nextState, 1, 1, false) * R::dbinom(success, trials, numvalue, false);
    double denumerator = R::dbeta(chain[x-1], 1, 1, false) * R::dbinom(success, trials, chain[x-1], false);
    double likelihood_ratio = numerator / denumerator;
    
    double acceptance_prob = std::min(1.0, likelihood_ratio);
    if( acceptance_prob == 1.0){
      chain[x] = nextState;
    }else if(acceptance_prob > Rcpp::runif(1)[0]){
      chain[x] = nextState;
    }else{
      chain[x] = chain[x-1];
    }
  }
  return chain;
}

// rm(list=ls())
// library(raster)
// n <- 100000
// chain <- sapply(rep(NA, n), c)
// chain[1] <- 0.5
// for(x in 2:n){
//   nextState <- rnorm(1, chain[x-1], sd = 0.16)
//   
//   likelihood_ratio <- 
//     (dbeta(nextState,1,1) * dbinom(x=4,size=10,prob=clamp(nextState,0,1)))/
//       (dbeta(chain[x-1],1,1) * dbinom(x=4,size=10,prob=clamp(chain[x-1],0,1)))
//     
//     acceptance_prob <- min(likelihood_ratio, 1)
//     if(acceptance_prob == 1){
//       chain[x] = nextState
//     } else if(acceptance_prob > runif(1)){
//       chain[x] = nextState
//     } else {
//       chain[x] = chain[x-1]
//     }
// }