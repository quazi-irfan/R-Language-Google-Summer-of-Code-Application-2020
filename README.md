# MCMC for Forensic Science - Solutions

### Easy: download data from https://lftdi.camden.rutgers.edu/provedit/files/ then use seqinr::read.abif to read and seqinr::plotabif to plot an fsa file.

```R
library("seqinr")
fsa_data <- read.abif("A04_RD14-0003-23d3d-0.0078IP-Q2.5_001.5sec.fsa")
plotabif(fsa_data)
```

### Medium: make a similar multi-panel ggplot with facet_grid.

```R
rm(list=ls())
library("seqinr")
fsa_data <- read.abif("A04_RD14-0003-23d3d-0.0078IP-Q2.5_001.5sec.fsa")
tscale = 1000
tmin = 1/tscale
tmax = fsa_data$Data[["SCAN.1"]]/tscale
irange = (tmin*tscale):(tmax*tscale)
x = irange/tscale

chanel.names = c(1:4,105)

temp <- data.frame()
for(i in 1:5)
{
  chanel = i
  DATA = paste("DATA", chanel.names[chanel], sep = ".")
  yscale = 1000
  y = fsa_data$Data[[DATA]][irange]/yscale
  
  data_label = sapply(rep(DATA, length(fsa_data$Data[[DATA]])), c) # 9335
  
  temp <- rbind(temp, data.frame(x, y, data_label))
}

library(ggplot2)
ggplot(temp, aes(x, y)) + geom_line() + facet_grid(data_label~.)
```

![Medium](Medium/output.png)

### Hard: Demonstrate your capability in one of the "Bayesian packages for general model fitting" listed here: https://cloud.r-project.org/web/views/Bayesian.html, or in writing an R package with C++ code.

Here we are performing multivariate regression using the WaffleDivorce dataset from rethinking package.

```R
rm(list=ls())
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d <- d[ , c("Divorce","MedianAgeMarriage","Marriage") ] 

summary(lm(Divorce ~ MedianAgeMarriage + Marriage, data=d))
```
Here we can see the estimation of the regression coefficients, and we can see the Median Age of Marriage is an important predictor in predicting Divorce rate.

```console
Coefficients:
                  Estimate Std. Error t value Pr(>|t|)    
(Intercept)       36.87665    7.66104   4.814 1.58e-05 ***
MedianAgeMarriage -0.99965    0.24593  -4.065 0.000182 ***
Marriage          -0.05686    0.08053  -0.706 0.483594    
```

Now, if we have prior knowledge about the parameters we are estimating, we have to use Bayesian statistics. We perform another multivariate regression using priors over all parameters. We will use the **rStan** package. We will write our model, and rStan will convert it to a model that Stan model.

```R
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
```
Here we can see that parameter estimates are very similar to what we found with our lm function because our priors has large standard devaition.
```console
       mean   sd  5.5% 94.5% n_eff Rhat4
a     37.02 7.53 25.12 48.40   226     1
bR    -0.06 0.08 -0.19  0.07   262     1
bA    -1.00 0.24 -1.38 -0.62   236     1
sigma  1.52 0.16  1.29  1.77   374     1
```

Internally, Stan is using MCMC method to sample parameter from each of the posterior. We need to investigate the trace plot of each parameter to check if the Markov Chain was staionary.


