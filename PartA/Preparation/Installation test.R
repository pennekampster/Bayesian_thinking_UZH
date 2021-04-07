
# Installation test for JAGS, NIMBLE, Stan and TMB

# Bayesian Intro course, 
# UZH, Switzerland,
# 12-14 April 2021

# This is an abbreviated version of 'Code Module Practical model fitting.R'

# When you execute all code in this file and 
# these model fitting engines (and your R) work well, 
# then in the end you should get a table such as this one
# (actual numbers may vary very slightly):

> print(tmp, 4)
     truth     lm   JAGS NIMBLE   Stan    DIY    TMB
mean   600 599.75 599.75 599.71 599.74 599.75 599.75
sd      30  29.65  29.69  29.69  29.66  29.65  29.64


#### Start of the test code right here:

### Some utility functions we need
getMLE <- function(sol, dig = 3){
  MLE <- sol$par
  VC <- solve(sol$hessian)
  ASE <- sqrt(diag(VC))
  print(cbind(MLE, ASE), dig)
}

nimble_summary <- function(samples, params=NULL, digits=3){
  if(!is.null(params)){
    samples <- jagsUI:::order.params(samples, params, FALSE, FALSE)
  }
  mat <- as.matrix(samples)
  nchain <- length(samples)
  niter <- nrow(samples[[1]])
  rhat <- sapply(1:ncol(samples[[1]]), function(i){
    coda::gelman.diag(samples[,i], autoburnin=FALSE)$psrf[1,1]
  })
  stats <- t(apply(mat, 2, function(x){
    x <- na.omit(x)
    c(mean=mean(x), sd=sd(x), quantile(x, c(0.025,0.5,0.975)))
  }))
  out <- data.frame(stats, rhat=rhat, check.names=FALSE)
  cat("Estimates based on",nchain,"chains of",niter,"iterations\n")
  round(out, digits=digits)
}

tmb_summary <- function(tmb_obj){
  npar <- length(tmb_obj$par)
  pnames_fixed <- names(tmb_obj$par)
  out <- summary(sdreport(tmb_obj))
  pnames <- rownames(out)
  pnames[1:npar] <- pnames_fixed
  pcount <- sapply(unique(pnames), function(x) sum(pnames==x))
  idx <- unlist(sapply(pcount, function(i){
    if(i == 1) return("")
    paste0("[",1:i,"]")
    }))
  pnames <- paste0(pnames, idx)
  rownames(out) <- pnames
  out
}


### 2 Data generation
set.seed(39)
y10 <- rnorm(n = 10, mean = 600, sd = 30)       # Sample of 10 birds
y1000 <- rnorm(n = 1000, mean = 600, sd = 30)   # Sample of 1000 birds
truth <- c(mean=600, sd=30)


### 3 MLEs using canned functions in R
out3 <- lm(y1000 ~ 1)
lm_est <- c(coef(out3), sigma=sigma(out3))
print(lm_est)

### 4 Bayesian analysis with JAGS
# -------------------------------
library(jagsUI)
dataList <- list(mass = y1000, n = length(y1000))
cat(file="model4.txt", "        # This code line is R
model {                         # Starting here, we have BUGS code
# Priors
pop.mean ~ dunif(0, 5000)       # Normal parameterized by precision
precision <- 1 / pop.var        # Precision = 1/variance
pop.var <- pop.sd * pop.sd
pop.sd ~ dunif(0, 100)

# Likelihood
for(i in 1:n){
  mass[i] ~ dnorm(pop.mean, precision)
}
}                               # This is the last line of BUGS code
")                              # ... and this is R again
inits <- function()
  list (pop.mean = rnorm(1, 600), pop.sd = runif(1, 1, 30))
params <- c("pop.mean", "pop.sd", "pop.var")
na <- 1000        # Number of iterations in the adaptive phase
ni <- 12000       # Number of draws from the posterior (in each chain)
nb <- 2000        # Number of draws to discard as burn-in
nc <- 3           # Number of chains
nt <- 1           # Thinning rate (nt = 1 means we do not thin)
out4 <- jags(data = dataList, inits = inits, parameters.to.save = params, model.file = "model5.4.txt", n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jags_est <- out4$summary[1:2,1]
print(jags_est)

### 5 Bayesian analysis with NIMBLE
library(nimble)
dataList <- list(mass = y1000, n = length(y1000))
model5 <- nimbleCode( {
# Priors and linear models
pop.mean ~ dunif(0, 5000)       # Normal parameterized by precision
precision <- 1 / pop.variance   # Precision = 1/variance
pop.variance <- pop.sd * pop.sd
pop.sd ~ dunif(0, 100)

# Likelihood
for(i in 1:n){
  mass[i] ~ dnorm(pop.mean, precision)
}
} )
inits <- function()
  list (pop.mean = rnorm(1, 600), pop.sd = runif(1, 1, 30))
params <- c("pop.mean", "pop.sd", "pop.variance")
ni <- 3000   ;   nb <- 1000   ;  nc <- 3   ;  nt <- 1
out5 <- nimbleMCMC(code = model5, constants = dataList, inits = inits, 
monitors = params, nburnin = nb, niter = ni, thin = nt, nchains = nc,
samplesAsCodaMCMC = TRUE)
nsum <- nimble_summary(out5, params)
nimble_est <- nsum[1:2,1]
print(nimble_est)

### 6 Bayesian analysis with Stan
library(rstan)
dataList <- list(n = length(y1000), y = y1000)
cat(file = "model6.stan",  # This line is still R code
"data {                    // This is the first line of Stan code
  int n;                   // Define the format of all data
  vector[n] y;             // ... including the dimension of vectors
}
parameters {               // Same for the parameters
  real mu;
  real<lower=0> sigma;
}
model {
  // Priors
  mu ~ normal(0, 1000);
  sigma ~ cauchy(0, 10);
  // Likelihood
  y ~ normal(mu, sigma);
}                          // This is the last line of Stan code
" )
ni <- 1200   ;   nb <- 200   ;  nc <- 3   ;  nt <- 1
system.time(
  out6 <- stan(file = "model6.stan", data=dataList, 
    chains=nc, iter=ni, warmup=nb, thin=nt)  )     # Ignore the warnings
stan_est <- summary(out6)$summary[1:2,1]
print(stan_est)

### 7 Do-it-yourself maximum likelihood estimates (DIY-MLEs)
NLL <- function(param, y) {
  mu <- param[1]           # Define the first element to be the mean ..
  sigma <- param[2]        # ... and the second the standard deviation
  L <- dnorm(y, mu, sigma) # Likelihood contribution for 1 observation
  LL <- log(L)             # Loglikelihood for 1 observation
  NLL <- -sum(LL)          # NLL for all observations in vector y
  return(NLL)
}
# Minimize NLL and use quadratic approximation of SEs using the Hessian
inits <- c('mu' = 500, 'sigma' = 10)
out7 <- optim(inits, NLL, y=y1000, hessian=TRUE)  # Large data set
getMLE(out7, 5)
diy_est <- out7$par

### 8 Get the MLEs with TMB
library(TMB)
dataList <- list(y = y1000, n = length(y1000))
cat(file="model8.cpp",
"#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_INTEGER(n);      //Sample size
  DATA_VECTOR(y);       //response
  //Describe parameters
  PARAMETER(mu);        //Mean
  PARAMETER(log_sigma); //log(standard deviation)
  Type sigma = exp(log_sigma);  //Type = match type of function output (double)
  Type loglik = 0.0;    //Initialize total log likelihood at 0
  for (int i=0; i<n; i++){ //Note index starts at 0 instead of 1!
    //Calculate log-likelihood of observation
    //value of true in final argument indicates function should return log(lik)
    loglik += dnorm(y(i), mu, sigma, true); //Add log-lik of obs i
  }
  return -loglik; //Return negative of total log likelihood
}
")
compile("model8.cpp")
dyn.load(dynlib("model8"))
params <- list(mu=0, log_sigma=0)
out8 <- MakeADFun(data = dataList,
            parameters = params, random=NULL,
            DLL = "model8", silent=TRUE)
starts <- rep(0, length(unlist(params)))
sol <- optim(starts, fn=out8$fn, gr=out8$gr, method="BFGS", hessian = TRUE)
tsum <- tmb_summary(out8)
tmb_est <- c(sol$par[1], exp(sol$par[2])) # backtransform SD estimate!
tmp <- cbind(cbind(truth=truth, lm=lm_est, JAGS=jags_est, 
   NIMBLE=nimble_est, Stan=stan_est, DIY=diy_est, TMB=tmb_est))
print(tmp, 4)