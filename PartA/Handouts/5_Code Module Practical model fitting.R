
# Code from:
# A tutorial on the use of the model fitting engines 
# JAGS, NIMBLE, Stan and TMB and on 
# Do-it-yourself maximum likelihood (DIY-MLE)

# Bayesian Intro course, 
# UZH, Switzerland,
# April 2021

# Written by:
# Marc Kéry (Swiss Ornithological Institute, Sempach, Switzerland) &
# Ken Kellner (College of Environmental Science and Forestry, SUNY Syracuse, NY, USA)
# 2020-2021

# This file contains material from the book "Introduction to WinBUGS for Ecologists",
# Academic Press, 2010; see https://www.mbr-pwrc.usgs.gov/software/kerybook/ 

# For an overview of our research programme on hierarchical modeling in ecology,
# see https://www.mbr-pwrc.usgs.gov/pubanalysis/roylebook/ 



### 2 Data generation
# -------------------
set.seed(39)

# Generate two samples of body mass measurements of male peregrines
y10 <- rnorm(n = 10, mean = 600, sd = 30)       # Sample of 10 birds
y1000 <- rnorm(n = 1000, mean = 600, sd = 30)   # Sample of 1000 birds

# Save the data-generating values of the parameters for later comparisons
truth <- c(mean=600, sd=30)

# Plot data (Fig. 2)
xlim = c(450, 750)
par(mfrow = c(1, 2), mar = c(6, 6, 6, 3), cex.lab = 1.5, cex.axis = 1.5)
hist(y10, col = 'grey ', xlim = xlim, main = 'Body mass (g) of 10 male peregrines')
hist(y1000, col = 'grey', xlim = xlim, main = ' Body mass (g) of 1000 male peregrines')




### 3 MLEs using canned functions in R
# ------------------------------------
summary(out3 <- lm(y10 ~ 1))         # small data set: not shown
summary(out3 <- lm(y1000 ~ 1))

# Compare estimates with truth
lm_est <- c(coef(out3), sigma=sigma(out3))
cbind(truth=truth, lm=lm_est)




### 4 Bayesian analysis with JAGS
# -------------------------------
library(jagsUI)

# Bundle and summarize data (here the larger variant)
str(dataList <- list(mass = y1000, n = length(y1000)))

# Write JAGS model file
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

# Function to generate starting values
inits <- function()
  list (pop.mean = rnorm(1, 600), pop.sd = runif(1, 1, 30))

# Parameters monitored (= to be estimated)
params <- c("pop.mean", "pop.sd", "pop.var")

# MCMC settings
na <- 1000        # Number of iterations in the adaptive phase
ni <- 12000       # Number of draws from the posterior (in each chain)
nb <- 2000        # Number of draws to discard as burn-in
nc <- 3           # Number of chains
nt <- 1           # Thinning rate (nt = 1 means we do not thin)

# Call JAGS (ART 1 min) and marvel at JAGS' progress bar
out4 <- jags(data = dataList, inits = inits, parameters.to.save = params, model.file = "model4.txt", n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = FALSE)

# Call JAGS in parallel (ART <1 min) and check convergence
out4 <- jags(data = dataList, inits = inits, parameters.to.save = params, model.file = "model4.txt", n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
par(mfrow=c(2, 2)); jagsUI::traceplot(out4)        # Produce Fig. 3

# Produce a summary of the fitted model object 
print(out4, 3)
print(out4$summary, 3)    # Same, but only main body of table

# What's in the object created ?
names(out4)
str(out4)          # Full resolution of overview (not shown)
str(out4, 1)       # Reduced resolution of overview of out4


# Here is the basic result from running the MCMC algorithm
str(out4$sims.list)

# And here for one parameter
print(out4$sims.list$pop.mean)

# Check for convergence
hist(out4$summary[,8]) # Rhat values in the eighth column of the summary
which(out4$summary[,8] > 1.1) # None in this case

# For trace plots "by hand" for the entire chains do:
par(mfrow = c(3,1), mar = c(5,5,4,2))
matplot(cbind(out4$samples[[1]][,1], out4$samples[[2]][,1], out4$sample[[3]][,1]), type = 'l', lty = 1, col = c('red', 'blue', 'green'), xlab = 'Iteration', ylab = 'Value', main = 'Population mean', frame = FALSE)
matplot(cbind(out4$samples[[1]][,2], out4$samples[[2]][,2], out4$sample[[3]][,2]), type = 'l', lty = 1, col = c('red', 'blue', 'green'), xlab = 'Iteration', ylab = 'Value', main = 'Population sd', frame = FALSE)
matplot(cbind(out4$samples[[1]][,3], out4$samples[[2]][,3], out4$sample[[3]][,3]), type = 'l', lty = 1, col = c('red', 'blue', 'green'), xlab = 'Iteration', ylab = 'Value', main = 'Population variance', frame = FALSE)

#  or just for the start of the chains, to see how rapidly they converge ...
n <- 10         # Choose last iteration in plot
par(mfrow = c(3,1), mar = c(5,5,4,2))
matplot(cbind(out4$samples[[1]][1:n,1], out4$samples[[2]][ 1:n,1], out4$sample[[3]][1:n,1]), type = 'l', lty = 1, col = c('red', 'blue', 'green'), xlab = 'Iteration', ylab = 'Value', main = 'Population mean', frame = FALSE)
matplot(cbind(out4$samples[[1]][1:n,2], out4$samples[[2]][1:n,2], out4$sample[[3]][1:n,2]), type = 'l', lty = 1, col = c('red', 'blue', 'green'), xlab = 'Iteration', ylab = 'Value', main = 'Population sd', frame = FALSE)
matplot(cbind(out4$samples[[1]][1:n,3], out4$samples[[2]][1:n,3], out4$sample[[3]][1:n,3]), type = 'l', lty = 1, col = c('red', 'blue', 'green'), xlab = 'Iteration', ylab = 'Value', main = 'Population variance', frame = FALSE)

# Produce graphical summaries

par(mfrow = c(3, 2), mar = c(5,5,4,2))    # Fig. 4
hist(out4$sims.list$pop.mean, col = "grey", breaks = 100, xlab = "pop.mean")
plot(density(out4$sims.list$pop.mean), type = 'l', lwd = 3, col = rgb(0,0,0,0.6), main = "", frame = FALSE)
hist(out4$sims.list$pop.sd, col = "grey", breaks = 100, xlab = "pop.sd")
plot(density(out4$sims.list$pop.sd), type = 'l', lwd = 3, col = rgb(0,0,0,0.6), main = "", frame = FALSE)
hist(out4$sims.list$pop.var, col = "grey", breaks = 100, xlab = "pop.var")
plot(density(out4$sims.list$pop.var), type = 'l', lwd = 3, col = rgb(0,0,0,0.6), main = "", frame = FALSE)

# Can play probability games based on produced draws from posteriors
sims <- out4$sims.list
sims$pop.mean < 599       # MANY logical tests ! (not shown)
mean(sims$pop.mean < 599) # Prob(mu < 599)

# Same in two dimensions
# Fig. 5: compute 2d probability game and plot bivariate posterior
test.true <- sims$pop.mean < 600 & sims$pop.sd > 30
mean(test.true)         # [1] 0.1851333

# Plot that same probability
plot(sims$pop.mean, sims$pop.sd, pch = 16, col = rgb(0,0,0,0.3), cex = 0.8, frame = FALSE)
points(sims$pop.mean[test.true], sims$pop.sd[test.true], pch = 16, col = rgb(1,0,0,0.6), cex = 0.8)
abline(h = 30, col = 'red', lwd = 1)
abline(v = 600, col = 'red', lwd = 1)

# alternatively do this to see bivariate posterior plots
pairs(cbind(out4$sims.list$pop.mean, out4$sims.list$pop.sd, out4$sims.list$pop.var))

# Numerical summaries "by hand"
summary(out4$sims.list$pop.mean)
summary(out4$sims.list$pop.sd)
sd(out4$sims.list$pop.mean)
sd(out4$sims.list$pop.sd)

# Compare likelihood with Bayesian estimates and with truth
jags_est <- out4$summary[1:2,1]
tmp <- cbind(truth=truth, lm=lm_est, JAGS=jags_est)
print(tmp, 3)




### 5 Bayesian analysis with NIMBLE
# ---------------------------------

# Load NIMBLE
library(nimble)

# Bundle and summarize data (same as for JAGS)
str(dataList <- list(mass = y1000, n = length(y1000)))

# Write Nimble model file
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

# Can use same function to generate starting values as for JAGS
inits <- function()
  list (pop.mean = rnorm(1, 600), pop.sd = runif(1, 1, 30))

# Parameters monitored: same as before
params <- c("pop.mean", "pop.sd", "pop.variance")

# MCMC settings
ni <- 3000   ;   nb <- 1000   ;  nc <- 3   ;  nt <- 1

## Run Nimble in a single step using function nimbleMCMC
# Call NIMBLE (ART 20 sec), check convergence and summarize posteriors
system.time(
  out5 <- nimbleMCMC(code = model5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )

# Overview of output
str(out5, 1)

# Use own function nimble_summary for producing a summary table
(nsum <- nimble_summary(out5, params))   # Produce posterior summary table

# We can produce trace plots using the function in R package coda.
par(mfrow=c(1, 3)); coda::traceplot(out5)   # not shown


## Excursion: Running NIMBLE with all different steps 
#  of the analysis dissected will (later) allow you 
#  MUCH more user control

# Create a NIMBLE model from BUGS code
rawModel <- nimbleModel(code = model5, constants=dataList, inits=inits())

# Configure the MCMC algorithm: create a default MCMC configuration
# and monitor our selected list of parameters (See Section 7 in Manual)
mcmcConfig <- configureMCMC(rawModel, monitors=params)

# Build the MCMC algorithm function
rawMCMC <- buildMCMC(mcmcConfig)

# Can run the algorithm now, but slow 
system.time(rawMCMC$run(30)) # Take 30 samples
as.matrix(rawMCMC$mvSamples)

# Convert (i.e., compile) our slow R code into fast C++ code
compModel <- compileNimble(rawModel)

# Also compile our MCMC function
compMCMC <- compileNimble(rawMCMC, project = rawModel)

# Run the compiled MCMC function, which will be significantly faster
system.time(compMCMC$run(30))

# Use compMCMC and the runMCMC function to generate posterior samples from 3 chains
system.time(
samples <- runMCMC(compMCMC, niter=ni, nburnin=nb, thin=nt, nchains=nc,
                   samplesAsCodaMCMC=TRUE) )  # Takes almost no time !
# Peak at samples
lapply(samples, head)

# Produce marginal posterior summaries using our utility function
(nsum <- nimble_summary(samples))

# Traceplots
par(mfrow=c(1,3)); coda::traceplot(samples)

# Make our usual comparison with truth for all estimates up to now.

# Compare estimates with truth
nimble_est <- nsum[1:2,1]
tmp <- cbind(truth=truth, lm=lm_est, JAGS=jags_est, NIMBLE=nimble_est)
print(tmp, 4)




### 6 Bayesian analysis with Stan
# -------------------------------

# Load Stan R package
library(rstan)

# Bundle and summarize data
str(dataList <- list(n = length(y1000), y = y1000)) # same as before, not shown

# Write text file with model description in Stan language
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

# HMC settings
ni <- 1200   ;   nb <- 200   ;  nc <- 3   ;  nt <- 1

# Call STAN (ART 28 / 2 sec)
system.time(
  out6 <- stan(file = "model6.stan", data=dataList, 
    chains=nc, iter=ni, warmup=nb, thin=nt)  )     # Ignore the warnings

# Much faster on second execution than first

# Check convergence and print posterior summaries
par(mfrow = c(2, 2))                # not sure why this is ignored !
rstan::traceplot(out6, c("mu", "sigma"))   # Wilted-flower color plots ...
print(out6, c("mu", "sigma"))
print(out6, dig = 1)

# Overview of what Stan produced
str(out6)           # not shown

# Visual output using CODA
library(coda)
sims <- extract(out6)
str(sims)
plot(mcmc(cbind(sims[[1]], sims[[2]], sims[[3]])))

# Compare estimates with truth
stan_est <- summary(out6)$summary[1:2,1]
tmp <- cbind(truth=truth, lm=lm_est, JAGS=jags_est, NIMBLE=nimble_est, Stan=stan_est)
print(tmp, 4)




### 7 Do-it-yourself maximum likelihood estimates (DIY-MLEs)
# ----------------------------------------------------------

# Definition of NLL for a normal linear model with an intercept only
NLL <- function(param, y) {
  mu <- param[1]           # Define the first element to be the mean ..
  sigma <- param[2]        # ... and the second the standard deviation
  L <- dnorm(y, mu, sigma) # Likelihood contribution for 1 observation
  LL <- log(L)             # Loglikelihood for 1 observation
  NLL <- -sum(LL)          # NLL for all observations in vector y
  return(NLL)
}

# An alternative NLL, which is more numerically accurate
NLL <- function(param, y) {
  mu <- param[1]
  sigma <- param[2]
  # log-likelihood contribution for 1 observation
  LL <- dnorm(y, mu, sigma, log = TRUE) 
  NLL <- -sum(LL)          # NLL for all observations in vector y
  return(NLL)
}

# Minimize NLL and use quadratic approximation of SEs using the Hessian
inits <- c('mu' = 500, 'sigma' = 10)
sol <- optim(inits, NLL, y=y10, hessian=TRUE)  # Small data set
MLE <- sol$par                 # Grab MLEs
VC <- solve(sol$hessian)       # Get variance-covariance matrix
ASE <- sqrt(diag(VC))          # Extract asymptotic SEs
print(cbind(MLE, ASE), 3)      # Print MLEs and SE's


# Look at output from optimizer
sol

# Compare with least-squares solution
summary(lm(y10~1))

# Define function to get asymptotic SEs and print out MLEs and SEs
getMLE <- function(sol, dig = 3){
  MLE <- sol$par
  VC <- solve(sol$hessian)
  ASE <- sqrt(diag(VC))
  print(cbind(MLE, ASE), dig)
}

out7 <- optim(inits, NLL, y=y1000, hessian=TRUE)  # Large data set
getMLE(out7, 5)


# Compare again with least-squares solution for the big data set
summary(lm(y1000~1))

# Compare estimates with truth and previous estimates
diy_est <- out7$par
tmp <- cbind(truth=truth, lm=lm_est, JAGS=jags_est,
   NIMBLE=nimble_est, Stan=stan_est, DIY=diy_est)
print(tmp, 4)




### 8 Get the MLEs with TMB instead
# ---------------------------------

library(TMB)

# Bundle and summarize data
str(dataList <- list(y = y1000, n = length(y1000)))

# Write TMB model file
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

# Compile and load TMB function
compile("model8.cpp") # Produces a lot of gibberish ... don't panick !
dyn.load(dynlib("model8"))

?TMB::compile               # Check out that function if you wish

# Provide dimensions and starting values for parameters
params <- list(mu=0, log_sigma=0)

# Create TMB object ready to be optimized
out8 <- MakeADFun(data = dataList,
            parameters = params, random=NULL,
            DLL = "model8", silent=TRUE)

?TMB::MakeADFun               # Check out that function

# Optimize TMB object and print results
starts <- rep(0, length(unlist(params)))
sol <- optim(starts, fn=out8$fn, gr=out8$gr, method="BFGS", hessian = TRUE)

# Use our utility function tmb_summary to get parameter estimates and standard errors from an optimized TMB object
(tsum <- tmb_summary(out8))  # Use our summary function




#### 9 Comparison of the all the parameter estimates for the model of the mean
# Compare results with truth and previous estimates
tmb_est <- c(sol$par[1], exp(sol$par[2])) # backtransform SD estimate!
tmp <- cbind(cbind(truth=truth, lm=lm_est, JAGS=jags_est, 
   NIMBLE=nimble_est, Stan=stan_est, DIY=diy_est, TMB=tmb_est))
print(tmp, 4)

### We would argue that you can consider 
#   these estimates as numerically identical 
#   for virtually all practical purposes

