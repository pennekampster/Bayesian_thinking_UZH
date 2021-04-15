
# Bayes course UZH April 2021

# Marc Kéry, Swiss Ornithological Institute

# Exercise: in Chapter 11, turn the random-effects model into the analogous 
# fixed-effects model. That is, make both the intercepts and the slopes fixed-effects
# (but still estimate 56 intercepts and 56 slopes !)

# We work with this data set:
str(dataList <- list(mass = as.numeric(mass), pop = as.numeric(pop), length = length, nPops = nPops, n = n) )

# and it is best if you copy into a new file the model in Section 11.4.2 and then modify that

# Solution: 
# ---------
# Random effects are sets of parameters that are given a common distribution
# with shared (hyper-)parameters that are estimated.
# Thus, we "simply" have to replace the "dnorm(XXX.mean, XXX.prec)" statement
# in the definition of the priors for intercept and slope with this: "dnorm(0, 'some small value')".
# Then, there are also a couple other small changes in the rest of the code on which we comment
# on below.

# This may be a little confusing at first, since at first sight, the sets of
# intercepts and slopes are still given a Normal distribution at that place.

# However, the key thing is that when we estimate them as fixed effects, then in these Normal 
# distributions there no longer are any shared hyperparameters that we estimate,
# but instead we choose simply constants for the mean and the precision of these two Normals.

# The first solution that was proposed by a participant was in fact correct, but it had one
# unexpected snag: the chosen prior for the intercepts (estimated as fixed effects) was
# not vague, inadvertently.

# We first with this proposed model and then compare it with another one that has wider
# Normal priors for the (fixed) intercepts. This is an interesting exercise in prior sensitivity
# analysis. It emphasizes the fact that whether a certain prior is vague or informative
# must be ascertained anew in every analysis !


# So here is the model that was originally suggested as the solution. We denote it by adding
# one 'F' in the name (F for 'fixed').

# Write JAGS model file
cat(file="model11.4.2F.txt", "
model {
# Priors
for (i in 1:nPops){

# Here, hashed out, is the original where intercepts and slopes are estimated
# as random effect
#  intercept[i] ~ dnorm(intercept.mean, intercept.prec)  # Random intercepts
#  slope[i] ~ dnorm(slope.mean, slope.prec)              # Random slopes

# And here is the solution for when we want to estimate them as fixed effects instead
# Note that the dispersion parameter in the Normal of JAGS is the precision, i.e., 1 over
# the variance, hence 1/1000 represents a variance of 1000, which is a SD of about 32.
  intercept[i] ~ dnorm(0, 1/ 1000)           # Fixed intercepts
  slope[i] ~ dnorm(0, 1 / 1000)              # Fixed slopes
}

intercept.mean ~ dnorm(0, 0.001)  # Mean hyperparameter for random intercepts
intercept.prec <- pow(intercept.sd, -2)
intercept.sd ~ dunif(0, 100)      # SD hyperparameter for random intercepts

slope.mean ~ dnorm(0, 0.001)      # Mean hyperparameter for random slopes
slope.prec <- pow(slope.sd, -2)
slope.sd ~ dunif(0, 100)          # SD hyperparameter for slopes

prec <- pow(residual.sd, -2)      # Residual precision
residual.sd ~ dunif(0, 100)       # Residual standard deviation

# 'Likelihood'
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], prec)
  mu[i] <- intercept[pop[i]] + slope[pop[i]] * length[i]
}
}
")

# In the inits, we also have to take out the hyperparameters, but instead now give 56 inits
# for each of the intercept and the slope vectors
# Function to generate starting values
inits <- function(){ list(intercept = rnorm(56, 0, 1), 
  slope = rnorm(56, 0, 1), residual.sd = rlnorm(1)) }

# Here we can drop the names of the hyperparameters that are no longer part of the model
# Parameters to estimate
params <- c("intercept", "slope", "residual.sd")

# MCMC settings
na <- 1000  ;  ni <- 3000  ;  nb <- 1000  ; nc <- 3  ; nt <- 1

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out11.4.2F <- jags(dataList, inits, params, "model11.4.2F.txt", 
  n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na,
  parallel = TRUE)
par(mfrow=c(2, 2)); jagsUI::traceplot(out11.4.2F)     # not shown
print(out11.4.2F, 2)

# This works all fine, convergence is fine etc but ......
# Look at the scale of the estimated intercepts, they are in the range of about 100.
# but our prior for them was a Normal with SD of just about 30 .... so this prior
# is likely to be inadvertently highly informative 
# Actually, even for the slopes the chosen prior may be informative.

# Hence, we repeat the analysis with a much wider variance for both the intercepts and 
# the slopes and then compare the estimates. We name this new model 'FF'

# All the rest of the code is identical though.

# Write JAGS model file
cat(file="model11.4.2FF.txt", "
model {
# Priors
for (i in 1:nPops){

# Here now we choose a variance of 1 Million
  intercept[i] ~ dnorm(0, 1/ 1000000)           # Fixed intercepts
  slope[i] ~ dnorm(0, 1 / 1000000)              # Fixed slopes
}

intercept.mean ~ dnorm(0, 0.001)  # Mean hyperparameter for random intercepts
intercept.prec <- pow(intercept.sd, -2)
intercept.sd ~ dunif(0, 100)      # SD hyperparameter for random intercepts

slope.mean ~ dnorm(0, 0.001)      # Mean hyperparameter for random slopes
slope.prec <- pow(slope.sd, -2)
slope.sd ~ dunif(0, 100)          # SD hyperparameter for slopes

prec <- pow(residual.sd, -2)      # Residual precision
residual.sd ~ dunif(0, 100)       # Residual standard deviation

# 'Likelihood'
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], prec)
  mu[i] <- intercept[pop[i]] + slope[pop[i]] * length[i]
}
}
")

inits <- function(){ list(intercept = rnorm(56, 0, 1), 
  slope = rnorm(56, 0, 1), residual.sd = rlnorm(1)) }
params <- c("intercept", "slope", "residual.sd")
na <- 1000  ;  ni <- 3000  ;  nb <- 1000  ; nc <- 3  ; nt <- 1

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out11.4.2FF <- jags(dataList, inits, params, "model11.4.2FF.txt", 
  n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na,
  parallel = TRUE)
par(mfrow=c(2, 2)); jagsUI::traceplot(out11.4.2FF)     # not shown
print(out11.4.2FF, 2)

# This looks good and we compare the estimates from the two versions of the model with the 
# truth in the data simulation.
lim <- c(100, 300)
par(mfrow = c(2, 2) )
plot(intercept.effects, out11.4.2F$mean$intercept, pch = 16, xlab = "True values in data simulation", ylab = "Estimates", main = "Intercepts\n(prior variance = 10^3)", frame = FALSE,
xlim = lim, ylim = lim)
abline(0, 1)

plot(intercept.effects, out11.4.2FF$mean$intercept, pch = 16, xlab = "True values in data simulation", ylab = "Estimates", main = "Intercepts\n(prior variance = 10^6)", frame = FALSE,
xlim = lim, ylim = lim)
abline(0, 1)

lim <- c(-100, 200)
plot(slope.effects, out11.4.2F$mean$slope, pch = 16, xlab = "True values in data simulation", ylab = "Estimates", main = "Slopes\n(prior variance = 10^3)", frame = FALSE,
xlim = lim, ylim = lim)
abline(0, 1)

plot(slope.effects, out11.4.2FF$mean$slope, pch = 16, xlab = "True values in data simulation", ylab = "Estimates", main = "Slopes\n(prior variance = 10^6)", frame = FALSE,
xlim = lim, ylim = lim)
abline(0, 1)

# We see that with a Variance of 10^6 we have chosen a really vague prior for both sets
# of parameters, but not with 10^3


# Then, we had the question of whether we could have also scaled the response
# and we recognized that we could just express the response in kilograms instead of 
# in grams. That should then also bring down the magnitude of the intercept and
# thus render the originally chosen prior with variance 1000 vague. Let us try that.


# Scale mass in kgs
# Bundle and summarize data
str(dataList <- list(massKG = as.numeric(mass)/1000, pop = as.numeric(pop), length = length, nPops = nPops, n = n) )

# Write JAGS model file
cat(file="model11.4.2FFF.txt", "
model {
# Priors
for (i in 1:nPops){
  intercept[i] ~ dnorm(0, 1/ 1000)           # Fixed intercepts
  slope[i] ~ dnorm(0, 1 / 1000)              # Fixed slopes
}

intercept.mean ~ dnorm(0, 0.001)  # Mean hyperparameter for random intercepts
intercept.prec <- pow(intercept.sd, -2)
intercept.sd ~ dunif(0, 100)      # SD hyperparameter for random intercepts

slope.mean ~ dnorm(0, 0.001)      # Mean hyperparameter for random slopes
slope.prec <- pow(slope.sd, -2)
slope.sd ~ dunif(0, 100)          # SD hyperparameter for slopes

prec <- pow(residual.sd, -2)      # Residual precision
residual.sd ~ dunif(0, 100)       # Residual standard deviation

# 'Likelihood'
for (i in 1:n) {
  massKG[i] ~ dnorm(mu[i], prec)
  mu[i] <- intercept[pop[i]] + slope[pop[i]] * length[i]
}
}
")

# Function to generate starting values
inits <- function(){ list(intercept = rnorm(56, 0, 1), 
  slope = rnorm(56, 0, 1), residual.sd = rlnorm(1)) }

# Parameters to estimate
params <- c("intercept", "slope", "residual.sd")

# MCMC settings
na <- 1000  ;  ni <- 3000  ;  nb <- 1000  ; nc <- 3  ; nt <- 1

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out11.4.2FFF <- jags(dataList, inits, params, "model11.4.2FFF.txt", 
  n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na,
  parallel = TRUE)

# We compare the three sets of estimates in a table
# Intercepts
cbind(out11.4.2F$mean$intercept, out11.4.2FF$mean$intercept, out11.4.2FFF$mean$intercept)
# Slopes
cbind(out11.4.2F$mean$slope, out11.4.2FF$mean$slope, out11.4.2FFF$mean$slope)

# Comparing columns 2 and 3 (models FF and FFF) we see that this solution 
# is also fine
