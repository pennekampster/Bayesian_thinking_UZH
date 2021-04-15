

# Exercise 2 on Wednesday
# Fit a Normal GLM to the same data instead of the Poisson GLM

# Note: this is not in general a good idea (to adopt a Normal for count data)
#  and this exercise is SOLELY meant to clarify to you the relationship between
#  these two GLMs and to show how we implement them in JAGS


# 14.4 Bayesian analysis with JAGS
# --------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, pop = as.numeric(pop), nPops = nPops, length = length, n = n) )

# Write JAGS model file
cat(file="model14.4N.txt", "
model {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(0, 0.01)      # Intercepts
  beta[i] ~ dnorm(0, 0.01)       # Slopes
}
tau <- pow(sd, -2)     # = 1 / (sd * sd)
sd ~ dunif(0, 100)

# Likelihood
for (i in 1:n) {
  C[i] ~ dnorm(lambda[i], tau)  # The random variable
  lambda[i] <- alpha[pop[i]] + beta[pop[i]]* length[i]
}

# Derived quantities
# Recover effects relative to baseline level (no. 1)
a.effe2 <- alpha[2] - alpha[1]   # Intercept Massif Central vs. Pyr.
a.effe3 <- alpha[3] - alpha[1]   # Intercept Jura vs. Pyr.
b.effe2 <- beta[2] - beta[1]     # Slope Massif Central vs. Pyr.
b.effe3 <- beta[3] - beta[1]     # Slope Jura vs. Pyr.

# Custom test
test1 <- beta[3] - beta[2]       # Slope Jura vs. Massif Central
}
")

# Function to generate starting values
inits <- function(){list(alpha=rlnorm(nPops, 3, 1), beta=rlnorm(nPops, 2, 1))}

# Parameters to estimate
params <- c("alpha", "beta", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1", "sd")

# MCMC settings
na <- 2000  ;  ni <- 25000  ;  nb <- 5000  ; nc <- 3  ; nt <- 5

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out14.4N <- jags(dataList, inits, params, "model14.4N.txt", n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
par(mfrow=c(2, 2)); jagsUI::traceplot(out14.4N)      # not shown
print(out14.4N, 3)                           # not shown

# Seems to look good :)
