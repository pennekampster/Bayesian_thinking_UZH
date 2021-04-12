##############################################################
###
### Frequentist and Bayesian analyses of a binomial proportion
###
##############################################################

### Code examples taken from the books of Royle & Dorazio (2008), 
#   Ntzoufras (2009), and Schaub and Kéry (2021)

# Last changes: 7 April 2021

# The tadpole example: we released 50 pollywogs in the bucket and later counted 20
# Define the data
r <- 20
N <- 50

# Main idea of maximum likelihood: 
# Use probability to make a guess at likely values of detection probability
# Parameter values that are more likely to produce the observed data are 
# a better guess when you want to pick an estimate ...

# Example: Fig. 2.3 (IPMbook)
all.possible <- 0:50           # All possible values of the observed data (the count r)
pmf1 <- dbinom(all.possible, size=N, prob=0.1)  # pmf 1
pmf2 <- dbinom(all.possible, size=N, prob=0.5)  # pmf 2
pmf3 <- dbinom(all.possible, size=N, prob=0.9)  # pmf 3

par(mfrow=c(1, 3), mar=c(5, 5, 4, 1), cex.lab=1.5, cex.axis=1.5, cex.main=2, las=1)
plot(all.possible, pmf1, type="h", lend="butt", lwd=3, frame=FALSE, xlab="Counts (y)", ylab="Probability of y", main=expression(paste(theta, " = 0.1")))
abline(v=20, col="blue", lwd=2)
plot(all.possible, pmf2, type="h", lend="butt", lwd=3, frame=FALSE, xlab="Counts (y)", ylab="", main=expression(paste(theta, " = 0.5")))
abline(v=20, col="blue", lwd=2)
plot(all.possible, pmf3, type="h", lend="butt", lwd=3, frame=FALSE, xlab="Counts (y)", ylab="", main=expression(paste(theta, " = 0.9")))
abline(v=20, col="blue", lwd=2)


### (1) Maximum likelihood (brute force approach)
#################################################
# = "just try out what's gives the highest function value")

# From chapter 2 in IPM book

# Brute-force search for MLE of theta for the tadpole data set (Fig. 2.4)
(try.theta <- seq(0, 1, by=0.01) )                 # Values to try out
like <- dbinom(20, 50, try.theta, log=FALSE)       # Likelihood
loglike <- dbinom(20, 50, try.theta, log=TRUE)     # Log-Likelihood
negloglike <- -dbinom(20, 50, try.theta, log=TRUE) # Negative log-likelihood (NLL)

par(mfrow=c(1, 3), mar=c(5, 5, 4, 1), cex.lab=1.5, cex.axis=1.5)
plot(x=try.theta, y=like, xlab=expression(paste("Detection probability (", theta, ")")), ylab="Likelihood", frame=FALSE, type="p", pch=16, col="black")
abline(v=try.theta[which(like == max(like))], col="red", lwd=3)
plot(x=try.theta, y=loglike, xlab=expression(paste("Detection probability (", theta, ")")), ylab="Log-Likelihood", frame=FALSE, type="p", pch=16, col='black')
abline(v=try.theta[which(loglike == max(loglike))], col="red", lwd=3)
plot(x=try.theta, y=negloglike, xlab=expression(paste("Detection probability (", theta, ")")), ylab="Negative log-Likelihood", frame=FALSE, type="p", pch=16, col="black")
abline(v=try.theta[which(negloglike == min(negloglike))], col="red", lwd=3)


# Check integral: likelihood function DOESN't sum to 1
# IS NOT a pdf !
sum(like)


### (2) Maximum likelihood (numerical optimization of function)
###############################################################

# Define negative log-likelihood function
nll1 <- function(p) -dbinom(r, size = N, prob = p, log = TRUE)
# Minimize function for observed data and return MLE
(fit1 <- optim(par = 0.5, fn = nll1, method = "BFGS"))


## Try the same, do things 'by hand' to be more explicit
# Note the binomial likelihood:
L(p| r, N) = (N! / (r! (N-r)! ) * p^r * (1-p)^(N-r)

# Can ignore combinatorial term, since constant, hence
L(p| r, N) = p^r * (1-p)^(N-r)

# Take log for log-likelihood
LL(p | r, N) = log(p^r * (1-p)^(N-r)) =
               r * log(p) + (N-r) * log(1-p)
			   
# Negate for negative log-likelihood
NLL(p | r, N) = -( log(p^r * (1-p)^(N-r))) =
               - (r * log(p) + (N-r) * log(1-p) )

nll2 <- function(p) -(r * log(p) + (N-r)*log(1-p))
(fit2 <- optim(par = 0.5, fn = nll2, method = "BFGS"))

cat("Maximum likelihood estimate of p (using dbinom): ", fit1$par, "\n\n")
cat("Maximum likelihood estimate of p (by hand): ", fit2$par, "\n\n")



### (3) Maximum likelihood (using canned function glm())
########################################################

# Estimate parameter on link scale
fm <- glm(cbind(r,N-r) ~ 1, family = binomial)
summary(fm)
predict(fm, type = "response", se = TRUE)








### (4) Bayesian inference using a random-walk Metropolis algorithm
###################################################################

# From chapter 2 in the IPM book and using 
# function 'demoMCMC' in the IPMbook package (Thanks to Mike Meredith)

# Choose initial value for logit(theta) and tuning parameters
ltheta1 <- 1        # Initial value for tadpole detection prob.
sigma_prop <- 1     # SD of Gaussian proposal distribution

# Array to hold the MCMC samples for logit(theta)
ltheta <- numeric()

# Initial value becomes first (and 'current') value in the chain
ltheta[1] <- ltheta1
ltheta              # Our posterior sample up to now (not shown)

# Randomly perturb the current value
set.seed(1)         # To initalize your RNGs identically to ours
( ltheta_star <- rnorm(1, ltheta[1], sigma_prop) )  #  [1] 0.3735462

# Now we have two values for logit(theta) to compare and we compute 
# the ratio of their posterior densities, ignoring the denominator in Bayes rule.

# Compute likelihood times prior evaluated for the proposed new value of ltheta
( pd_star <- dbinom(20, 50, plogis(ltheta_star)) * dbeta(plogis(ltheta_star), 1, 1) )
# [1] 0.002716919

# Compute likelihood times prior evaluated for the current value of ltheta
( pd_1 <- dbinom(20, 50, plogis(ltheta[1])) * dbeta(plogis(ltheta[1]), 1, 1) )
# [1] 6.951277e-07

# Compute posterior density ratio R
( R <- pd_star / pd_1 )     # [1] 3908.518

# The proposed new value for ltheta yields a much larger value of the product 
# of likelihood times prior (R >> 1), and so we accept it into our MCMC sample right away.

# Add theta_star into MCMC sample to become ltheta[2]
ltheta[2] <- ltheta_star
ltheta              # Our posterior sample up to now (not shown)

# Our posterior sample of ltheta now has two values, 1 and 0.3735. 
# We repeat these steps by perturbing ltheta[2] and drawing a new proposal for logit(theta) :
( ltheta_star <- rnorm(1, ltheta[2], sigma_prop) )       # [1] 0.5571895

# We compare again the posterior densities for ltheta[2], which is the current value, 
# and the new proposal, by building their ratio
pd_star <- dbinom(20, 50, plogis(ltheta_star)) * dbeta(plogis(ltheta_star), 1, 1)
pd_t <- dbinom(20, 50, plogis(ltheta[2])) * dbeta(plogis(ltheta[2]), 1, 1)
( R <- pd_star / pd_t )      #  [1] 0.1398872

# Now R < 1, so the new proposal for ltheta has a much lower posterior density than 
# the current value. Therefore, we only accept the proposed new value with a probability 
# of about 14% and with a probability of about 86%, we duplicate the current value and 
# add it to the chain. We flip a coin to decide whether we add the new ltheta_star to 
# our sample or whether we push forward the current value to become the new value. 

( keep.ltheta_star <- rbinom(1, 1, R) )        # [1] 0

# Bad luck for this theta_star: it is rejected and hence the ‘old’ value of theta in 
# the chain is plugged in as the new, current value in the chain instead.
ltheta[3] <- ltheta[2]
ltheta              # Our posterior sample up to now (not shown)

# Now we would continue by perturbing theta[3] according to the Normal proposal distribution
# and so forth. Fig. 2.7 (left) shows our nascent Markov chain that gives our posterior 
# samples for detection probability (on the logit scale) in the tadpole example. 
# The following loop allows us to extend this code of our RW-MH algorithm for an arbitrary 
# number of steps (here, we choose 100,000).

# Iteration 4 to T of RW-MH algorithm
T <- 100000            # Choose chain length
for (t in 4:T){       # Continue where we left off
  if(t %% 10000 == 0) {print(t)}
  ltheta_star <- rnorm(1, ltheta[t-1], sigma_prop) 
  pd_star <- dbinom(20, 50, plogis(ltheta_star)) * dbeta(plogis(ltheta_star), 1, 1) 
  pd_t <- dbinom(20, 50, plogis(ltheta[t-1])) * dbeta(plogis(ltheta[t-1]), 1, 1)
  R <- min(1, pd_star / pd_t)  # Note more general solution here
  keep.ltheta_star <- rbinom(1, 1, R)
  ltheta[t] <- ifelse(keep.ltheta_star == 1, ltheta_star, ltheta[t-1])
  # ltheta              # Our posterior sample up to now (not shown)
}

# This takes almost no time at all ! 
# We plot the results for our Bayesian estimate of tadpole detection probability (Fig. 2.7).

# Fig. 2.7 left
par(mfrow=c(1, 3), mar=c(6, 7, 6, 3), cex.lab=2, cex.axis=2, cex.main=2, las=1)
plot(1:10, plogis(ltheta[1:10]), xlab='Iteration', ylab=expression(theta), type='l', frame=FALSE, lwd=3, main='First ten iterations')
abline(h=0.4, col='red', lwd=2) # The maximum likelihood estimate
abline(h=mean(plogis(ltheta[1:10])), lty=2, col='blue', lwd=2)  # posterior mean

# Update trace-plot of time series of posterior draws (Fig. 2-7 middle)
plot(1:T, plogis(ltheta), xlab='Iteration', ylab= expression(theta), type='l', frame=FALSE, lwd=1, main='All iterations')
abline(h=0.4, col='red', lwd=3) # The maximum likelihood estimate
abline(h=mean(plogis(ltheta)), lty=2, col='blue', lwd=3)  # Posterior mean

# Plot histogram of posterior samples of tadpole detection probability
# Fig. 2-7 right
hist(plogis(ltheta), breaks=50, col='lightgray', xlab=expression(theta), main=expression(bold(paste('Posterior distribution of ', theta))), border=NA)
abline(v=0.4, col='red', lwd=3)    # The maximum likelihood estimate
abline(v=mean(plogis(ltheta)), lty=2, col='blue', lwd=3) # Posterior mean



# Do some experiments with function demoMCMC in the IPMbook package
library(IPMbook)
?demoMCMC               # Check out the function (also check out Chapter 2 from the IPM book)
out <- demoMCMC(y=20, N=50, niter=25000, mu.ltheta=0, sd.ltheta=100, prop.sd=1, init=0)

# Show convergence
# ----------------
tmp <- demoMCMC(y=20, N=50, niter=2500, mu.ltheta=0, sd.ltheta=100, prop.sd=0.1, init=10)

# No convergence within 2500 iterations
tmp <- demoMCMC(y=20, N=50, niter=2500, mu.ltheta=0, sd.ltheta=100, prop.sd=0.1, init=100)

# But convergence is reached after about 3k iterations
tmp <- demoMCMC(y=20, N=50, niter=25000, mu.ltheta=0, sd.ltheta=100, prop.sd=0.1, init=100)

# ... and you get convergence within 2500 iters with longer step length
tmp <- demoMCMC(y=20, N=50, niter=2500, mu.ltheta=0, sd.ltheta=100, prop.sd=1, init=100)


# Look at effects of the tuning parameter (SD of Gaussian proposal density)
# -------------------------------------------------------------------------
# Very, very small step size: very inefficient MCMC
str(out <- demoMCMC(prop.s = 0.01))

# Very small step size: fairly inefficient
str(out <- demoMCMC(prop.s = 0.1))

# Larger than default step size: efficiency goes down
str(out <- demoMCMC(prop.s = 10))

# Much larger step size..... brrrrr !
str(out <- demoMCMC(prop.s = 100))

# Brutally large step size..... ACH !
str(out <- demoMCMC(prop.s = 1000))

# Default step size: pretty good for this case
str(out <- demoMCMC(prop.s = 1))



##### (5) The four analyses of the tadpole data
###############################################
# From Chapter 2 in the IPM book
theta.vals <- seq(0, 1, 0.001)

# Define likelihood function (same for all four analyses)
like <- dbinom(20, 50, theta.vals)
sc.like <- like * (50 + 1)          # Scale likelihood. Multiplied by n + 1 because the area under the curve for trial size 1 is 1/(n + 1). 

# Define four prior distributions
prior0 <- dbeta(theta.vals, 1, 1) 
prior1 <- dbeta(theta.vals, 4, 6)
prior2 <- dbeta(theta.vals, 40, 60)
prior3 <- dbeta(theta.vals, 60, 40)

# Derive four posterior distributions
post0 <- dbeta(theta.vals, 20 + 1, 30 + 1)
post1 <- dbeta(theta.vals, 20 + 4, 30 + 6)
post2 <- dbeta(theta.vals, 20 + 40, 30 + 60)
post3 <- dbeta(theta.vals, 20 + 60, 30 + 40)

library(scales)
co <- viridis_pal(option='E')(20)[c(18, 11, 2)]
lwd <- 3; cx <- 1.5
par(mfrow=c(2, 2), mar=c(5, 5, 4, 2), cex.axis=cx, cex.lab=cx, cex.main=cx)

# Analysis 1 with vague prior
plot(theta.vals, post0, type ="l", col=co[3], xlab="", ylab="Scaled likelihood or density", las=1, frame=FALSE, lwd=2, ylim=c(0, 10))
mtext("Vague prior", side=3, line=0.5, font=2)
lines(theta.vals, sc.like, lty=2, lwd=2, col=co[2])
lines(theta.vals, prior0, lwd=2, col=co[1])
legend(0.5, 10, c("Prior dist.", "Likelihood function", "Posterior dist."), col=co, lty=1, lwd=2, bty="n")

# Analysis 2 with informative prior 1
plot(theta.vals, post1, type="l", lwd=2, col=co[3], xlab="", ylab="", las=1, frame=FALSE, ylim=c(0, 10))
mtext("Informative prior 1", side=3, line=0.5, font=2)
lines(theta.vals, sc.like, lwd=2, col=co[2])
lines(theta.vals, prior1, lty=1, lwd=2, col=co[1])

# Analysis 3 with informative prior 2
plot(theta.vals, post2, type="l", lwd=2, col=co[3], xlab= expression(theta), ylab="Scaled likelihood or density", las=1, frame=FALSE, ylim=c(0, 10))
mtext("Informative prior 2", side=3, line=0.5, font=2)
lines(theta.vals, sc.like, lwd=2, col=co[2])
lines(theta.vals, prior2, lty=1, lwd=2, col=co[1])

# Analysis 4 with informative prior 3
plot(theta.vals, post3, type="l", lwd=2, col=co[3], xlab= expression(theta), ylab='', las=1, frame=FALSE, ylim=c(0, 10))
mtext("Informative prior 3", side=3, line=0.5, font=2)
lines(theta.vals, sc.like, lwd=2, col=co[2])
lines(theta.vals, prior3, lty=1, lwd= 2, col=co[1])


# (6) How information in the data changes our state of knowledge
################################################################
# From Chapter 2 in the IPM book

theta.vals <- seq(0, 1, 0.001)
sc.prior1 <- dbeta(theta.vals, 1, 1)
post1 <- dbeta(theta.vals, 20 + 1, 30 + 1)
sc.post1 <- post1 / max(post1)

co <- c("grey92", "grey60", "black")
plot(theta.vals, sc.prior1, type="l", lwd=5, col=co[1],  xlab=expression(theta), ylab="", xlim=c(0,1), ylim=c(0,1), las=1, frame=FALSE, axes=FALSE)
axis(1)
lines(theta.vals, sc.post1, col=co[2], lwd=5)
abline(v=0.4, col=co[3], lwd=5)
legend(0.6, 0.9, c("Complete ignorance", "Improved state of knowledge", "Certainty (fixed parameter)"), col=co, lty=1, lwd=5, bty="n")























######################################################
# Random walk MCMC for binomial proportion
# (older version, taken from p. 48 of Ntzoufras (2009)
######################################################


### Define function that does RW-MCMC
mcmc.fn <- function(y = 399, N = 845, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 0.1, initial.theta = 0){

# Function that does random-walk MCMC for a binomial proportion (theta = logit(p))
# Code taken from p. 48 in the excellent book by Ntzoufras (2009)

# Initalize calculations
  start.time = Sys.time()		 # Set timer
  y <- y ; N <- N			     # Observed data
  niters <- niters			     # Number of MCMC iterations		
  mu.theta <- mu.theta; s.theta <- s.theta   # "Flat normal" prior
  prop.s <- prop.s			     # Proposal parameter (determines step length of random walk)
  theta <- numeric(niters)		 # Set up vector for posterior draws of theta
  acc.counter <- 0			     # Initialize counter for acceptance
  current.theta <- initial.theta # Initial value for theta

# Start MCMC algorithm
for (t in 1:niters){			# Repeat niters times
   prop.theta <- rnorm(1, current.theta, prop.s)		# Propose a value prop.theta
   loga <- (( prop.theta * y - N * log(1 + exp(prop.theta)))	# Compare likelihood times prior (which is proportional to posterior density)
   - (current.theta * y - N * log(1 + exp(current.theta)))	# between the new (proposed) and the old (current value)
   + dnorm(prop.theta, mu.theta, s.theta, log = TRUE)
   - dnorm(current.theta, mu.theta, s.theta, log = TRUE) )

   u <- log(runif(1))
   if (u < loga){ current.theta <- prop.theta		# If new (proposed) theta leads to a higher product of likelihood times prior,
   							# then take it as new current.theta with prob a
							# Otherwise keep the old value of theta
							
   acc.counter <- acc.counter + 1		#  Counts the number of acceptances
	        }
#   browser()						# If unhashed, allows to inspect values of logu and loga at each iteration
   theta[t] <- current.theta
}
p <- plogis(theta)			# Compute p from logit(p)
acc.prob <- acc.counter/niters		# Acceptance probability
cat("Acceptance probability:", round(acc.prob, 2), "\n")
end.time = Sys.time()			# Stop time
elapsed.time = round(difftime(end.time, start.time, units='secs'), 2)  # Compute elapsed time
cat(paste(paste('Posterior samples drawn in ', elapsed.time, sep=''), ' seconds\n\n', sep='')) # Output run time

par(mfrow = c(2,2))			# Plots of theta=logit(p) and of p
plot(theta, type = "l", ylab = "theta (=logit(p))")		# Plot 1: time-series plot of theta = logit(p)
plot(p, type = "l", ylim = c(0,1))				# Plot 2: time-series plot of p
abline(h = y/N, col = "red")		# Add maximum likelihood estimate
abline(h = mean(p), col = "blue")	# Add posterior mean
hist(p, breaks = 100, col = "grey", main = "", freq = FALSE)	# plot 3: Histogram of posterior with smoothed line
smooth <- density(p, adjust = 2)
lines(smooth$x, smooth$y, type = 'l', lwd = 2, col = "blue")
plot(acf(p, plot = FALSE), main = "", lwd = 3)			# Plot 4: Autocorrelation function plot
} # end of function



### Execute function
####################
mcmc.fn(y = 20, N = 50, niters = 25000, mu.theta = 0, s.theta = 100, 
prop.s = 0.1, initial.theta = 0)

mcmc.fn(y = 20, N = 50, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 1, initial.theta = 0)

mcmc.fn(y = 20, N = 50, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 2, initial.theta = 0)

mcmc.fn(y = 20, N = 50, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 100, initial.theta = 0)

mcmc.fn(y = 20, N = 50, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 0.1, initial.theta = 10)
