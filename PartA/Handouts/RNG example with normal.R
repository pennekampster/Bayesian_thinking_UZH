

# This is an anaology for how statistical inference via MCMC works:
# We draw random numbers from a distribution which we typically do not even know
# ... and from those random numbers we describe features of the unknown distribution
# ... such as its mean, variance, percentiles.
# Here we illustrate this with the RNG for normal numbers in R
# Can imagine that we didn't know the mean and the sd of the distribution, but
# ... want to use a sample of values produced by the RNG for that distribution
# ... to estimate it.  


par(mfrow = c(2,2), mar = c(4,3,2,1))

# 10 draws from distribution
out <- rnorm(10, mean = 5, sd = 2)
hist(out, breaks = 50, col = "grey", main = 'Cannot believe it is a normal !')
mean(out)
sd(out)


# 100 draws from distribution
out <- rnorm(100, mean = 5, sd = 2)
hist(out, breaks = 50, col = "grey", main = 'Well, perhaps ...')
mean(out)
sd(out)


# 1,000 draws from distribution
out <- rnorm(1000, mean = 5, sd = 2)
hist(out, breaks = 50, col = "grey", main = 'Oh !')
mean(out)
sd(out)


# 1,000,000 draws from distribution
#out <- rnorm(10^6, mean = 5, sd = 2)
#hist(out, breaks = 100, col = "grey", main = 'Now THAT is pretty good !')
#mean(out)
#sd(out)


# 100,000,000 draws from distribution
out <- rnorm(10^7, mean = 5, sd = 2)
hist(out, breaks = 1000, col = "grey", main = 'Darn good!')
mean(out)
sd(out)


# 100,000,000 draws from distribution
out <- rnorm(10^8, mean = 5, sd = 2)
hist(out, breaks = 1000, col = "grey", main = 'Darn good!')
mean(out)
sd(out)




norm.dens <- function(y, mu = 750, sigma = 187){
#   (1 / sqrt(2 * pi * sigma)) * exp(-((y-mu)^2) / (2 * sigma^2))
   CONST <- 1 / sqrt(2 * pi * sigma)
   EXP <- -((y-mu)^2) / (2 * sigma^2))
   all <- CONST * exp(EXP)
   all
}



