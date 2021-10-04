###########################################################
# Simplified approach to test the latent variable approach#
###########################################################

#######################################################
# Create a total number of observations of Adult trees#
#######################################################

# Define total number of species and observations
nSpecies     = 100
totalObs     = 100000
# Create abundance distribution (Lisa&Florian's function)

### Function begin
createAbundance <- function(species, total = 130000){
  
  abundance = rep(1, species)
  
  decay = 20 /(species + 50)
  prob = dexp(1:species, rate = decay)
  abundance = abundance + as.vector(rmultinom(1, total-species, prob = prob))
  
  names(abundance) = 1:species
  return(abundance)
}
### Function end

# Use the function to create the abundance distribution
Abundances  <- createAbundance(nSpecies,total = totalObs)

# Create one large vector with observations (for sampling)
AbundancesV <- vector(mode="double",length=sum(Abundances))
k            = 1
for(i in 1:length(Abundances)){
  AbundancesV[k:(k+Abundances[i])] = i
  k = k+Abundances[i]+1
}    

#######################################################
# Create set of communities that influence focal plots#
#######################################################

# Define total number of communities and individuals in each
nComms       = 100
nCommObs     = 100
noReps       = 10

# Creates matrix with abundance of each species(row) per plot (column)
CommsAT <- data.frame(matrix(0, nrow = nSpecies,ncol = nComms*noReps))
CommsAP <- data.frame(matrix(0, nrow = nCommObs,ncol = nComms))
for(i in 1:nComms){
  CommsAP[,i]  <- sample(AbundancesV, nCommObs, replace = FALSE, prob = NULL)
  for(ii in 1:nSpecies) {
    CommsAT[ii,(noReps*(i-1)+1):(noReps*(i-1)+noReps)] = sum(CommsAP[,i]==ii)
  }
} 

#######################################################
# Calculate no.of seedlings from community composition#
#######################################################

# Define strength of density dependence and 
CNDD         = -0.5 # Simple case, only CNDD and same for all
SeedlPerTree = 1000 # Realized in focal plot case of no CNDD

# Define total number of plots and individuals per plot
nPlotObs     = 20

# Calculate surviving seedling densities in focal plot
CommsS <- round(SeedlPerTree*exp(CNDD*CommsAT))

#######################################################
# Create sample of observed adults in each focal plot #
#######################################################

# Creates matrix with abundance of each species(row) per plot (column)
PlotA   <- data.frame(matrix(0, nrow = nSpecies,ncol = nComms*noReps))
for(i in 1:nComms){
  for(ii in 1:noReps){
    PlotASample <- sample(CommsAP[,i],nPlotObs,replace = FALSE, prob = NULL)
    for(iii in 1:nSpecies) {
      PlotA[iii,(noReps*(i-1)+ii)] = sum(PlotASample==iii)
    }
  }
} 

#######################################################
#Assess error when focal plot adults are the predictor#
#######################################################

# Regress observed seedlings against observed adults
FittedCNDDsLM <- vector(mode="double",length=nSpecies)
for(i in 1:nSpecies) {
  # Correcting for number of trees in community and plot (discussion point)
  PlotADens <- (nCommObs/nPlotObs)*c(t(PlotA[,(noReps*(i-1)+1):(noReps*(i-1)+noReps)]))
  
  # Check to make sure we have the true signal in CommsAT
  #PlotADens <- c(t(CommsAT[,(noReps*(i-1)+1):(noReps*(i-1)+noReps)]))
  
  fit<-glm(c(t(CommsS[,(noReps*(i-1)+1):(noReps*(i-1)+noReps)])) ~ 1 + PlotADens, family = poisson(link = "log"))
  FittedCNDDsLM[i] <- fit$coefficients[2]
}
plot(Abundances/sum(Abundances),FittedCNDDsLM)

#######################################################
# Use data to fit error-in-variable model with Stan   #
#######################################################

# Code mostly from mbjoseph.github.io, shared by Florian
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(1234)

n.reps <- noReps
n.repeated <- nComms
n <- nComms

# run for one species
i = 1

# true covariate values
x <- c(t(CommsAT[i,seq(from = 1, to = n.reps*n, by = n.reps)]))
y <- c(t(CommsS[i,seq(from = 1, to = n.reps*n, by = n.reps)]))

# indx assigns measurements to sample units
indx <- c(1:n, rep(1:n, each = n.reps - 1))
indx <- sort(indx)
nobs <- length(indx)
xobs <- (nCommObs/nPlotObs)*c(t(PlotA[i,]))

# Check that the procedure still retains true signal(not going well yet!)
#xobs <- c(t(CommsAT[i,]))

# Calculate means for error estimation in predictor
xmeans <- vector(mode="double",length=nobs)
for(i in 1:n) {xmeans[((n.reps*(i-1))+1):((n.reps*(i-1))+n.reps)] <- mean(xobs[(1+n.reps*floor((n.reps*(i-1))/n.reps)):(n.reps+n.reps*floor((n.reps*(i-1))/n.reps))])}
xmeans <- c(t(xmeans))
plot(x[indx], xobs,
     xlab = "True covariate value",
     ylab = "Observed covariate value")
abline(0, 1, lty = 2)

# write the .stan file
cat("
data{
  int n;
  int nobs;
  int xobs[nobs];
  int y[n];
  int indx[nobs];
  int noReps;
  real xmeans[nobs];
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigmay;
  real alpha2;
  real beta2;
  real<lower=0> sigmay2;
  real<lower=0> sigmax;
  real <lower=0> x[n];
}
model {
  // priors
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  alpha2 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  sigmay ~ normal(0, 1);
  sigmay2 ~ normal(0, 1);
  sigmax ~ normal(0, 1);
  // model structure  
  for (i in 1:nobs){
    xobs[i] ~ normal(x[indx[i]], sigmax);
  }
  for (i in 1:n){
    y[i] ~ normal(exp(alpha + beta*x[indx[1+(i-1)*noReps]]), sigmay);
  }
  
  for (i in 1:n){
    y[i] ~ normal(exp(alpha2 + beta2*xobs[i]), sigmay2);
  }
}
  ",
file = "latent_x.stan")


stan_d <- c("y", "xobs", "nobs", "n", "indx","xmeans","noReps")
chains <- 3
iter <- 3000
thin <- 1
mod1 <- stan(file = "latent_x.stan", data = stan_d,
             chains = chains, iter = iter,
             thin = thin)

posteriors <- rstan::extract(mod1)

# Show posterior distribution (now for beta = CNDD strength)
# Both the latent variable model and naive model are shown
hist(posteriors$beta, breaks = 30,
     main = "Posterior for CNDD",
     xlab = "Strength of CNDD",
     ylim=c(0,5000))
hist(posteriors$beta2, breaks = 30,
     main = "Posterior for CNDD", add=T, col="#DD000030")

# Show estimates vs true values
dim(posteriors$x)
plot(colMeans(posteriors$x), x)
