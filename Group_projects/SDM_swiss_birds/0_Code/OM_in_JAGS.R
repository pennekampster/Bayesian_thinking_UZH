# Course Bayesian thinking - UZH 2021
# Author: Gesa von Hirschheydt


#  Add working path  ####
#-----------------------#

# This code needs to be adjusted

# set working directory
setwd("H:/GitHub/Git_repositories_to_pull_from/Bayesian_thinking_UZH/Group_projects/SDM_swiss_birds/0_Code")

# define path to data folder
folder <- "H:/GitHub/Git_repositories_to_pull_from/Bayesian_thinking_UZH/Group_projects/SDM_swiss_birds/1_Data/"



#  Load data  ####
#----------------#

kite <- read.table(paste(folder,"kite.csv",sep=""), header=TRUE, sep =";")
kite_for_jags <- read.table(paste(folder,"kite_for_jags.csv",sep=""), header=TRUE, sep =";")


# load packages
library(tidyverse)  # for summarise() and group_by()
library(jagsUI)  # for jags()

nvisits <- summarise(group_by(kite, site), visit=max(visit, na.rm=T))$visit
cbind(kite_for_jags, nvisits) # for a quick check whether sites are in correct order



#  MODEL0 in JAGS  ####
#---------------------#

# Bundle data and summarize data bundle
str(data <- list(y=kite_for_jags,
                 nsites=dim(kite_for_jags)[1],
                 nvisits=nvisits))

# Specify model in BUGS language
sink("model0.txt")
cat("model{
    
    # Priors
    psi ~ dunif(0,1)
    p ~ dunif(0,1)
    
    # Likelihood
    for(i in 1:nsites){
      z[i] ~ dbern(psi)
    
      for(k in 1:nvisits[i]){
        y[i,k] ~ dbern(z[i]*p)
    }
   }
 }
", fill=TRUE)
sink()


# Initial values
zst <- summarise(group_by(kite, site), zst=max(det, na.rm=T))$zst   # to avoid data/model/inits conflict
inits <- function(){list(z=zst)}

# Parameters monitored
params <- c("psi","p")

# MCMC settings
nc <- 4; nb <- 1000; ni <- 5000; nt <- 1

# Call JAGS and summarize posteriors
output <- jags(data, inits, params, "model0.txt",
            n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb)
# save(output, file="output.model0.Rdata")
load(file="output.model0.Rdata") # it will have the name "output"
print(output, dig=2)


# Visualize output
out <- data.frame(output$sims.list)
params.titles <- colnames(out)

# Traceplots
par(mfrow=c(3,1), mar=c(3,3,2,1))
chain.length <- (ni-nb)/nt
for(i in 1:length(params.titles)){
  plot(x=1:chain.length, y=out[1:chain.length,i], type="l", col=2,
       main=params.titles[i], xlab="", ylab="", las=1, frame.plot=FALSE)
  for(j in 2:nc){
    lines(x=1:chain.length, y=out[(chain.length*(j-1)+1):(chain.length*j),i], col=(j+1))
  }
  abline(h=mean(out[,i]), lty=2, col=1)
}


# Overview posterior densities
par(mfrow=c(2,2), mar=c(5,1,3,3), yaxt="n")
for(i in 1:4){   # length(params.titles)
  dens <- density(out[,i])
  CIlower <- quantile(out[,i],probs=c(0.025))
  CIupper <- quantile(out[,i],probs=c(0.975))
  CIrange.y <- dens[[2]][dens[[1]]>CIlower & dens[[1]]<CIupper]
  CIrange.x <- seq(CIlower, CIupper, length.out=length(CIrange.y))
  plot(dens, bty="n", lwd=1, main=paste(params.titles[i]),
       xlab="Parameter estimate", ylab="", # yaxt="n", xaxt="n",
       xlim=c(min(dens$x), max(dens$x)),
       ylim=c(0, max(dens$y)))
  library(scales)
  polygon(x=c(CIrange.x, rev(CIrange.x)), y=c(CIrange.y, rep(0, length(CIrange.y))),
          border=NA, col=alpha("grey", 0.5))
  points(x=mean(out[,i]), y=0, pch="|")
  abline(v=0, col="blue", lty=2)
  text(x=0, y=max(dens$y), col="blue",
       labels=paste(round(1-ecdf(out[,i])(0),2)), pos=4) # right tail; left tail: round(ecdf(...)), pos=2
}
par(mfrow=c(1,1), mar=c(5,5,4,2), yaxt="s")




