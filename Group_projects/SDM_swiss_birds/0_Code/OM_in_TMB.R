# Course Bayesian thinking - UZH 2021
# Authors: Gesa von Hirschheydt & Andrea Vallejo Vargas


#  Add working path  ####
#-----------------------#

# This code needs to be adjusted

# set working directory
setwd("C:/Users/Gesa/Documents/GitHub/GitHub_cloned_repositories_to_pull_from/Bayesian_thinking_UZH/Group_projects/SDM_swiss_birds/0_Code/")

# define path to data folder
folder <- "C:/Users/Gesa/Documents/GitHub/GitHub_cloned_repositories_to_pull_from/Bayesian_thinking_UZH/Group_projects/SDM_swiss_birds/1_Data/"



#  Load data  ####
#----------------#

kite <- read.table(paste(folder,"kite.csv",sep=""), header=TRUE, sep =";")
sites <- read.table(paste(folder,"sites.csv",sep=""), header=TRUE, sep =";")


# load packages
library(TMB)
library(tidyverse)  # for summarise() and group_by() function
library(robustHD) # for standardize()



#  occupancy template by Andrea  ####
#-----------------------------------#

## Compile and load the model
# dyn.unload(dynlib("occupancy"))
compile("occupancy.cpp")
dyn.load(dynlib("occupancy"))

## Data and parameters
nvisits <- summarise(group_by(kite, site), visit=max(visit, na.rm=T))$visit
cumdet <- summarise(group_by(kite, site), zst=max(det, na.rm=T))$zst # cumulative detection

elev.st <- standardize(sites$elevation)
# elev.st <- (sites$elevation-mean(sites$elevation))/sd(sites$elevation)
X <- as.matrix(cbind(rep(1,nrow(sites)), elev.st, sites$farmland)); dim(X); head(X)

time.st <- standardize(kite$time)
Z <- as.matrix(cbind(rep(1,nrow(kite)), time.st)); dim(Z); head(Z)

str(data <- list(y=kite$det,
             J=nvisits,
             X=X,
             Z=Z,
             nd=(-1)*(cumdet-1))) # nd is now exactly the opposite (1/0 reversed) of cumdet
(parameters <- list(beta=rep(0,ncol(X)),
                    alpha=rep(0, ncol(Z))))

## Make a function object
obj <- MakeADFun(data, parameters, DLL="occupancy")

## Call function minimizer
(opt <- nlminb(obj$par, obj$fn, obj$gr))

## Get parameter uncertainties and convergence diagnostics
(sdr <- sdreport(obj))
summary(sdr, "fixed", p.value=TRUE)
(AIC <- 2*opt$objective + 2*length(opt$par))
(BIC <- 2*opt$objective + log(nrow(kite))*length(opt$par))

(report <- obj$report())
