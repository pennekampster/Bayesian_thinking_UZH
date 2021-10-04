# Course Bayesian thinking - UZH 2021
# Authors: Gesa von Hirschheydt & Andrea Vallejo Vargas


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


# load packages
library(TMB)
library(tidyverse)  # for summarise() and group_by() function



#  occupancy template by Andrea  ####
#-----------------------------------#

## Compile and load the model
compile("occupancy.cpp")
dyn.load(dynlib("occupancy"))

## Data and parameters
nvisits <- summarise(group_by(kite, site), visit=max(visit, na.rm=T))$visit
cumdet <- summarise(group_by(kite, site), zst=max(det, na.rm=T))$zst

str(data <- list(y=kite$det,
             J=nvisits,
             knownOcc=cumdet,
             nd=(-1)*(cumdet-1)   # nd is now exactly the opposite (1/0 reversed) of knownOcc
             ))
str(parameters <- list(logit_psi=rep(0,times=length(nvisits)),
                       logit_p=rep(0,times=sum(nvisits))))

## Make a function object
obj <- MakeADFun(data, parameters, DLL="occupancy")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)  # takes a long while...
opt # --> no convergence?

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj) # this seems to never end...
sdr
