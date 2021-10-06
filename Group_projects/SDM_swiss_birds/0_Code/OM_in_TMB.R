# Course Bayesian thinking - UZH 2021
# Authors: Gesa von Hirschheydt & Andrea Vallejo Vargas
# Supervisors: Andrea Havron, Leila Schuh


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
sites.sdiv <- read.table(paste(folder,"sites.sdiv.csv",sep=""), header=TRUE, sep =";")
colnames(sites.sdiv)[25:27] <- c("ent","ent1","ent2")

# round(cor(sites.sdiv[,c(3:5,7,25:27)]),2)
# plot(forest ~ elevation, data=sites.sdi, pch=16)
cov_names <- c("forest","ent","ent1","ent2","farmland","elevation")
covariates <- sites.sdiv[,cov_names]; head(covariates)
par(mfrow=c(1,2))
for(i in c(2:4,6)){
  covariates[,i] <- standardize(sites.sdiv[,cov_names[i]])
}; head(covariates)

# load packages
library(TMB)
library(tidyverse)  # for summarise() and group_by() function
library(robustHD) # for standardize()

# create table to store model AIC
comparison <- data.frame(model=c("-","random site effect",
                                 "forest","forest + ent","forest + ent1","forest + ent2",
                                 "best1 + farmland","best2 + elevation","best3 + time[p]"),
                         AIC=NA, BIC=NA,
                         signif=NA)



#  ACTIVE MODEL  ####
#-------------------#

## Compile and load the model
# dyn.unload(dynlib("occupancy"))
compile("occupancy.cpp")
dyn.load(dynlib("occupancy"))

## Data and parameters
nvisits <- summarise(group_by(kite, site), visit=max(visit, na.rm=T))$visit
cumdet <- summarise(group_by(kite, site), zst=max(det, na.rm=T))$zst # cumulative detection
# site_y <- as.numeric(as.factor(kite$site)); site_y <- site_y-1; head(site_y)

cov_names <- c("forest","ent","ent1","ent2","farmland","elevation")
X <- as.matrix(cbind(rep(1,nrow(sites)), covariates[,cov_names])); dim(X); head(X)

time.st <- standardize(kite$time)
Z <- as.matrix(cbind(rep(1,nrow(kite)), time.st)); dim(Z); head(Z)

str(data <- list(y=kite$det,
             J=nvisits,
             X=X,
             Z=Z,
             nd=(-1)*(cumdet-1))) # nd is now exactly the opposite (1/0 reversed) of cumdet
(parameters <- list(beta=rep(0,ncol(X)),
                    alpha=rep(0, ncol(Z)),
                    log_sig_u=0,
                    u=rep(0,nrow(X))))

## Make a function object
obj <- MakeADFun(data, parameters, random="u", DLL="occupancy") # if fitting a random effect --> Laplace approximation
# obj <- MakeADFun(data, parameters, DLL="occupancy")

## Call function minimizer
(opt <- nlminb(obj$par, obj$fn, obj$gr))
# (sim_y <- obj$simulate())

## Get parameter uncertainties and convergence diagnostics
(sdr <- sdreport(obj))
summary(sdr, "fixed", p.value=TRUE)
summary(sdr, "random", p.value=TRUE)
(AIC <- 2*opt$objective + 2*length(opt$par))
(BIC <- 2*opt$objective + log(nrow(kite))*length(opt$par))

(report <- obj$report())



#  empty occupancy model  ####
#----------------------------#

## Compile and load the model
# dyn.unload(dynlib("occupancy_empty"))
compile("occupancy_empty.cpp")
dyn.load(dynlib("occupancy_empty"))

## Data and parameters
nvisits <- summarise(group_by(kite, site), visit=max(visit, na.rm=T))$visit
cumdet <- summarise(group_by(kite, site), zst=max(det, na.rm=T))$zst # cumulative detection

str(data <- list(y=kite$det,
                 J=nvisits,
                 nd=(-1)*(cumdet-1))) # nd is now exactly the opposite (1/0 reversed) of cumdet
(parameters <- list(logit_psi=0,
                    logit_p=0))

## Make a function object
obj <- MakeADFun(data, parameters, DLL="occupancy_empty")

## Call function minimizer
(opt <- nlminb(obj$par, obj$fn, obj$gr))

## Get parameter uncertainties and convergence diagnostics
(sdr <- sdreport(obj))
summary(sdr, "fixed", p.value=TRUE)
(report <- obj$report())
(AIC <- 2*opt$objective + 2*length(opt$par))
(BIC <- 2*opt$objective + log(nrow(kite))*length(opt$par))

# save AIC/BIC for model selection
index <- which(comparison$model=="-")
comparison[index,c("AIC","BIC")] <- round(c(AIC, BIC),0)
comparison[index,"signif"] <- "-"


#  random site effect on occ  ####
#--------------------------------#

## Compile and load the model
# dyn.unload(dynlib("occupancy_random1"))
compile("occupancy_random1.cpp")
dyn.load(dynlib("occupancy_random1"))

## Data and parameters
nvisits <- summarise(group_by(kite, site), visit=max(visit, na.rm=T))$visit
cumdet <- summarise(group_by(kite, site), zst=max(det, na.rm=T))$zst # cumulative detection

str(data <- list(y=kite$det,
                 J=nvisits,
                 X=as.matrix(rep(1,nrow(sites))),
                 nd=(-1)*(cumdet-1))) # nd is now exactly the opposite (1/0 reversed) of cumdet
(parameters <- list(logit_p=0,
                    log_sig_u=0,
                    u=rep(0,nrow(sites))))

## Make a function object
obj <- MakeADFun(data, parameters, random="u", DLL="occupancy_random1")

## Call function minimizer
(opt <- nlminb(obj$par, obj$fn, obj$gr))

## Get parameter uncertainties and convergence diagnostics
(sdr <- sdreport(obj))
summary(sdr, "fixed", p.value=TRUE)
summary(sdr, "random", p.value=TRUE) # gives strange output... why only 2 values?
(report <- obj$report())
(AIC <- 2*opt$objective + 2*length(opt$par))
(BIC <- 2*opt$objective + log(nrow(kite))*length(opt$par))

# save AIC/BIC for model selection
index <- which(comparison$model=="random site effect")
comparison[index,c("AIC","BIC")] <- round(c(AIC, BIC),0)
# comparison[index,"signif"] <- "-"


#  fixed effects on occ only  ####
#--------------------------------#

## Compile and load the model
# dyn.unload(dynlib("occupancy_fixed1"))
compile("occupancy_fixed1.cpp")
dyn.load(dynlib("occupancy_fixed1"))

## Data and parameters
nvisits <- summarise(group_by(kite, site), visit=max(visit, na.rm=T))$visit
cumdet <- summarise(group_by(kite, site), zst=max(det, na.rm=T))$zst # cumulative detection

#  forest  ##
#-----------#
cov_names <- c("forest") # ,"ent","ent1","ent2","farmland","elevation"
X <- as.matrix(cbind(rep(1,nrow(sites)), covariates[,cov_names])); dim(X); head(X)

str(data <- list(y=kite$det,
                 J=nvisits,
                 X=X,
                 nd=(-1)*(cumdet-1))) # nd is now exactly the opposite (1/0 reversed) of cumdet
(parameters <- list(logit_p=0,
                    beta=rep(0,dim(X)[2])))

## Make a function object
obj <- MakeADFun(data, parameters, DLL="occupancy_fixed1")

## Call function minimizer
(opt <- nlminb(obj$par, obj$fn, obj$gr))

## Get parameter uncertainties and convergence diagnostics
(sdr <- sdreport(obj))
summary(sdr, "fixed", p.value=TRUE)
(report <- obj$report())
(AIC <- 2*opt$objective + 2*length(opt$par))
(BIC <- 2*opt$objective + log(nrow(kite))*length(opt$par))

# save AIC/BIC for model selection
index <- which(comparison$model=="forest")
comparison[index,c("AIC","BIC")] <- round(c(AIC, BIC),0)
# comparison[index,"signif"] <- "none"

