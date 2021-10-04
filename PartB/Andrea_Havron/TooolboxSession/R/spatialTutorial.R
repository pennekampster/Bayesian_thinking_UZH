############################################################################
## Preliminaries
############################################################################

#install packages from CRAN
install.packages('TMB','mvtnorm','lattice', 'tmbstan')

#install INLA
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

#if using RStudio, first setup RStudio for TMB. This will help with debugging
TMB:::setupRStudio()

#setup TMB
#TMB model is a C++ file and therefore needs to be compiled
#In order to compile TMB models on Windows, you first need to install Rtools. 
#See https://cran.r-project.org/bin/windows/Rtools/
TMB::compile('src/spatial.cpp')
#Compiling the C++ file will create a .dll file. 
#This .dll file is the model in binary, executable form and 
# needs to be loaded into the R environment for the TMB model to run
dyn.load(TMB::dynlib('src/spatial'))

#This .dll file needs to be unloaded before compiling the code again
# dyn.unload(TMB::dynlib('src/spatial'))

#When compiling, if you get:
#error: ld returned 1 exit status
#it means you need to unload the .dll and try again
#but make sure to alwyas load the .dll again before running the model!

#To see loaded .dll files:
getLoadedDLLs()

######################################################################
## Example - Spatial Poisson model with and without a covariate
######################################################################
library(TMB)
library(INLA)
library(mvtnorm)
library(lattice)
library(tmbstan)
#read in simulation functions
source('R/simdata.R')

#No covariate example
#Spatial extent is: x=c(0,100); y=c(0,100)
#Spatial range (distance at which correlation is neglible)
#Spatial variance (variability of spatial structure)

#Simulate data: - this will take a few minutes to run!
#True parameters
Range <- 50
sp.var <- 1
Beta <- 2

#If this takes too long, modify the spatial extent of the grid (Line 57 in simdata)
dat <- sim.spatial.dat(cov.mod = 'nonspatial', beta = Beta, Range = Range,
               sig2 = sp.var, seed = 5)

#View data using levelplot from the lattice package
lattice::levelplot(dat$Y~dat$loc.x*dat$loc.y)

#Sample data from the grid
samp <- samp.dat(200, dat)

#Build INLA mesh. For details on this step, see 2.6 and 2.7 in https://becarioprecario.bitbucket.io/spde-gitbook/index.html
Loc <- samp[,1:2]
mesh <- inla.mesh.2d(Loc, max.edge = c(10,20), offset = c(5,25))
plot(mesh)
#calculate sparse distance matrix components
spde <- inla.spde2.matern(mesh)

#create data list and initial parameter list
Dat <- list(y = samp$Y,
            x = matrix(samp$X,ncol=length(Beta)),
            v_i = mesh$idx$loc-1, #maps mesh location to obs. C++ start indexing at 0, so need to -1
            M0 = spde$param.inla$M0, #sparse distance matrix
            M1 = spde$param.inla$M1, #sparse distance matrix
            M2 = spde$param.inla$M2) #sparse distance matrix
Par <- list(beta = rep(0,length(Beta)),
            ln_kappa = 0,
            ln_tau = 0,
            omega = rep(0,mesh$n)) #spatial random effects estimated for each mesh vertex


#Fit TMB model
#load .dll file
dyn.load(TMB::dynlib('src/spatial'))
obj <- MakeADFun(Dat,Par,random='omega',DLL = 'spatial')#add silent = TRUE to stop printing to console
## if you get the following error:
#    Error in .Call("getParameterOrder", data, parameters, new.env(), PACKAGE = DLL) : 
#      "getParameterOrder" not available for .Call() for package "spatial"
##  the .dll is not loaded into the R environment, use dyn.load(dynlib('spatial')) to load the .dll library
opt <- nlminb(obj$par, obj$fn, obj$gr)

#Reports
sdr <- sdreport(obj)
summary(sdr, 'fixed', p.value = TRUE)
summary(sdr, "report")
report <- obj$report(obj$env$last.par.best)
#report marginal negative log likelihood
opt$objective

#Check model convergence
opt$convergence #0: converged, 1:not converged
sdr$pdHess #True means Hessian is positive definite


#Spatial predictions
#INLA projection matrix
P <- inla.mesh.project(mesh, as.matrix(dat[,1:2]))
omega.proj <- P$A %*% report$omega
eta.proj <- opt$par['beta'] + omega.proj
lambda.proj <- exp(eta.proj)
levelplot(omega.proj~dat$loc.x*dat$loc.y)
levelplot(lambda.proj~dat$loc.x*dat$loc.y)
#compare with true data
levelplot(dat$Y~dat$loc.x*dat$loc.y)

a <- Sys.time()
#Run model in tmbstan
stan.mod <- tmbstan(obj)
b <- Sys.time()