
load("spider_ Alopacce.Rdata") #Data borrowed from mvabund package
citation("mvabund")

## Load TMB TRUE
library(TMB)

## Compile and load the model
compile("Alopacce.cpp")
dyn.load(dynlib("Alopacce"))

## Data and parameters

dat <- list(abund = Alopacce$abund,
						x = Alopacce$fallen.leaves)

with(Alopacce, plot(fallen.leaves, abund))

pars <- list(coefs=c(0,0), log_k=0)

## Make a function object
obj <- MakeADFun(dat, pars, DLL="Alopacce")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

# Compare with a packaged GLM
library(glmmTMB)

mod <- glmmTMB(abund ~ fallen.leaves, family="nbinom2", data =Alopacce)
mod

# What if we want to do multiple regression? i.e. include more covariates than fallen.leaves
#Use model.matrix

?model.matrix

Xmat <- model.matrix(~fallen.leaves + moss , data=Alopacce)
head(Xmat)

dat2 <- list(abund = Alopacce$abund,
						Xmat = Xmat)

pars2 <- list(coefs = rep(0, ncol(Xmat)), log_k=0)

compile("Alopacce_mm.cpp")
dyn.load(dynlib("Alopacce_mm"))

## Make a function object
obj <- MakeADFun(dat2, pars2, DLL="Alopacce_mm")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

# Compare with a packaged GLM

mod <- glmmTMB(abund ~ fallen.leaves +moss, family="nbinom2", data =Alopacce)
mod
