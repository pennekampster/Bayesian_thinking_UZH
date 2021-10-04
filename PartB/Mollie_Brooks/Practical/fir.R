## Load TMB TRUE
library(TMB)
load("fir.Rdata")
## Compile and load the model
compile("fir.cpp")
dyn.load(dynlib("fir"))

## Data and parameters
Ldat = list(dbh=firdat$DBH,
            cones= firdat$cones,
            wave_non=model.matrix(~wave_non, data=firdat)
)

Lpin = list(a=c(1,0), b=c(1,0), log_k=1)

## Make a function object
obj = MakeADFun(Ldat, Lpin, DLL="fir")

## Call function minimizer
opt = nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr = sdreport(obj)
summary(sdr)
