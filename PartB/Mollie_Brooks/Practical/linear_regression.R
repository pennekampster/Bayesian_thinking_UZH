library(TMB)

#Create data
x <- -5:5
atrue <- 1
btrue <- .5
set.seed(111)
y <- rnorm(n = length(x), mean = atrue  + btrue *x)

#Organize data for TMB
dat <- list(x=x, y=y)

#Set initial values for parameters to estimate
par <- list(a=0, b=0, log_sigma=0)

#Compile the model
compile("linear_regression.cpp")
dyn.load(dynlib("linear_regression"))
  
#Combine data, parameters, and model to make objective function (i.e. nll)
obj <- MakeADFun(dat, par, DLL="linear_regression")

#Optimize the objective function (i.e. minimize the nll)
opt <- nlminb(obj$par, obj$fn, obj$gr)

#Output results
summary(sdreport(obj))

yfit <- obj$report()$yfit

plot(x,y)
points(x, yfit, type="l")







