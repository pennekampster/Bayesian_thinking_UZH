library(TMB)

#Simulate data
set.seed(111)
trap_sd <- 8 # std dev among traps
date_sd <- 3 # std dev among dates
resid_sd <- 1 # residual normal std dev
ndates <- 30 # number of dates
ntraps <- 4 # number of traps
mu <- 5 # mean of control group 
trt_eff <- 10 #treatment effect
dev_trap <- rnorm(ntraps, 0, trap_sd) # random effect of trap
dev_date <- rnorm(ndates, 0, date_sd) # random effect of date

dat <- expand.grid(trap=factor(1:ntraps), date=factor(1:ndates), trt=c(0,1))
dat0 <- transform(dat, y=mu+trt_eff*trt + dev_trap[trap] + dev_date[date])
dat0$y <- dat0$y+rnorm(nrow(dat0),0,resid_sd)

library(ggplot2)
dat2 = transform(dat0, treatment=factor(trt, labels=c("control", "manipulated")))
p1 <- ggplot(dat2, aes(treatment, y, group=trap:date))+geom_line(aes(colour=factor(date)))
p1
p1 + facet_grid(~trap, labeller = label_both)+ theme(axis.text.x = element_text(angle = 45))
ggsave("LMM.png", height=4, width=6)

X <- model.matrix(~-1+trap+trt, data=dat0) # fixed effects w/o an intercept
Z <- model.matrix(~-1+date, data=dat0) # random effects w/o an intercept

dat <- list(y=dat0$y, X=X, Z=Z, ntraps=ntraps)
pars <- list(log_re_sd=0, log_resid_sd=1, beta=rep(0, ncol(X)), dev=rep(0, ncol(Z)))

compile("LMM.cpp")
dyn.load(dynlib("LMM"))
obj <- MakeADFun(data = dat,
							parameters=pars,
							random=c("dev"), DLL="LMM")
opt  <-  nlminb(obj $par, obj $fn, obj $gr)
sdr <- sdreport(obj) 
sdr

as.list(sdr,"Est")

# Compare with a package
library(glmmTMB)

mod <- glmmTMB(y ~ -1 + trt + trap + (1|date), data =dat0) #different data format
fixef(mod)
ranef(mod)

#Compare random effect estimates from both versions
plot(c(ranef(mod)$cond$date)[[1]], unname(sdr$par.random))
