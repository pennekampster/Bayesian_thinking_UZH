# ----------------------------------------------------------------------------------------------------------------------------------
# Functional response estimation
# ----------------------------------------------------------------------------------------------------------------------------------
#
# The following R-code of the source functions is based on the method of:
#
# Rosenbaum, B. & Rall, B.C. (2018). Fitting functional responses: Direct parameter estimation by
# simulating differential equations. Methods in Ecology and Evolution.
#
# This method estimates the functional response parameters directly by simulating differential equations
#
# For this, the package odeintr is used.
#
# Here, the method is extended by including the numerical response estimation together with the functional response estimation
#
# Note that for this reason, the starting values for c (conversion efficiency) needs to be provided as well.
#
# For more details about the method, see: https://www.biorxiv.org/content/10.1101/498030v2
#
# For a manual about the use of the method and a detailed description of the functions, 
# see the paper of Rosenbaum and Rall (2018) mentioned above
#
# ----------------------------------------------------------------------------------------------------------------------------------
#
# Uriah Daugaard
#
# Department of Evolutionary Biology and Environmental Studies, University of Zurich, Switzerland
#
# May 2019
#
# ----------------------------------------------------------------------------------------------------------------------------------

library(odeintr)

FR.gen.pred = '
dxdt[0] = -b * pow(x[0],1+q) / (1 + b * h * pow(x[0],1+q)) * x[1] + r * x[0] * (1.0 - x[0]/K);
dxdt[1] = x[1] * c * b * pow(x[0],1+q) / (1 + b * h * pow(x[0],1+q));
'
compile_sys("FR_gen_pred", FR.gen.pred, pars = c("b","h","q","r","K","c"), method = "rk54")

eq.odeint.general.pred = function(start.v, b, h, q, r, K, c, Tt, timesteplength){
  FR_gen_pred_set_params(b = b, h = h, q = q, r=r, K=K, c=c)
  FR_gen_pred(start.v,Tt,timesteplength)
}

eaten.odeint.general.pred = function(N0, b, h, q, r, K, c, Tt, P, steps=100){
  Neaten.est = vector()
  predrep.est = vector()
  for(i.eaten in 1:length(N0)){
    temp = eq.odeint.general.pred(start.v = c(N0[i.eaten], P[i.eaten]),
                                  b = b,
                                  h = h,
                                  q = q,
                                  r = r,
                                  K = K,
                                  c = c,
                                  Tt= Tt[i.eaten],
                                  timesteplength=Tt[i.eaten]/steps)
    Neaten.est[i.eaten] = N0[i.eaten] - temp[steps+1,2] # how many at the beginning - how many at the end = how many eaten
    predrep.est[i.eaten] = temp[steps+1,3] - P[i.eaten] # how many at the end - how many at the beginning = how many new
    
  }
  ret = list(Neaten.est, predrep.est, temp)
  return(ret)}


nll.odeint.general.pred = function(Ndead, N0, b.log, h, q, r, K.log, c, Tt, P, P.end, steps=100, sigma){
  if(sigma <= 0) return(Inf)
  temp2 = eaten.odeint.general.pred(N0=N0,
                                    b=exp(b.log),
                                    h=h,
                                    q=q,
                                    r=r,
                                    K=exp(K.log),
                                    c=c,
                                    Tt=Tt,
                                    P=P,
                                    steps=steps)
  
  y = temp2[[1]]
  z = temp2[[2]]
  
  nll = -1*sum(dnorm(x = log(N0-Ndead),
                     mean = log(N0-y),
                     sd = sigma,
                     log=T))

  nll2 = -1*sum(dpois(x = P.end, 
                      lambda = P+z, 
                      log=T))
  nll = nll + nll2
  return(nll)
}
