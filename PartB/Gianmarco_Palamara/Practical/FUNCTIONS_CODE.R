#########################################################################
###functions used in Test_Daphnia.R
### Workshop Bayesian Thinking in Ecology   UZH 06-10-2021

### Code written by Gian Marco Palamara
#########################################################################

#########################################################################
# fecundity and survival function 
#########################################################################

# fecundity function
# takes vector of states (at a given time step), parameters, number of adult stages
# returns sum of the number of eggs per per adult stage (at the next time step)

Fe  <- function(FECUNDITY,ZERO_INFLATION,y,Dt,parms,Ne,Na,DET){ ## note rate*Dt
  
  sum_a <- 0
  eggs  <- 0
  
  for(i in  (length(y)-Na):(length(y))) eggs  <- eggs + parms["f"]*y[i]*Dt/Ne
  for(i in  (length(y)-Na):(length(y))) sum_a <- sum_a + y[i]
  #if(eggs==(sum_a*parms["f"])) print("eggs_ok")
  if(FECUNDITY=="P"){
    if(ZERO_INFLATION=="F")output <- ifelse(DET==TRUE,eggs,rpois(1,eggs))
    if(ZERO_INFLATION=="T"){
      z <- rbinom(1,sum_a,1-parms["pi"])
      output <- ifelse(DET==TRUE,eggs,rpois(1,parms["f"]*z*Dt/Ne))
    }
  }
  if(FECUNDITY=="N"){
    if(ZERO_INFLATION=="F"){
      p      <- parms["r"]/(parms["r"]+parms["f"])
      rr     <- (sum_a*parms["r"])/Ne
      lambda <- rgamma(1,rr,p/(1-p))
      output <- ifelse(DET==TRUE,eggs,rpois(1,lambda))
    }
    if(ZERO_INFLATION=="T"){
      z      <- rbinom(1,sum_a,1-parms["pi"])
      p      <- parms["r"]/(parms["r"]+parms["f"])
      rr     <- (z*parms["r"])/Ne
      lambda <- rgamma(1,rr,p/(1-p))
      output <- ifelse(DET==TRUE,eggs,rpois(1,lambda))
      
    }
  }
  return(output)
}


# survival function
# takes vector of states (at a given time step), parameters, number of adult stages
# returns vector of states without eggs (at the next time step)

Su  <- function(TIME_DEP,CROWDING,y,Dt,parms,time,Ne,Na,DET){
  
  y_new <- numeric(length(y))
  
  for(i in 1:(Ne-1))y_new[i+1] <- y[i]
  
  for(i in Ne:(length(y)-1)) {
    x0 <- 0.001
    mort <- parms["H"]*Dt+ifelse(TIME_DEP=="T",parms["HH"]*Dt*time,0)+ifelse(CROWDING=="T",parms["HK"],1/1000000000)*sum(y[-1])*Dt^2
    mort <- ifelse(mort>=x0,mort,x0*exp(mort/x0-1))
    ss <- exp(-mort)
    #ss <- exp(-max(parms["H"]*Dt+ifelse(TIME_DEP=="T",parms["HH"]*Dt*time,0)+ifelse(CROWDING=="T",parms["HK"],1/10000)*sum(y[-1])*Dt^2,0))
    y_new[i+1] <- ifelse(DET==TRUE,y[i]*ss,rbinom(1,y[i],ss))
  }
  return(y_new[-1])
}

#########################################################################
### define model
#########################################################################

model <- function(FECUNDITY,ZERO_INFLATION,TIME_DEP,CROWDING,y,Dt,parms,time,Ne,Na,DET){
  return(c(Fe(FECUNDITY,ZERO_INFLATION,y,Dt,parms,Ne,Na,DET),
           Su(TIME_DEP,CROWDING,y,Dt,parms,time,Ne,Na,DET)))
}

#########################################################################
## function to extract data
#########################################################################

prepare.LTE.data.list <- function(filename,separator=","){
  data.LTE <- read.csv(filename,sep=separator,stringsAsFactors = FALSE)
  data.LTE$date.n <- as.Date(data.LTE$date,format="%d.%m.%y")
  
  data.LTE$doy <- as.numeric(strftime(data.LTE$date.n, format = "%j"))
  
  # data.LTE$TotalP
  
  rind <- grep("<",data.LTE$TotalP)
  
  data.LTE[,"TotalP_num"] <- data.LTE$TotalP
  for(i in rind)
  {
    temp <- unlist(strsplit(data.LTE$TotalP[i],split="<"))[2]
    
    data.LTE$TotalP_num[i] <- as.numeric(temp)/2 # set values below detection limit to 50% of the detection limit
  }
  
  # sort data
  trts <- unique(data.LTE$trt)#[-5]
  clones <- unique(data.LTE$clone)#[-5]
  
  ldata.LTE <- list()
  
  for (t in trts)
  {
    ind.t <- which(data.LTE$trt==t)
    data.t <- data.LTE[ind.t,]
    
    ldata.LTE[[t]] <- list()
    
    for(cl in clones)
    {
      
      Dcl <- ifelse(cl=="g100","D1",
                    ifelse(cl=="GAR24","D2",
                           ifelse(cl=="VARE249","D3",
                                  ifelse(cl=="mixed","Ds",
                                         ifelse(cl=="algae_only","none",break)))))
      
      ldata.LTE[[t]][[Dcl]] <- list()
      
      ind.cl <- which(data.t$clone == cl) 
      data.t.cl <- data.t[ind.cl,]
      
      jars <- unique(data.t.cl$jar)
      
      for(j in jars)
      {
        ind.j <- which(data.t.cl$jar==j)
        dataj <- data.t.cl[ind.j,]
        tind <- order(dataj$doy) # only works if experiment is run within one calendar year
        dataj <- dataj[tind,]
        
        ### babies that died
        rind <- which(dataj$reset_Y_N==1)
        if(length(rind)>0) dataj <- dataj[max(rind):nrow(dataj),]
        
        t0 <- min(dataj$doy)
        times <- dataj$doy - t0  
        
        if(cl=="g100")
        {
          C_D1e <- as.numeric(dataj$egg_sum)
          C_D1a <- as.numeric(dataj$adult_non)+as.numeric(dataj$adult_egg)  
          C_D1j <- as.numeric(dataj$babies)+as.numeric(dataj$juv)  
          C_D2e <- rep(0,length(times))
          C_D2a <- rep(0,length(times))
          C_D2j <- rep(0,length(times))
          C_D3e <- rep(0,length(times))
          C_D3a <- rep(0,length(times))
          C_D3j <- rep(0,length(times))
          C_Dse <- C_D1e+C_D2e+C_D3e  
          C_Dsa <- C_D1a+C_D2a+C_D3a  
          C_Dsj <- C_D1j+C_D2j+C_D3j
          C_P   <- as.numeric(dataj$TotalP_num)
          C_ALG <- as.numeric(dataj$mgC)
        }
        
        if(cl=="GAR24")
        {
          C_D1e <- rep(0,length(times))  
          C_D1a <- rep(0,length(times))
          C_D1j <- rep(0,length(times))
          C_D2e <- as.numeric(dataj$egg_sum)
          C_D2a <- as.numeric(dataj$adult_non)+as.numeric(dataj$adult_egg)  
          C_D2j <- as.numeric(dataj$babies)+as.numeric(dataj$juv)  
          C_D3e <- rep(0,length(times))  
          C_D3a <- rep(0,length(times))
          C_D3j <- rep(0,length(times))
          C_Dse <- C_D1e+C_D2e+C_D3e  
          C_Dsa <- C_D1a+C_D2a+C_D3a  
          C_Dsj <- C_D1j+C_D2j+C_D3j
          C_P   <- as.numeric(dataj$TotalP_num)
          C_ALG <- as.numeric(dataj$mgC)
        }
        
        if(cl=="VARE249")
        {
          C_D1e <- rep(0,length(times))  
          C_D1a <- rep(0,length(times))
          C_D1j <- rep(0,length(times))
          C_D2e <- rep(0,length(times))  
          C_D2a <- rep(0,length(times))
          C_D2j <- rep(0,length(times))
          C_D3e <- as.numeric(dataj$egg_sum)
          C_D3a <- as.numeric(dataj$adult_non)+as.numeric(dataj$adult_egg)  
          C_D3j <- as.numeric(dataj$babies)+as.numeric(dataj$juv)  
          C_Dse <- C_D1e+C_D2e+C_D3e  
          C_Dsa <- C_D1a+C_D2a+C_D3a  
          C_Dsj <- C_D1j+C_D2j+C_D3j
          C_P   <- as.numeric(dataj$TotalP_num)
          C_ALG <- as.numeric(dataj$mgC)
        }
        
        if(cl=="mixed")
        {
          C_D1e <- rep(NA,length(times))
          C_D1a <- rep(NA,length(times))
          C_D1j <- rep(NA,length(times))
          C_D2e <- rep(NA,length(times))
          C_D2a <- rep(NA,length(times))
          C_D2j <- rep(NA,length(times))
          C_D3e <- rep(NA,length(times))
          C_D3a <- rep(NA,length(times))
          C_D3j <- rep(NA,length(times))
          C_Dse <- as.numeric(dataj$egg_sum)  
          C_Dsa <- as.numeric(dataj$adult_non)+as.numeric(dataj$adult_egg)  
          C_Dsj <- as.numeric(dataj$babies)+as.numeric(dataj$juv)
          C_P   <- as.numeric(dataj$TotalP_num)
          C_ALG <- as.numeric(dataj$mgC)
        }
        
        if(cl=="algae_only")
        {
          C_D1e <- rep(0,length(times))
          C_D1a <- rep(0,length(times))
          C_D1j <- rep(0,length(times))
          C_D2e <- rep(0,length(times))
          C_D2a <- rep(0,length(times))
          C_D2j <- rep(0,length(times))
          C_D3e <- rep(0,length(times))
          C_D3a <- rep(0,length(times))
          C_D3j <- rep(0,length(times))
          C_Dse <- rep(0,length(times))
          C_Dsa <- rep(0,length(times))
          C_Dsj <- rep(0,length(times))
          C_P   <- as.numeric(dataj$TotalP_num)
          C_ALG <- as.numeric(dataj$mgC)
        }
        
        data.jl <- cbind(t=times,C_D1e,C_D1a,C_D1j,C_D2e,C_D2a,C_D2j,C_D3e,C_D3a,C_D3j,C_Dse,C_Dsa,C_Dsj,C_P,C_ALG)
        ldata.LTE[[t]][[Dcl]][[paste0("jar_",j)]] <- data.jl
      }
    }
  }
  return(ldata.LTE)
}

#########################################################################
### SET PRIORS 
#########################################################################
### NON HIERARCHICAL MODELS

INSTAR_TIME <- 3

mf <- ifelse(ZERO_INFLATION=="F",INSTAR_TIME,INSTAR_TIME/2);
mr <- 1; mpi <- 0.5
mH <- 0.1; mHH <- 0;  mHK <- 0.0001

### proportion of the means to get standar deviation
PF <- 1; PR <- 0.25; PPI <- 0.3
PH <- 1; PHH <- 1; PHK <- 1

sdf <- PF*mf; sdr <- PR*mr; sdpi <- PPI*mpi
sdH <- PH*mH; sdHH <- ifelse(mHH==0,0.01,PHH*mHH); sdHK <- PHK*mHK

#########################################################################
## conversion functions that give the mean and sd of associated normal (to put into R)
mfun <- function(m,s){log(m)-log(1+s^2/m^2)/2}
sdfun <- function(m,s){sqrt(log(1+(s/m)^2))}

## mean sd and precision tau (to put into Jags) of associated normals to put into lognormal priors
mlf  <- mfun(mf,sdf)
sdlf <- sdfun(mf,sdf)
tauf <- (1/sdlf)^2

mlH  <- mfun(mH,sdH)
sdlH <- sdfun(mH,sdH)
tauH <- (1/sdlH)^2

mlHH  <- mfun(mHH,sdHH)
sdlHH <- sdfun(mHH,sdHH)
tauHH <- (1/sdlHH)^2

mlHK  <- mfun(mHK,sdHK)
sdlHK <- sdfun(mHK,sdHK)
tauHK <- (1/sdlHK)^2

mlr  <- mfun(mr,sdr)
sdlr <- sdfun(mr,sdr)
taur <- (1/sdlr)^2

mlpi  <- mfun(mpi,sdpi)
sdlpi <- sdfun(mpi,sdpi)
taupi <- (1/sdlpi)^2


#########################################################################
### HIERARCHICAL MODELS
### define priors for mean and variance of the hyperparameters distributions

### Fraction ratio between mean and variance in hyperparameters
FRAC <- 2

### means of means and standard deviations 
m_mf <- mf; m_mr <- mr; m_mpi <- mpi
m_mH <- mH; m_mHH <- mHH; m_mHK <- mHK

m_sdf <- sdf; m_sdr <- sdr; m_sdpi <- sdpi
m_sdH <- sdH; m_sdHH <- sdHH; m_sdHK <- sdHK

### standard deviations of means and standard deviations
sd_mf <- sdf; sd_mr <- sdr; sd_mpi <- sdpi
sd_mH <- sdH; sd_mHH <- sdHH; sd_mHK <- sdHK

#sd_mf <- m_mf/FRAC; sd_mr <- m_mr/FRAC; sd_mpi <- m_mpi/FRAC
#sd_mH <- m_mH/FRAC; sd_mHH <- m_mH/FRAC; sd_mHK <- m_mHK/FRAC

sd_sdf <- m_sdf/FRAC; sd_sdr <- m_sdr/FRAC; sd_sdpi <- m_sdpi/FRAC
sd_sdH <- m_sdH/FRAC; sd_sdHH <- m_sdHH/FRAC; sd_sdHK <- m_sdHK/FRAC

## mean sd and precision tau (to put into Jags) of associated normals to put into lognormal priors
ml_mf   <- mfun(m_mf,sd_mf)
sdl_mf  <- sdfun(m_mf,sd_mf)
tau_mf  <- (1/sdl_mf)^2

ml_mr   <- mfun(m_mr,sd_mr)
sdl_mr  <- sdfun(m_mr,sd_mr)
tau_mr  <- (1/sdl_mr)^2

ml_mpi   <- mfun(m_mpi,sd_mpi)
sdl_mpi  <- sdfun(m_mpi,sd_mpi)
tau_mpi  <- (1/sdl_mpi)^2

ml_mH   <- mfun(m_mH,sd_mH)
sdl_mH  <- sdfun(m_mH,sd_mH)
tau_mH  <- (1/sdl_mH)^2

ml_mHH   <- mfun(m_mHH,sd_mHH)
sdl_mHH  <- sdfun(m_mHH,sd_mHH)
tau_mHH  <- (1/sdl_mHH)^2

ml_mHK   <- mfun(m_mHK,sd_mHK)
sdl_mHK  <- sdfun(m_mHK,sd_mHK)
tau_mHK  <- (1/sdl_mHK)^2

ml_sdf   <- mfun(m_sdf,sd_sdf)
sdl_sdf  <- sdfun(m_sdf,sd_sdf)
tau_sdf  <- (1/sdl_sdf)^2

ml_sdr   <- mfun(m_sdr,sd_sdr)
sdl_sdr  <- sdfun(m_sdr,sd_sdr)
tau_sdr  <- (1/sdl_sdr)^2

ml_sdpi   <- mfun(m_sdpi,sd_sdpi)
sdl_sdpi  <- sdfun(m_sdpi,sd_sdpi)
tau_sdpi  <- (1/sdl_sdpi)^2

ml_sdH   <- mfun(m_sdH,sd_sdH)
sdl_sdH  <- sdfun(m_sdH,sd_sdH)
tau_sdH  <- (1/sdl_sdH)^2

ml_sdHH   <- mfun(m_sdHH,sd_sdHH)
sdl_sdHH  <- sdfun(m_sdHH,sd_sdHH)
tau_sdHH  <- (1/sdl_sdHH)^2

ml_sdHK   <- mfun(m_sdHK,sd_sdHK)
sdl_sdHK  <- sdfun(m_sdHK,sd_sdHK)
tau_sdHK  <- (1/sdl_sdHK)^2

#########################################################################
### LogPosterior 
#########################################################################
# 
## debug
# FILE_ACRONYM <- "Opt_Chains_"
# MODEL_NAME <- "PH"
# treat <- "dz"
# clone <- "D2"
# Priors<- 1
# Ne <- 3
# Nj <- 5
# Na <- 45
# Dt <- 1
# lik <- F
# replication<-F
# DEBUG<-F

LogPosterior <- function(FILE_ACRONYM,MODEL_NAME,treat,clone,Priors,Ne,Nj,Na,Dt,lik,replication,DEBUG){
  
  ##########################################################################################
  ## get info about model structure and parameters
  ##########################################################################################
  
  Npars <- 2
  
  INF_METH  <- substr(MODEL_NAME, nchar(MODEL_NAME), nchar(MODEL_NAME))
  FECUNDITY <- substr(MODEL_NAME, 1, 1)
  if(FECUNDITY=="N") Npars <- Npars + 1
  
  ZERO_INFLATION <- F
  if(substr(MODEL_NAME, 2, 2)=="Z")ZERO_INFLATION<-T
  
  if(ZERO_INFLATION==F){
    if(nchar(MODEL_NAME)==4) {
      TIME_DEP <- "T"
      CROWDING <- "T"
      Npars <- Npars+2
    }
    if(nchar(MODEL_NAME)==3) {
      if(substr(MODEL_NAME, 2, 2)=="T"){
        TIME_DEP <- "T"
        CROWDING <- "F"
      }
      if(substr(MODEL_NAME, 2, 2)=="C"){
        TIME_DEP <- "F"
        CROWDING <- "T"
      }
      Npars <- Npars+1
    }
    if(nchar(MODEL_NAME)==2){
      TIME_DEP <- "F"
      CROWDING <- "F"
    }
  }
  if(ZERO_INFLATION==T){
    Npars <- Npars + 1
    if(nchar(MODEL_NAME)==5) {
      TIME_DEP <- "T"
      CROWDING <- "T"
      Npars <- Npars+2
    }
    if(nchar(MODEL_NAME)==4) {
      if(substr(MODEL_NAME, 3, 3)=="T"){
        TIME_DEP <- "T"
        CROWDING <- "F"
      }
      if(substr(MODEL_NAME, 3, 3)=="C"){
        TIME_DEP <- "F"
        CROWDING <- "T"
      }
      Npars <- Npars+1
    }
    if(nchar(MODEL_NAME)==3){
      TIME_DEP <- "F"
      CROWDING <- "F"
    }
  }
  
  ##########################################################################################
  ### get real data
  ##########################################################################################
  tt<-treat
  cc<-clone
  
  dd.list <- prepare.LTE.data.list(filename="Daphnia_LTE_master_20171013.csv",separator=";")
  dd.to.fit.list <- list()
  
  ## get all the replicates per treatment per clone
  TREATMENT <- treat
  CLONE     <- clone
  dd.eval   <- dd.list[[TREATMENT]][[CLONE]]
  
  ## count the number of replicates
  ## and prepare arrays 
  R <- length(dd.eval)
  n_obs <- numeric(R)
  t_end <- numeric(R)
  
  t_obs <- numeric()
  y_e_obs <- numeric()
  y_j_obs <- numeric()
  y_a_obs <- numeric()
  
  ## replicates loop
  for(k in 1:R) {
    if(CLONE=="D1")pp <- c(2,3,4)
    if(CLONE=="D2")pp <- c(5,6,7)
    if(CLONE=="D3")pp <- c(8,9,10)
    if(CLONE=="Ds")pp <- c(11,12,13)
    
    ## get specific replicate as a dataframe
    dd <- data.frame(dd.eval[[k]])[,c(1,pp)]
    
    ## remove NAs (measurments where daphnia has not been counted)
    nn <- which(is.na(rowSums(dd[,2:4])==T))
    if(length(nn)>0)dd<-dd[-nn,]
    
    ## count how many obsrevations per replciate
    n_obs[k] <- length(dd[,"t"])
    
    ## get observation times and individual counts
    t_obs   <- c(t_obs,dd[,"t"])
    y_e_obs <- c(y_e_obs,dd[,2])
    y_j_obs <- c(y_j_obs,dd[,4])
    y_a_obs <- c(y_a_obs,dd[,3])
    
    ## get time of last measurment
    t_end[k] <- t_obs[sum(n_obs)] 
    
    y_obs <- cbind(t_obs,y_e_obs,y_j_obs,y_a_obs)
    
  }###replicates cycle
  
  ## classifiy the measurments
  class_obs <- rep(NA,length(t_obs))
  ll <- which(rowSums(y_obs[,2:4])==0)
  
  ## EXT if there is a final extinction LAZ if there is lazarus effect
  class_obs[ll] <- ifelse((t_obs[ll] %in% t_end)==TRUE,"EXT","LAZ")
  
  ## put everything into the list
  y_obs <- data.frame(cbind(y_obs,class_obs))
  dd.to.fit.list[[tt]][[cc]] <- list(n_obs,t_end,y_obs)
  
  tt <- TREATMENT
  cc <- CLONE
  
  ## extract from the list
  n_obs_to_fit <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][1])))
  x_end        <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][2])))
  R            <- length(n_obs_to_fit)
  
  x_obs_to_fit   <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,1])))
  y_e_obs_to_fit <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,2])))
  y_j_obs_to_fit <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,3])))
  y_a_obs_to_fit <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,4])))
  
  ## remove initial conditions
  n_obs_to_fit <- n_obs_to_fit-1
  rr           <- which(x_obs_to_fit==0)
  
  x_obs_to_fit   <- x_obs_to_fit[-rr]
  y_e_obs_to_fit <- y_e_obs_to_fit[-rr]
  y_j_obs_to_fit <- y_j_obs_to_fit[-rr]
  y_a_obs_to_fit <- y_a_obs_to_fit[-rr]
  
  ## get bounds of nobs and t_obs per each replicate
  nn_s<- numeric(R)
  nn_e<- numeric(R)
  xx_s<-numeric()
  xx_e<-numeric()
  for(k in 1:R){
    nn_s[k] <- ifelse(k==1,1,sum(n_obs_to_fit[1:(k-1)])+1)
    xx_s[k] <- ifelse(k==1,1,sum(x_end[1:(k-1)])+1)
    nn_e[k] <- sum(n_obs_to_fit[1:k])
    xx_e[k] <- sum(x_end[1:k])
  }
  
  
  ##########################################################################################
  #### compute log prior, log likelihood and log posterior
  ##########################################################################################
  
  LogPri    <- 0
  LogLik    <- 0
  LogLikInt <- 0
  
  ## get states and parameters and put them into the function
  ##########################################################################################
  
  if(INF_METH=="J"){
    
    ##########################################################################################
    ### prepare structure, extract and check chains
    ##########################################################################################
    
    filename <- paste(FILE_ACRONYM,MODEL_NAME,"_",treat,"_",clone,"_",Priors,".csv",sep="")
    Chains  <- read.csv(filename)
    samples <- dim(Chains)[1]
    
    LogPri    <- numeric(samples)
    LogLik    <- numeric(samples)
    LogLikInt <- numeric(samples)
    LogLikMean <- 0
    
    y  <- array(0,dim=c(Na+Nj+Ne,sum(t_end),samples))
    y_e<- array(0,dim=c(Ne,sum(x_end),samples))
    y_j<- array(0,dim=c(Nj,sum(x_end),samples))
    y_a<- array(0,dim=c(Na,sum(x_end),samples))
    
    if(FECUNDITY=="N")     lambda <- array(0,dim=c(sum(t_end),samples))
    if(ZERO_INFLATION)     z      <- array(0,dim=c(sum(t_end),samples))
    
    f  <- Chains[,"f"]
    H  <- Chains[,"H"]
    if(ZERO_INFLATION)    pi <- Chains[,"pi"]
    if(FECUNDITY=="N")    r  <- Chains[,"r"]
    if(TIME_DEP=="T")     HH <- Chains[,"HH"] 
    if(CROWDING=="T")     HK <- Chains[,"HK"] 
    
    for(ii in 1:(Ne+Nj+Na))for(jj in 1:sum(t_end))    y[ii,jj,]  <-Chains[,paste("y.",ii,".",jj,".",sep="")]
    
    if(FECUNDITY=="N") {
      for(k in 1:R ){
        for( i in xx_s[k]:(xx_e[k]) ) {
          lambda[i,]<-Chains[,paste("lambda.",i,".",sep="")]
        }}
    }
    
    if(ZERO_INFLATION){
      for(k in 1:R){
        for( i in xx_s[k]:(xx_e[k]) ) {
          z[i,] <- Chains[,paste("z.",i,".",sep="")]
        }}
    }
    Corr_z <- ifelse(ZERO_INFLATION,sum(t_end),0)
    Corr <- ifelse(FECUNDITY=="N",(sum(t_end)),0)
    if(Npars+2+sum(t_end)*(Nj+Na+Ne)+Corr+Corr_z==dim(Chains)[2])print("Chains extracted successfully. Calculating Posteriors Density")
    
    ##########################################################################################
    ## assign states
    ##########################################################################################
    
    for ( k in 1:R ){ 
      for (i in xx_s[k]:xx_e[k]){
        for(j in 1:Ne )y_e[j,i,] <- y[j,i,]
        for(j in 1:Nj )y_j[j,i,] <- y[Ne+j,i,]
        for(j in 1:Na )y_a[j,i,] <- y[Ne+Nj+j,i,]
      }
    }
    ##########################################################################################
    
    if(lik==F){
      ##########################################################################################
      ### logprior
      ##########################################################################################
      
      if(Priors==1){
        LogPri <- LogPri + log(dlnormTrunc(f,mlf,sdlf,0,30)) + dlnorm(H,mlH,sdlH,log=TRUE)
        if(ZERO_INFLATION)LogPri <- LogPri + log(dnormTrunc(pi,mpi,sdpi,0,1))
        if(FECUNDITY=="N")LogPri <- LogPri + dlnorm(r,mlr,sdlr,log=TRUE)  
        if(TIME_DEP=="T") LogPri <- LogPri + dnorm(HH,mHH,sdHH,log=TRUE)
        if(CROWDING=="T") LogPri <- LogPri + dlnorm(HK,mlHK,sdlHK,log=TRUE)
      }
      if(Priors==2){
        LogPri <- LogPri + dunif(f,0,30,log=TRUE) + dunif(H,0,1,log=TRUE)
        if(ZERO_INFLATION)LogPri <- LogPri + dunif(pi,0,1,log=TRUE)
        if(FECUNDITY=="N")LogPri <- LogPri + dunif(r,0,100,log=TRUE)
        if(TIME_DEP=="T") LogPri <- LogPri + dunif(HH,0,1,log=TRUE)
        if(CROWDING=="T") LogPri <- LogPri + dunif(HK,0,10,log=TRUE)
      }
      
      ##########################################################################################
      ### loglikelihood (intial conditions)
      ##########################################################################################
      
      ## (should be all 0)
      # for ( k in 1:R ) {
      #   for ( j in 1:Ne )LogLikInt <- LogLikInt + dpois(y[j,xx_s[k],],0,log=TRUE)
      #   LogLikInt <- LogLikInt + dbinom(y[Ne+1,xx_s[k],],3,1,log=TRUE)
      #   for ( j in (Ne+2):(Nj+Na+Ne) )LogLikInt <- LogLikInt + dpois(y[j,xx_s[k],],0,log=TRUE)
      # }
      
      ##########################################################################################
      ### loglikelihood (Model component)
      ##########################################################################################
      
      for(k in 1:R ) {
        for ( i in xx_s[k]:(xx_e[k]-1) ) {
          
          
          if(FECUNDITY=="P"){
            if(ZERO_INFLATION){
              incr_z      <- dbinom(z[i+1,],colSums(y_a[,i,]),1-pi,log=T)
              incr_pois   <- dpois(y[1,i+1,],f*z[i+1,]*Dt/Ne,log=T)
              incr_f      <- incr_pois + incr_z
            }
            if(ZERO_INFLATION==F)incr_f   <- dpois(y[1,i+1,],f*colSums(y_a[,i,])*Dt/Ne,log=T)
          }
          if(FECUNDITY=="N") {
            if(ZERO_INFLATION){
              incr_z        <- dbinom(z[i+1,],colSums(y_a[,i,]),1-pi,log=T)
              incr_lambda   <- dgamma(lambda[i+1,],pmax(r*(ifelse(z[i+1,]<0.5,1,z[i+1,]))/Ne,0.001),pmax(r/(f*Dt),0.001),log=T)
              incr_pois     <- dpois(y[1,i+1,],ifelse(z[i+1,]==0,0,lambda[i+1,]),log = T) 
              incr_f        <- incr_pois + incr_lambda + incr_z
              
              #incr_nb       <- dnbinom(y[1,i+1,],(r*colSums(y_a[,i,]))*Dt/Ne,r/(r+f*Dt),log=T)
            }
            if(ZERO_INFLATION==F){
              incr_lambda   <- dgamma(lambda[i+1,],pmax(r*(ifelse(colSums(y_a[,i,])<0.5,1,colSums(y_a[,i,])))/Ne,0.001),pmax(r/(f*Dt),0.001),log=T)
              incr_pois     <- dpois(y[1,i+1,],ifelse(colSums(y_a[,i,])==0,0,lambda[i+1,]),log = T) 
              incr_f        <- incr_pois + incr_lambda
              
              #incr_nb       <- dnbinom(y[1,i+1,],(r*colSums(y_a[,i,]))*Dt/Ne,r/(r+f*Dt),log=T)
            }
          }
          
          ## debug
          if(DEBUG==T){
            ind.samp<-10
            print(paste("***********************************************",sep=""))
            print(paste("DEBUG replica            ",k,sep=""))
            print(paste("DEBUG time step          ",i,sep=""))
            #print(paste("MEAN (m)          ",m[ind.samp],sep=""))
            if(FECUNDITY=="P"){
              if(ZERO_INFLATION){
                
                print(paste("fecundity component z",signif(incr_z[ind.samp],digits=5),sep=""))
                print(paste("fecundity component pois",signif(incr_pois[ind.samp],digits=5),sep=""))
                
                if(is.finite(incr_f[ind.samp])==F){
                  print(paste("y[1,",i+1,"] = ",y[1,i+1,ind.samp],sep=""))
                  print(paste("z[",i+1,"] = ",z[i+1,ind.samp],sep=""))
                }
              }
              if(ZERO_INFLATION==F){
                print(paste("fecundity component pois",signif(incr_f[ind.samp],digits=5),sep=""))
                if(is.finite(incr_f[ind.samp])==F){
                  print(paste("y[1,",i+1,"] = ",y[1,i+1,ind.samp],sep=""))
                  print(paste("z[",i+1,"] = ",z[i+1,ind.samp],sep=""))
                }
              }
            }
            if(FECUNDITY=="N"){
              print(paste("fecundity component  lambda ",signif(incr_lambda[ind.samp],digits=5),sep=""))
              print(paste("fecundity component  pois   ",signif(incr_pois[ind.samp],digits=5),sep=""))
              print(paste("fecundity component  nebin   ",signif(incr_nb[ind.samp],digits=5),sep=""))
              if(is.finite(incr_f[ind.samp])==F){
                print(paste("y[1,",i+1,"] = ",y[1,i+1,ind.samp],sep=""))
                print(paste("r/(3*f*Dt) = ",signif(r[ind.samp],digits=3),"/",signif(3*f[ind.samp],digits=3),sep=""))
                print(paste("(r*colSums(y_a[,",i,",])+0.01)*Dt/3 = (",signif(r[ind.samp],digits=3),"*",
                            colSums(y_a[,i,])[ind.samp]," + 0.01)/3",sep=""))
                print(paste(",log=T)",sep=""))
                print(paste("is NaN",sep=""))
              }
              #break
            }
          }
          
          LogLikInt <- LogLikInt + incr_f
          j<-1
          
          #for(j in 1:(Ne-1))LogLikInt <- LogLikInt + dbinom(y[j+1,i+1,],y[j,i,],1,log=T)
          for(j in 1:(Ne-1)){
            increment <- dbinom(y[j+1,i+1,],y[j,i,],.99,log=T)
            iind <- !(is.finite(increment))
            if(any(iind)){
              for(ii in which(iind)){
                
                print(paste("time ",i,"stage ",j,"sample ",ii,sep=""))
                print(paste("state ",y[j,i,ii],"to state ",y[j+1,i+1,ii],sep=""))
              }}
            
            LogLikInt <- LogLikInt + increment
          }
          for(j in Ne:(Ne+Nj+Na-1)){
            x0    <- 0.001
            surv  <- rep(0,length(H))
            mort  <- H*Dt+ ifelse(TIME_DEP=="T",HH*(i-xx_s[k])*Dt,0) + ifelse(CROWDING=="T",HK,1/1000000000)*(colSums(y_j[,i,])+colSums(y_a[,i,]))*Dt^2
            mort  <- ifelse(mort>=x0,mort,x0*exp(mort/x0-1))
            surv  <- exp(-mort)
            #surv  <- exp(-pmax(mort,surv))
            
            incr_m  <- dbinom(y[j+1,i+1,],y[j,i,],surv,log=T)
            
            LogLikInt<- LogLikInt + incr_m
            if(DEBUG==T){
              if(incr_m[ind.samp]!=0){
                #print(paste("DEBUG replica            ",k,sep=""))
                #print(paste("DEBUG time step          ",i,sep=""))
                print(paste("DEBUG age class          ",j,sep=""))
                #print(paste("survival 1                 ",signif(surv_1[ind.samp],digits=5),sep=""))
                #print(paste("survival 2                 ",signif(surv_2[ind.samp],digits=5),sep=""))
                print(paste("mortality                 ",signif(mort[ind.samp],digits=5),sep=""))
                
                print(paste("survival                 ",signif(surv[ind.samp],digits=5),sep=""))
                print(paste("mortality component      ",signif(incr_m[ind.samp],digits=5),sep=""))
              }
            }
            
          }## close loop in age classes
          
          if(DEBUG==T){
            print(paste("***********************************************",sep=""))
            print(paste("Likelihood      ",signif(LogLikInt[ind.samp],digits=5),sep=""))
            print(paste("mean Likelihood ",signif(mean(LogLikInt),digits=5),sep=""))
            
            nns <- length(which(is.finite(LogLikInt)==F))
            if(nns>0){print(paste(nns," NaNS appear in the likelihood at replica ",k," time ",i,sep=""))
              print(which(is.finite(LogLikInt)==F))
            }
          }
          
        }# close loop in time steps
      }# close loop in replicates
      
      ##########################################################################################
      
    }### compute only if lik is F
    
    ##########################################################################################
    ### Observation part
    ##########################################################################################
    
    for ( k in 1:R ) {
      i <- nn_s[k]
      for ( i in nn_s[k]:nn_e[k] ) {
        
        if(DEBUG==T){
          # cat("\nk=",k," i=",i," t=",x_obs_to_fit[i],
          #     "\n  e_obs=",y_e_obs_to_fit[i]," y_e=",y_e[,xx_s[k]+x_obs_to_fit[i]-1,ind.samp]," loglik=",dnorm(y_e_obs[i],y_e[,xx_s[k]+x_obs_to_fit[i]-1,],1,log=TRUE)[ind.samp],
          #     "\n  j_obs=",y_j_obs_to_fit[i]," y_j=",colSums(y_j[,xx_s[k]+x_obs_to_fit[i]-1,])[ind.samp]," loglik=",dnorm(y_j_obs[i],colSums(y_j[,xx_s[k]+x_obs_to_fit[i]-1,]),1,log=TRUE)[ind.samp],
          #     "\n  a_obs=",y_a_obs_to_fit[i]," y_a=",colSums(y_a[,xx_s[k]+x_obs_to_fit[i]-1,])[ind.samp]," loglik=",dnorm(y_a_obs[i],colSums(y_a[,xx_s[k]+x_obs_to_fit[i]-1,]),1,log=TRUE)[ind.samp],
          #     sep="")
        }
        
        LogLik     <-  LogLik     + dnorm(y_e_obs_to_fit[i],colSums(y_e[,xx_s[k]+x_obs_to_fit[i]-1,]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_e_obs_to_fit[i],sum(apply(y_e[,xx_s[k]+x_obs_to_fit[i]-1,],1,mean)),1,log=TRUE)
        
        LogLik     <-  LogLik     + dnorm(y_j_obs_to_fit[i],colSums(y_j[,xx_s[k]+x_obs_to_fit[i]-1,]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_j_obs_to_fit[i],sum(apply(y_j[,xx_s[k]+x_obs_to_fit[i]-1,],1,mean)),1,log=TRUE)
        
        LogLik     <-  LogLik     + dnorm(y_a_obs_to_fit[i],colSums(y_a[,xx_s[k]+x_obs_to_fit[i]-1,]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_a_obs_to_fit[i],sum(apply(y_a[,xx_s[k]+x_obs_to_fit[i]-1,],1,mean)),1,log=TRUE)
      }
    }
    
    
    ##########################################################################################
    
  }### close INF_METH "J"
  
  if(INF_METH=="R"){
    
    ##########################################################################################
    ### initialization and file extraction
    ##########################################################################################
    
    filename <- paste(FILE_ACRONYM,MODEL_NAME,"_",treat,"_",clone,"_REP_",1,"_",Priors,".csv",sep="")
    Chains  <- read.csv(filename)
    samples <- dim(Chains)[1]
    
    ##########################################################################################
    ### replication loop: compute log posterior for all replicates
    ##########################################################################################
    
    if(is.numeric(replication)==F){
      
      ##########################################################################################
      ### containers
      ##########################################################################################
      
      LogPri              <- array(0,dim = c(samples,R))
      colnames(LogPri)    <- paste("LogPri_",seq(1,R,1),sep="")
      LogLik    <- array(0,dim = c(samples,R))
      colnames(LogLik)    <- paste("LogLik_",seq(1,R,1),sep="")
      LogLikInt <- array(0,dim = c(samples,R))
      colnames(LogLikInt)    <- paste("LogLikInt_",seq(1,R,1),sep="")
      LogLikMean <- numeric(R)
      
      ##########################################################################################
      
      k<-1
      for(k in 1:R){
        ##########################################################################################
        ### prepare structure, extract and check chains
        ##########################################################################################
        
        ### select replica specific observation
        x_obs_to_fit_rep   <- x_obs_to_fit[nn_s[k]:nn_e[k]] 
        y_e_obs_to_fit_rep <- y_e_obs_to_fit[nn_s[k]:nn_e[k]] 
        y_j_obs_to_fit_rep <- y_j_obs_to_fit[nn_s[k]:nn_e[k]] 
        y_a_obs_to_fit_rep <- y_a_obs_to_fit[nn_s[k]:nn_e[k]] 
        
        if(k>1){
          filename <- paste(FILE_ACRONYM,MODEL_NAME,"_",treat,"_",clone,"_REP_",k,"_",Priors,".csv",sep="")
          Chains  <- read.csv(filename)
        }
        if(dim(Chains)[1]==samples){
          print(paste("Getting Replica ",k,sep=""))
          samples <- dim(Chains)[1]
        }
        
        y  <- array(0,dim=c(Na+Nj+Ne,t_end[k],samples))
        y_e<- array(0,dim=c(Ne,x_end[k],samples))
        y_j<- array(0,dim=c(Nj,x_end[k],samples))
        y_a<- array(0,dim=c(Na,x_end[k],samples))
        
        if(FECUNDITY=="N")     lambda <- array(0,dim=c(t_end[k],samples))
        if(ZERO_INFLATION)     z <- array(0,dim=c(t_end[k],samples))
        
        f <- Chains[,"f"]
        H <- Chains[,"H"]
        if(ZERO_INFLATION)      pi <- Chains[,"pi"]
        if(FECUNDITY=="N")      r  <- Chains[,"r"]
        if(TIME_DEP=="T")       HH <- Chains[,"HH"] 
        if(CROWDING=="T")       HK <- Chains[,"HK"] 
        
        for(ii in 1:(Ne+Na+Nj)) for(jj in 1:t_end[k]) y[ii,jj,] <-Chains[,paste("y.",ii,".",jj,".",sep="")]
        if(FECUNDITY=="N")     for( i in 1:t_end[k]) lambda[i,]<-Chains[,paste("lambda.",i,".",sep="")]
        if(ZERO_INFLATION)     for( i in 1:t_end[k]) z[i,]     <- Chains[,paste("z.",i,".",sep="")]
        
        Corr_z <- ifelse(ZERO_INFLATION,t_end[k],0)
        Corr   <- ifelse(FECUNDITY=="N",t_end[k],0)
        if(Npars+2+t_end[k]*(Nj+Na+Ne)+Corr+Corr_z==dim(Chains)[2])print(paste("Chains for replicate ",k," extracted successfully. Calculating Posteriors Density",sep=""))
        
        ##########################################################################################
        ## assign states
        ##########################################################################################
        
        for (i in 1:x_end[k]){
          for(j in 1:Ne )y_e[j,i,] <- y[j,i,]
          for(j in 1:Nj )y_j[j,i,] <- y[Ne+j,i,]
          for(j in 1:Na )y_a[j,i,] <- y[Ne+Nj+j,i,]
        }
        
        ##########################################################################################
        
        if(lik==F){
          ##########################################################################################
          ### logprior
          ##########################################################################################
          
          if(Priors==1){
            LogPri[,k] <- LogPri[,k] + log(dlnormTrunc(f,mlf,sdlf,0,30)) + dlnorm(H,mlH,sdlH,log=TRUE)
            if(ZERO_INFLATION)LogPri[,k] <- LogPri[,k] + log(dnormTrunc(pi,mpi,sdpi,0,1))
            if(FECUNDITY=="N")LogPri[,k] <- LogPri[,k] + dlnorm(r,mlr,sdlr,log=TRUE)  
            if(TIME_DEP=="T") LogPri[,k] <- LogPri[,k] + dnorm(HH,mHH,sdHH,log=TRUE)
            if(CROWDING=="T") LogPri[,k] <- LogPri[,k] + dlnorm(HK,mlHK,sdlHK,log=TRUE)
          }
          if(Priors==2){
            LogPri[,k] <- LogPri[,k] + dunif(f,0,30,log=TRUE) + dunif(H,0,1,log=TRUE)
            if(ZERO_INFLATION)LogPri[,k] <- LogPri[,k] + dunif(pi,0,1,log=TRUE)
            if(FECUNDITY=="N")LogPri[,k] <- LogPri[,k] + dunif(r,0,100,log=TRUE)
            if(TIME_DEP=="T") LogPri[,k] <- LogPri[,k] + dunif(HH,0,1,log=TRUE)
            if(CROWDING=="T") LogPri[,k] <- LogPri[,k] + dunif(HK,0,1,log=TRUE)
          } 
          
          ##########################################################################################
          ### loglikelihood (intial conditions)
          ##########################################################################################
          # 
          # ## (should be all 0)
          # for ( j in 1:(Ne) )LogLikInt[,k] <- LogLikInt[,k] + dpois(y[j,1,],0,log=TRUE)
          # 
          # LogLikInt[,k] <- LogLikInt[,k] + dbinom(y[Ne+1,1,],3,1,log=T)
          # for ( j in (Ne+2):(Ne+Nj+Na) )LogLikInt[,k] <- LogLikInt[,k] + dpois(y[j,1,],0,log=TRUE)
          # 
          ##########################################################################################
          ### loglikelihood (Model component)
          ##########################################################################################
          i<-1
          for ( i in 1:(x_end[k]-1) ) {
            
            
            if(FECUNDITY=="P"){
              if(ZERO_INFLATION){
                incr_z      <- dbinom(z[i+1,],colSums(y_a[,i,]),1-pi,log=T)
                incr_pois   <- dpois(y[1,i+1,],f*z[i+1,]*Dt/Ne,log=T)
                incr_f      <- incr_pois + incr_z
              }
              if(ZERO_INFLATION==F)incr_f   <- dpois(y[1,i+1,],f*colSums(y_a[,i,])*Dt/Ne,log=T)
            }
            if(FECUNDITY=="N") {
              if(ZERO_INFLATION){
                incr_z        <- dbinom(z[i+1,],colSums(y_a[,i,]),1-pi,log=T)
                incr_lambda   <- dgamma(lambda[i+1,],pmax(r*(ifelse(z[i+1,]<0.5,1,z[i+1,]))/Ne,0.001),pmax(r/(f*Dt),0.001),log=T)
                incr_pois     <- dpois(y[1,i+1,],ifelse(z[i+1,]==0,0,lambda[i+1,]),log = T) 
                incr_f        <- incr_pois + incr_lambda + incr_z
                
                #incr_nb       <- dnbinom(y[1,i+1,],(r*colSums(y_a[,i,]))*Dt/Ne,r/(r+f*Dt),log=T)
              }
              if(ZERO_INFLATION==F){
                incr_lambda   <- dgamma(lambda[i+1,],pmax(r*(ifelse(colSums(y_a[,i,])<0.5,1,colSums(y_a[,i,])))/Ne,0.001),pmax(r/(f*Dt),0.001),log=T)
                incr_pois     <- dpois(y[1,i+1,],ifelse(colSums(y_a[,i,])==0,0,lambda[i+1,]),log = T) 
                incr_f        <- incr_pois + incr_lambda
                
                #incr_nb       <- dnbinom(y[1,i+1,],(r*colSums(y_a[,i,]))*Dt/Ne,r/(r+f*Dt),log=T)
              }
            }
            
            ## debug
            if(DEBUG==T){
              ind.samp<-5
              print(paste("***********************************************",sep=""))
              print(paste("DEBUG replica            ",k,sep=""))
              print(paste("DEBUG time step          ",i,sep=""))
              print(paste("MEAN (m)          ",colSums(y_a[,i,ind.samp]),sep=""))
              print(paste("fecundity component  lambda ",signif(incr_lambda[ind.samp],digits=5),sep=""))
              print(paste("fecundity component  pois   ",signif(incr_pois[ind.samp],digits=5),sep=""))
              print(paste("fecundity component  nebin   ",signif(incr_nb[ind.samp],digits=5),sep=""))
              if(is.finite(incr_f[ind.samp])==F){
                print(paste("y[1,",i+1,"] = ",y[1,i+1,ind.samp],sep=""))
                print(paste("r/(3*f*Dt) = ",signif(r[ind.samp],digits=3),"/",signif(3*f[ind.samp],digits=3),sep=""))
                print(paste("(r*colSums(y_a[,",i,",])+0.01)*Dt/3 = (",signif(r[ind.samp],digits=3),"*",
                            colSums(y_a[,i,])[ind.samp]," + 0.01)/3",sep=""))
                print(paste(",log=T)",sep=""))
                print(paste("is NaN",sep=""))
                #break
              }
            }
            
            LogLikInt[,k] <- LogLikInt[,k] + incr_f
            
            
            #for(j in 1:(Ne-1))LogLikInt[,k] <- LogLikInt[,k] + dbinom(y[j+1,i+1,],y[j,i,],1,log=T)
            for(j in 1:(Ne-1)){
              increment <- dbinom(y[j+1,i+1,],y[j,i,],.99,log=T)
              iind <- !(is.finite(increment))
              if(any(iind)){
                for(ii in which(iind)){
                  
                  print(paste("time ",i,"stage ",j,"sample ",ii,sep=""))
                  print(paste("state ",y[j,i,ii],"to state ",y[j+1,i+1,ii],sep=""))
                }}
              
              LogLikInt <- LogLikInt + increment
            }
            
            for(j in Ne:(Ne+Nj+Na-1)){
              x0 <- 0.001
              surv <- rep(0,length(H))
              mort <- H*Dt + ifelse(TIME_DEP=="T",HH*i*Dt,0) + ifelse(CROWDING=="T",HK,1/1000000000)*(colSums(y_j[,i,])+colSums(y_a[,i,]))*Dt^2
              mort <- ifelse(mort>=x0,mort,x0*exp(mort/x0-1))
              surv <- exp(-mort)
              #surv <- exp(-pmax(mort,surv))
              
              incr_m <- dbinom(y[j+1,i+1,],y[j,i,],surv,log=T)
              LogLikInt[,k] <- LogLikInt[,k] + incr_m
              
              if(DEBUG==T){
                if(incr_m[ind.samp]!=0){
                  #print(paste("DEBUG replica            ",k,sep=""))
                  #print(paste("DEBUG time step          ",i,sep=""))
                  print(paste("DEBUG age class          ",j,sep=""))
                  print(paste("mortality                ",signif(mort[ind.samp],digits=5),sep=""))
                  print(paste("survival                 ",signif(surv[ind.samp],digits=5),sep=""))
                  print(paste("mortality component      ",signif(incr_m[ind.samp],digits=5),sep=""))
                }
              }
              
            }
            if(DEBUG==T){
              print(paste("***********************************************",sep=""))
              print(paste("Likelihood      ",signif(LogLikInt[ind.samp],digits=5),sep=""))
              print(paste("mean Likelihood ",signif(mean(LogLikInt),digits=5),sep=""))
              
              nns <- length(which(is.finite(LogLikInt)==F))
              if(nns>0){print(paste(nns," NaNS appear in the likelihood at replica ",k," time ",i,sep=""))
                print(which(is.finite(LogLikInt)==F))
              }
              if(length(which(is.finite(LogLikInt[,k])==F))>0)print(paste("Problem in Replica_",k,"_time_",i,sep=""))
            }
            
            if(length(which(is.finite(LogLikInt)==F))>0){
              ind.samp <- 10
              print(paste("infinite in Replica_",k,"_time_",i,sep=""))
              print(paste("***********************************************",sep=""))
              print(paste("Likelihood      ",signif(LogLikInt[ind.samp],digits=5),sep=""))
              print(paste("mean Likelihood ",signif(mean(LogLikInt),digits=5),sep=""))
              print(paste("fertility      ",signif(incr_f[ind.samp],digits=5),sep=""))
              print(paste("mortality      ",signif(incr_m[ind.samp],digits=5),sep=""))
            }  
            
            
          }
          
        }### compute only if lik is false
        
        ##########################################################################################
        ### Observation part
        ##########################################################################################
        
        i<-10
        for ( i in 1:n_obs_to_fit[k] ) {
          LogLik[,k]    <-  LogLik[,k] + dnorm(y_e_obs_to_fit_rep[i],colSums(y_e[,x_obs_to_fit_rep[i],]),1,log=TRUE)
          LogLikMean[k] <-  LogLikMean[k] + dnorm(y_e_obs_to_fit_rep[i],sum(apply(y_e[,x_obs_to_fit_rep[i],],1,mean)),1,log=TRUE)
          
          LogLik[,k]    <-  LogLik[,k] + dnorm(y_j_obs_to_fit_rep[i],colSums(y_j[,x_obs_to_fit_rep[i],]),1,log=TRUE)
          LogLikMean[k] <-  LogLikMean[k] + dnorm(y_j_obs_to_fit_rep[i],sum(apply(y_j[,x_obs_to_fit_rep[i],],1,mean)),1,log=TRUE)
          
          LogLik[,k] <-  LogLik[,k] + dnorm(y_a_obs_to_fit_rep[i],colSums(y_a[,x_obs_to_fit_rep[i],]),1,log=TRUE)
          LogLikMean[k] <-  LogLikMean[k] + dnorm(y_a_obs_to_fit_rep[i],sum(apply(y_a[,x_obs_to_fit_rep[i],],1,mean)),1,log=TRUE)
        }
        
      }### close loop in replicates
    }### close replication loop
    
    ##########################################################################################
    ### compute log posterior for one replicate
    ##########################################################################################
    
    if(is.numeric(replication)==T){
      
      ##########################################################################################
      ### containers
      ##########################################################################################
      
      LogPri <- array(0,dim = c(samples))
      LogLik <- array(0,dim = c(samples))
      LogLikInt <- array(0,dim = c(samples))
      LogLikMean <- 0
      
      k <- replication
      ##########################################################################################
      ### prepare structure, extract and check chains
      ##########################################################################################
      
      ### select replica specific observation
      x_obs_to_fit_rep   <- x_obs_to_fit[nn_s[k]:nn_e[k]] 
      y_e_obs_to_fit_rep <- y_e_obs_to_fit[nn_s[k]:nn_e[k]] 
      y_j_obs_to_fit_rep <- y_j_obs_to_fit[nn_s[k]:nn_e[k]] 
      y_a_obs_to_fit_rep <- y_a_obs_to_fit[nn_s[k]:nn_e[k]] 
      
      if(k>1){
        filename <- paste(FILE_ACRONYM,MODEL_NAME,"_",treat,"_",clone,"_REP_",k,"_",Priors,".csv",sep="")
        Chains  <- read.csv(filename)
      }
      if(dim(Chains)[1]==samples){
        print(paste("Getting Replica ",k,sep=""))
        samples <- dim(Chains)[1]
      }
      
      y  <- array(0,dim=c(Na+Nj+Ne,t_end[k],samples))
      y_e<- array(0,dim=c(Ne,x_end[k],samples))
      y_j<- array(0,dim=c(Nj,x_end[k],samples))
      y_a<- array(0,dim=c(Na,x_end[k],samples))
      
      if(FECUNDITY=="N") lambda <- array(0,dim=c(t_end[k],samples))
      if(ZERO_INFLATION) z <- array(0,dim=c(t_end[k],samples))
      
      f <- Chains[,"f"]
      H <- Chains[,"H"]
      if(ZERO_INFLATION)      pi <- Chains[,"pi"]
      if(FECUNDITY=="N")      r  <- Chains[,"r"]
      if(TIME_DEP=="T")       HH <- Chains[,"HH"] 
      if(CROWDING=="T")       HK <- Chains[,"HK"] 
      
      for(ii in 1:(Ne+Na+Nj))for(jj in 1:t_end[k]) y[ii,jj,]<-Chains[,paste("y.",ii,".",jj,".",sep="")]
      if(FECUNDITY=="N")     for( i in 1:t_end[k]) lambda[i,]<-Chains[,paste("lambda.",i,".",sep="")]
      if(ZERO_INFLATION)     for( i in 1:t_end[k]) z[i,]     <-Chains[,paste("z.",i,".",sep="")]
      
      Corr   <- ifelse(FECUNDITY=="N",t_end[k],0)
      Corr_z <- ifelse(ZERO_INFLATION,t_end[k],0)
      if(Npars+2+t_end[k]*(Nj+Na+Ne)+Corr+Corr_z==dim(Chains)[2])print(paste("Chains for replicate ",k," extracted successfully. Calculating Posteriors Density",sep=""))
      
      ## assign states
      for (i in 1:x_end[k]){
        for(j in 1:Ne )y_e[j,i,] <- y[j,i,]
        for(j in 1:Nj )y_j[j,i,] <- y[Ne+j,i,]
        for(j in 1:Na )y_a[j,i,] <- y[Ne+Nj+j,i,]
      }
      
      ##########################################################################################
      
      if(lik==F){
        ##########################################################################################
        ### logprior
        ##########################################################################################
        
        if(Priors==1){
          LogPri <- LogPri + log(dlnormTrunc(f,mlf,sdlf,0,30)) + dlnorm(H,mlH,sdlH,log=TRUE)
          if(ZERO_INFLATION)LogPri <- LogPri + log(dnormTrunc(pi,mpi,sdpi,0,1))
          if(FECUNDITY=="N")LogPri <- LogPri + dlnorm(r,mlr,sdlr,log=TRUE)  
          if(TIME_DEP=="T") LogPri <- LogPri + dnorm(HH,mHH,sdHH,log=TRUE)
          if(CROWDING=="T") LogPri <- LogPri + dlnorm(HK,mlHK,sdlHK,log=TRUE)
        }
        if(Priors==2){
          LogPri <- LogPri + dunif(f,0,30,log=TRUE) + dunif(H,0,1,log=TRUE)
          if(ZERO_INFLATION)LogPri <- LogPri + dunif(pi,0,1,log=TRUE)
          if(FECUNDITY=="N")LogPri <- LogPri + dunif(r,0,100,log=TRUE)
          if(TIME_DEP=="T") LogPri <- LogPri + dunif(HH,0,1,log=TRUE)
          if(CROWDING=="T") LogPri <- LogPri + dunif(HK,0,1,log=TRUE)
        } 
        
        ##########################################################################################
        ### loglikelihood (intial conditions)
        ##########################################################################################
        # 
        # ## (should be all 0)
        # for ( j in 1:Ne)LogLikInt <- LogLikInt + dpois(y[j,1,],0,log=TRUE)
        # LogLikInt <- LogLikInt + dbinom(y[Ne+1,1,],3,1,log=T)
        # for ( j in 3:(Nj+Na+1) )LogLikInt <- LogLikInt + dpois(y[j,1,],0,log=TRUE)
        # 
        ##########################################################################################
        ### loglikelihood (Model component)
        ##########################################################################################
        i<-1
        LogLikInt <- rep(0,samples)
        incr_f    <- rep(0,samples) 
        incr_f    <- rep(0,samples) 
        for ( i in 1:(x_end[k]-1) ) {
          
          if(FECUNDITY=="P"){
            if(ZERO_INFLATION){
              incr_z      <- dbinom(z[i+1,],colSums(y_a[,i,]),1-pi,log=T)
              incr_pois   <- dpois(y[1,i+1,],(f/Ne)*z[i+1,]*Dt,log=T)
              incr_f      <- incr_pois + incr_z
            }
            if(ZERO_INFLATION==F)incr_f   <- dpois(y[1,i+1,],(f/Ne)*colSums(y_a[,i,])*Dt,log=T)
          }
          if(FECUNDITY=="N") {
            if(ZERO_INFLATION){
              incr_z        <- dbinom(z[i+1,],colSums(y_a[,i,]),1-pi,log=T)
              incr_lambda   <- dgamma(lambda[i+1,],pmax(r*(ifelse(z[i+1,]<0.5,1,z[i+1,]))/Ne,0.001),pmax(r/(f*Dt),0.001),log=T)
              incr_pois     <- dpois(y[1,i+1,],ifelse(z[i+1,]==0,0,lambda[i+1,]),log = T) 
              incr_f        <- incr_pois + incr_lambda + incr_z
              
              #incr_nb       <- dnbinom(y[1,i+1,],(r*colSums(y_a[,i,]))*Dt/Ne,r/(r+f*Dt),log=T)
            }
            if(ZERO_INFLATION==F){
              incr_lambda   <- dgamma(lambda[i+1,],pmax(r*(ifelse(colSums(y_a[,i,])<0.5,1,colSums(y_a[,i,])))/Ne,0.001),pmax(r/(f*Dt),0.001),log=T)
              incr_pois     <- dpois(y[1,i+1,],ifelse(colSums(y_a[,i,])==0,0,lambda[i+1,]),log = T) 
              incr_f        <- incr_pois + incr_lambda
              
              #incr_nb       <- dnbinom(y[1,i+1,],(r*colSums(y_a[,i,]))*Dt/Ne,r/(r+f*Dt),log=T)
            }
          }
          
          ## debug
          if(DEBUG==T){
            ind.samp<-5
            print(paste("***********************************************",sep=""))
            print(paste("DEBUG replica            ",k,sep=""))
            print(paste("DEBUG time step          ",i,sep=""))
            print(paste("MEAN (m)          ",colSums(y_a[,i,])[ind.samp],sep=""))
            if(FECUNDITY=="NEGBIN"){
              print(paste("fecundity component  lambda ",signif(incr_lambda[ind.samp],digits=5),sep=""))
              print(paste("fecundity component  pois   ",signif(incr_pois[ind.samp],digits=5),sep=""))
              print(paste("fecundity component  nebin   ",signif(incr_nb[ind.samp],digits=5),sep=""))
            }
            if(is.finite(incr_f[ind.samp])==F){
              print(paste("dnbinom(y[1,",i+1,"] = ",y[1,i+1,ind.samp],sep=""))
              print(paste(",r/(3*f*Dt) = ",signif(r[ind.samp],digits=3),"/",signif(3*f[ind.samp],digits=3),sep=""))
              print(paste(",(r*colSums(y_a[,",i,",])+0.01)*Dt/3 = (",signif(r[ind.samp],digits=3),"*",
                          colSums(y_a[,i,])[ind.samp]," + 0.01)/3",sep=""))
              print(paste(",log=T)",sep=""))
              print(paste("is NaN",sep=""))
              #break
            }
          }
          
          LogLikInt <- LogLikInt + incr_f
          #for(j in 1:(Ne-1))LogLikInt <- LogLikInt + dbinom(y[j+1,i+1,],y[j,i,],1,log=T)
          for(j in 1:(Ne-1)){
            increment <- dbinom(y[j+1,i+1,],y[j,i,],.99,log=T)
            iind <- !(is.finite(increment))
            if(any(iind)){
              for(ii in which(iind)){
                
                print(paste("time ",i,"stage ",j,"sample ",ii,sep=""))
                print(paste("state ",y[j,i,ii],"to state ",y[j+1,i+1,ii],sep=""))
              }}
            
            LogLikInt <- LogLikInt + increment
          }
          
          for(j in Ne:(Ne+Nj+Na-1)){
            x0 <- 0.001
            surv <- rep(0,length(H))
            mort <- H*Dt + ifelse(TIME_DEP=="T",HH*i*Dt,0) + ifelse(CROWDING=="T",HK,1/1000000000)*(colSums(y_j[,i,])+colSums(y_a[,i,]))*Dt^2
            mort <- ifelse(mort>=x0,mort,x0*exp(mort/x0-1))
            surv  <- exp(-mort)
            #surv  <- exp(-pmax(mort,surv))
            
            incr_m <- dbinom(y[j+1,i+1,],y[j,i,],surv,log=T)
            LogLikInt <- LogLikInt + incr_m
            
            if(DEBUG==T){
              if(incr_m[ind.samp]!=0){
                #print(paste("DEBUG replica            ",k,sep=""))
                #print(paste("DEBUG time step          ",i,sep=""))
                print(paste("DEBUG age class          ",j,sep=""))
                print(paste("mortality                ",signif(mort[ind.samp],digits=5),sep=""))
                print(paste("survival                 ",signif(surv[ind.samp],digits=5),sep=""))
                print(paste("mortality component      ",signif(incr_m[ind.samp],digits=5),sep=""))
              }
            }
          }
          ## debug
          if(DEBUG==T){
            print(paste("***********************************************",sep=""))
            print(paste("Likelihood      ",signif(LogLikInt[ind.samp],digits=5),sep=""))
            print(paste("mean Likelihood ",signif(mean(LogLikInt),digits=5),sep=""))
            
            nns <- length(which(is.finite(LogLikInt)==F))
            if(nns>0){print(paste(nns," NaNS appear in the likelihood at replica ",k," time ",i,sep=""))
              print(which(is.finite(LogLikInt)==F))
            }
            if(length(which(is.finite(LogLikInt)==F))>0)print(paste("Problem in Replica_",k,"_time_",i,sep=""))
          }
          
          if(length(which(is.finite(LogLikInt)==F))>0){
            ind.samp <- 10
            print(paste("infinite in Replica_",k,"_time_",i,sep=""))
            print(paste("***********************************************",sep=""))
            print(paste("Likelihood      ",signif(LogLikInt[ind.samp],digits=5),sep=""))
            print(paste("mean Likelihood ",signif(mean(LogLikInt),digits=5),sep=""))
            print(paste("fertility      ",signif(incr_f[ind.samp],digits=5),sep=""))
            print(paste("mortality      ",signif(incr_m[ind.samp],digits=5),sep=""))
          }  
          
        }
        ##########################################################################################
      }## close lik = F
      
      ##########################################################################################
      ### Observation part
      ##########################################################################################
      
      for ( i in 1:n_obs_to_fit[k] ) {
        LogLik     <-  LogLik + dnorm(y_e_obs_to_fit_rep[i],colSums(y_e[,x_obs_to_fit_rep[i],]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_e_obs_to_fit_rep[i],sum(apply(y_e[,x_obs_to_fit_rep[i],],1,mean)),1,log=TRUE)
        
        LogLik     <-  LogLik + dnorm(y_j_obs_to_fit_rep[i],colSums(y_j[,x_obs_to_fit_rep[i],]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_j_obs_to_fit_rep[i],sum(apply(y_j[,x_obs_to_fit_rep[i],],1,mean)),1,log=TRUE)
        
        LogLik     <-  LogLik + dnorm(y_a_obs_to_fit_rep[i],colSums(y_a[,x_obs_to_fit_rep[i],]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_a_obs_to_fit_rep[i],sum(apply(y_a[,x_obs_to_fit_rep[i],],1,mean)),1,log=TRUE)
      }
      ##LogLik + LogPri
      ##########################################################################################
    }### close single replicate estimation
    
  }### close INF_METH "R"
  
  if(INF_METH=="H"){
    
    ##########################################################################################
    ### prepare structure, extract and check chains
    ##########################################################################################
    
    filename <- paste(FILE_ACRONYM,MODEL_NAME,"_",treat,"_",clone,"_",Priors,".csv",sep="")
    Chains  <- read.csv(filename)
    samples <- dim(Chains)[1]
    
    LogPri    <- numeric(samples)
    LogLik    <- numeric(samples)
    LogLikInt <- numeric(samples)
    LogLikMean <- 0
    
    y  <- array(0,dim=c(Na+Nj+Ne,sum(t_end),samples))
    y_e<- array(0,dim=c(Ne,sum(x_end),samples))
    y_j<- array(0,dim=c(Nj,sum(x_end),samples))
    y_a<- array(0,dim=c(Na,sum(x_end),samples))
    
    if(ZERO_INFLATION) z <- array(0,dim=c(sum(t_end),samples))
    if(FECUNDITY=="N") lambda <- array(0,dim=c(sum(t_end),samples))
    
    f_f <- array(0,dim=c(R,samples))
    H_H <- array(0,dim=c(R,samples))
    
    m_f  <- Chains[,"m_f"]
    sd_f <- Chains[,"sd_f"]
    for(kk in 1:R) f_f[kk,] <- Chains[,paste("f.",kk,".",sep="")]
    
    m_H <- Chains[,"m_H"]
    sd_H <- Chains[,"sd_H"]
    for(kk in 1:R) H_H[kk,] <- Chains[,paste("H.",kk,".",sep="")]
    
    if(ZERO_INFLATION){
      pi_pi  <- array(0,dim=c(R,samples))
      m_pi   <- Chains[,"m_pi"]
      sd_pi  <- Chains[,"sd_pi"]
      for(kk in 1:R) pi_pi[kk,] <- Chains[,paste("pi.",kk,".",sep="")]
    }
    if(FECUNDITY=="N"){
      r_r <- array(0,dim=c(R,samples))
      m_r  <- Chains[,"m_r"]
      sd_r  <- Chains[,"sd_r"]
      for(kk in 1:R) r_r[kk,] <- Chains[,paste("r.",kk,".",sep="")]
    }
    if(TIME_DEP=="T"){
      HH_HH <- array(0,dim=c(R,samples))
      m_HH <- Chains[,"m_HH"]
      sd_HH <- Chains[,"sd_HH"]
      for(kk in 1:R) HH_HH[kk,] <- Chains[,paste("HH.",kk,".",sep="")]
    }
    if(CROWDING=="T"){
      HK_HK <- array(0,dim=c(R,samples))
      m_HK <- Chains[,"m_HK"] 
      sd_HK <- Chains[,"sd_HK"] 
      for(kk in 1:R) HK_HK[kk,] <- Chains[,paste("HK.",kk,".",sep="")]
    }  
    
    if(FECUNDITY=="N") {
      for(k in 1:R ){
        for( i in xx_s[k]:(xx_e[k]) ) {
          lambda[i,]<-Chains[,paste("lambda.",i,".",sep="")]
        }}
    }
    
    if(ZERO_INFLATION){
      for(k in 1:R){
        for( i in xx_s[k]:(xx_e[k]) ) {
          z[i,] <- Chains[,paste("z.",i,".",sep="")]
        }}
    }
    
    for(ii in 1:(Ne+Na+Nj)) for(jj in 1:sum(t_end)) y[ii,jj,]<-Chains[,paste("y.",ii,".",jj,".",sep="")]
    
    Corr   <- ifelse(FECUNDITY=="N",(sum(t_end)),0)
    Corr_z <- ifelse(ZERO_INFLATION,(sum(t_end)),0)
    if(7*Npars+2+sum(t_end)*(Nj+Na+Ne)+Corr+Corr_z==dim(Chains)[2])print("Chains extracted successfully. Calculating Posteriors Density")
    
    ##########################################################################################
    ## assign states
    ##########################################################################################
    
    for ( k in 1:R ){ 
      for (i in xx_s[k]:xx_e[k]){
        for(j in 1:Ne )y_e[j,i,] <- y[j,i,]
        for(j in 1:Nj )y_j[j,i,] <- y[Ne+j,i,]
        for(j in 1:Na )y_a[j,i,] <- y[Ne+Nj+j,i,]
      }
    }
    ##########################################################################################
    
    if(lik==F){
      
      ##########################################################################################
      ### logprior
      ##########################################################################################
      
      if(Priors==1){
        LogPri <- LogPri + log(dlnormTrunc(m_f,ml_mf,sdl_mf,0,30)) + dlnorm(sd_f,ml_sdf,sdl_sdf,log=TRUE)
        for(kk in 1:R) LogPri <- LogPri + log(dlnormTrunc(f_f[kk,],log(m_f)-log(1+sd_f^2/m_f^2)/2, sqrt(log(1+(sd_f/m_f)^2)),0,30))
        
        LogPri <- LogPri + dlnorm(m_H,ml_mH,sdl_mH,log=TRUE) + dlnorm(sd_H,ml_sdH,sdl_sdH,log=TRUE)
        for(kk in 1:R) LogPri <- LogPri + dlnorm(H_H[kk,],log(m_H)-log(1+sd_H^2/m_H^2)/2, sqrt(log(1+(sd_H/m_H)^2)),log=TRUE)
        
        if(ZERO_INFLATION){
          LogPri <- LogPri + log(dnormTrunc(m_pi,m_mpi,sd_mpi,0,1)) + dlnorm(sd_pi,ml_sdpi,sdl_sdpi,log=TRUE)
          for(kk in 1:R) LogPri <- LogPri + log(dnormTrunc(pi_pi[kk,],m_pi,sd_pi,0,1))
        }
        if(FECUNDITY=="N"){
          LogPri <- LogPri + dlnorm(m_r,ml_mr,sdl_mr,log=TRUE) + dlnorm(sd_r,ml_sdr,sdl_sdr,log=TRUE)
          for(kk in 1:R) LogPri <- LogPri + dlnorm(r_r[kk,],log(m_r)-log(1+sd_r^2/m_r^2)/2, sqrt(log(1+(sd_r/m_r)^2)),log=TRUE)
        }
        if(TIME_DEP=="T"){
          LogPri <- LogPri + dnorm(m_HH,m_mHH,sd_mHH,log=TRUE) + dlnorm(sd_HH,ml_sdHH,sdl_sdHH,log=TRUE)
          for(kk in 1:R) LogPri <- LogPri + dnorm(HH_HH[kk,],m_HH,sd_HH,log=TRUE)
        }
        if(CROWDING=="T"){
          LogPri <- LogPri + dlnorm(m_HK,ml_mHK,sdl_mHK,log=TRUE) + dlnorm(sd_HK,ml_sdHK,sdl_sdHK,log=TRUE)
          for(kk in 1:R) LogPri <- LogPri + dlnorm(HK_HK[kk,],log(m_HK)-log(1+sd_HK^2/m_HK^2)/2, sqrt(log(1+(sd_HK/m_HK)^2)),log=TRUE)
        }
        
      }  
      if(Priors==2){
        for(kk in 1:R) LogPri <- LogPri + dunif(f_f[kk,],0,30,log=TRUE)
        for(kk in 1:R) LogPri <- LogPri + dunif(H_H[kk,],0,1,log=TRUE)
        if(ZERO_INFLATION){
          for(kk in 1:R) LogPri <- LogPri + dunif(pi_pi[kk,],0,1,log=TRUE)
        }
        if(FECUNDITY=="N"){
          for(kk in 1:R) LogPri <- LogPri + dunif(r_r[kk,],0,100,log=TRUE)
        }
        if(TIME_DEP=="T"){
          for(kk in 1:R) LogPri <- LogPri + dunif(HH_HH[kk,],0,1,log=TRUE)
        }
        if(CROWDING=="T"){
          for(kk in 1:R) LogPri <- LogPri + dunif(HK_HK[kk,],0,1,log=TRUE)
        }
      }
      
      ##########################################################################################
      ### loglikelihood (intial conditions)
      ##########################################################################################
      
      ## (should be all 0)
      # for ( k in 1:R ) {
      #   for ( j in 1:Ne )LogLikInt <- LogLikInt + dpois(y[j,xx_s[k],],0,log=TRUE)
      #   LogLikInt <- LogLikInt + dbinom(y[Ne+1,xx_s[k],],3,1,log=TRUE)
      #   for ( j in (Ne+2):(Nj+Na+Ne) )LogLikInt <- LogLikInt + dpois(y[j,xx_s[k],],0,log=TRUE)
      # }
      
      ##########################################################################################
      ### loglikelihood (Model component)
      ##########################################################################################
      
      k<-1
      for(k in 1:R ) {
        
        i<-1
        for ( i in xx_s[k]:(xx_e[k]-1) ) {
          
          if(FECUNDITY=="P"){
            if(ZERO_INFLATION){
              incr_z      <- dbinom(z[i+1,],colSums(y_a[,i,]),1-pi_pi[k,],log=T)
              incr_pois   <- dpois(y[1,i+1,],f_f[k,]*z[i+1,]*Dt/Ne,log=T)
              incr_f      <- incr_pois + incr_z
            }
            if(ZERO_INFLATION==F)incr_f   <- dpois(y[1,i+1,],f_f[k,]*colSums(y_a[,i,])*Dt/Ne,log=T)
          }
          if(FECUNDITY=="N") {
            if(ZERO_INFLATION){
              incr_z        <- dbinom(z[i+1,],colSums(y_a[,i,]),1-pi_pi[k,],log=T)
              incr_lambda   <- dgamma(lambda[i+1,],pmax(r_r[k,]*(ifelse(z[i+1,]<0.5,1,z[i+1,]))/Ne,0.001),pmax(r_r[k,]/(f_f[k,]*Dt),0.001),log=T)
              incr_pois     <- dpois(y[1,i+1,],ifelse(z[i+1,]==0,0,lambda[i+1,]),log = T) 
              incr_f        <- incr_pois + incr_lambda + incr_z
              
              #incr_nb       <- dnbinom(y[1,i+1,],(r_r[k,]*colSums(y_a[,i,]))*Dt/3,r_r[k,]/(r_r[k,]+3*f_f[k,]*Dt),log=T)
            }
            if(ZERO_INFLATION==F){
              incr_lambda   <- dgamma(lambda[i+1,],pmax(r_r[k,]*(ifelse(colSums(y_a[,i,])<0.5,1,colSums(y_a[,i,])))/Ne,0.001),pmax(r_r[k,]/(f_f[k,]*Dt),0.001),log=T)
              incr_pois     <- dpois(y[1,i+1,],ifelse(colSums(y_a[,i,])==0,0,lambda[i+1,]),log = T) 
              incr_f        <- incr_pois + incr_lambda
              
              #incr_nb       <- dnbinom(y[1,i+1,],(r_r[k,]*colSums(y_a[,i,]))*Dt/3,r_r[k,]/(r_r[k,]+3*f_f[k,]*Dt),log=T)
            }
          }
          
          ## debug
          if(DEBUG==T){
            ind.samp<-10
            print(paste("***********************************************",sep=""))
            print(paste("DEBUG replica            ",k,sep=""))
            print(paste("DEBUG time step          ",i,sep=""))
            print(paste("MEAN (m)          ",colSums(y_a[i,ind.samp]),sep=""))
            print(paste("fecundity component  lambda ",signif(incr_lambda[ind.samp],digits=5),sep=""))
            print(paste("fecundity component  pois   ",signif(incr_pois[ind.samp],digits=5),sep=""))
            print(paste("fecundity component  nebin   ",signif(incr_nb[ind.samp],digits=5),sep=""))
            
            if(is.finite(incr_f[ind.samp])==F){
              print(paste("dnbinom(y[1,",i+1,"] = ",y[1,i+1,ind.samp],sep=""))
              print(paste(",r_r[",k,"]/(r_r[",k,"]+3*f_f[",k,"]*Dt) = ",signif(r_r[k,ind.samp],digits=3),"/",signif(3*f_f[k,ind.samp],digits=3),sep=""))
              #print(paste(",(r*colSums(y_a[,",i,",])+0.01)*Dt/3 = (",signif(r[ind.samp],digits=3),"*",
              #           colSums(y_a[,i,])[ind.samp]," + 0.01)/3",sep=""))
              #print(paste(",log=T)",sep=""))
              #print(paste("is NaN",sep=""))
              #break
            }
          }
          LogLikInt <- LogLikInt + incr_f
          
          for(j in 1:(Ne-1)){
            increment <- dbinom(y[j+1,i+1,],y[j,i,],.99,log=T)
            iind <- !(is.finite(increment))
            if(any(iind)){
              for(ii in which(iind)){
                
                print(paste("time ",i,"stage ",j,"sample ",ii,sep=""))
                print(paste("state ",y[j,i,ii],"to state ",y[j+1,i+1,ii],sep=""))
              }}
            
            LogLikInt <- LogLikInt + increment
          }
          
          
          j<-1
          for(j in Ne:(Ne+Nj+Na-1)){
            x0    <- 0.001
            surv  <- rep(0,length(H_H[k,]))
            mort  <- H_H[k,]*Dt+ ifelse(TIME_DEP=="T",HH_HH[k,]*(i-xx_s[k])*Dt,0) + ifelse(CROWDING=="T",HK_HK[k,],1/1000000000)*(colSums(y_j[,i,])+colSums(y_a[,i,]))*Dt^2
            mort  <- ifelse(mort>=x0,mort,x0*exp(mort/x0-1))
            surv  <- exp(-mort)
            #surv  <- exp(-pmax(mort,surv))
            
            incr_m  <- dbinom(y[j+1,i+1,],y[j,i,],surv,log=T)
            
            LogLikInt<- LogLikInt + incr_m
            
            
            
            # surv  <- exp(-pmax(H_H[k,]*Dt + ifelse(TIME_DEP=="T",HH_HH[k,]*(i-xx_s[k])*Dt,0) + ifelse(CROWDING=="T",HK_HK[k,],1/10000)*(colSums(y_j[,i,])+colSums(y_a[,i,]))*Dt^2,0))
            # LogLikInt<- LogLikInt + dbinom(y[j+1,i+1,],y[j,i,],surv,log=T)
          } ##close loop in age classes
          
          
          if(length(which(is.finite(LogLikInt)==F))>0)print(paste("NaNs after Replica_",k,"_time_",i,sep=""))
        }# close loop in time steps
      }# close loop in replicates
      
    }### close lik==F
    
    ##########################################################################################
    ### Observation part (deviance)
    ##########################################################################################
    
    k<-1;i<-1
    for ( k in 1:R ) {
      for ( i in nn_s[k]:nn_e[k] ) {
        LogLik     <-  LogLik     + dnorm(y_e_obs_to_fit[i],colSums(y_e[,xx_s[k]+x_obs_to_fit[i]-1,]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_e_obs_to_fit[i],sum(apply(y_e[,xx_s[k]+x_obs_to_fit[i]-1,],1,mean)),1,log=TRUE)
        
        LogLik     <-  LogLik + dnorm(y_j_obs_to_fit[i],colSums(y_j[,xx_s[k]+x_obs_to_fit[i]-1,]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_j_obs_to_fit[i],sum(apply(y_j[,xx_s[k]+x_obs_to_fit[i]-1,],1,mean)),1,log=TRUE)
        
        LogLik     <-  LogLik + dnorm(y_a_obs_to_fit[i],colSums(y_a[,xx_s[k]+x_obs_to_fit[i]-1,]),1,log=TRUE)
        LogLikMean <-  LogLikMean + dnorm(y_a_obs_to_fit[i],sum(apply(y_a[,xx_s[k]+x_obs_to_fit[i]-1,],1,mean)),1,log=TRUE)
      }
    }
    ## LogLik + LogPri  
    
    
  }## close INF_METH "H"
  
  ## return likelihood oder posterior
  if(lik==F)return(cbind(LogPri,LogLikInt,LogLik))
  if(lik==T)return(LogLikMean)
}

#########################################################################
### consistency function (work in progress) 
#########################################################################

make_values_consistent <- function(xx_s,xx_e,y){
  for ( k in 1:length(xx_e) ) {
    for ( j in (xx_s[k]+1):xx_e[k] ) {
      if(sum(y[(Nj+2):(Nj+Na+1),j-1])==0){
        if(y[1,j]>0) {
          print("inconsistency in egg numbers removed")
          print(j)
          y[1,j]<-0
        }
      }
    }}
  
  #ddii$y[1,]     <- 0
  y[,xx_s]  <- 0
  y[2,xx_s] <- 3
  for ( i in 2:nrow(y) ) {
    for ( k in 1:length(xx_e) ) {
      for ( j in (xx_s[k]+1):xx_e[k] ) {
        if(y[i,j]>y[i-1,j-1]){ 
          print("inconsistency in age structure removed")
          print(j)
        }
        
        y[i,j] <- min(y[i,j],y[i-1,j-1])
      }
    }
  }
  return(y)
}


#########################################################################
### functions for plotting population time series
#########################################################################

plot_shade <- function(Predictor,Med, U95, L95,YMin,YMax,title){
  gray<-rgb(0.8,0.8,0.8,1)
  black<-rgb(0,0,0,1)
  plot(Predictor,Med,type="l",cex=1,col=black,bg=black,ylim=c(YMin,YMax),lwd=2, 
       xlab="time (days)",ylab="individuals",main=title,
       cex.lab=1.5,cex.axis=1.5,cex.main=1.5,font.main=1)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(U95,rev(Med))
  polygon(xs,ys, col=gray, border=gray)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(Med,rev(L95))
  polygon(xs,ys, col=gray, border=gray)
  lines(Predictor, Med,type="l",cex=1,col=black,bg=black,ylim=c(YMin,YMax), xlab="", ylab="",xaxt='n',yaxt='n',lwd=2)
}

Plot_two_models <- function(Predictor, Med1, U951, L951, Med2, U952, L952,YMin,YMax,title){
  
  gray<-rgb(0.8,0.8,0.8,1)
  black<-rgb(0,0,0,1)
  plot(Predictor, Med1,type="l",cex=1,col=black,bg=black,ylim=c(YMin,YMax), 
       xlab="time (days)",ylab="individuals",main=title,
       cex.lab=1.5,cex.axis=1.5,cex.main=1.5,font.main=1)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(U951,rev(Med1))
  polygon(xs,ys, col=gray, border=gray)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(Med1,rev(L951))
  polygon(xs,ys, col=gray, border=gray)
  lines(Predictor, Med1,type="l",cex=1,col=black,bg=black,ylim=c(YMin,YMax), xlab="", ylab="",xaxt='n',yaxt='n',lwd=2)
  
  red<-rgb(1,0,0,1)
  pinkt<-rgb(1,0.6,0.6,0.4)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(U952,rev(Med2))
  polygon(xs,ys, col=pinkt, border=pinkt)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(Med2,rev(L952))
  polygon(xs,ys, col=pinkt, border=pinkt)
  lines(Predictor, Med2,type="l",cex=1,col=red,bg=red,ylim=c(YMin,YMax), xlab="", ylab="",xaxt='n',yaxt='n',lwd=2)
}

Plot_three_models <- function(Predictor, Med1, U951, L951,Med2, U952, L952,Med3, U953, L953,YMin,YMax,title){
  
  gray<-rgb(0.8,0.8,0.8,1)
  black<-rgb(0,0,0,1)
  plot(Predictor, Med1,type="l",cex=1,col=black,bg=black,ylim=c(YMin,YMax), 
       xlab="time (days)",ylab="individuals",main=title,
       cex.lab=1.5,cex.axis=1.5,cex.main=1.5,font.main=1)
  
  xs<-append(Predictor,rev(Predictor))
  ys<-append(U951,rev(Med1))
  polygon(xs,ys, col=gray, border=gray)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(Med1,rev(L951))
  polygon(xs,ys, col=gray, border=gray)
  lines(Predictor, Med1,type="l",cex=1,col=black,bg=black,
        ylim=c(YMin,YMax), xlab="", ylab="",xaxt='n',yaxt='n',lwd=2)
  
  red<-rgb(0,0,1,1)
  pinkt<-rgb(0.6,0.6,1,0.4)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(U952,rev(Med2))
  polygon(xs,ys, col=pinkt, border=pinkt)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(Med2,rev(L952))
  polygon(xs,ys, col=pinkt, border=pinkt)
  lines(Predictor, Med2,type="l",cex=1,col=red,bg=red,ylim=c(YMin,YMax), xlab="", ylab="",xaxt='n',yaxt='n',lwd=2)
  
  green  <-rgb(1,0,0,1)
  greenkt<-rgb(1,0.6,0.6,0.4)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(U953,rev(Med3))
  polygon(xs,ys, col=greenkt, border=greenkt)
  xs<-append(Predictor,rev(Predictor))
  ys<-append(Med3,rev(L953))
  polygon(xs,ys, col=greenkt, border=greenkt)
  lines(Predictor, Med3,type="l",cex=1,col=green,bg=green,
        ylim=c(YMin,YMax), xlab="", ylab="",xaxt='n',yaxt='n',lwd=2)
  
}

## quantile functions
lb <- function(x) quantile(x,prob=0.025,na.rm=TRUE)
ub <- function(x) quantile(x,prob=0.975,na.rm=TRUE)
#fn <- function(x) 3*sd(x)

