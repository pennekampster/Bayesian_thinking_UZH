############################################################################################################
### prototype to play and make inference with Daphnia data
### Workshop Bayesian Thinking in Ecology   UZH 06-10-2021

### Code written by Gian Marco Palamara
############################################################################################################

## set folder if not working on server
## setwd("Path/to/your/folder")

### set language
Sys.setenv(LANG = "en")

### clear memory and set seed
rm(list=ls())
dev.off()
set.seed(54321)

### plotting of output
### if TRUE everythin is plotted into the termial
### if FALSE a pdf with all the outputs is createed
PLOT_IN_TERMINAL      <- T


## select clone 
## CLONE_TEST D1  --> Taken From Dutch lake 
## CLONE_TEST D2  --> Taken from Italian lake (very polluted)
## CLONE_TEST D3  --> Taken from swiss Lake
## CLONE_TEST Ds  --> mixture of the three clones
CLONE_TEST     <- "D3"

## select treatment 
## TREAT_TEST cont --> Control treatment
## TREAT_TEST dr   --> Herbicide treatment (diuron)
## TREAT_TEST dz   --> pesticide treatment (diazinon)
## TREAT_TEST dr_dz--> booth herbicide and pesticide treatment
TREAT_TEST     <- "cont"

## select Priors 
## PRIORS <- 1 informative priors
## PRIORS <- 2 uniform priors
PRIORS         <- 1

## Choice of FECUNDITY distribution
## P Poisson 
## N Negative Binomial
FECUNDITY      <- "P"
## Zero inflation TRUE or FALSE
ZERO_INFLATION <- "F"

## Choice of MORTALITY distribution
## Time dependnet mortality "TRUE" or "FALSE"
TIME_DEP       <- "F"
## Desntiy dependence "TRUE" or "FALSE"
CROWDING       <- "F"

## Inference method
## INF_METH <- J Joint fit;
## INF_METH <- R rep by rep fit;
## INF_METH <- H hierarchical fit;
INF_METH       <- "H"

############################################################################################################
### MACROS FOR THE SCRIPT
############################################################################################################

### getting and plotting data (they need to be always TRUE)
GET_REAL_DATA  <- T
PLOT_REAL_DATA <- T

### inference
INFERENCE      <- T

### Set chain_length, number of chains and thinning
CHAIN_LENGTH <- 1000
N_CHAINS     <- 1
THIN         <- 1
THIN_PLOT    <- 1

### plotting posterior
PLOT_POSTERIOR        <- T
### keep it FALSE
STATE_RESOLUTION_PLOT <- F
## plot of correlaion plots
PLOT_CORR             <- T

### More technical details
INITIALIZE_CHAINS     <- F
INTERPOLATE           <- F
MAKE_CONSISTENT       <- F

### If TRUE computes DIC using my own functions
### and compares it to JAGS outputs
COMPUTE_DIC           <- T

## if check chain T plots time series and predicted fit
CHECK_CHAIN <- T

## plots one (only inferred states) or two (inferred states and inferred trajectories)
## or three trajectories (if only with hierarchical models)
PLOT_SHADE   <- F
PLOT_SHADE_2 <- T
PLOT_SHADE_3 <- T
## plot of all replica all together with joint fit
PLOT_ALL     <- T


############################################################################################################
### START OF THE PROGRAM, LOAD LIBRARIES
############################################################################################################

## install inference packages (coda and Rjags)
if ( !require("coda") )  { install.packages("coda");  library("coda")  }
if ( !require("rjags") ) { install.packages("rjags"); library("rjags") }
if ( !require("EnvStats") ) { install.packages("EnvStats"); library("EnvStats") }
load.module("dic")
# note that jags must be installed separately from http://mcmc-jags.sourceforge.net/

## package for correlation plots
if ( ! require(IDPmisc) )   { install.packages("IDPmisc");   library(IDPmisc) }

## load functions
source("FUNCTIONS_CODE.R")

############################################################################################################
## argoments passed from terminal
############################################################################################################

com <- commandArgs()

## Macros for treatment and model
## select treatment and clone
if(length(com)>2){
  TREAT_TEST     <- com[3]
  CLONE_TEST     <- com[4]
  ## select confiuration of priors (1 informative, 2 uniform) 
  PRIORS         <- as.numeric(com[5])
  ## select model structure 
  FECUNDITY      <- substr(com[6], 1, 1) ## POISSON or NEGBIN
  ZERO_INFLATION <- com[7] ## True or False
  TIME_DEP       <- com[8] ## True or False
  CROWDING       <- com[9] ## True or False
  ## select inference method  
  INF_METH       <- com[10] ## (J) Joint fit, (R) Rep by Rep fit, (H) HIerarchical fit 
  #OBS_PROC      <- com[11]
}## potential to add flexible observation process

MODEL_NAME <- paste(ifelse(FECUNDITY=="P","P","N"),
                    ifelse(ZERO_INFLATION,"Z",""),
                    ifelse(TIME_DEP=="T","T",""),
                    ifelse(CROWDING=="T","C",""),
                    INF_METH,             sep="")

### set file names
file.jags.model  <- paste("rjags_test_model_",MODEL_NAME,"_",TREAT_TEST,"_",CLONE_TEST,"_",PRIORS,".jags",sep="")
file.plot.output <- paste("OUTPUT_",MODEL_NAME,"_",TREAT_TEST,"_",CLONE_TEST,"_",PRIORS,".pdf",sep="")

if(PLOT_IN_TERMINAL==F) pdf(file.plot.output)

############################################################################################################
### GET AND PLOT REAL DATA
############################################################################################################
## get real data

if(GET_REAL_DATA==T){
  
  treatments <- c("cont","dr","dz","dr_dz")
  clones     <- c("D1","D2","D3","Ds")
  
  ## load dataframe and create list to put data for fitting
  dd.list <- prepare.LTE.data.list(filename="Daphnia_LTE_master_20171013.csv",separator=";")
  
  dd.to.fit.list <- list()
  
  ## treatments and clones loop
  for(tt in treatments){
    for(cc in clones){
      
      ## get all the replicates per treatment per clone
      TREATMENT <- tt
      CLONE <- cc
      dd.eval <- dd.list[[TREATMENT]][[CLONE]]
      
      ## count the number of replicates
      ## and prepare arrays 
      R <- length(dd.eval)
      n_obs <- numeric(R)
      t_end <- numeric(R)
      
      t_obs <- numeric()
      y_e_obs <- numeric()
      y_j_obs <- numeric()
      y_a_obs <- numeric()
      
      
      k<-5
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
      
    }###close clones cycle
  }#close treatments cycle
}
## data ready to be fitted to the model    

if(PLOT_REAL_DATA==T){
  
  #########################################################################
  ## plot real time series  
  # treatments <- c("cont","dr","dz","dr_dz")
  # clones <- c("D1","D2","D3","Ds")
  
  treatments <- TREAT_TEST
  clones     <- CLONE_TEST
  
  
  ## treatments and clones loop
  for(tt in treatments){
    for(cc in clones){
      
      ## extract from the list
      n_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][1])))
      t_end <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][2])))
      
      t_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,1])))
      y_e_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,2])))
      y_j_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,3])))
      y_a_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,4])))
      # get plotting parameters
      R <- length(n_obs)
      
      nn_s<- numeric(R)
      nn_e<- numeric(R)
      for(k in 1:R){
        nn_s[k] <- ifelse(k==1,1,sum(n_obs[1:(k-1)])+1)
        nn_e[k] <- sum(n_obs[1:k])
      }
      x_max <- max(t_obs)
      e_max <- max(y_e_obs)
      j_max <- max(y_j_obs)
      a_max <- max(y_a_obs)
      
      #check
      if(max(t_end)==x_max) print("ok")
      
      ## graphical setting
      #X11()
      par(mfcol=c(3,2))
      
      
      #########################################################################
      ## plot 3 classes
      
      ## eggs
      x <- t_obs[1:n_obs[1]]
      y <- y_e_obs[1:n_obs[1]]
      plot(x,y,main=paste("Eggs_",tt,"_",cc),type="b",xlab="time",ylab="eggs",
           xlim=c(0,x_max),ylim=c(0,e_max),lwd=2)
      if(R>1) for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        y <- y_e_obs[nn_s[k]:nn_e[k]]
        lines(x,y,lwd=2,type="b",col=k)
      }
      ## juveniles
      x <- t_obs[1:n_obs[1]]
      y <- y_j_obs[1:n_obs[1]]
      plot(x,y,main=paste("Juveniles_",tt,"_",cc),type="b",xlab="time",ylab="juveniles",
           xlim=c(0,x_max),ylim=c(0,j_max),lwd=2)
      if(R>1) for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        y <- y_j_obs[nn_s[k]:nn_e[k]]
        lines(x,y,lwd=2,type="b",col=k)
      }
      ## adults
      x <- t_obs[1:n_obs[1]]
      y <- y_a_obs[1:n_obs[1]]
      plot(x,y,main=paste("Adults_",tt,"_",cc),type="b",xlab="time",ylab="adults",
           xlim=c(0,x_max),ylim=c(0,a_max),lwd=2)
      if(R>1) for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        y <- y_a_obs[nn_s[k]:nn_e[k]]
        lines(x,y,lwd=2,type="b",col=k)
      }
      
      #########################################################################
      ## plot phase space diagrams
      
      ## Juveniles vs Adults
      x <- y_j_obs[1:n_obs[1]]
      y <- y_a_obs[1:n_obs[1]]
      plot(x,y,main=paste(tt,"_",cc),type="b",xlab="Juveniles",ylab="adults",
           xlim=c(0,j_max),ylim=c(0,a_max),lwd=2)
      if(R>1) for(k in 2:R){
        x <- y_j_obs[nn_s[k]:nn_e[k]]
        y <- y_a_obs[nn_s[k]:nn_e[k]]
        lines(x,y,lwd=2,type="b",col=k)
      }
      
      ## Eggs vs Adults
      x <- y_e_obs[1:n_obs[1]]
      y <- y_a_obs[1:n_obs[1]]
      plot(x,y,main=paste(tt,"_",cc),type="b",xlab="eggs",ylab="adults",
           xlim=c(0,e_max),ylim=c(0,a_max),lwd=2)
      if(R>1) for(k in 2:R){
        x <- y_e_obs[nn_s[k]:nn_e[k]]
        y <- y_a_obs[nn_s[k]:nn_e[k]]
        lines(x,y,lwd=2,type="b",col=k)
      }
      
      ## Eggs vs Juveniles
      x <- y_e_obs[1:n_obs[1]]
      y <- y_j_obs[1:n_obs[1]]
      plot(x,y,main=paste(tt,"_",cc),type="b",xlab="eggs",ylab="juveniles",
           xlim=c(0,e_max),ylim=c(0,j_max),lwd=2)
      if(R>1) for(k in 2:R){
        x <- y_e_obs[nn_s[k]:nn_e[k]]
        y <- y_j_obs[nn_s[k]:nn_e[k]]
        lines(x,y,lwd=2,type="b",col=k)
      }
      legend("bottomright",title="REP",legend=c("1","2","3","4","5"),col=c(1,2,3,4,5),lwd=2)
      
    }}## close treatments clones cycles
}

############################################################################################################
### SET PRIORS 
############################################################################################################
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

############################################################################################################
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


############################################################################################################
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


############################################################################################################
### INFERENCE AND PLOTTING OF INFERENCE RESULTS
############################################################################################################

if(INFERENCE==T){
  
  ############################################################################################################
  ### prepare data for fitting
  ############################################################################################################
  
  ## time step for model and observations and time of simulation (days)
  Dt    <- 1  
  
  ## maturation age and maximum age (days)
  age.max.j <-  5
  age.max.a <-  50 
  
  ## instar time (number of embryo classes)
  Ne <- INSTAR_TIME
  ## maturation time (number of juvenile classes) use an even number
  Nj <- round(age.max.j/Dt)
  ## adult time (number of adult classes) use an even number
  Na <- round((age.max.a-age.max.j)/Dt)
  
  ## all treatments, to put in loops 
  tt <- TREAT_TEST
  cc <- CLONE_TEST
  
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
  
  ## some check
  if(length(x_obs_to_fit)==sum(n_obs_to_fit))print("Data ready to be fitted")
  
  ############################################################################################################
  #### inference methods
  ############################################################################################################
  
  if(INF_METH=="J"){
    
    ######################################################################
    ## jags function
    ######################################################################
    
    model.code <- paste("var y_e[Ne,sum(x_end)],y_j[Nj,sum(x_end)],y_a[Na,sum(x_end)]",
                        ifelse(ZERO_INFLATION=="T",",z[sum(x_end)]",""),
                        ifelse(FECUNDITY=="N",",lambda[sum(x_end)]\n","\n"),
                        "model {\n",
                        "  f ~ ",ifelse(PRIORS==1,"dlnorm(mlf,tauf)I(0,30)\n","dunif(0,30)\n"),                        
                        "  H ~ ",ifelse(PRIORS==1,"dlnorm(mlH,tauH)\n","dunif(0,1)\n"),sep="")
    if(ZERO_INFLATION=="T") model.code <- paste(model.code," pi ~ ",ifelse(PRIORS==1,"dnorm(mpi,(1/sdpi)^2)I(0,1)\n","dunif(0,1)\n"),sep="")
    if(FECUNDITY=="N")      model.code <- paste(model.code," r  ~ ",ifelse(PRIORS==1,"dlnorm(mlr,taur)\n","dunif(0,100)\n"),sep="")
    if(TIME_DEP=="T")       model.code <- paste(model.code," HH ~ ",ifelse(PRIORS==1,"dnorm(mHH,(1/sdHH)^2)\n","dunif(0,1)\n"),sep="")
    if(CROWDING=="T")       model.code <- paste(model.code," HK ~ ",ifelse(PRIORS==1,"dlnorm(mlHK,tauHK)\n","dunif(0,1)\n"),sep="")
    model.code <- paste(model.code,
                        "  for ( k in 1:R ) {\n",
                        "    for ( j in 1:Ne ) {\n",
                        "      y[j,xx_s[k]] ~ dpois(0)\n",
                        "    }\n",
                        "    y[Ne+1,xx_s[k]] ~ dbin(1,3)\n",
                        "    for ( j in (Ne+2):(Ne+Nj+Na) ) {\n",
                        "      y[j,xx_s[k]] ~ dpois(0)\n",
                        "    }\n",
                        "  }\n",
                        sep="")
    
    if(FECUNDITY=="N"){
      model.code<-paste(model.code,
                        "  for ( k in 1:R ) {\n",
                        "    lambda[xx_s[k]] ~ dpois(0)\n",
                        "  }\n",sep="")
    }
    if(ZERO_INFLATION=="T"){
      model.code<-paste(model.code,
                        "  for ( k in 1:R ) {\n",
                        "    z[xx_s[k]] ~ dpois(0)\n",
                        "  }\n",sep="")
    }
    model.code<-paste(model.code,
                      "  for ( k in 1:R ) {\n",
                      "    for ( i in xx_s[k]:(xx_e[k]-1) ) {\n",
                      "      for( j in 1:Ne){\n",
                      "        y_e[j,i] <- y[j,i]\n",
                      "      }\n",
                      "      for ( j in 1:Nj ) {\n",
                      "        y_j[j,i] <- y[Ne+j,i]\n",
                      "      }\n",
                      "      for ( j in 1:Na ) {\n",
                      "        y_a[j,i] <- y[Ne+Nj+j,i]\n",
                      "      }\n",
                      sep="")
    
    if(FECUNDITY=="P"){
      if(ZERO_INFLATION=="T"){
        model.code<-paste(model.code,
                          "      z[i+1]   ~ dbin(1-pi,sum(y_a[,i]))\n",
                          "      y[1,i+1] ~ dpois((f/Ne)*z[i+1]*Dt)\n",sep="")
      }
      if(ZERO_INFLATION=="F"){
        model.code<-paste(model.code,
                          "      y[1,i+1] ~ dpois((f/Ne)*sum(y_a[,i])*Dt)\n",sep="")
      }
    }
    if(FECUNDITY=="N"){
      if(ZERO_INFLATION=="T"){
        model.code<-paste(model.code,
                          "      z[i+1]      ~ dbin(1-pi,sum(y_a[,i]))\n",
                          "      lambda[i+1] ~ dgamma(max(r*(ifelse(z[i+1]<0.5,1,z[i+1]))/Ne,0.001),max(r/(f*Dt),0.001))T(1e-6,)\n",
                          "      y[1,i+1]    ~ dpois(ifelse(z[i+1]==0,0,lambda[i+1]))\n",sep="")
      }
      if(ZERO_INFLATION=="F"){
        model.code<-paste(model.code,
                          "      lambda[i+1]  ~ dgamma(max(r*(ifelse(sum(y_a[,i])<0.5,1,sum(y_a[,i])))/Ne,0.001),max(r/(f*Dt),0.001)) T(1e-6,)\n",
                          "      y[1,i+1]     ~ dpois(ifelse(sum(y_a[,i])==0,0,lambda[i+1]))\n",sep="")
      }
    }
    model.code <- paste(model.code,
                        "      for(j in 1:(Ne-1)){\n",
                        "        y[j+1,i+1] ~ dbin(.99,y[j,i])\n",
                        "      }\n",
                        "      for(j in Ne:(Ne+Nj+Na-1)){\n",
                        "        y[j+1,i+1] ~ dbin( exp(-ifelse(\n",
                        "(H*Dt",ifelse(TIME_DEP=="T","+HH*(i-xx_s[k])*Dt",""),ifelse(CROWDING=="T","+HK*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),")>0.001,\n",
                        "(H*Dt",ifelse(TIME_DEP=="T","+HH*(i-xx_s[k])*Dt",""),ifelse(CROWDING=="T","+HK*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),"),\n",
                        "0.001*exp((H*Dt",ifelse(TIME_DEP=="T","+HH*(i-xx_s[k])*Dt",""),ifelse(CROWDING=="T","+HK*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),")/0.001-1))),\n",
                        "y[j,i])  \n",
                        "      }\n",
                        "    }\n",
                        
                        "    for( j in 1:Ne){\n",
                        "      y_e[j,xx_e[k]] <- y[j,xx_e[k]]\n",
                        "    }\n",
                        "    for ( j in 1:Nj ) {\n",
                        "      y_j[j,xx_e[k]] <- y[Ne+j,xx_e[k]]\n",
                        "    }\n",
                        "    for ( j in 1:Na ) {\n",
                        "      y_a[j,xx_e[k]] <- y[Ne+Nj+j,xx_e[k]]\n",
                        "    }\n",
                        "  }\n",
                        
                        "  for ( k in 1:R ) {\n",
                        "    for ( i in nn_s[k]:nn_e[k] ) {\n",
                        "      y_e_obs[i] ~ dnorm(sum(y_e[,xx_s[k]+x_obs[i]-1]),1)\n",
                        "      y_j_obs[i] ~ dnorm(sum(y_j[,xx_s[k]+x_obs[i]-1])*Dt,1)\n",
                        "      y_a_obs[i] ~ dnorm(sum(y_a[,xx_s[k]+x_obs[i]-1])*Dt,1)\n",
                        "    }\n",
                        "  }\n",
                        "}\n",
                        sep="")
    
    cat(model.code,file=file.jags.model)
    
    ######################################################################
    ### prepare lists containers and labels and set initial test
    ######################################################################
    
    filename <- paste("Test_Chains_",MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,".csv",sep="")
    treat    <- paste(MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,sep="")
    
    ### lists for jags initial test to optimize
    
    ddll <- list(x_end=x_end,R=R,Ne=Ne,Nj=Nj,Na=Na,Dt=Dt,x_obs=x_obs_to_fit,
                 y_e_obs=y_e_obs_to_fit,y_j_obs=y_j_obs_to_fit,y_a_obs=y_a_obs_to_fit,
                 nn_s=nn_s,nn_e=nn_e,xx_s=xx_s,xx_e=xx_e,
                 mlf=mlf,tauf=tauf,mlH=mlH,tauH=tauH)
    ddii <- list(f=mf,H=mH)
    ccss <- c("f","H")
    
    if(ZERO_INFLATION=="T"){
      ddll <- c(ddll,mpi=mpi,sdpi=sdpi)
      ddii <- c(ddii,pi=mpi)
      ccss <- c(ccss,"pi")
    }
    if(FECUNDITY=="N"){
      ddll <- c(ddll,mlr=mlr,taur=taur)
      ddii <- c(ddii,r=mr)
      ccss <- c(ccss,"r")
    }
    if(TIME_DEP =="T"){
      ddll <- c(ddll,mHH=mHH,sdHH=sdHH)
      ddii <- c(ddii,HH=mHH)
      ccss <- c(ccss,"HH")
    }
    if(CROWDING =="T"){
      ddll <- c(ddll,mlHK=mlHK,tauHK=tauHK)
      ddii <- c(ddii,HK=mHK)
      ccss <- c(ccss,"HK")
    }
    
    N_p <- length(ddii)
    
    ### interpolated initial conditions
    y_start <- array(rep(0,(Ne+Nj+Na)*sum(x_end)),dim=c(Ne+Nj+Na,sum(x_end)))
    ddii$y=matrix(0,nrow=Ne+Nj+Na,ncol=sum(x_end))
    
    if(INTERPOLATE == T){    
      
      parms <- c(
        #fecunditity
        f  = mf,
        r  = mr,
        pi = mpi,
        #death rate
        H  = mH,
        HH = mHH,
        HK = mHK
      )
      n <- max(x_end)
      #state vector (eggs + individuals per (juvenile and adult) age class) X timesteps X number of replicates
      Iy    <- array(0,dim=c(Na + Nj + Ne,n,R))
      
      ## set initial conditions for every replicate
      ICIN <- c(0,3,rep(0,Nj-1),rep(0,Na))
      for(k in 1:R) Iy[,1,k] <- model(FECUNDITY,TIME_DEP,CROWDING,ICIN,Dt,parms,1,Ne,Na,F)
      
      ## run the model
      for(k in 1:R){
        for(i in 1:(n-1)){
          Iy[,i+1,k] <- model(FECUNDITY,TIME_DEP,CROWDING,Iy[,i,k],Dt,parms,i,Ne,Na,F)
        }
      }
      
      IIy <- cbind(ICIN,Iy[,1:(x_end[1]-1),1])
      for(k in 2:R) IIy <- cbind(IIy,cbind(ICIN,Iy[,1:(x_end[k]-1),k]))
      
      y_start <- IIy
    }      
    if(MAKE_CONSISTENT==T){
      for ( k in 1:length(xx_e) ) {
        for ( j in (xx_s[k]+1):xx_e[k] ) {
          if(sum(ddii$y[(Nj+2):(Nj+Na+1),j-1])==0){
            if(ddii$y[1,j]>0) {
              print("inconsistency in egg numbers removed")
              print(j)
              ddii$y[1,j]<-0
            }}
        }}
      
      #ddii$y[1,]     <- 0
      ddii$y[,xx_s]  <- 0
      ddii$y[2,xx_s] <- 3
      for ( i in 2:nrow(ddii$y) ) {
        for ( k in 1:length(xx_e) ) {
          for ( j in (xx_s[k]+1):xx_e[k] ) {
            if(ddii$y[i,j]>ddii$y[i-1,j-1]){ 
              print("inconsistency in age structure removed")
              print(j)
            }
            ddii$y[i,j] <- min(ddii$y[i,j],ddii$y[i-1,j-1])
          }}}
    }
    
    if(FECUNDITY=="P")ccss_y <- c(ccss,"y")
    if(FECUNDITY=="N"){
      ccss_y <- c(ccss,"lambda","y")
      #ddii$lambda <- rep(1,sum(x_end))
    }
    if(ZERO_INFLATION=="T")ccss_y <- c(ccss_y,"z")
    ccss_y <- c(ccss_y,"deviance")
    
    ######################################################################
    ### Get a chain to optimize initial values
    ######################################################################
    if(INITIALIZE_CHAINS){
      jags.res <- jags.model(file.jags.model,data=ddll,inits=ddii,n.chains=1)
      
      ## message on the terminal
      print(paste("FITTING PARAMETERS for model ",MODEL_NAME," treatment ",tt,"_",cc,"_Test_chain",sep=""))
      
      ### inference
      jags.res2 <- coda.samples(jags.res,ccss_y,CHAIN_LENGTH,thin=THIN)
      
      write.csv(as.matrix(jags.res2),file=filename)
      
      ll <- LogPosterior("Test_Chains_",MODEL_NAME,tt,cc,PRIORS,Nj,Na,1,F,F,DEBUG=F)
      ind.max <- which(rowSums(ll)==max(rowSums(ll)))[1]
      
      print(paste("GETTING INITIAL PARAMETERS for model ",MODEL_NAME," treatment ",tt,"_",cc,sep=""))
      dd_t <- read.csv(paste("Test_Chains_",MODEL_NAME,"_",TREAT_TEST,"_",CLONE_TEST,"_",PRIORS,".csv",sep=""))
      
      ddii <- list(f=dd_t[,"f"][ind.max],H=dd_t[,"H"][ind.max])
      if(ZERO_INFLATION=="T")     ddii <- c(ddii,pi=dd_t[,"pi"][ind.max])
      if(FECUNDITY=="N")          ddii <- c(ddii,r =dd_t[,"r"][ind.max])
      if(TIME_DEP =="T")          ddii <- c(ddii,HH=dd_t[,"HH"][ind.max])
      if(CROWDING =="T")          ddii <- c(ddii,HK=dd_t[,"HK"][ind.max])
      
      
      ddii$y <- matrix(NA,ncol=sum(x_end),nrow=Na+Nj+1)
      for(j in 1:sum(x_end)){for(i in 1:(Na+Nj+1)){
        ddii$y[i,j] <-dd_t[ind.max,paste("y.",i,".",j,".",sep="")]
      }}
      
      if(FECUNDITY=="N"){
        ddii$lambda  <- rep(NA,sum(x_end))
        #ddii$m      <- rep(1,sum(x_end))
        for(i in 1:sum(x_end))ddii$lambda[i]  <-dd_t[ind.max,paste("lambda.",i,".",sep="")]
        #for(i in 1:sum(x_end))ddii$m[i]      <-dd_t[ind.max,paste("m.",i,".",sep="")]
      }
      if(ZERO_INFLATION=="T"){
        ddii$z <-rep(NA,sum(x_end))
        for(i in 1:sum(x_end))ddii$z[i] <-dd_t[ind.max,paste("z.",i,".",sep="")]
      }
    }
    ######################################################################
    ### initialization, fit of the model and plot of results
    ######################################################################
    
    filename <- paste("Opt_Chains_",MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,".csv",sep="")
    
    jags.res <- jags.model(file.jags.model,data=ddll,inits=ddii,
                           #n.adapt=0,
                           n.chains=N_CHAINS)
    
    ## message on the terminal
    print(paste("FITTING PARAMETERS for model ",MODEL_NAME," treatment ",tt,"_",cc,"_",N_CHAINS,"_chains",sep=""))
    
    ### inference
    jags.res2 <- coda.samples(jags.res,ccss_y,CHAIN_LENGTH,thin=THIN)
    
    
    ### DEBUG pdf("TEST_1.pdf")
    ### plotting
    plot(jags.res2[,ccss],ylab=treat)
    
    dev <- jags.res2[[1]][,"deviance"]
    if(N_CHAINS>1) for(nc in 2:N_CHAINS) dev <- c(dev,jags.res2[[nc]][,"deviance"])
    
    PV_1  <- var(dev)/2
    DIC_1 <- mean(dev)+PV_1
    ### use the function with lik=T to get the LogLikelihood computed at the mean parameters (states) 
    ### PV_2  <- -2*LogPosterior(MODEL_NAME,tt,cc,pp,5,45,1,T,T,F)
    ### DIC_2 <- mean(dev)+PV_2
    
    plot(jags.res2[,"deviance"],ylab=treat,
         main=paste("Deviance (DIC ",signif(DIC_1,digits=4),")",sep=""))
    
    ## compute and plot correlation between parameters
    if(PLOT_CORR==T){
      chains <- as.matrix(jags.res2[,ccss])
      thin <- THIN_PLOT
      sample <- seq(length(chains[,1])/2,length(chains[,1]),by=thin)
      par(mfrow=c(1,1))
      this <- paste(MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,sep="")
      ipairs(chains[sample,], main=this)
      cor(chains)
    }
    
    ## save chains
    write.csv(as.matrix(jags.res2),file=filename)
    
    if(PLOT_POSTERIOR){
      llpp <- LogPosterior("Opt_Chains_",MODEL_NAME,tt,cc,1,Ne,Nj,Na,Dt,F,F,F)
      par(mfrow=c(2,1))
      plot(rowSums(llpp),type="l",ylab="logPosterior")
      plot(-2*(llpp[,"LogLik"]),type="l",ylab="Deviance")
    }
    ######################################################################
  }### close INF_METH "J"
  
  if(INF_METH=="R"){
    ############################################################################################################
    ### containers
    ############################################################################################################
    
    dev <- array(0,dim = c(ceiling(N_CHAINS*CHAIN_LENGTH/THIN),R))
    PV_1  <- array(0,dim = c(R))
    DIC_1 <- array(0,dim = c(R))
    PV_2  <- array(0,dim = c(R))
    DIC_2 <- array(0,dim = c(R))
    
    kk<-1
    ############################################################################################################
    ### replicate loop
    ############################################################################################################
    
    for(kk in 1:R){
      x_obs_to_fit_rep   <- x_obs_to_fit[nn_s[kk]:nn_e[kk]] 
      y_e_obs_to_fit_rep <- y_e_obs_to_fit[nn_s[kk]:nn_e[kk]] 
      y_j_obs_to_fit_rep <- y_j_obs_to_fit[nn_s[kk]:nn_e[kk]] 
      y_a_obs_to_fit_rep <- y_a_obs_to_fit[nn_s[kk]:nn_e[kk]] 
      
      #####################################################################
      # jags function
      #####################################################################
      
      model.code <- paste("var y_e[Ne,x_end[k]],y_j[Nj,x_end[k]],y_a[Na,x_end[k]]",
                          ifelse(ZERO_INFLATION=="T",",z[sum(x_end[k])]",""),
                          ifelse(FECUNDITY=="N",",lambda[sum(x_end[k])]\n","\n"),
                          "model {\n",
                          "  f  ~ ",ifelse(PRIORS==1,"dlnorm(mlf,tauf)I(0,30)\n","dunif(0,30)\n"),
                          "  H  ~ ",ifelse(PRIORS==1,"dlnorm(mlH,tauH)\n","dunif(0,1)\n"),sep="")
      if(ZERO_INFLATION=="T") model.code <- paste(model.code,"  pi ~ ",ifelse(PRIORS==1,"dnorm(mpi,(1/sdpi)^2)I(0,1)\n","dunif(0,1)\n"),sep="")
      if(FECUNDITY=="N")      model.code <- paste(model.code,"  r  ~ ",ifelse(PRIORS==1,"dlnorm(mlr,taur)\n","dunif(0,100)\n"),sep="")
      if(TIME_DEP=="T")       model.code <- paste(model.code,"  HH ~ ",ifelse(PRIORS==1,"dnorm(mHH,(1/sdHH)^2)\n","dunif(0,1)\n"),sep="")
      if(CROWDING=="T")       model.code <- paste(model.code,"  HK ~ ",ifelse(PRIORS==1,"dlnorm(mlHK,tauHK)\n","dunif(0,1)\n"),sep="")
      model.code <- paste(model.code,
                          "    for ( j in 1:Ne ) {\n",
                          "      y[j,1] ~ dpois(0)\n",
                          "    }\n",
                          "    y[Ne+1,1] ~ dbin(.99,3)\n",
                          "    for ( j in (Ne+2):(Ne+Nj+Na) ) {\n",
                          "      y[j,1] ~ dpois(0)\n",
                          "    }\n",sep="")
      
      if(FECUNDITY=="N"){
        model.code<-paste(model.code,
                          "    lambda[1] ~ dpois(0)\n",sep="")
      }
      if(ZERO_INFLATION=="T"){
        model.code<-paste(model.code,
                          "    z[1] ~ dpois(0)\n",sep="")
      }
      
      model.code<-paste(model.code, 
                        "    for ( i in 1:(x_end[k]-1) ) {\n",
                        "      for( j in 1:Ne){\n",
                        "        y_e[j,i] <- y[j,i]\n",
                        "      }\n",
                        "      for ( j in 1:Nj ) {\n",
                        "        y_j[j,i] <- y[Ne+j,i]\n",
                        "      }\n",
                        "      for ( j in 1:Na ) {\n",
                        "        y_a[j,i] <- y[Ne+Nj+j,i]\n",
                        "      }\n",sep="")
      if(FECUNDITY=="P"){
        if(ZERO_INFLATION=="T"){
          model.code<-paste(model.code,
                            "z[i+1]   ~ dbin(1-pi,sum(y_a[,i]))\n",
                            "y[1,i+1] ~ dpois((f/Ne)*z[i+1]*Dt)\n",sep="")
        }
        if(ZERO_INFLATION=="F"){
          model.code<-paste(model.code,
                            "y[1,i+1] ~ dpois((f/Ne)*sum(y_a[,i])*Dt)\n",sep="")
        }
      }
      if(FECUNDITY=="N"){
        if(ZERO_INFLATION=="T"){
          model.code<-paste(model.code,
                            "z[i+1]      ~ dbin(1-pi,sum(y_a[,i]))\n",
                            "lambda[i+1] ~ dgamma(max(r*(ifelse(z[i+1]<0.5,1,z[i+1]))/Ne,0.001),max(r/(f*Dt),0.001))T(1e-6,)\n",
                            "y[1,i+1]    ~ dpois(ifelse(z[i+1]==0,0,lambda[i+1]))\n",sep="")
        }
        if(ZERO_INFLATION=="F"){
          model.code<-paste(model.code,
                            "lambda[i+1]  ~ dgamma(max(r*(ifelse(sum(y_a[,i])<0.5,1,sum(y_a[,i])))/Ne,0.001),max(r/(f*Dt),0.001)) T(1e-6,)\n",
                            "y[1,i+1]     ~ dpois(ifelse(sum(y_a[,i])==0,0,lambda[i+1]))\n",sep="")
        }
      }
      model.code <- paste(model.code,
                          "      for(j in 1:(Ne-1)){\n",
                          "        y[j+1,i+1] ~ dbin(.99,y[j,i])\n",
                          "      }\n",
                          
                          "      for(j in Ne:(Ne+Nj+Na-1)){\n",
                          "        y[j+1,i+1] ~ dbin( exp(-ifelse(\n",
                          
                          "(H*Dt",ifelse(TIME_DEP=="T","+HH*i*Dt",""),ifelse(CROWDING=="T","+HK*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),")>0.001,\n",
                          "(H*Dt",ifelse(TIME_DEP=="T","+HH*i*Dt",""),ifelse(CROWDING=="T","+HK*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),"),\n",
                          "0.001*exp((H*Dt",ifelse(TIME_DEP=="T","+HH*i*Dt",""),ifelse(CROWDING=="T","+HK*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),")/0.001-1))),\n",
                          "    y[j,i])\n",
                          "      }\n",
                          "    }\n",
                          "    for( j in 1:Ne){\n",
                          "      y_e[j,x_end[k]] <- y[j,x_end[k]]\n",
                          "    }\n",
                          "    for ( j in 1:Nj ) {\n",
                          "      y_j[j,x_end[k]] <- y[Ne+j,x_end[k]]\n",
                          "    }\n",
                          "    for ( j in 1:Na ) {\n",
                          "      y_a[j,x_end[k]] <- y[Ne+Nj+j,x_end[k]]\n",
                          "    }\n",
                          "    for ( i in 1:n_obs ) {\n",
                          "      y_e_obs[i] ~ dnorm(sum(y_e[,x_obs[i]]),1)\n",
                          "      y_j_obs[i] ~ dnorm(sum(y_j[,x_obs[i]])*Dt,1)\n",
                          "      y_a_obs[i] ~ dnorm(sum(y_a[,x_obs[i]])*Dt,1)\n",
                          "    }\n",
                          "}\n",
                          sep="")
      
      cat(model.code,file=file.jags.model)
      
      ######################################################################
      ### prepare lists containers and labels
      ######################################################################
      
      filename <- paste("Test_Chains_",MODEL_NAME,"_",tt,"_",cc,"_REP_",kk,"_",PRIORS,".csv",sep="")
      treat    <- paste(MODEL_NAME,"_",tt,"_",cc,"_REP_",kk,"_",PRIORS,sep="")
      
      ### lists for jags
      
      ddll <- list(x_end=x_end,k=kk,Ne=Ne,Nj=Nj,Na=Na,Dt=Dt,x_obs=x_obs_to_fit_rep,
                   y_e_obs=y_e_obs_to_fit_rep,y_j_obs=y_j_obs_to_fit_rep,y_a_obs=y_a_obs_to_fit_rep,
                   n_obs=n_obs_to_fit[kk],
                   mlf=mlf,tauf=tauf,mlH=mlH,tauH=tauH)
      
      ddii <- list(f=mf,H=mH)
      ccss <- c("f","H")
      if(ZERO_INFLATION=="T"){
        ddll <- c(ddll,mpi=mpi,sdpi=sdpi)
        ddii <- c(ddii,pi =mpi)  
        ccss <- c(ccss,"pi")  
      }
      
      if(FECUNDITY=="N"){
        ddll <- c(ddll,mlr=mlr,taur=taur)
        ddii <- c(ddii,r=mr)
        ccss <- c(ccss,"r")
      }
      if(TIME_DEP =="T"){
        ddll <- c(ddll,mHH=mHH,sdHH=sdHH)
        ddii <- c(ddii,HH=mHH)
        ccss <- c(ccss,"HH")
      }
      if(CROWDING =="T"){
        ddll <- c(ddll,mlHK=mlHK,tauHK=tauHK)
        ddii <- c(ddii,HK=mHK)
        ccss <- c(ccss,"HK")
      }
      
      y_start <- array(rep(0,(Ne+Nj+Na)*x_end[kk]),dim=c(Ne+Nj+Na,x_end[kk]))
      ddii$y=matrix(0,nrow=Ne+Nj+Na,ncol=x_end[kk])
      
      if(INTERPOLATE == T){    
        
        parms <- c(
          #fecunditity
          f  = mf,
          pi = mpi,
          r  = mr,
          
          #death rate
          H = mH,
          HH = mHH,
          HK = mHK
        )
        n <- x_end[kk]
        #state vector (eggs + individuals per (juvenile and adult) age class) X timesteps X number of replicates
        Iy    <- array(0,dim=c(Na + Nj + 1,n))
        
        ## set initial conditions for every replicate
        ICIN <- c(0,3,rep(0,Nj-1),rep(0,Na))
        Iy[,1] <- model(FECUNDITY,TIME_DEP,CROWDING,ICIN,Dt,parms,1,Ne,Na,F)
        
        ## run the model
        i <- 1
        for(i in 1:(n-1)){
          #Iy[,i+1] <- model(FECUNDITY,TIME_DEP,CROWDING,Iy[,i],Dt,i,parms,Na,F)
          Iy[,i+1] <- c(Fe(FECUNDITY,Iy[,i],Dt,parms,Ne,Na,F),Su(TIME_DEP,CROWDING,Iy[,i],Dt,parms,i,Ne,Na,F))
          
        }
        
        IIy <- cbind(ICIN,Iy[,1:(n-1)])
        y_start <- IIy
        
      }      
      if(MAKE_CONSISTENT==T){
        ## first loop
        for ( j in 2:x_end[kk] ) {
          if(sum(ddii$y[(Nj+2):(Nj+Na+1),j-1])==0){
            if(ddii$y[1,j]>0) {
              print("inconsistency in egg numbers removed")
              print(j)
              ddii$y[1,j]<-0
            }}}
        
        # make y_start consistent with parents:
        #ddii$y[1,]     <- 0
        ddii$y[,1]  <- 0
        ddii$y[2,1] <- 3
        for ( i in 2:nrow(ddii$y) ) {
          for ( j in 2:x_end[kk] ) {
            if(ddii$y[i,j]>ddii$y[i-1,j-1]){ 
              print("inconsistency in age structure removed")
              print(j)
            }
            ddii$y[i,j] <- min(ddii$y[i,j],ddii$y[i-1,j-1])
          }}
      }
      
      if(FECUNDITY=="P")ccss_y <- c(ccss,"y")
      if(FECUNDITY=="N"){
        ccss_y <- c(ccss,"lambda","y")
        #ddii$lambda <- rep(1,x_end[kk])
      }
      if(ZERO_INFLATION=="T")ccss_y <- c(ccss_y,"z")
      
      ccss_y <- c(ccss_y,"deviance")
      
      ######################################################################
      ### Get a chain to optimize initial values
      ######################################################################
      if(INITIALIZE_CHAINS){
        jags.res <- jags.model(file.jags.model,data=ddll,inits=ddii,n.chains=1)
        
        ## message on the terminal
        print(paste("FITTING PARAMETERS for model ",MODEL_NAME," treatment ",tt,"_",cc,"_Replicate_",kk,"_Test_chain",sep=""))
        
        ### inference
        jags.res2 <- coda.samples(jags.res,ccss_y,CHAIN_LENGTH,thin=THIN)
        
        write.csv(as.matrix(jags.res2),file=filename)
        ll <- LogPosterior("Test_Chains_",MODEL_NAME,tt,cc,PRIORS,5,45,1,F,kk,F)
        ind.max <- which(rowSums(ll)==max(rowSums(ll)))[1]
        
        dd_t <- read.csv(paste("Test_Chains_",MODEL_NAME,"_",TREAT_TEST,"_",CLONE_TEST,"_REP_",kk,"_",PRIORS,".csv",sep=""))
        
        print(paste("GETTING INITIAL PARAMETERS for model ",MODEL_NAME," treatment ",tt,"_",cc,sep=""))
        
        ddii <- list()
        ddii <- list(f=dd_t[,"f"][ind.max],H=dd_t[,"H"][ind.max])
        if(ZERO_INFLATION)     ddii <- c(ddii,pi=dd_t[,"pi"][ind.max])
        if(FECUNDITY=="N")     ddii <- c(ddii,r =dd_t[,"r"][ind.max])
        if(TIME_DEP =="T")     ddii <- c(ddii,HH=dd_t[,"HH"][ind.max])
        if(CROWDING =="T")     ddii <- c(ddii,HK=dd_t[,"HK"][ind.max])
        
        ddii$y <- matrix(NA,ncol=x_end[kk],nrow=Na+Nj+1)
        for(j in 1:x_end[kk]){for(i in 1:(Na+Nj+1)){
          ddii$y[i,j] <-dd_t[ind.max,paste("y.",i,".",j,".",sep="")]
        }}
        
        if(FECUNDITY=="N"){
          ddii$lambda <- rep(1,x_end[kk])
          #ddii$m      <- rep(1,x_end[kk])
          for(i in 1:x_end[kk])ddii$lambda[i] <-dd_t[ind.max,paste("lambda.",i,".",sep="")]
          #for(i in 1:x_end[kk])ddii$m[i]      <-dd_t[ind.max,paste("m.",i,".",sep="")]
        }
        if(ZERO_INFLATION=="T"){
          ddii$z <-rep(NA,x_end[k])
          for(i in 1:x_end[k])ddii$z[i] <-dd_t[ind.max,paste("z.",i,".",sep="")]
        }
        
      }
      ######################################################################
      ### initialization, fit of the model and plot of results
      ######################################################################
      
      filename <- paste("Opt_Chains_",MODEL_NAME,"_",tt,"_",cc,"_REP_",kk,"_",PRIORS,".csv",sep="")
      
      jags.res <- jags.model(file.jags.model,
                             data=ddll,
                             inits=ddii,
                             #n.adapt=0,
                             n.chains=N_CHAINS)
      
      ## message on the terminal
      print(paste("FITTING PARAMETERS for model_",MODEL_NAME,"_treatment ",tt,"_",cc,"_Replicate_",kk,"_",N_CHAINS,"chains",sep=""))
      
      ### inference
      jags.res2 <- coda.samples(jags.res,ccss_y,CHAIN_LENGTH,thin=THIN)
      
      ### plotting
      plot(jags.res2[,ccss],ylab=treat)
      
      dev_t <- jags.res2[[1]][,"deviance"]
      if(N_CHAINS>1) for(nc in 2:N_CHAINS)  dev_t <- c(dev_t,jags.res2[[nc]][,"deviance"])
      dev[,kk] <- dev_t     
      
      PV_1[kk]  <- var(dev[,kk])/2
      DIC_1[kk] <- mean(dev[,kk])+PV_1[kk]
      # PV_2[kk]  <- -2*LogPosterior(MODEL_NAME,tt,cc,pp,Ne,Nj,Na,1,T,kk,F)
      # DIC_2[kk] <- mean(dev[,kk])+PV_2[kk]
      
      plot(jags.res2[,"deviance"],ylab=treat,
           main=paste("Deviance (DIC ",signif(DIC_1[kk],digits=4),")",sep=""))
      
      if(PLOT_CORR==T){
        chains <- as.matrix(jags.res2[,ccss])
        thin <- THIN_PLOT
        sample <- seq(length(chains[,1])/2,length(chains[,1]),by=thin)
        par(mfrow=c(1,1))
        this <- paste(MODEL_NAME,"_",tt,"_",cc,"_REP_",kk,"_",PRIORS,sep="")
        ipairs(chains[sample,], main=this)
        cor(chains)
      }
      
      ## save chain
      write.csv(as.matrix(jags.res2),file= filename)
      
      if(PLOT_POSTERIOR){
        llpp <- LogPosterior("Opt_Chains_",MODEL_NAME,tt,cc,1,Ne,Nj,Na,Dt,F,kk,F)
        par(mfrow=c(2,1))
        plot(rowSums(llpp),type="l",ylab="logPosterior",main=paste("REP_",kk,sep=""))
        plot(-2*(llpp[,"LogLik"]),type="l",ylab="Deviance",main=paste("REP_",kk,sep=""))
      }
    }# close loop in replicates
    
    ######################################################################
    ### MODEL SELECTION for all replicates
    ######################################################################
    if(PLOT_POSTERIOR){
      print("Summing up DICs and log posteriors")
      ### sum deviances across replicates
      treat    <- paste(MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,sep="")
      llpp <- LogPosterior("Opt_Chains_",MODEL_NAME,tt,cc,1,Ne,Nj,Na,Dt,F,T,F)
      
      dev_tot <- rowSums(dev[,1:R])
      PV_tot  <- var(dev_tot)/2
      DIC_tot <- mean(dev_tot)+PV_tot
      DIC_tot <- signif(DIC_tot,digits=4)
      
      par(mfrow=c(2,1))
      plot(rowSums(llpp),type="l",ylab="logPosterior",
           main=paste("Sum of logposteriors Across Replicates",sep=""))
      plot(dev_tot,type="l",ylab="Deviance", 
           main=paste("Sum of Deviances Across Replicates",sep=""))
    }
  }### close "R" METHOD
  
  if(INF_METH=="H"){
    
    #####################################################################
    # jags function
    #####################################################################
    
    model.code <- paste("var y_e[Ne,sum(x_end)],y_j[Nj,sum(x_end)],y_a[Na,sum(x_end)]",
                        ifelse(ZERO_INFLATION=="T",",z[sum(x_end)]\n","\n"),
                        ifelse(FECUNDITY=="N",",lambda[sum(x_end)]\n","\n"),
                        "model {\n",
                        "  m_f   ~ dlnorm(ml_mf,tau_mf)I(0,30)\n",
                        "  sd_f  ~ dlnorm(ml_sdf,tau_sdf)\n",
                        "  m_H   ~ dlnorm(ml_mH,tau_mH)\n",
                        "  sd_H  ~ dlnorm(ml_sdH,tau_sdH)\n",
                        sep="")
    if(ZERO_INFLATION=="T"){
      model.code <- paste(model.code,                                    # pi ~ N(0.5,0.15)
                          "  m_pi   ~ dnorm(m_mpi,(1/sd_mpi)^2)I(0,1)\n",# N(0.5,0.15)
                          "  sd_pi  ~ dlnorm(ml_sdpi,tau_sdpi)\n",       #LN(0.15,0.075)
                          sep="")
    }
    if(FECUNDITY=="N"){
      model.code <- paste(model.code,                          #r ~ LN(1,0.25)
                          "  m_r   ~ dlnorm(ml_mr,tau_mr)\n",  #LN(1,0.25)
                          "  sd_r  ~ dlnorm(ml_sdr,tau_sdr)\n",#LN(0.25,0.125)
                          sep="")
    }
    if(TIME_DEP=="T"){
      model.code <- paste(model.code,                              #HH ~ N(0,1)
                          "  m_HH   ~ dnorm(m_mHH,(1/sd_mHH)^2)\n",#N(0,1)
                          "  sd_HH  ~ dlnorm(ml_sdHH,tau_sdHH)\n", #LN(1,0.5)
                          sep="")
    }
    if(CROWDING=="T"){
      model.code <- paste(model.code,                             ##HK ~LN(0.0001,0.0001)
                          "  m_HK   ~ dlnorm(ml_mHK,tau_mHK)\n",  #LN(0.0001,0.0001)
                          "  sd_HK  ~ dlnorm(ml_sdHK,tau_sdHK)\n",#LN(0.0001,0.0001/2)
                          sep="")
    }
    
    model.code <- paste(model.code,          
                        "  for ( ii in 1:R ) {\n",
                        "    f[ii]  ~ ",ifelse(PRIORS==1,"dlnorm( log(m_f)-log(1+sd_f^2/m_f^2)/2, 1/log(1+(sd_f/m_f)^2))T(0,30)\n","dunif(0,10)\n"),
                        "    H[ii]  ~ ",ifelse(PRIORS==1,"dlnorm( log(m_H)-log(1+sd_H^2/m_H^2)/2, 1/log(1+(sd_H/m_H)^2))\n","dunif(0,1)\n"),
                        sep="")
    if(ZERO_INFLATION=="T"){
      model.code <- paste(model.code,
                          "    pi[ii] ~ ",ifelse(PRIORS==1,"dnorm( m_pi, (1/sd_pi)^2)T(0,1)\n","dunif(0,1)\n"),
                          sep="")
    }
    if(FECUNDITY=="N"){
      model.code <- paste(model.code,
                          "    r[ii] ~ ",ifelse(PRIORS==1,"dlnorm( log(m_r)-log(1+sd_r^2/m_r^2)/2, 1/log(1+(sd_r/m_r)^2))\n","dunif(0,100)\n"),
                          sep="")                        
    }
    if(TIME_DEP=="T"){
      model.code <- paste(model.code,
                          "    HH[ii] ~ ",ifelse(PRIORS==1,"dnorm( m_HH, (1/sd_HH)^2)\n","dunif(0,1)\n"),
                          sep="")
    }
    if(CROWDING=="T"){
      model.code <- paste(model.code,
                          "    HK[ii] ~ ",ifelse(PRIORS==1,"dlnorm( log(m_HK)-log(1+sd_HK^2/m_HK^2)/2, 1/log(1+(sd_HK/m_HK)^2))\n","dunif(0,1)\n"),
                          sep="")
    }
    model.code <- paste(model.code,
                        "  }\n",
                        "  for ( k in 1:R ) {\n",
                        "    for ( j in 1:Ne ) {\n",
                        "      y[j,xx_s[k]] ~ dpois(0)\n",
                        "    }\n",
                        "    y[Ne+1,xx_s[k]] ~ dbin(1,3)\n",
                        "    for ( j in (Ne+2):(Ne+Nj+Na) ) {\n",
                        "      y[j,xx_s[k]] ~ dpois(0)\n",
                        "    }\n",
                        "  }\n",
                        sep="")
                        if(FECUNDITY=="N"){
                          model.code<-paste(model.code,
                                            "  for ( k in 1:R ) {\n",
                                            "    lambda[xx_s[k]] ~ dpois(0)\n",
                                            "  }\n",sep="")
                        }
                        if(ZERO_INFLATION=="T"){
                          model.code<-paste(model.code,
                                            "  for ( k in 1:R ) {\n",
                                            "    z[xx_s[k]] ~ dpois(0)\n",
                                            "  }\n",sep="")
                        }
    
                        model.code<-paste(model.code,
                        "  for ( k in 1:R ) {\n",
                        "    for ( i in xx_s[k]:(xx_e[k]-1) ) {\n",
                        "      for( j in 1:Ne){\n",
                        "        y_e[j,i] <- y[j,i]\n",
                        "      }\n",
                        "      for ( j in 1:Nj ) {\n",
                        "        y_j[j,i] <- y[Ne+j,i]\n",
                        "      }\n",
                        "      for ( j in 1:Na ) {\n",
                        "        y_a[j,i] <- y[Ne+Nj+j,i]\n",
                        "      }\n",
                        sep="")

    if(FECUNDITY=="P"){
      if(ZERO_INFLATION=="T"){
        model.code<-paste(model.code,
                          "z[i+1]   ~ dbin(1-pi[k],sum(y_a[,i]))\n",
                          "y[1,i+1] ~ dpois((f[k]/Ne)*z[i+1]*Dt)\n",sep="")
      }
      if(ZERO_INFLATION=="F"){
        model.code<-paste(model.code,
                          "y[1,i+1] ~ dpois((f[k]/Ne)*sum(y_a[,i])*Dt)\n",sep="")
      }
    }
    if(FECUNDITY=="N"){
      if(ZERO_INFLATION=="T"){
        model.code<-paste(model.code,
                          "z[i+1]      ~ dbin(1-pi[k],sum(y_a[,i]))\n",
                          "lambda[i+1] ~ dgamma(max(r[k]*(ifelse(z[i+1]<0.5,1,z[i+1]))/Ne,0.001),max(r[k]/(f[k]*Dt),0.001))T(1e-6,)\n",
                          "y[1,i+1]    ~ dpois(ifelse(z[i+1]==0,0,lambda[i+1]))\n",sep="")
      }
      if(ZERO_INFLATION=="F"){
        model.code<-paste(model.code,
                          "lambda[i+1]  ~ dgamma(max(r[k]*(ifelse(sum(y_a[,i])<0.5,1,sum(y_a[,i])))/Ne,0.001),max(r[k]/(f[k]*Dt),0.001)) T(1e-6,)\n",
                          "y[1,i+1]     ~ dpois(ifelse(sum(y_a[,i])==0,0,lambda[i+1]))\n",sep="")
      }
    }
                        
    model.code <- paste(model.code,
                        "      for(j in 1:(Ne-1)){\n",
                        "        y[j+1,i+1] ~ dbin(.99,y[j,i])\n",
                        "      }\n",
                        "      for(j in Ne:(Ne+Nj+Na-1)){\n",
                        "        y[j+1,i+1] ~ dbin( exp(-ifelse(\n",
                        "(H[k]*Dt",ifelse(TIME_DEP=="T","+HH[k]*(i-xx_s[k])*Dt",""),ifelse(CROWDING=="T","+HK[k]*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),")>0.001,\n",
                        "(H[k]*Dt",ifelse(TIME_DEP=="T","+HH[k]*(i-xx_s[k])*Dt",""),ifelse(CROWDING=="T","+HK[k]*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),"),\n",
                        "0.001*exp((H[k]*Dt",ifelse(TIME_DEP=="T","+HH[k]*(i-xx_s[k])*Dt",""),ifelse(CROWDING=="T","+HK[k]*(sum(y_j[,i])+sum(y_a[,i]))*Dt^2","+(sum(y_j[,i])+sum(y_a[,i]))*Dt^2/1000000000"),")/0.001-1))),\n",
                        "y[j,i])  \n",
                        "      }\n",
                        "    }\n",
                        
                        "    for( j in 1:Ne){\n",
                        "      y_e[j,xx_e[k]] <- y[j,xx_e[k]]\n",
                        "    }\n",
                        "    for ( j in 1:Nj ) {\n",
                        "      y_j[j,xx_e[k]] <- y[Ne+j,xx_e[k]]\n",
                        "    }\n",
                        "    for ( j in 1:Na ) {\n",
                        "      y_a[j,xx_e[k]] <- y[Ne+Nj+j,xx_e[k]]\n",
                        "    }\n",
                        "  }\n",
                        
                        "  for ( k in 1:R ) {\n",
                        "    for ( i in nn_s[k]:nn_e[k] ) {\n",
                        "      y_e_obs[i] ~ dnorm(sum(y_e[,xx_s[k]+x_obs[i]-1]),1)\n",
                        "      y_j_obs[i] ~ dnorm(sum(y_j[,xx_s[k]+x_obs[i]-1])*Dt,1)\n",
                        "      y_a_obs[i] ~ dnorm(sum(y_a[,xx_s[k]+x_obs[i]-1])*Dt,1)\n",
                        "    }\n",
                        "  }\n",
                        "}\n",
                        sep="")
    
    cat(model.code,file=file.jags.model)
    
    ######################################################################
    ### prepare containers and labels
    ######################################################################
    
    filename <- paste("Test_Chains_",MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,".csv",sep="")
    treat    <- paste(MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,sep="")
    
    ### lists for jags
    ddll <- list(x_end=x_end,R=R,Ne=Ne,Nj=Nj,Na=Na,Dt=Dt,x_obs=x_obs_to_fit,
                 y_e_obs=y_e_obs_to_fit,y_j_obs=y_j_obs_to_fit,y_a_obs=y_a_obs_to_fit,
                 nn_s=nn_s,nn_e=nn_e,xx_s=xx_s,xx_e=xx_e,
                 
                 ml_mf=ml_mf,tau_mf=tau_mf,
                 ml_sdf=ml_sdf,tau_sdf=tau_sdf,
                 ml_mH=ml_mH,tau_mH=tau_mH,
                 ml_sdH=ml_sdH,tau_sdH=tau_sdH)
    
    ddii    <- list(m_f=m_mf,sd_f=m_sdf,m_H=m_mH,sd_H=m_sdH)
    ccss_H  <- c("m_f","sd_f","m_H","sd_H")
    ccss    <- c("f","H")
    if(ZERO_INFLATION=="T"){
      ddll   <- c(ddll,m_mpi=m_mpi,sd_mpi=sd_mpi,
                  ml_sdpi=ml_sdpi,tau_sdpi=tau_sdpi)
      ddii   <- c(ddii,m_pi=m_mpi,sd_pi=m_sdpi)
      ccss_H <- c(ccss_H,"m_pi","sd_pi")
      ccss   <- c(ccss,"pi")
    }
    if(FECUNDITY=="N"){
      ddll   <- c(ddll,ml_mr=ml_mr,tau_mr=tau_mr,
                  ml_sdr=ml_sdr,tau_sdr=tau_sdr)
      ddii   <- c(ddii,m_r=m_mr,sd_r=m_sdr)
      ccss_H <- c(ccss_H,"m_r","sd_r")
      ccss   <- c(ccss,"r")
    }
    if(TIME_DEP =="T"){
      ddll   <- c(ddll,m_mHH=m_mHH,sd_mHH=sd_mHH,
                  ml_sdHH=ml_sdHH,tau_sdHH=tau_sdHH)
      ddii   <- c(ddii,m_HH=m_mHH,sd_HH=m_sdHH)
      ccss_H <- c(ccss_H,"m_HH","sd_HH")
      ccss   <- c(ccss,"HH")
    }
    if(CROWDING =="T"){
      ddll   <- c(ddll,ml_mHK=ml_mHK,tau_mHK=tau_mHK,
                  ml_sdHK=ml_sdHK,tau_sdHK=tau_sdHK)
      ddii   <- c(ddii,m_HK=m_mHK,sd_HK=m_sdHK)
      ccss_H <- c(ccss_H,"m_HK","sd_HK")    
      ccss   <- c(ccss,"HK")    
    }
    
    ### interpolated initial conditions
    y_start <- array(rep(0,(Ne+Nj+Na)*sum(x_end)),dim=c(Ne+Nj+Na,sum(x_end)))
    ddii$y=matrix(0,nrow=Ne+Nj+Na,ncol=sum(x_end))
    
    if(INTERPOLATE == T){    
      
      parms <- c(
        #death rate
        H  = m_mH,
        HH = m_mHH,
        #fecunditity
        f  = m_mf,
        r  = m_mr,
        pi = m_pi,
        HK = m_mHK
      )
      n <- max(x_end)
      #state vector (eggs + individuals per (juvenile and adult) age class) X timesteps X number of replicates
      Iy    <- array(0,dim=c(Na + Nj + Ne,n,R))
      
      ## set initial conditions for every replicate
      ICIN <- c(0,3,rep(0,Nj-1),rep(0,Na))
      for(k in 1:R) Iy[,1,k] <- model(FECUNDITY,TIME_DEP,CROWDING,ICIN,Dt,parms,1,Ne,Na,F)
      
      ## run the model
      for(k in 1:R){
        for(i in 1:(n-1)){
          Iy[,i+1,k] <- model(FECUNDITY,TIME_DEP,CROWDING,Iy[,i,k],Dt,parms,i,Ne,Na,F)
        }
      }
      
      IIy <- cbind(ICIN,Iy[,1:(x_end[1]-1),1])
      for(k in 2:R) IIy <- cbind(IIy,cbind(ICIN,Iy[,1:(x_end[k]-1),k]))
      
      y_start <- IIy
    }      
    if(MAKE_CONSISTENT==T){
      ## first loop
      for ( k in 1:length(xx_e) ) {
        for ( j in (xx_s[k]+1):xx_e[k] ) {
          if(sum(ddii$y[(Nj+2):(Nj+Na+1),j-1])==0){
            if(ddii$y[1,j]>0) {
              print("inconsistency in egg numbers removed")
              print(j)
              ddii$y[1,j]<-0
            }}}}
      
      #ddii$y[1,]     <- 0
      ddii$y[,xx_s]  <- 0
      ddii$y[2,xx_s] <- 3
      for ( i in 2:nrow(ddii$y) ) {
        for ( k in 1:length(xx_e) ) {
          for ( j in (xx_s[k]+1):xx_e[k] ) {
            if(ddii$y[i,j]>ddii$y[i-1,j-1]){ 
              print("inconsistency in age structure removed")
              print(j)
            }
            ddii$y[i,j] <- min(ddii$y[i,j],ddii$y[i-1,j-1])
          }}}
    }
    
    if(FECUNDITY=="P")ccss_y <- c(ccss_H,ccss,"y")
    if(FECUNDITY=="N")ccss_y <- c(ccss_H,ccss,"lambda","y")
    if(ZERO_INFLATION=="T")ccss_y <- c(ccss_y,"z")
    
    ccss_y <- c(ccss_y,"deviance")
    
    ######################################################################
    ### Get a chain to optimize initial values
    ######################################################################
    if(INITIALIZE_CHAINS){
      jags.res <- jags.model(file.jags.model,data=ddll,inits=ddii,n.chains=1)
      
      ## message on the terminal
      print(paste("FITTING PARAMETERS for model ",MODEL_NAME," treatment ",tt,"_",cc,"_TEST_chain",sep=""))
      
      ### inference
      jags.res2 <- coda.samples(jags.res,ccss_y,CHAIN_LENGTH,thin=THIN)
      
      write.csv(as.matrix(jags.res2),file=filename)
      ll <- LogPosterior("Test_Chains_",MODEL_NAME,tt,cc,PRIORS,5,45,1,F,T,F)
      ind.max <- which(rowSums(ll)==max(rowSums(ll)))[1]
      
      print(paste("GETTING INITIAL PARAMETERS for model ",MODEL_NAME," treatment ",tt,"_",cc,sep=""))
      dd_t <- read.csv(paste("Test_Chains_",MODEL_NAME,"_",TREAT_TEST,"_",CLONE_TEST,"_",PRIORS,".csv",sep=""))
      
      ddii    <- list(m_f=dd_t[,"m_f"][ind.max],sd_f=dd_t[,"sd_f"][ind.max],
                      m_H=dd_t[,"m_H"][ind.max],sd_H=dd_t[,"sd_H"][ind.max])
      if(ZERO_INFLATION=="T")ddii <- c(ddii,m_pi=dd_t[,"m_pi"][ind.max],sd_pi=dd_t[,"sd_pi"][ind.max])
      if(FECUNDITY=="N")     ddii <- c(ddii,m_r=dd_t[,"m_r"][ind.max],sd_r=dd_t[,"sd_r"][ind.max])
      if(TIME_DEP =="T")     ddii <- c(ddii,m_HH=dd_t[,"m_HH"][ind.max],sd_HH=dd_t[,"sd_HH"][ind.max])
      if(CROWDING =="T")     ddii <- c(ddii,m_HK=dd_t[,"m_HK"][ind.max],sd_HK=dd_t[,"sd_HK"][ind.max])
      
      
      ddii$y <- matrix(NA,ncol=sum(x_end),nrow=Na+Nj+1)
      for(j in 1:sum(x_end)){for(i in 1:(Na+Nj+1)){
        ddii$y[i,j] <-dd_t[ind.max,paste("y.",i,".",j,".",sep="")]
      }}
      
      if(FECUNDITY == "N"){
        ddii$lambda <- rep(1,sum(x_end))
        for(i in 1:sum(x_end))ddii$lambda[i] <- dd_t[ind.max,paste("lambda.",i,".",sep="")]
      }
      if(ZERO_INFLATION=="T"){
        ddii$z <-rep(NA,sum(x_end))
        for(i in 1:sum(x_end))ddii$z[i] <-dd_t[ind.max,paste("z.",i,".",sep="")]
      }
      
    }
    ######################################################################
    ### initialization, fit of the model and plot of results
    ######################################################################
    
    filename <- paste("Opt_Chains_",MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,".csv",sep="")
    
    jags.res <- jags.model(file.jags.model,
                           data=ddll,
                           inits=ddii,
                           #n.adapt=0,
                           n.chains=N_CHAINS)
    
    ## message on the terminal
    print(paste("FITTING PARAMETERS for model ",MODEL_NAME," treatment ",tt,"_",cc,"_",N_CHAINS,"_chains",sep=""))
    
    ### inference
    jags.res2 <- coda.samples(jags.res,ccss_y,CHAIN_LENGTH,thin=THIN)
    
    ### plotting
    plot(jags.res2[,ccss_H],ylab=treat)
    
    plot(jags.res2[,c("f[1]","f[2]","f[3]","f[4]","f[5]")],ylab=treat)
    
    if(ZERO_INFLATION=="T")plot(jags.res2[,c("pi[1]","pi[2]","pi[3]","pi[4]","pi[5]")],ylab=treat)
    
    if(FECUNDITY=="N") plot(jags.res2[,c("r[1]","r[2]","r[3]","r[4]","r[5]")],ylab=treat)
    
    plot(jags.res2[,c("H[1]","H[2]","H[3]","H[4]","H[5]")],ylab=treat)
    
    if(TIME_DEP=="T") plot(jags.res2[,c("HH[1]","HH[2]","HH[3]","HH[4]","HH[5]")],ylab=treat)
    
    if(CROWDING=="T") plot(jags.res2[,c("HK[1]","HK[2]","HK[3]","HK[4]","HK[5]")],ylab=treat)
    
    dev <- jags.res2[[1]][,"deviance"]
    if(N_CHAINS>1) for(nc in 2:N_CHAINS) dev <- c(dev,jags.res2[[nc]][,"deviance"])
    
    PV_1  <- var(dev)/2
    DIC_1 <- mean(dev)+PV_1
    ### use the function with lik=T to get the LogLikelihood computed at the mean parameters (states) 
    ##PV_2  <- -2*LogPosterior(MODEL_NAME,tt,cc,pp,5,45,1,T,T,F)
    ##DIC_2 <- mean(dev)+PV_2
    
    plot(jags.res2[,"deviance"],ylab=treat,
         main=paste("Deviance (DIC ",signif(DIC_1,digits=4),")",sep=""))
    
    if(PLOT_CORR==T){
      chains <- as.matrix(jags.res2[,ccss_H])
      par(mfrow=c(1,1))
      this <- paste(MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,sep="")
      ipairs(chains, main=this)
      cor(chains)
      kk<-1
      for(kk in 1:R){
        parms_kk <- c(paste("f[",kk,"]",sep=""),paste("H[",kk,"]",sep=""))
        if(FECUNDITY=="N")      parms_kk <- c(parms_kk,paste("r[",kk,"]",sep=""))
        if(ZERO_INFLATION=="T") parms_kk <- c(parms_kk,paste("pi[",kk,"]",sep=""))
        if(TIME_DEP=="T")       parms_kk <- c(parms_kk,paste("HH[",kk,"]",sep=""))
        if(CROWDING=="T")       parms_kk <- c(parms_kk,paste("HK[",kk,"]",sep=""))
        chains <- as.matrix(jags.res2[,parms_kk])
        ipairs(chains, main=this)
        cor(chains)
      }
    }
    
    ## save chain
    write.csv(as.matrix(jags.res2),file= filename)
    
    llpp <- LogPosterior("Opt_Chains_",MODEL_NAME,tt,cc,1,Ne,Nj,Na,Dt,F,F,F)
    par(mfrow=c(2,1))
    plot(rowSums(llpp),type="l",ylab="logPosterior")
    plot(-2*(llpp[,"LogLik"]),type="l",ylab="Deviance")
    
  }### close "H" METHOD
  
}### close inference

############################################################################################################
### PLOTTING OF FITTED MODELS RESULTS and COMPUTATION OF MAXIMUM POSTERIOR TRAJECTORY
############################################################################################################

#### take the fitted states and the fitted parameters from the chains and plot them
if(CHECK_CHAIN==T){
  
  ############################################################################################################
  ### get labels
  ############################################################################################################
  
  ## time step for model and observations and time of simulation (days)
  Dt    <- 1  
  
  ## maturation age and maximum age (years)
  age.max.j <-  5
  age.max.a <-  50 
  
  Ne <- INSTAR_TIME
  ## maturation time (number of juvenile classes) use an even number
  Nj <- round(age.max.j/Dt)
  ## adult time (number of adult classes) use an even number
  Na <- round((age.max.a-age.max.j)/Dt)
  
  ## all treatments, to put in loops 
  tt <- TREAT_TEST
  cc <- CLONE_TEST
  
  parnames <- c("f_","H_")
  if(ZERO_INFLATION =="T")parnames <- c(parnames,"pi_")
  if(FECUNDITY=="N")      parnames <- c(parnames,"r_")
  if(TIME_DEP=="T")       parnames <- c(parnames,"HH_")
  if(CROWDING=="T")       parnames <- c(parnames,"HK_")
  
  ############################################################################################################
  ### inference methods
  ############################################################################################################
  
  if(INF_METH=="J"){
    ############################################################################################################    
    ### get file
    ############################################################################################################
    
    filename <- paste("Opt_Chains_",MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,".csv",sep="")
    Chains     <- read.csv(filename)
    
    dim(Chains)
    
    ## get t_end and R from real data
    n_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][1])))
    t_end <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][2])))
    R <- length(n_obs)
    ## check
    if(dim(Chains)[2] == 2+length(parnames)+sum(t_end)*(Na+Nj+Ne))print("PLOTTING TIME SERIES. DATA AND CHAINS EXTRACTED CORRECTLY")
    
    ## set thinning and sampling 
    thin <- THIN_PLOT
    #sample <- seq(length(Chains[,1])/2,length(Chains[,1]),by=thin)
    sample <- seq(1,dim(Chains)[1],by=THIN_PLOT)
    
    ############################################################################################################################
    ## get fitted states
    ############################################################################################################################
    
    ## define multidimensional arrays with chains for every fitted state
    y_c   <- array(0,dim=c(length(sample),Ne+Na+Nj,sum(t_end)))
    y_c_o <- array(0,dim=c(length(sample),3,sum(t_end)))
    
    for(i in 1:(Ne+Na+Nj)){
      for(j in 1:sum(t_end)){
        y_c[,i,j] <- Chains[sample,paste("y.",i,".",j,".",sep="")]
      }}
    
    for(j in 1:sum(t_end)) y_c_o[,1,j] <- rowSums(y_c[,1:Ne,j],dims=1)
    for(j in 1:sum(t_end)) y_c_o[,2,j] <- rowSums(y_c[,(Ne+1):(Ne+Nj),j],dims=1)
    for(j in 1:sum(t_end)) y_c_o[,3,j] <- rowSums(y_c[,(Ne+Nj+2):(Ne+Nj+Na),j],dims=1)
    
    ############################################################################################################################
    ## get maxposterior states
    ############################################################################################################################
    if(PLOT_POSTERIOR){
      y_max_post   <- array(0,dim = c(Na+Nj+Ne,sum(t_end)))
      y_max_post_o <- array(0,dim = c(3,sum(t_end)))
      
      ll <- LogPosterior("Opt_Chains_",MODEL_NAME,tt,cc,PRIORS,Ne,Nj,Na,1,F,T,F)
      ind.max <- which(rowSums(ll)==max(rowSums(ll)))
      
      for(j in 1:sum(x_end))for(i in 1:(Na+Nj+Ne)) y_max_post[i,j]<-Chains[ind.max,paste("y.",i,".",j,".",sep="")]
      
      for(j in 1:sum(t_end)) y_max_post_o[1,j] <- sum(y_max_post[1:Ne,j])
      for(j in 1:sum(t_end)) y_max_post_o[2,j] <- sum(y_max_post[(Ne+1):(Nj+Ne),j])
      for(j in 1:sum(t_end)) y_max_post_o[3,j] <- sum(y_max_post[(Ne+Nj+1):(Nj+Na+Ne),j])
    }    
    ############################################################################################################################
    ## get states from fitted parameters
    ############################################################################################################################
    
    #state vector (eggs + individuals per (juvenile and adult) age class) X timesteps X chain element
    max_tend <- max(t_end)
    y_s <- array(0,dim=c(Na+Nj+Ne,max_tend,length(sample)))
    ## sample from the chains of the parameters and simulate the model from such samples
    ii<-1
    for(ii in 1:length(sample)){
      
      fitted.parms <- c(
        #survival
        H =   Chains[sample[ii],"H"],
        HH =  Chains[sample[ii],"HH"],
        HK =  Chains[sample[ii],"HK"],
        #fecunditity
        r  =   Chains[sample[ii],"r"],
        f  =   Chains[sample[ii],"f"],
        pi =   ifelse(ZERO_INFLATION,Chains[sample[ii],"pi"],0)
        
      )
      
      ## set initial conditions for every replicate
      CIN <- c(rep(0,Ne),3,rep(0,Nj-1),rep(0,Na))
      y_s[,1,ii] <- model(FECUNDITY, ZERO_INFLATION, TIME_DEP, CROWDING,CIN,Dt,fitted.parms,0,Ne,Na,F)
      
      ## run the model
      for(i in 1:(max_tend-1)){
        y_s[,i+1,ii] <- model(FECUNDITY,ZERO_INFLATION,TIME_DEP,CROWDING,y_s[,i,ii],Dt,fitted.parms,i,Ne,Na,F)
      }
    }
    
    ## get mean and 95% confidence interval
    
    yy_e_s <- colSums(y_s[1:Ne,,],dim=1)
    yy_j_s <- colSums(y_s[(Ne+1):(Ne+Nj),,],dim=1)
    yy_a_s <- colSums(y_s[(Ne+Nj+1):(Ne+Na+Nj),,],dim=1)
    
    Med_e_s <- apply(yy_e_s,1,mean)
    Med_j_s <- apply(yy_j_s,1,mean)
    Med_a_s <- apply(yy_a_s,1,mean)
    
    U95_e_s <- apply(yy_e_s,1,ub)
    U95_j_s <- apply(yy_j_s,1,ub)
    U95_a_s <- apply(yy_a_s,1,ub)
    
    L95_e_s <- apply(yy_e_s,1,lb)
    L95_j_s <- apply(yy_j_s,1,lb)
    L95_a_s <- apply(yy_a_s,1,lb)
    
    
    ##############################################################################################
    ## get real data
    ##############################################################################################
    
    ## extract real data from the list
    t_obs   <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,1])))
    y_e_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,2])))
    y_j_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,3])))
    y_a_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,4])))
    
    # get plotting parameters
    nn_s<- numeric(R)
    nn_e<- numeric(R)
    xx_s<- numeric(R)
    xx_e<- numeric(R)
    
    for(k in 1:R){
      nn_s[k] <- ifelse(k==1,1,sum(n_obs[1:(k-1)])+1)
      xx_s[k] <- ifelse(k==1,1,sum(t_end[1:(k-1)])+1)
      nn_e[k] <- sum(n_obs[1:k])
      xx_e[k] <- sum(t_end[1:k])
    }
    
    x_max <- max(t_obs)
    e_max <- max(y_e_obs)
    j_max <- max(y_j_obs)
    a_max <- max(y_a_obs)
    
    ##############################################################################################
    # check and plot
    ##############################################################################################
    if(max(t_end)==x_max) print("second check, time series ready to be plotted")
    
    k<-1
    ## loop in replicates to get and plot the fitted states
    for(k in 1:R){
      ## graphical setting
      par(mfcol=c(3,1))
      
      treatment <- paste(tt,"_",cc,"_","REP_",k,sep="")
      ## real data
      x <- t_obs[nn_s[k]:nn_e[k]]  
      e <- y_e_obs[nn_s[k]:nn_e[k]]
      j <- y_j_obs[nn_s[k]:nn_e[k]]
      a <- y_a_obs[nn_s[k]:nn_e[k]]
      
      ## fitted states mean and upper and lower 95
      predictor <- 1:t_end[k]
      yy_e <- y_c_o[,1,xx_s[k]:xx_e[k]]
      yy_j <- y_c_o[,2,xx_s[k]:xx_e[k]]
      yy_a <- y_c_o[,3,xx_s[k]:xx_e[k]]
      
      Med_e<-numeric(t_end[k]);Med_j<-numeric(t_end[k]);Med_a<-numeric(t_end[k])
      L95_e<-numeric(t_end[k]);L95_j<-numeric(t_end[k]);L95_a<-numeric(t_end[k])
      U95_e<-numeric(t_end[k]);U95_j<-numeric(t_end[k]);U95_a<-numeric(t_end[k])
      
      for(ii in 1:t_end[k]){
        Med_e[ii] <- mean(yy_e[,ii]);L95_e[ii]<-lb(yy_e[,ii]);U95_e[ii]<-ub(yy_e[,ii])
        Med_a[ii] <- mean(yy_a[,ii]);L95_a[ii]<-lb(yy_a[,ii]);U95_a[ii]<-ub(yy_a[,ii])
        Med_j[ii] <- mean(yy_j[,ii]);L95_j[ii]<-lb(yy_j[,ii]);U95_j[ii]<-ub(yy_j[,ii])
      }
      
      ### max posterior states
      if(PLOT_POSTERIOR){
        yy_e_mp <- y_max_post_o[1,xx_s[k]:xx_e[k]]
        yy_j_mp <- y_max_post_o[2,xx_s[k]:xx_e[k]]
        yy_a_mp <- y_max_post_o[3,xx_s[k]:xx_e[k]]
      }
      
      ##############################################################################################
      # plot
      ##############################################################################################
      
      if(PLOT_SHADE==T){
        plot_shade(predictor,Med_e,U95_e,L95_e,0,max(e,max(U95_e)),paste("Eggs",treatment))
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
        plot_shade(predictor,Med_j,U95_j,L95_j,0,max(j,max(U95_j)),paste("Juveniles",treatment))
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
        plot_shade(predictor,Med_a,U95_a,L95_a,0,max(a,max(U95_a)),paste("Adults",treatment))
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        legend("topright",legend=c("Fitted states","Real Data"),col=c("black","red"),lwd=2)
      }#Close PLOT_SHADE
      
      ## plot shade for fitted states and states from fitted paramaters
      if(PLOT_SHADE_2==T){
        Plot_two_models(predictor, Med_e_s[1:length(predictor)],U95_e_s[1:length(predictor)], L95_e_s[1:length(predictor)],Med_e, U95_e, L95_e, 
                        0,max(e,max(U95_e),max(U95_e_s)),paste("Eggs",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_e_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        Plot_two_models(predictor, Med_j_s[1:length(predictor)],U95_j_s[1:length(predictor)], L95_j_s[1:length(predictor)],Med_j, U95_j, L95_j,
                        0,max(j,max(U95_j),max(U95_j_s)),paste("Juveniles",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_j_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        Plot_two_models(predictor, Med_a_s[1:length(predictor)],U95_a_s[1:length(predictor)], L95_a_s[1:length(predictor)],Med_a, U95_a, L95_a,
                        0,max(a,max(U95_a),max(U95_a_s)),paste("Adults",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_a_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        legend("topright",legend=c("Fitted parameters","Fitted states","max posterior","real data"),col=c("black","red","red","red"),
               pch=c(-1,-1,-1,1),lty=c(1,1,2,1),lwd=2)
        if(STATE_RESOLUTION_PLOT){
          Plot_two_models(predictor, Med_e_s[1:length(predictor)],U95_e_s[1:length(predictor)], L95_e_s[1:length(predictor)],Med_e, U95_e, L95_e, 
                          0,1.2*max(e,max(U95_e)),paste("Eggs",treatment))
          if(PLOT_POSTERIOR)lines(predictor,yy_e_mp,col="red",type="l",lty=2,pch=1,cex=2)
          lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
          
          Plot_two_models(predictor, Med_j_s[1:length(predictor)],U95_j_s[1:length(predictor)], L95_j_s[1:length(predictor)],Med_j, U95_j, L95_j,
                          0,1.2*max(j,max(U95_j)),paste("Juveniles",treatment))
          if(PLOT_POSTERIOR)lines(predictor,yy_j_mp,col="red",type="l",lty=2,pch=1,cex=2)
          lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
          
          Plot_two_models(predictor, Med_a_s[1:length(predictor)],U95_a_s[1:length(predictor)], L95_a_s[1:length(predictor)],Med_a, U95_a, L95_a,
                          0,1.2*max(a,max(U95_a)),paste("Adults",treatment))
          if(PLOT_POSTERIOR)lines(predictor,yy_a_mp,col="red",type="l",lty=2,pch=1,cex=2)
          lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
          
          legend("topright",legend=c("Fitted parameters","Fitted states","max posterior","real data"),col=c("black","red","red","red"),
                 pch=c(-1,-1,-1,1),lty=c(1,1,2,1),lwd=2)
          
        }
      } #Close PLOT_SHADE_2
      
    }## close loops in replicates
    
    if(PLOT_ALL==T){
      #########################################################################
      ## plot 3 classes
      ## graphical setting
      par(mfcol=c(3,1))
      treatment <- paste(tt,"_",cc,sep="")
      ## eggs
      predictor <- 1:max_tend
      x <- t_obs[1:n_obs[1]]
      
      e <- y_e_obs[1:n_obs[1]]
      plot_shade(predictor,Med_e_s,U95_e_s,L95_e_s,0,max(e_max,max(U95_e_s)),paste("Eggs",treatment))
      lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
      for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        e <- y_e_obs[nn_s[k]:nn_e[k]]
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
      }
      ## juveniles
      x <- t_obs[1:n_obs[1]]
      j <- y_j_obs[1:n_obs[1]]
      plot_shade(predictor,Med_j_s,U95_j_s,L95_j_s,0,max(j_max,max(U95_j_s)),paste("Juveniles",treatment))
      lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
      for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        j <- y_j_obs[nn_s[k]:nn_e[k]]
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
      }
      
      ## adults
      x <- t_obs[1:n_obs[1]]
      a <- y_a_obs[1:n_obs[1]]
      plot_shade(predictor,Med_a_s,U95_a_s,L95_a_s,0,max(a_max,max(U95_a_s)),paste("Adults",treatment))
      lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
      for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        a <- y_a_obs[nn_s[k]:nn_e[k]]
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
      }
      
      legend("topright",legend=c("Fitted parameters","real data"),col=c("black","red"),
             pch=c(-1,1),lty=c(1,1),lwd=2)
      
    } ## close PLOT_ALL
  }## close INF_METH "J"
  
  if(INF_METH=="R"){
    kk<-1
    for(kk in 1:R){
      ############################################################################################################################
      ### get file
      ############################################################################################################################
      
      filename <- paste("Opt_Chains_",MODEL_NAME,"_",tt,"_",cc,"_REP_",kk,"_",PRIORS,".csv",sep="")
      Chains <- read.csv(filename)
      
      ## get t_end and R from real data
      n_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][1])))
      t_end <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][2])))
      R <- length(n_obs)
      
      ## set thinning and sampling 
      thin <- THIN_PLOT
      #sample <- seq(length(Chains[,1])/2,length(Chains[,1]),by=thin)
      sample <- seq(1,dim(Chains)[1],by=THIN_PLOT)
      
      ############################################################################################################################
      ## get fitted states
      ############################################################################################################################
      
      ## define multidimensional arrays with chains for every fitted state
      y_c   <- array(0,dim=c(length(sample),Ne+Na+Nj,t_end[kk]))
      y_c_o <- array(0,dim=c(length(sample),3,t_end[kk]))
      
      for(i in 1:(Ne+Na+Nj)){
        for(j in 1:t_end[kk]){
          y_c[,i,j] <- Chains[sample,paste("y.",i,".",j,".",sep="")]
        }}
      
      for(j in 1:t_end[kk]) y_c_o[,1,j] <- rowSums(y_c[,1:Ne,j],dims=1)
      for(j in 1:t_end[kk]) y_c_o[,2,j] <- rowSums(y_c[,2:(Nj+1),j],dims=1)
      for(j in 1:t_end[kk]) y_c_o[,3,j] <- rowSums(y_c[,(Nj+2):(Nj+Na+1),j],dims=1)
      ############################################################################################################################
      ## get maxposterior states
      ############################################################################################################################
      if(PLOT_POSTERIOR){
        y_max_post   <- array(0,dim = c(Na+Nj+Ne,t_end[kk]))
        y_max_post_o <- array(0,dim = c(3,t_end[kk]))
        
        ll <- LogPosterior("Opt_Chains_",MODEL_NAME,tt,cc,PRIORS,Ne,Nj,Na,1,F,kk,F)
        ind.max <- which(rowSums(ll)==max(rowSums(ll)))
        
        for(j in 1:t_end[kk])for(i in 1:(Na+Nj+Ne)) y_max_post[i,j]<-Chains[ind.max,paste("y.",i,".",j,".",sep="")]
        
        for(j in 1:t_end[kk]) y_max_post_o[1,j] <- sum(y_max_post[1:Ne,j])
        for(j in 1:t_end[kk]) y_max_post_o[2,j] <- sum(y_max_post[(Ne+1):(Ne+Nj),j])
        for(j in 1:t_end[kk]) y_max_post_o[3,j] <- sum(y_max_post[(Ne+Nj+1):(Ne+Nj+Na),j])
      }      
      ############################################################################################################################
      ## get states from fitted parameters
      ############################################################################################################################
      
      #state vector (eggs + individuals per (juvenile and adult) age class) X timesteps X chain element
      max_tend <- max(t_end[kk])
      y_s <- array(0,dim=c(Na + Nj + Ne,t_end[kk],length(sample)))
      ## sample from the chains of the parameters and simulate the model from such samples
      ii<-1
      for(ii in 1:length(sample)){
        
        fitted.parms <- c(
          #survival
          H =   Chains[sample[ii],"H"],
          HH =   Chains[sample[ii],"HH"],
          #fecunditity
          f =   Chains[sample[ii],"f"],
          r =   Chains[sample[ii],"r"],
          pi =  ifelse(ZERO_INFLATION,Chains[sample[ii],"pi"],1),
          HK =  Chains[sample[ii],"HK"]
        )
        
        ## set initial conditions for every replicate
        CIN <- c(rep(0,Ne),3,rep(0,Nj-1),rep(0,Na))
        y_s[,1,ii] <- model(FECUNDITY,ZERO_INFLATION,TIME_DEP,CROWDING,CIN,Dt,fitted.parms,0,Ne,Na,F)
        
        ## run the model
        for(i in 1:(max_tend-1)){
          y_s[,i+1,ii] <- model(FECUNDITY,ZERO_INFLATION,TIME_DEP,CROWDING,y_s[,i,ii],Dt,fitted.parms,i,Ne,Na,F)
        }
      }
      
      ## get mean and 95% confidence interval
      
      yy_e_s <- colSums(y_s[1:Ne,,])
      yy_j_s <- colSums(y_s[(Ne+1):(Ne+Nj),,],dim=1)
      yy_a_s <- colSums(y_s[(Ne+Nj+1):(Ne+Na+Nj),,],dim=1)
      
      Med_e_s <- apply(yy_e_s,1,mean)
      Med_j_s <- apply(yy_j_s,1,mean)
      Med_a_s <- apply(yy_a_s,1,mean)
      
      U95_e_s <- apply(yy_e_s,1,ub)
      U95_j_s <- apply(yy_j_s,1,ub)
      U95_a_s <- apply(yy_a_s,1,ub)
      
      L95_e_s <- apply(yy_e_s,1,lb)
      L95_j_s <- apply(yy_j_s,1,lb)
      L95_a_s <- apply(yy_a_s,1,lb)
      
      ##############################################################################################
      ## get real data
      ##############################################################################################      
      
      ## extract real data from the list
      t_obs   <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,1])))
      y_e_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,2])))
      y_j_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,3])))
      y_a_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,4])))
      
      # get plotting parameters
      nn_s<- numeric(R)
      nn_e<- numeric(R)
      xx_s<- numeric(R)
      xx_e<- numeric(R)
      
      for(k in 1:R){
        nn_s[k] <- ifelse(k==1,1,sum(n_obs[1:(k-1)])+1)
        xx_s[k] <- ifelse(k==1,1,sum(t_end[1:(k-1)])+1)
        nn_e[k] <- sum(n_obs[1:k])
        xx_e[k] <- sum(t_end[1:k])
      }
      
      x_max <- max(t_obs)
      e_max <- max(y_e_obs)
      j_max <- max(y_j_obs)
      a_max <- max(y_a_obs)
      
      ##############################################################################################
      # check
      ##############################################################################################
      
      #check
      if(max(t_end)==x_max) print(paste("ok_REP_",kk,sep=""))
      
      ## graphical setting
      par(mfcol=c(3,1))
      
      treatment <- paste(tt,"_",cc,"_","REP_",kk,sep="")
      ## real data
      x <- t_obs[nn_s[kk]:nn_e[kk]]  
      e <- y_e_obs[nn_s[kk]:nn_e[kk]]
      j <- y_j_obs[nn_s[kk]:nn_e[kk]]
      a <- y_a_obs[nn_s[kk]:nn_e[kk]]
      
      ## fitted states mean and upper and lower 95
      predictor <- 1:t_end[kk]
      yy_e <- y_c_o[,1,1:t_end[kk]]
      yy_j <- y_c_o[,2,1:t_end[kk]]
      yy_a <- y_c_o[,3,1:t_end[kk]]
      
      Med_e<-numeric(t_end[kk]);Med_j<-numeric(t_end[kk]);Med_a<-numeric(t_end[kk])
      L95_e<-numeric(t_end[kk]);L95_j<-numeric(t_end[kk]);L95_a<-numeric(t_end[kk])
      U95_e<-numeric(t_end[kk]);U95_j<-numeric(t_end[kk]);U95_a<-numeric(t_end[kk])
      
      for(ii in 1:t_end[kk]){
        Med_e[ii] <- mean(yy_e[,ii]);L95_e[ii]<-lb(yy_e[,ii]);U95_e[ii]<-ub(yy_e[,ii])
        Med_a[ii] <- mean(yy_a[,ii]);L95_a[ii]<-lb(yy_a[,ii]);U95_a[ii]<-ub(yy_a[,ii])
        Med_j[ii] <- mean(yy_j[,ii]);L95_j[ii]<-lb(yy_j[,ii]);U95_j[ii]<-ub(yy_j[,ii])
      }
      
      ### max posterior states
      if(PLOT_POSTERIOR){
        yy_e_mp <- y_max_post_o[1,1:t_end[kk]]
        yy_j_mp <- y_max_post_o[2,1:t_end[kk]]
        yy_a_mp <- y_max_post_o[3,1:t_end[kk]]
      }
      
      ##############################################################################################
      # plot
      ##############################################################################################
      
      if(PLOT_SHADE==T){
        plot_shade(predictor,Med_e,U95_e,L95_e,0,max(e,max(U95_e)),paste("Eggs",treatment))
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
        plot_shade(predictor,Med_j,U95_j,L95_j,0,max(j,max(U95_j)),paste("Juveniles",treatment))
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
        plot_shade(predictor,Med_a,U95_a,L95_a,0,max(a,max(U95_a)),paste("Adults",treatment))
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        legend("topright",legend=c("Fitted states","Real Data"),col=c("black","red"),lwd=2)
      }#Close PLOT_SHADE
      
      ## plot shade for fitted states and fitted paramaters
      if(PLOT_SHADE_2==T){
        Plot_two_models(predictor, Med_e_s[1:length(predictor)],U95_e_s[1:length(predictor)], L95_e_s[1:length(predictor)],Med_e, U95_e, L95_e, 
                        0,max(e,max(U95_e),max(U95_e_s)),paste("Eggs",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_e_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        Plot_two_models(predictor, Med_j_s[1:length(predictor)],U95_j_s[1:length(predictor)], L95_j_s[1:length(predictor)],Med_j, U95_j, L95_j,
                        0,max(j,max(U95_j),max(U95_j_s)),paste("Juveniles",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_j_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        Plot_two_models(predictor, Med_a_s[1:length(predictor)],U95_a_s[1:length(predictor)], L95_a_s[1:length(predictor)],Med_a, U95_a, L95_a,
                        0,max(a,max(U95_a),max(U95_a_s)),paste("Adults",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_a_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        legend("topright",legend=c("Fitted parameters","Fitted states","max posterior","real data"),col=c("black","red","red","red"),
               pch=c(-1,-1,-1,1),lty=c(1,1,2,1),lwd=2)
        if(STATE_RESOLUTION_PLOT){
          Plot_two_models(predictor, Med_e_s[1:length(predictor)],U95_e_s[1:length(predictor)], L95_e_s[1:length(predictor)],Med_e, U95_e, L95_e, 
                          0,1.2*max(e,max(U95_e)),paste("Eggs",treatment))
          if(PLOT_POSTERIOR)lines(predictor,yy_e_mp,col="red",type="l",lty=2,pch=1,cex=2)
          lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
          
          Plot_two_models(predictor, Med_j_s[1:length(predictor)],U95_j_s[1:length(predictor)], L95_j_s[1:length(predictor)],Med_j, U95_j, L95_j,
                          0,1.2*max(j,max(U95_j)),paste("Juveniles",treatment))
          if(PLOT_POSTERIOR)lines(predictor,yy_j_mp,col="red",type="l",lty=2,pch=1,cex=2)
          lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
          
          Plot_two_models(predictor, Med_a_s[1:length(predictor)],U95_a_s[1:length(predictor)], L95_a_s[1:length(predictor)],Med_a, U95_a, L95_a,
                          0,1.2*max(a,max(U95_a)),paste("Adults",treatment))
          if(PLOT_POSTERIOR)lines(predictor,yy_a_mp,col="red",type="l",lty=2,pch=1,cex=2)
          lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
          
          legend("topright",legend=c("Fitted parameters","Fitted states","max posterior","real data"),col=c("black","red","red","red"),
                 pch=c(-1,-1,-1,1),lty=c(1,1,2,1),lwd=2)
        }
      } #Close PLOT_SHADE_2
      
    }## close replicate loop
    
  }## close INF_METH "R"
  
  if(INF_METH=="H"){
    ############################################################################################################################
    ### get file
    ############################################################################################################################
    
    filename <- paste("Opt_Chains_",MODEL_NAME,"_",tt,"_",cc,"_",PRIORS,".csv",sep="")
    
    Chains <- read.csv(filename)
    
    ## get t_end and R from real data
    n_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][1])))
    t_end <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][2])))
    R <- length(n_obs)
    
    ## check
    if(dim(Chains)[2] == 2+length(parnames)*(R+2)+sum(t_end)*(Na+Nj+Ne)){
      print("PLOTTING TIME SERIES. DATA AND CHAINS EXTRACTED CORRECTLY")
    }
    
    ## set thinning and sampling 
    thin <- THIN_PLOT
    #sample <- seq(length(Chains[,1])/2,length(Chains[,1]),by=thin)
    sample <- seq(1,dim(Chains)[1],by=THIN_PLOT)
    
    ############################################################################################################################
    ## get fitted states 
    ############################################################################################################################
    
    ## define multidimensional arrays with chains for every fitted state
    y_c   <- array(0,dim=c(length(sample),Ne+Na+Nj,sum(t_end)))
    y_c_o <- array(0,dim=c(length(sample),3,sum(t_end)))
    
    for(i in 1:(Ne+Na+Nj)){
      for(j in 1:sum(t_end)){
        y_c[,i,j] <- Chains[sample,paste("y.",i,".",j,".",sep="")]
      }}
    
    for(j in 1:sum(t_end)) y_c_o[,1,j] <- rowSums(y_c[,1:Ne,j],dims=1)
    for(j in 1:sum(t_end)) y_c_o[,2,j] <- rowSums(y_c[,(Ne+1):(Ne+Nj),j],dims=1)
    for(j in 1:sum(t_end)) y_c_o[,3,j] <- rowSums(y_c[,(Ne+Nj+2):(Ne+Nj+Na),j],dims=1)
    
    ############################################################################################################################
    ## get maxposterior states
    ############################################################################################################################
    
    if(PLOT_POSTERIOR){
      y_max_post   <- array(0,dim = c(Na+Nj+Ne,sum(t_end)))
      y_max_post_o <- array(0,dim = c(3,sum(t_end)))
      
      ll <- LogPosterior("Opt_Chains_",MODEL_NAME,tt,cc,PRIORS,Ne,Nj,Na,1,F,T,F)
      ind.max <- which(rowSums(ll)==max(rowSums(ll)))
      
      for(j in 1:sum(x_end))for(i in 1:(Na+Nj+Ne)) y_max_post[i,j]<-Chains[ind.max,paste("y.",i,".",j,".",sep="")]
      
      for(j in 1:sum(t_end)) y_max_post_o[1,j] <- sum(y_max_post[1:Ne,j])
      for(j in 1:sum(t_end)) y_max_post_o[2,j] <- sum(y_max_post[(Ne+1):(Nj+Ne),j])
      for(j in 1:sum(t_end)) y_max_post_o[3,j] <- sum(y_max_post[(Ne+Nj+1):(Nj+Na+Ne),j])
    }    
    
    ############################################################################################################################
    ## get states from replica specific parameters
    ############################################################################################################################
    
    max_tend <- max(t_end)
    Med_rs <-array(0,dim=c(3,max_tend,R)); 
    L95_rs  <-array(0,dim=c(3,max_tend,R));
    U95_rs  <-array(0,dim=c(3,max_tend,R));
    
    kk <- 1
    for(kk in 1:R){
      print(paste("Replica specific parameter fitting, REP ",kk,sep=""))
      y_rs <- array(0,dim=c(Na + Nj + Ne,max_tend,length(sample)))
      
      ii <- 1
      for(ii in 1:length(sample)){
        
        fitted.parms <- c(
          #survival
          H  =   Chains[sample[ii],paste("H.",kk,".",sep="")],
          pi =   ifelse(ZERO_INFLATION=="T",Chains[sample[ii],paste("pi.",kk,".",sep="")],0),
          
          #pi =   Chains[sample[ii],paste("pi.",kk,".",sep="")],
          f  =   Chains[sample[ii],paste("f.",kk,".",sep="")]
        )  
        if(FECUNDITY=="N")  fitted.parms <- c(fitted.parms, r  =  Chains[sample[ii],paste("r.",kk,".",sep="")])
        if(TIME_DEP=="T")       fitted.parms <- c(fitted.parms, HH =  Chains[sample[ii],paste("HH.",kk,".",sep="")])
        if(CROWDING=="T")       fitted.parms <- c(fitted.parms, HK =  Chains[sample[ii],paste("HK.",kk,".",sep="")])
        
        ## set initial conditions for every replicate
        CIN <- c(rep(0,Ne),3,rep(0,Nj-1),rep(0,Na))
        y_rs[,1,ii] <- model(FECUNDITY,ZERO_INFLATION,TIME_DEP,CROWDING,CIN,Dt,fitted.parms,0,Ne,Na,F)
        
        ## run the model
        for(i in 1:(max_tend-1)){
          y_rs[,i+1,ii] <- model(FECUNDITY,ZERO_INFLATION,TIME_DEP,CROWDING,y_rs[,i,ii],Dt,fitted.parms,i,Ne,Na,F)
        }
      }
      
      ## get mean and 95% confidence interval
      
      yy_e_rs <- colSums(y_rs[1:Ne,,],dim=1)
      yy_j_rs <- colSums(y_rs[(Ne+1):(Ne+Nj),,],dim=1)
      yy_a_rs <- colSums(y_rs[(Ne+Nj+1):(Ne+Na+Nj),,],dim=1)
      
      Med_rs[1,,kk] <- apply(yy_e_rs,1,mean,na.rm=TRUE)
      Med_rs[2,,kk] <- apply(yy_j_rs,1,mean,na.rm=TRUE)
      Med_rs[3,,kk] <- apply(yy_a_rs,1,mean,na.rm=TRUE)
      
      U95_rs[1,,kk] <- apply(yy_e_rs,1,ub)
      U95_rs[2,,kk] <- apply(yy_j_rs,1,ub)
      U95_rs[3,,kk] <- apply(yy_a_rs,1,ub)
      
      L95_rs[1,,kk] <- apply(yy_e_rs,1,lb)
      L95_rs[2,,kk] <- apply(yy_j_rs,1,lb)
      L95_rs[3,,kk] <- apply(yy_a_rs,1,lb)  
    }  
    
    ############################################################################################################################
    ## get states from fitted parameters
    ############################################################################################################################
    
    ############################################################################################################################
    ### first way
    
    #state vector (eggs + individuals per (juvenile and adult) age class) X timesteps X chain element
    max_tend <- max(t_end)
    y_s <- array(0,dim=c(Na + Nj + Ne,max_tend,length(sample)))
    ## sample from the chains of the parameters and simulate the model from such samples
    
    print("getting fitted states from hyperparameters")
    for(ii in 1:length(sample)){
      
      ## conversion functions that give the mean and sd of associated normal (to put into R)
      mfun  <- function(m,s){log(m)-log(1+s^2/m^2)/2}
      sdfun <- function(m,s){sqrt(log(1+(s/m)^2))}
      
      ### means and standard deviations of the lognormals
      m_f   <- Chains[sample[ii],"m_f"]
      sd_f  <- Chains[sample[ii],"sd_f"]
      mlf   <- mfun(m_f,sd_f)
      sdlf  <- sdfun(m_f,sd_f)
      
      if(ZERO_INFLATION=="T"){
      m_pi   <- Chains[sample[ii],"m_pi"]
      sd_pi  <- Chains[sample[ii],"sd_pi"]
      mlpi   <- mfun(m_pi,sd_pi)
      sdlpi  <- sdfun(m_pi,sd_pi)
      }
      
      m_H   <- Chains[sample[ii],"m_H"]
      sd_H  <- Chains[sample[ii],"sd_H"]
      mlH   <- mfun(m_H,sd_H)
      sdlH  <- sdfun(m_H,sd_H)
      
      if(FECUNDITY=="N"){
        m_r  <- Chains[sample[ii],"m_r"]
        sd_r  <- Chains[sample[ii],"sd_r"]  
        mlr   <- mfun(m_r,sd_r)
        sdlr  <- sdfun(m_r,sd_r)
      }
      
      if(TIME_DEP=="T"){
        m_HH   <- Chains[sample[ii],"m_HH"]
        sd_HH  <- Chains[sample[ii],"sd_HH"]
        mlHH   <- mfun(m_HH,sd_HH)
        sdlHH  <- sdfun(m_HH,sd_HH)
      }
      
      if(CROWDING=="T"){
        m_HK  <- Chains[sample[ii],"m_HK"]
        sd_HK  <- Chains[sample[ii],"sd_HK"]
        mlHK   <- mfun(m_HK,sd_HK)
        sdlHK  <- sdfun(m_HK,sd_HK)
      }
      
      fitted.parms <- c(
        #fecunditity
        f  =   ifelse(PRIORS==1,rlnormTrunc(1,mlf,sdlf,0,30),runif(1,0,10)),
        #survival
        H =   ifelse(PRIORS==1,rlnorm(1,mlH,sdlH),runif(1,0,1))
      )
      
      if(ZERO_INFLATION=="T")fitted.parms<-c(fitted.parms,pi =ifelse(PRIORS==1,rnormTrunc(1,m_pi,sd_pi,0,1),runif(1,0,1)))
      
      if(FECUNDITY=="N")     fitted.parms<-c(fitted.parms,r  =ifelse(PRIORS==1,rlnorm(1,mlr,sdlr),runif(1,0,100)))
      if(TIME_DEP=="T")      fitted.parms<-c(fitted.parms,HH =ifelse(PRIORS==1,rnorm(1,m_HH,sd_HH),runif(1,0,1)))
      if(CROWDING=="T")      fitted.parms<-c(fitted.parms,HK =ifelse(PRIORS==1,rlnorm(1,mlHK,sdlHK),runif(1,0,1)))
      
      ## set initial conditions for every replicate
      CIN <- c(rep(0,Ne),3,rep(0,Nj-1),rep(0,Na))
      y_s[,1,ii] <- model(FECUNDITY,ZERO_INFLATION,TIME_DEP,CROWDING,CIN,Dt,fitted.parms,0,Ne,Na,F)
      
      ## run the model
      for(i in 1:(max_tend-1)){
        y_s[,i+1,ii] <- model(FECUNDITY,ZERO_INFLATION,TIME_DEP,CROWDING,y_s[,i,ii],Dt,fitted.parms,i,Ne,Na,F)
      }
      # print(ii)
      # if(sum(y_s[1,,ii])>10000){
      #   print(fitted.parms)
      #   print(y_s[1,,ii])
      # }
    }
    
    ############################################################################################################################
    ### second way
    # 
    # #state vector (eggs + individuals per (juvenile and adult) age class) X timesteps X chain element
    # max_tend <- max(t_end)
    # y_s <- array(0,dim=c(Na + Nj + 1,max_tend,length(sample)))
    # ## sample from the chains of the parameters and simulate the model from such samples
    # 
    # ## conversion functions that give the mean and sd of associated normal (to put into R)
    # mfun  <- function(m,s){log(m)-log(1+s^2/m^2)/2}
    # sdfun <- function(m,s){sqrt(log(1+(s/m)^2))}
    # 
    # ### means and standard deviations of the lognormals
    # m_f   <- Chains[,"m_f"]
    # sd_f  <- Chains[,"sd_f"]
    # mlf   <- mfun(m_f,sd_f)
    # sdlf  <- sdfun(m_f,sd_f)
    # 
    # m_H   <- Chains[,"m_H"]
    # sd_H  <- Chains[,"sd_H"]
    # mlH   <- mfun(m_H,sd_H)
    # sdlH  <- sdfun(m_H,sd_H)
    # 
    # if(FECUNDITY=="N"){
    #   m_r  <- Chains[,"m_r"]
    #   sd_r  <- Chains[,"sd_r"]  
    #   mlr   <- mfun(m_r,sd_r)
    #   sdlr  <- sdfun(m_r,sd_r)
    # }
    # 
    # if(TIME_DEP=="T"){
    #   m_HH   <- Chains[,"m_HH"]
    #   sd_HH  <- Chains[,"sd_HH"]
    #   mlHH   <- mfun(m_HH,sd_HH)
    #   sdlHH  <- sdfun(m_HH,sd_HH)
    # }
    # 
    # if(CROWDING=="T"){
    #   m_HK  <- Chains[,"m_HK"]
    #   sd_HK  <- Chains[,"sd_HK"]
    #   mlHK   <- mfun(m_HK,sd_HK)
    #   sdlHK  <- sdfun(m_HK,sd_HK)
    # }
    # 
    # for(ii in 1:length(sample)){
    # fitted.parms <- c(
    #   #fecunditity
    #   f =   ifelse(PRIORS==1,rlnorm(1,mlf[ii],sdlf[ii]),runif(1,0,10)),
    #   #survival
    #   H =   ifelse(PRIORS==1,rlnorm(1,mlH[ii],sdlH[ii]),runif(1,0,1))
    # )
    # 
    # if(FECUNDITY=="N")fitted.parms<-c(fitted.parms,r = ifelse(PRIORS==1,rlnorm(1,mlr[ii],sdlr[ii]),runif(1,0,100)))
    # if(TIME_DEP=="T")      fitted.parms<-c(fitted.parms,HH= ifelse(PRIORS==1,rlnorm(1,mlHH[ii],sdlHH[ii]),runif(1,0,1)))
    # if(CROWDING=="T")      fitted.parms<-c(fitted.parms,HK =  ifelse(PRIORS==1,rlnorm(1,mlHK[ii],sdlHK[ii]),runif(1,0,1)))
    # 
    # ## set initial conditions for every replicate
    # CIN <- c(0,3,rep(0,Nj-1),rep(0,Na))
    # y_s[,1,ii] <- model(FECUNDITY,TIME_DEP,CROWDING,CIN,Dt,fitted.parms,0,Na,F)
    # 
    # ## run the model
    # for(i in 1:(max_tend-1)){
    #   y_s[,i+1,ii] <- model(FECUNDITY,TIME_DEP,CROWDING,y_s[,i,ii],Dt,fitted.parms,i,Na,F)
    # }
    # print(ii)
    # if(sum(y_s[1,,ii])>10000){
    #   print(fitted.parms)
    #   print(y_s[1,,ii])
    # }
    # }
    
    ############################################################################################################################
    ## get mean and 95% confidence interval
    
    yy_e_s <- colSums(y_s[1:Ne,,],dim=1)
    yy_j_s <- colSums(y_s[(Ne+1):(Ne+Nj),,],dim=1)
    yy_a_s <- colSums(y_s[(Ne+Nj+1):(Ne+Na+Nj),,],dim=1)
    
    Med_e_s <- apply(yy_e_s,1,mean)
    Med_j_s <- apply(yy_j_s,1,mean)
    Med_a_s <- apply(yy_a_s,1,mean)
    
    U95_e_s <- apply(yy_e_s,1,ub)
    U95_j_s <- apply(yy_j_s,1,ub)
    U95_a_s <- apply(yy_a_s,1,ub)
    
    L95_e_s <- apply(yy_e_s,1,lb)
    L95_j_s <- apply(yy_j_s,1,lb)
    L95_a_s <- apply(yy_a_s,1,lb)
    
    
    ##############################################################################################
    ## get real data
    ##############################################################################################
    ## extract real data from the list
    t_obs   <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,1])))
    y_e_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,2])))
    y_j_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,3])))
    y_a_obs <- as.numeric(as.character(unlist(dd.to.fit.list[[tt]][[cc]][3][[1]][,4])))
    
    # get plotting parameters
    nn_s<- numeric(R)
    nn_e<- numeric(R)
    xx_s<- numeric(R)
    xx_e<- numeric(R)
    
    for(k in 1:R){
      nn_s[k] <- ifelse(k==1,1,sum(n_obs[1:(k-1)])+1)
      xx_s[k] <- ifelse(k==1,1,sum(t_end[1:(k-1)])+1)
      nn_e[k] <- sum(n_obs[1:k])
      xx_e[k] <- sum(t_end[1:k])
    }
    
    x_max <- max(t_obs)
    e_max <- max(y_e_obs)
    j_max <- max(y_j_obs)
    a_max <- max(y_a_obs)
    
    #check
    if(max(t_end)==x_max) print("ok ready to plot\n")
    
    
    ##############################################################################################
    ### plotting
    ##############################################################################################
    ## loop in replicates to get and plot the fitted states with or without the states 
    ## obtained through the fitted paramaters
    k<-1
    for(k in 1:R){
      ## graphical setting
      par(mfcol=c(3,1))
      
      treatment <- paste(tt,"_",cc,"_","REP_",k,sep="")
      predictor <- 1:t_end[k]
      
      ## real data
      x <- t_obs[nn_s[k]:nn_e[k]]  
      e <- y_e_obs[nn_s[k]:nn_e[k]]
      j <- y_j_obs[nn_s[k]:nn_e[k]]
      a <- y_a_obs[nn_s[k]:nn_e[k]]
      
      ## fitted states mean and upper and lower 95
      yy_e <- y_c_o[,1,xx_s[k]:xx_e[k]]
      yy_j <- y_c_o[,2,xx_s[k]:xx_e[k]]
      yy_a <- y_c_o[,3,xx_s[k]:xx_e[k]]
      
      Med_e<-numeric(t_end[k]);Med_j<-numeric(t_end[k]);Med_a<-numeric(t_end[k])
      L95_e<-numeric(t_end[k]);L95_j<-numeric(t_end[k]);L95_a<-numeric(t_end[k])
      U95_e<-numeric(t_end[k]);U95_j<-numeric(t_end[k]);U95_a<-numeric(t_end[k])
      
      for(ii in 1:t_end[k]){
        Med_e[ii] <- mean(yy_e[,ii]);L95_e[ii]<-lb(yy_e[,ii]);U95_e[ii]<-ub(yy_e[,ii])
        Med_a[ii] <- mean(yy_a[,ii]);L95_a[ii]<-lb(yy_a[,ii]);U95_a[ii]<-ub(yy_a[,ii])
        Med_j[ii] <- mean(yy_j[,ii]);L95_j[ii]<-lb(yy_j[,ii]);U95_j[ii]<-ub(yy_j[,ii])
      }
      
      ### max posterior states
      if(PLOT_POSTERIOR){
      yy_e_mp <- y_max_post_o[1,xx_s[k]:xx_e[k]]
      yy_j_mp <- y_max_post_o[2,xx_s[k]:xx_e[k]]
      yy_a_mp <- y_max_post_o[3,xx_s[k]:xx_e[k]]
      }
      
      ############################################################################################################################
      ## plot shade for fitted states only
      ############################################################################################################################
      
      if(PLOT_SHADE==T){
        plot_shade(predictor,Med_e,U95_e,L95_e,0,max(e,max(U95_e)),paste("Eggs",treatment))
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
        plot_shade(predictor,Med_j,U95_j,L95_j,0,max(j,max(U95_j)),paste("Juveniles",treatment))
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
        plot_shade(predictor,Med_a,U95_a,L95_a,0,max(a,max(U95_a)),paste("Adults",treatment))
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        legend("topright",legend=c("Fitted states","Real Data"),col=c("black","red"),lwd=2)
      }#Close PLOT_SHADE
      
      if(PLOT_SHADE_2==T){
        
        ### fitted states and fitted parameters
        Plot_two_models(predictor, 
                        Med_e_s[1:length(predictor)],U95_e_s[1:length(predictor)],L95_e_s[1:length(predictor)],
                        Med_e, U95_e, L95_e, 
                        0,max(e,max(U95_e),max(U95_e_s)),
                        paste("Eggs",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_e_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)
        
        Plot_two_models(predictor, 
                        Med_j_s[1:length(predictor)],U95_j_s[1:length(predictor)], L95_j_s[1:length(predictor)],
                        Med_j, U95_j, L95_j,
                        0,max(j,max(U95_j),max(U95_j_s)),
                        paste("Juveniles",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_j_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)
        
        Plot_two_models(predictor, 
                        Med_a_s[1:length(predictor)],U95_a_s[1:length(predictor)], L95_a_s[1:length(predictor)],
                        Med_a, U95_a, L95_a,
                        0,max(a,max(U95_a),max(U95_a_s)),
                        paste("Adults",treatment))
        if(PLOT_POSTERIOR)lines(predictor,yy_a_mp,col="red",type="l",lty=2,pch=1,cex=2)
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        legend("topright",legend=c("Fitted parameters","Fitted states","max posterior","real data"),col=c("black","red","red","red"),
               pch=c(-1,-1,-1,1),lty=c(1,1,2,1),lwd=2)
        
        ## fitted states and replica specific parameteres
        Plot_two_models(predictor, 
                        Med_rs[1,1:length(predictor),k], U95_rs[1,1:length(predictor),k], L95_rs[1,1:length(predictor),k],
                        Med_e, U95_e, L95_e, 
                        0,max(e,max(U95_e),max(U95_rs[1,1:length(predictor),k])),paste("Eggs",treatment))
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
        if(PLOT_POSTERIOR)lines(predictor,yy_e_mp,col="red",type="l",lty=2,pch=1,cex=2)
        
        Plot_two_models(predictor, 
                        Med_rs[2,1:length(predictor),k], U95_rs[2,1:length(predictor),k], L95_rs[2,1:length(predictor),k], 
                        Med_j, U95_j, L95_j,
                        0,max(j,max(U95_j),max(U95_rs[2,1:length(predictor),k])),
                        paste("Juveniles",treatment))
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)
        if(PLOT_POSTERIOR)lines(predictor,yy_j_mp,col="red",type="l",lty=2,pch=1,cex=2)
        
        Plot_two_models(predictor, 
                        Med_rs[3,1:length(predictor),k], U95_rs[3,1:length(predictor),k], L95_rs[3,1:length(predictor),k],
                        Med_a, U95_a, L95_a,
                        0,max(a,max(U95_a),max(U95_rs[3,1:length(predictor),k])),
                        paste("Adults",treatment))
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
        if(PLOT_POSTERIOR)lines(predictor,yy_a_mp,col="red",type="l",lty=2,pch=1,cex=2)
        
        legend("topright",legend=c("rep spec parameters","Fitted states","max posterior","real data"),col=c("black","red","red","red"),
               pch=c(-1,-1,-1,1),lty=c(1,1,2,1),lwd=2)
      } #Close PLOT_SHADE_2
      
      ## plot shade for fitted states and fitted paramaters
      if(PLOT_SHADE_3==T){
        Plot_three_models(predictor, 
                          Med_e_s[1:length(predictor)],U95_e_s[1:length(predictor)], L95_e_s[1:length(predictor)],
                          Med_rs[1,1:length(predictor),k], U95_rs[1,1:length(predictor),k], L95_rs[1,1:length(predictor),k],
                          Med_e, U95_e, L95_e,
                          0,max(e,max(U95_rs[1,,k]),max(U95_e_s),U95_e),
                          paste("Eggs",treatment))
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        Plot_three_models(predictor, 
                          Med_j_s[1:length(predictor)],U95_j_s[1:length(predictor)], L95_j_s[1:length(predictor)],
                          Med_rs[2,1:length(predictor),k], U95_rs[2,1:length(predictor),k], L95_rs[2,1:length(predictor),k], 
                          Med_j, U95_j, L95_j,
                          0,max(j,max(U95_rs[2,,k]),max(U95_j_s),U95_j),
                          paste("Juveniles",treatment))
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        Plot_three_models(predictor, 
                          Med_a_s[1:length(predictor)],U95_a_s[1:length(predictor)], L95_a_s[1:length(predictor)],
                          Med_rs[3,1:length(predictor),k], U95_rs[3,1:length(predictor),k], L95_rs[3,1:length(predictor),k],
                          Med_a, U95_a, L95_a,
                          0,max(a,max(U95_rs[3,,k]),max(U95_a_s),U95_a),
                          paste("Adults",treatment))
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
        
        legend("topright",legend=c("Fitted mean parameters","Rep specific Fitted parms","fitted states","real data"),col=c("black","blue","red","red"),
               pch=c(-1,-1,-1,1),lty=c(1,1,1,1),lwd=2)
      } #Close PLOT_SHADE_3
      
    }## close loops in replicates
    
    if(PLOT_ALL==T){
      #########################################################################
      ## plot 3 classes
      ## graphical setting
      par(mfcol=c(3,1))
      treatment <- paste(tt,"_",cc,sep="")
      ## eggs
      predictor <- 1:max_tend
      x <- t_obs[1:n_obs[1]]
      
      e <- y_e_obs[1:n_obs[1]]
      plot_shade(predictor,Med_e_s,U95_e_s,L95_e_s,0,max(e_max,max(U95_e_s)),paste("Eggs",treatment))
      lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
      for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        e <- y_e_obs[nn_s[k]:nn_e[k]]
        lines(x,e,lwd=1,col="red",type="b",pch=16,cex=2)    
      }
      ## juveniles
      x <- t_obs[1:n_obs[1]]
      j <- y_j_obs[1:n_obs[1]]
      plot_shade(predictor,Med_j_s,U95_j_s,L95_j_s,0,max(j_max,max(U95_j_s)),paste("Juveniles",treatment))
      lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
      for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        j <- y_j_obs[nn_s[k]:nn_e[k]]
        lines(x,j,lwd=1,col="red",type="b",pch=16,cex=2)    
      }
      
      ## adults
      x <- t_obs[1:n_obs[1]]
      a <- y_a_obs[1:n_obs[1]]
      plot_shade(predictor,Med_a_s,U95_a_s,L95_a_s,0,max(a_max,max(U95_a_s)),paste("Adults",treatment))
      lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
      for(k in 2:R){
        x <- t_obs[nn_s[k]:nn_e[k]]
        a <- y_a_obs[nn_s[k]:nn_e[k]]
        lines(x,a,lwd=1,col="red",type="b",pch=16,cex=2)    
      }
      
      legend("topright",legend=c("Fitted parameters","real data"),col=c("black","red"),
             pch=c(-1,1),lty=c(1,1),lwd=2)
      
    } ## close PLOT_ALL
  }## close INF_METH "H"
  
}## close CHECK CHAIN

############################################################################################################
if(PLOT_IN_TERMINAL==F) dev.off()  
