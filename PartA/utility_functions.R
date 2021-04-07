# Three utility functions in the Bayesian intro course

getMLE <- function(sol, dig = 3){
  MLE <- sol$par
  VC <- solve(sol$hessian)
  ASE <- sqrt(diag(VC))
  print(cbind(MLE, ASE), dig)
}

nimble_summary <- function(samples, params=NULL, digits=3){
  if(!is.null(params)){
    samples <- jagsUI:::order.params(samples, params, FALSE, FALSE)
  }
  mat <- as.matrix(samples)
  nchain <- length(samples)
  niter <- nrow(samples[[1]])
  rhat <- sapply(1:ncol(samples[[1]]), function(i){
    coda::gelman.diag(samples[,i], autoburnin=FALSE)$psrf[1,1]
  })
  stats <- t(apply(mat, 2, function(x){
    x <- na.omit(x)
    c(mean=mean(x), sd=sd(x), quantile(x, c(0.025,0.5,0.975)))
  }))
  out <- data.frame(stats, rhat=rhat, check.names=FALSE)
  cat("Estimates based on",nchain,"chains of",niter,"iterations\n")
  round(out, digits=digits)
}

tmb_summary <- function(tmb_obj){
  npar <- length(tmb_obj$par)
  pnames_fixed <- names(tmb_obj$par)
  out <- summary(sdreport(tmb_obj))
  pnames <- rownames(out)
  pnames[1:npar] <- pnames_fixed
  pcount <- sapply(unique(pnames), function(x) sum(pnames==x))
  idx <- unlist(sapply(pcount, function(i){
    if(i == 1) return("")
    paste0("[",1:i,"]")
    }))
  pnames <- paste0(pnames, idx)
  rownames(out) <- pnames
  out
}