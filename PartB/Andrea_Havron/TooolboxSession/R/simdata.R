 
#' Simulate spatial MVN using a Matern covariance function
#'
#' @param Range distance at which spatial correlation is negligible
#' @param sig2 spatial variance
#' @param Dmat distance matrix
#' @param Nu smoothness parameter - usually fixed at 1
#'
#' @return vector of spatial effects
#' @export
#'
#' @examples
sim.omega <- function(Range, sig2, Dmat, Nu = 1){
  Kappa <- sqrt(8)/Range
  N <- dim(Dmat)[1]
  Sigma <- matrix(NA,N,N)
  #matern covariance
  Sigma <- sig2 * 2^(1-Nu)/gamma(Nu) * (Kappa * Dmat) * besselK(Kappa * Dmat, Nu)
  diag(Sigma) <- sig2
  omega <- t(rmvnorm(1, rep(0,N), Sigma, method = "chol"))
  
  return(omega)
}

#' Simulate covariate
#'
#' @param N sample size
#' @param d distance matrix
#' @param mod indicates whether or not to simulate a spatial covariate
#'
#' @return single column covariate matrix
#' @export
#'
#' @examples
sim.covariate <- function(N, d, Range = 20, sig2 = 1, mod){
  if(mod == 'spatial'){
    X <- as.matrix(sim.omega(Range,sig2,d))
  } 
  if(mod == 'nonspatial'){
    X <- as.matrix(runif(N,10,20))
  }
  return(X)
}

#' Simulate spatial count data from a binomial likelohood
#'
#' @param cov.mod 'spatial': covariate is spatially structured; 'nonspatial': covariate is iid Uniform(10,20)
#' @param beta Coefficients for covariates. If a single integer, the model is intercept only and does not perform covariate regression
#' @param Range Spatial range, distance at which spatial correlation is ~ 10%
#' @param sig2 Spatial variance
#' @param seed random seed
#'
#' @return data.frame with locations, covariates, and observations
#' @export
sim.spatial.dat <- function(cov.mod, beta, Range, sig2, seed){
  set.seed(seed)
  grid.xy <- expand.grid(x=seq(0,100,2), y=seq(0,100,2))
  d <- as.matrix(dist(grid.xy, upper = TRUE))
  Omega <- sim.omega(Range, sig2, d)
  X <- matrix(1,nrow(grid.xy),length(beta))
  if(length(beta) > 1){#if true, simulate covariate data based on cov.mod
    for(k in 2:length(beta)){
      X[,k] <- sim.covariate(N=nrow(grid.xy), d = d, mod = cov.mod)
    }
  }
  eta <- X%*%beta + matrix(Omega,ncol=1)
  Lambda <- exp(eta)
  sim <- rpois(nrow(grid.xy), Lambda)
  out.df <- data.frame(loc.x = grid.xy[,1], loc.y = grid.xy[,2], X = X, Y = sim, Omega = Omega)

  return(out.df)
}

#' Sample simulated data
#'
#' @param samp.size size of sample
#' @param dat dataframe or matrix with all data needed for TMB model
#'
#' @return sampled dataset
#' @export
samp.dat <- function(samp.size, dat){
  samp.idx <- sample(1:nrow(dat), samp.size)
  out <- data.frame(dat[samp.idx,])
  return(out)
}

