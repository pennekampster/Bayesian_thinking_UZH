//spatial Poisson model
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; //namespace with GMRF function
  using namespace Eigen; //namespace with SparseMatrix declaration

  DATA_VECTOR(y);
  DATA_MATRIX(x);
  DATA_IVECTOR(v_i); //location index
  DATA_SPARSE_MATRIX(M0); //sparse distance matrices from INLA
  DATA_SPARSE_MATRIX(M1); //sparse distance matrices from INLA
  DATA_SPARSE_MATRIX(M2); //sparse distance matrices from INLA

  PARAMETER_VECTOR(beta); //coefficients
  PARAMETER(ln_kappa); //spatial correlation decay
  PARAMETER(ln_tau); //spatial precision
  PARAMETER_VECTOR(omega); //spatial random effect
  
  int n = y.size(); //take the size of vector y and pass to integer n

  Type tau = exp(ln_tau);
  Type kappa = exp(ln_kappa);
  // marginal spatial variance
  Type marg_sp_sd = 1/(2*sqrt(M_PI)*kappa*tau);
  // distance at which spatial correlation approx. 10%
  Type Range= sqrt(8)/kappa;

  Type nll = 0.0;

  //Spatial Likelihood
  //Define sparse precision matrix
  SparseMatrix<Type> Q = pow(kappa,4) * M0 + 2 * pow(kappa,2) * M1 + M2; //Lindgren et al. 2011
  //Likelihood of spatial random effects
  nll += GMRF(Q)(omega); // +=: nll = nll + ...
  SIMULATE{ 
    omega = GMRF(Q).simulate();
  }
  //scale omega using precision parameter
  omega = omega/tau;

  vector<Type> eta = x*beta; //matrix multiplication
  vector<Type> lambda(n);
  //Data Likelihood
  for(int i=0; i<n; i++){ //TMB starts counting at 0    
    eta(i) += omega(v_i(i)); //v_i links observations to INLA mesh
    //log link
    lambda(i) = exp(eta(i));
    //Poisson likelihood, the last argument, true, tells TMB to take the log of the likelihood
    nll -= dpois(y(i), lambda(i), true); //sum up the negative log likelihoods
    SIMULATE{
      y = rpois(lambda(i));
    }
  } 

  SIMULATE{
    REPORT(y);
    REPORT(omega);
  }

  REPORT(Range);
  ADREPORT(Range); //ADREPORT() will return standard errors
  REPORT(marg_sp_sd);
  ADREPORT(marg_sp_sd);
  REPORT(lambda);
  REPORT(omega);
  REPORT(nll);

  return nll;

}
