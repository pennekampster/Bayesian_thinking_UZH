#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
 //Input section
  DATA_VECTOR(x);
  DATA_VECTOR(y);
    
  PARAMETER(a); //intercept
  PARAMETER(b); //slope
  PARAMETER(log_sigma); //log std dev of residuals
  
  Type nll = 0; //negative log likelihood
  
  Type sigma = exp(log_sigma); //back transform
  int n = y.size();
  vector<Type> yfit(n); //allocate storage vector of length n
  
  //Calculation section
  yfit = a + b * x;
  for(int i=0; i<n; i++)
  {
    nll -= dnorm(y[i], yfit[i], sigma, true);
  }
  //Output section
  ADREPORT(sigma); //use ADREPORT to output derived parameters
  REPORT(yfit);
  
  return nll; //always return the objective function last
}
