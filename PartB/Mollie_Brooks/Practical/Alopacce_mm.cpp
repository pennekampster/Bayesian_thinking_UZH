#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //INPUT SECTION
  DATA_VECTOR(abund);
  DATA_MATRIX(Xmat); //a single covariate, must be same length as response variable
  
  PARAMETER_VECTOR(coefs);
  PARAMETER( log_k ); //to keep k positive
  Type k= exp(log_k); //dispersion parameter
  
  Type nll = 0;
  
  vector<Type> mu(abund.size()); //storage for expectations
  vector<Type> var(abund.size()); //storage for variance
  
  mu = exp(Xmat*coefs); //calculate expectations
  
  //CALCULATION SECTION
  for(int i=0; i<abund.size(); i++)
  {
    var(i)=mu(i)+mu(i)*mu(i)/k;
    
    nll -=dnbinom2(abund(i), mu(i), var(i), true);
  }

  return nll;

}
