#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR( cones );
  DATA_VECTOR( dbh ); 
  DATA_MATRIX(wave_non ); 
  PARAMETER_VECTOR( a ); 
  PARAMETER_VECTOR( b ); 
  PARAMETER( log_k ); //to keep k positive
  Type k= exp(log_k);
  Type nll;
  
  vector<Type> a_vec=wave_non*a;
  vector<Type> b_vec=wave_non*b;
  vector<Type> mu=a_vec* pow(dbh, b_vec);
  vector<Type> var=mu+mu*mu/k;
  nll = -sum(dnbinom2(cones, mu, var, true));
  
  return( nll);
}
