#include <TMB.hpp>                                

template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(y);
	DATA_MATRIX(X);                             
	DATA_MATRIX(Z);
	DATA_INTEGER(ntraps);
	PARAMETER(log_re_sd);
	Type re_sd=exp(log_re_sd);
	ADREPORT(re_sd);
	PARAMETER(log_resid_sd);
	Type resid_sd=exp(log_resid_sd);
	ADREPORT(resid_sd);
	PARAMETER_VECTOR(beta);
	PARAMETER_VECTOR(dev);
	

	//CALCULATE NLL
	Type nll = 0;       
	vector<Type> Xbeta= X*beta;
	vector<Type> Zdev= Z*dev;
	for(int i=0; i<y.size(); i++)
	{
		nll-= dnorm(y(i), Xbeta(i)+Zdev(i), resid_sd, true);
	}
	for(int i=0; i<dev.size(); i++)
	{
		nll-= dnorm(dev(i), Type(0), re_sd, true);
	}           
	return nll;
}
