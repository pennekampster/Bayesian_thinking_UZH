#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() () {
    

    // data
    DATA_IVECTOR(y); //length=nobservations
    DATA_IVECTOR(nd); //length=nsites; nd=1=never detected, nd=0=detected at least once at this site
    DATA_IVECTOR(J); //length=nsites; number of visits per site

    PARAMETER(logit_psi); // length=1
	PARAMETER(logit_p); // length=1
	Type psi = invlogit(logit_psi);
	Type p = invlogit(logit_p);
    
    Type nll = 0;
	
	int N = nd.size();
    
    int k = 0; // k counts through all lines in y
    for(int i=0; i<N; i++){ // loop over sites
      Type cp = 1.0; // cp is a derived variable/help variable to calculate p
      for(int j=0; j<J(i); j++){ // loop over visits
        cp *= pow(p, y(k)) * pow(1-p, 1-y(k)); // calculates probability/product of detection history
		 //     when y(k)=1       when y(k)=0
        k++;
      }
      if(nd(i)==0){
        nll -= log(cp*psi);
      }
      if(nd(i)==1){
        nll -= log(cp*psi + (1-psi));
      }
    }
	
	REPORT(psi);
	REPORT(p);
    return nll;
  }
