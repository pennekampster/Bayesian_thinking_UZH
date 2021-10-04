#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() () {
    

    // data
    DATA_IVECTOR(y); //length: n.sites x n.reps
    DATA_IVECTOR(knownOcc); //length: n.sites, 1 if a site has a detection? or is known to be occupied?
    DATA_IVECTOR(nd); //length = n.sites, nd=0: at least one detection at site; nd=1: no detection at site?
    DATA_IVECTOR(J); //Number of reps per site length: n.sites
    
    PARAMETER_VECTOR(logit_psi); //length = n. sites   WHY does this have a length of nsites?
    PARAMETER_VECTOR(logit_p); //length = n.sites x n.reps  WHY does this have a length of nsites*nvisits
    
    vector<Type> psi = invlogit(logit_psi);
    vector<Type> p = invlogit(logit_p);
    
    Type nll = 0;
    int N = knownOcc.size(); //N=number of sites
	//int N = nd.size();
    
    int k = 0; // k counts through all lines in y
    for(int i=0; i<N; i++){
      Type cp = 1.0; // What is cp? Why can we not use p directly?
      for(int j=0; j<J(i); j++){
        cp *= pow(p(k), y(k)) * pow(1-p(k), 1-y(k)); // What does this do?
		 //     when y(k)=1       when y(k)=0
        k++;
      }
      if(knownOcc(i)==1){
        psi(i) = 1.0;
      }
      if(nd(i)==0){
		//psi(i) = 1.0;
        nll -= log(cp*psi(i) );
      }
      if(nd(i)==1){
        nll -= log(cp*psi(i) + (1-psi(i)));
      }
    }
    
    return nll;
  }
