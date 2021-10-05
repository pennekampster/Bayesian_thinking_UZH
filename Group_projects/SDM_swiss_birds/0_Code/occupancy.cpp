#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() () {
    

    // data
    DATA_IVECTOR(y); //length: n.sites x n.reps
    //DATA_IVECTOR(knownOcc); //length: n.sites, 1 if a site has a detection? or is known to be occupied?
    DATA_IVECTOR(nd); //length = n.sites, nd=0: at least one detection at site; nd=1: no detection at site?
    DATA_IVECTOR(J); //Number of reps per site length: n.sites
	DATA_MATRIX(X); //matrix of site covariates
    DATA_MATRIX(Z); //matrix of visit covariates

    //PARAMETER_VECTOR(logit_psi(i)); //length = n. sites   WHY does this have a length of nsites?
    //PARAMETER_VECTOR(logit_p); //length = n.sites x n.reps  WHY does this have a length of nsites*nvisits
    //PARAMETER(logit_psi(i)); // length=1
	//PARAMETER(logit_p); // length=1
	PARAMETER_VECTOR(beta); //parameter vector of length #covariates + 1
	PARAMETER_VECTOR(alpha);
	vector<Type> Xbeta = X*beta;
	vector<Type> Zalpha = Z*alpha;
	
    vector<Type> psi = invlogit(Xbeta);
	vector<Type> p = invlogit(Zalpha);
    //Type p = invlogit(logit_p);
    
    Type nll = 0;
    //int N = knownOcc.size(); //N=number of sites
	int N = nd.size();
    
    int k = 0; // k counts through all lines in y
    for(int i=0; i<N; i++){
      Type cp = 1.0; // cp is a derived variable/help variable to calculate p
      for(int j=0; j<J(i); j++){
        cp *= pow(p(k), y(k)) * pow(1-p(k), 1-y(k)); // What does this do?
		 //     when y(k)=1       when y(k)=0
        k++;
      }
      //if(knownOcc(i)==1){
      //  psi(i) = 1.0;
      //}
      if(nd(i)==0){
		//psi(i) = 1.0;
        nll -= log(cp*psi(i));
      }
      if(nd(i)==1){
        nll -= log(cp*psi(i) + (1-psi(i)));
      }
    }
	
	Type mean_psi = psi.sum()/N;
    Type mean_p = p.sum()/y.size();
	REPORT(mean_p);
	REPORT(mean_psi);
    return nll;
  }
