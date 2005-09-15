#include "header.h"

void rkappa(int *n, double *loc, double *scale, double *shape1,
	    double *shape2, double *sim){

  int i;

  RANDIN;
  if ( *scale <=0)
    error("Invalid scale parameter !");
  if ( *shape2 == -1){
    /*This is the Generalized Logistic distribution*/
    for (i = 0; i < *n; i++)
      sim[i] = *loc - *scale / *shape1 * 
	(1 - R_pow(1/UNIF - 1, - *shape1));
  }
  
  if ( *shape2 == 0){
    /*This is the Generalized Extreme Value Distribution*/
    if ( *shape1 == 0) 
      for (i = 0; i < *n; i++)
	sim[i] = *loc - *scale * log( EXP );
    else
      for (i = 0; i<*n; i++)
	sim[i] = *loc + *scale * (R_pow(EXP, - *shape1) - 1) / *shape1;
  }
  else if ( *shape2 == 1 ){
    /*This is the Generalized Pareto Distribution*/
    if ( *shape1 == 0)
      for (i = 0; i <*n; i++)
	sim[i] = *loc + *scale * EXP;
    else
      for (i = 0; i <*n; i++)
	sim[i] = *loc + *scale * (R_pow(UNIF, - *shape1) - 1) / *shape1;
    }
  
  else{
    /*This is the Four-parameter Kappa distribution*/
    for (i = 0; i <*n; i++)
      sim[i] = *loc - *scale / *shape1 * 
	(1 - R_pow( (1 - R_pow(UNIF, *shape2)) / *shape2, 
	 - *shape1) );
  }
  RANDOUT;
}
					    

  
      
