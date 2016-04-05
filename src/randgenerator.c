//** File Name 'randgenerator.c' **//

#include "randgenerator.h"
#include "mathematics.h"

#ifdef TEST
#endif

/* This generates the random number*/
double drand48 (void)
{
//   srand(22);    
//   return ((double)rand() + 0.5)/((double)RAND_MAX + 1.0); // return some uniform dist number
     return((double)unif_rand());     // calling from Rmath.h
}


// random number for uniform distribution
double myrunif(int constant)
{
       double U;
       U = (drand48 ());
       return U;
}       

void runif_val(int *n, int *constant, double *out)
{
       int i;
       for(i=0; i< *n; i++){
       out[i] = myrunif(*constant);
       }
       return;
}       


/* This defines the random exponential variate*/
double rexpon(double scale)
{
	double U;
	U = -log (drand48 ()) / scale;
	return U;
}

/* This includes the "rexpon" and works with it*/
void rexp_val(int *n, double *scale, double *result)
{
     int i;
     for (i=0; i< *n; i++) {
         result[i]=rexpon(*scale);
     }
     return;
}

/* This defines the random normal variate*/
double rnormal(double mu, double sigma)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
   
	if  (iset == 0) {
		do {
			v1=2.0*drand48()-1.0;
			v2=2.0*drand48()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=mu + sigma*v1*fac;
		iset=1;
		return mu + sigma*v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

/* This includes the "rnormal" and works with it*/
void rnorm_val(int *n, double *m, double *s, double *result)
{
     int i;
      for (i = 0; i < *n; i++) {
       result[i] = rnormal(*m, *s);
      }
     return;
}

/* This defines the random gamma variate*/
double rgammaa (double shape, double rate)
{
//   int i;
   double  b, h, r, g, f, x, r2, d, gamma1=10;

   b = shape - 1.0;
   
   if (b >=0) {
	   h = sqrt (3 * shape - 0.75);
       do {
          do {
             do {
                r = drand48 ();
                g = r - pow(r, 2);
             } while (g <= 0.0);
             f = (r - 0.5) * h /sqrt(g);
   			 x = b + f;
           } while (x <= 0.0);
   			r2 = drand48 ();
   			d = 64 * g * (pow(r2*g, 2));
   			if (d <= 0) {
               gamma1 = x;
               break;
            }
   			if (d*x < (x - 2*f*f)) {
                gamma1 = x;
                 break;
            }
       } while (log(d) > 2*(b*log(x/b) - f));
       gamma1 = x;
       gamma1 = (1 / rate) * gamma1;
   }
   else if (b < 0) {
	   x = rgammaa (shape+1, 1);
	   r = drand48 ();
	   x = x*pow(r, 1/shape);
	   gamma1 = x / rate;
   }
   return gamma1;
}     	

/* This includes the "rgammaa" and works with it*/
void rgamma_val(int *n, double *shape, double *rate, double *result)
{
     int i;
      for (i = 0; i < *n; i++) {
       result[i] = rgammaa(*shape, *rate);
      }
     return;
}

/* This defines the randomly generated inverse gamma distribution */
double rigammaa (double shape, double rate)
{
    double temp;
    	temp = rgammaa (shape, rate);
    	temp = 1 / temp;
    return temp;
} 
  
void rigamma_val(int *n, double *shape, double *rate, double *result)
{
     int i;
      for (i = 0; i < *n; i++) {
       result[i] = rigammaa(*shape, *rate);
      }
     return;
}

/*
// This defines the random gamma variate using GUI library
double rgammaa2 (double shape, double rate)
{
  // set up GSL RNG //
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  // end of GSL setup //
  out=gsl_ran_gamma(r, shape , rate);
  return out;
}     	
*/


/*****************************************/
    
/* Multivariate random normal */
double rnorm_for_multivariate(void)  
{ 
  double x,y,z;
  /* use the system dependent random number generator */
  do 
  { 
      x=2.0*drand48()-1.0;
      y=2.0*drand48()-1.0;
//    x=(random()/2147483648.0)*2.0-1.0;
//    y=(random()/2147483648.0)*2.0-1.0;
    z=x*x+y*y;
  }
  while(z>=1);
  z=x*sqrt(-2.0*log(z)/z);

  return(z);
}
// Cholesky Decompostion
void chol_for_multivariate(double *s, int *n, double *ltrg)
{ 
  int i,j,k, nn;
  double sum,eps;
  nn=*n;

  eps=1.e-5; // tolerance for near singularity
//  ltrg[0*nn+0]=sqrt(s[0*nn+0]);
  ltrg[0]=sqrt( *s );

  for(k=1;k<nn;k++){ 
    for(i=0;i<k;i++){ 
      sum=0.0;
      for(j=0;j<i;j++) {
          sum+=ltrg[i*nn+j]*ltrg[k*nn+j];
      }
      if(fabs(s[k*nn+i]-sum)>eps) {
         ltrg[k*nn+i]=(s[k*nn+i]-sum)/ltrg[i*nn+i];
      }
      else{
         ltrg[k*nn+i]=0.0;
      }
      ltrg[i*nn+k]=0.0;
    }
    for(j=0,sum=0.0;j<k;j++){
        sum+=ltrg[k*nn+j]*ltrg[k*nn+j];
    }
    if(s[k*nn+k]-sum<=0.0)
      { // Rprintf("Error! Not posi-definite matrix!");
      //;
      //exit(-1);
      }
    else {
       ltrg[k*nn+k]=sqrt(s[k*nn+k]-sum);
    }
  }

  return;
}

// Multivariate Normal RV;
// p is the row of mu or for s2 it is either row or col
// y is the output
void mvrnormal(int *n, double *mu, double *s2, int *p, double *y)
{ 
  double rnorm_for_multivariate(void);
  double *lowert,*z0,zz;
  int i,j,k, nn, pp;
  void chol_for_multivariate(double *, int *, double *);

  nn=*n;
  pp=*p;

//   z0=(double *) malloc(pp * sizeof(double));
   z0=(double *) malloc((size_t)((pp)*sizeof(double)));

//  lowert=(double *) malloc((pp*pp) * sizeof(double));
   lowert=(double *) malloc((size_t)((pp*pp) * sizeof(double)));

  chol_for_multivariate(s2, p, lowert);

  #if TEST
      Rprintf("Cholesky decomposition\n");
      for(i=0;i<pp;i++)
      { for(j=0;j<pp;j++) Rprintf("%f ", lowert[i*pp+j]); Rprintf("\n"); }
      Rprintf("\n");
  #endif

  for(i=0;i<nn;i++){ 
      for(j=0;j<pp;j++) {
          z0[j]=norm_rand(); //z0[j]=rnorm_for_multivariate();
      }
      for(j=0;j<pp;j++) { 
          for(zz=0.,k=0;k<=j;k++){
              zz+=lowert[j*pp+k]*z0[k];
          }
      y[i+j*nn]=zz+mu[j];
    }
  }

  free(lowert); free(z0);
}

/******************************************************************************/
