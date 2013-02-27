//** File Name 'prediction_xb_gp.c' **//

#include"main_gp.h"
#include "covariance.h"
#include "common.h"
#include <math.h>
#include "mathematics.h"
#include "randgenerator.h"


// Iteration and Prediction of O_lt for all sites and all time points "XB"
void z_pr_its_gp(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *N, int *valN, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *betap,   
     double *X, double *valX, double *op, int *constant, double *zpred)
{
     int its1, r1, rn, n1, T1, rT1, N1, col, i, j, k, ns, p1;
     its1 = *its;
     r1 = *r;
     n1 = *n;
     T1 = *T;
     rn = r1*n1;
     rT1 = *rT;
     N1 = *N;
     col = *constant;
     ns = *nsite;
     p1 = *p;
     
//     unsigned iseed = 44;
//     srand(iseed); 
     
     double *phi, *nu, *sig_e, *sig_eta, *beta, *o, *zpr;

     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     o = (double *) malloc((size_t)((N1*col)*sizeof(double)));       
     zpr = (double *) malloc((size_t)((rT1*col*ns)*sizeof(double)));       

     GetRNGstate();                                   
     for(i=0; i < its1; i++) {
         phi[0] = phip[i];
         if(cov[0]==4){
           nu[0] = nup[i];
         }
         else{
           nu[0] = 0.0;
         }
         sig_e[0] = sig_ep[i];     
         sig_eta[0] = sig_etap[i];
         for(j=0; j < p1; j++) {
         beta[j] = betap[j+i*p1];
         }              
         for(j=0; j < N1; j++) {
         o[j] = op[j+i*N1];
         }

         z_pr_gp(cov, nsite, n, r, rT, T, p, N, valN, d, d12, phi, nu, 
         sig_e, sig_eta, beta, X, valX, o, constant, zpr);
     
         for(k=0; k < ns; k++){       
         for(j=0; j < rT1; j ++){
             zpred[j+k*rT1+i*rT1*ns] = zpr[j+k*rT1];                                                
         }        
         }

       printR (i, its1); 

       } // end of iteration loop
       PutRNGstate();
            
       free(phi); free(nu); free(sig_e); 
       free(sig_eta); free(beta); free(o); free(zpr);
       return;
}       


// Prediction for all sites for all time point "XB"
void z_pr_gp(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *N, int *valN, double *d, double *d12, double *phip, double *nup, 
     double *sig_ep, double *sig_etap, double *betap, double *X, double *valX,
     double *op, int *constant, double *zpred)
{
	 int i, j, l, t, r1, T1, col, rT1, n1, nn, ns, nns, p1, N1, valN1;
 
	 col = *constant;
	 r1 = *r;
	 T1 = *T;
     rT1 = *rT;
     n1 = *n;
     nn = n1*n1;
     ns = *nsite;
     nns = n1*ns;
     p1 = *p;
     N1 = *N;
     valN1 = *valN;
      
    double *S_eta, *Si_eta, *S_eta12, *S_eta12c, *det; 
    S_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_eta12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_eta12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phip, nup, d, S_eta);
    MInv(S_eta, Si_eta, n, det);    
    covF(cov, n, nsite, phip, nup, d12, S_eta12);
    
/*        
      // exponential covariance
      if(cov[0] == 1){
        for(i = 0; i < (n1*n1); i++){
          S_eta[i] = exp(-1.0*phip[0]*d[i]);
        }
        MInv(S_eta, Si_eta, n, det);
        for(i=0; i < nns; i++){
          S_eta12[i] = exp(-phip[0]*d12[i]);        
        }
      }
      // gaussian covariance
      if(cov[0] == 2){
        for(i = 0; i < (n1*n1); i++){
          S_eta[i] = exp(-1.0*phip[0]*phip[0]*d[i]*d[i]);
        }
        MInv(S_eta, Si_eta, n, det);
        for(i=0; i < nns; i++){
          S_eta12[i] = exp(-1.0*phip[0]*phip[0]*d12[i]*d12[i]);        
        }
      }
      // spherical covariance
      if(cov[0] == 3){
        for(i = 0; i < (n1*n1); i++){
         if(d[i] > 0 && d[i] <= 1.0/phip[0]){
         S_eta[i] = 1.0-1.5*phip[0]*d[i]+0.5*(phip[0]*d[i])*(phip[0]*d[i])*(phip[0]*d[i]);
         }
         else if(d[i] >= 1.0/phip[0]){
         S_eta[i] = 0.0;
         }
         else{
         S_eta[i] = 1.0;        
         }        
        }
        MInv(S_eta, Si_eta, n, det);
        for(i=0; i < nns; i++){
         if(d12[i] > 0 && d12[i] <= 1.0/phip[0]){
         S_eta12[i] = 1.0-1.5*phip[0]*d12[i]+0.5*(phip[0]*d12[i])*(phip[0]*d12[i])*(phip[0]*d12[i]);
         }
         else if(d12[i] >= 1.0/phip[0]){
         S_eta12[i] = 0.0;
         }
         else{
         S_eta12[i] = 1.0;        
         }        
        }
      }
      // matern covariance, nu = 3/2
      if(cov[0] == 4){
        for(i = 0; i < (n1*n1); i++){
          S_eta[i] = (1.0+phip[0]*d[i])*exp(-1.0*phip[0]*d[i]);
        }
        MInv(S_eta, Si_eta, n, det);
        for(i=0; i < nns; i++){
          S_eta12[i] = (1.0+phip[0]*d12[i])*exp(-1.0*phip[0]*d12[i]);        
        }
      }
*/

    double *XB;
    XB = (double *) malloc((size_t)((N1*col)*sizeof(double))); 

    MProd(betap, constant, p, X, N, XB);


    double *s21, *del, *opp, *XB1, *valX1, *valXB1, *oox, *part2;
    s21 = (double *) malloc((size_t)((col)*sizeof(double))); 
    del = (double *) malloc((size_t)((ns*col)*sizeof(double))); 
    opp = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    valX1 = (double *) malloc((size_t)((p1*ns*col)*sizeof(double))); 
    valXB1 = (double *) malloc((size_t)((ns*col)*sizeof(double))); 
    oox = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    part2 = (double *) malloc((size_t)((col)*sizeof(double))); 

    double *opre, *mu, *sig, *out, *out1;
    opre = (double *) malloc((size_t)((col)*sizeof(double))); 
    mu = (double *) malloc((size_t)((col)*sizeof(double))); 
    sig = (double *) malloc((size_t)((col)*sizeof(double)));     
    out = (double *) malloc((size_t)((col)*sizeof(double))); 
    out1 = (double *) malloc((size_t)((col)*sizeof(double)));                                                   
    
    mu[0] = 0.0;
    for(j=0; j < ns; j++) {
         extn_12(j, n, S_eta12,S_eta12c);
         
         xTay(S_eta12c, Si_eta, S_eta12c, n, s21);
         if(s21[0] > 1.0){
             s21[0] = 1.0-pow(1,-320);
         }
         if(s21[0] == 1.0){
             s21[0] = 1.0-pow(1,-320);
         }
         sig[0] = sig_etap[0] * (1.0 - s21[0]);
         
        for(l=0; l < r1; l++) {
             for(t=0; t < T1; t++) {             
             extract_alt2(l, t, n, rT, T, op, opp); 
             extract_alt2(l, t, n, rT, T, XB, XB1);
             extract_X21(l, t, nsite, rT, T, p, valX, valX1);  // X1 = p x n
             MProd(valX1, nsite, p, betap, constant, valXB1);  // 1 x n 
                   for(i=0; i < n1; i++) {
                   oox[i] = opp[i]-XB1[i];
                   }         
                   xTay(S_eta12c, Si_eta, oox, n, part2);
    	     opre[0] = valXB1[j]+part2[0];		  
             mvrnormal(constant, mu, sig, constant, out);
             mvrnormal(constant, mu, sig_ep, constant, out1);           
             zpred[t+l*T1+j*(rT1)]=opre[0] + out[0]+ out1[0];
             }
	    }	  
    }

    free(S_eta); free(Si_eta); free(S_eta12); free(S_eta12c); free(det);
    free(XB); free(s21); free(del); free(opp); free(XB1); 
    free(valX1); free(valXB1); free(oox); free(part2); 
    free(opre); free(mu); free(sig); free(out); free(out1);
    return;
}




////////////////////////////// END /////////////////////////////////////////////
