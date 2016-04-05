//** File Name 'prediction_xb_ar.c' **//

#include"main_ar.h"
#include "covariance.h"
#include "common.h"
#include <math.h>
#include "mathematics.h"
#include "randgenerator.h"


// Iteration and Prediction of O_lt for all sites and all time points "XB"
void z_pr_its_ar(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *N, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_l0p, 
     double *rhop, double *betap, double *mu_lp,  
     double *X, double *valX, double *op, int *constant, double *zpred)
{
     int its1, r1, rT1, N1, col, i, j, k, ns, p1;
     its1 = *its;
     r1 = *r;
//     n1 = *n;
     rT1 = *rT;
     N1 = *N;
     col = *constant;
     ns = *nsite;
     p1 = *p;
     
//     unsigned iseed = 44;
//     srand(iseed); 
     
     double *phi, *nu, *sig_e, *sig_eta, *sig_l0, *rho, *beta, *mu_l, *o, *zpr;

     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_l0 = (double *) malloc((size_t)((r1*col)*sizeof(double)));       
     rho = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1)*sizeof(double)));       
     mu_l = (double *) malloc((size_t)((r1)*sizeof(double)));       
     o = (double *) malloc((size_t)((N1)*sizeof(double)));       
     zpr = (double *) malloc((size_t)((rT1*ns)*sizeof(double)));       

     GetRNGstate();                                   
     for(i=0; i < its1; i++) {
         phi[0] = phip[i];
         if(cov[0]==4){
           nu[0] = nup[i];
         }
         else{
           nu[0] = nup[0];
         }
         sig_e[0] = sig_ep[i];     
         sig_eta[0] = sig_etap[i];
         rho[0] = rhop[i];
         for(j=0; j < p1; j++) {
         beta[j] = betap[j+i*p1];
         }              
         for(j=0; j < r1; j++) {
         mu_l[j] = mu_lp[j+i*r1];
         }
         for(j=0; j < r1; j++) {
         sig_l0[j] = sig_l0p[j+i*r1];
         }
         for(j=0; j < N1; j++) {
         o[j] = op[j+i*N1];
         }

         z_pr_ar(cov, nsite, n, r, rT, T, p, N, d, d12, 
         phi, nu, sig_e, sig_eta, sig_l0, rho, beta, mu_l, X, valX,
         o, constant, zpr);
                                              
         for(k=0; k < ns; k++){       
         for(j=0; j < rT1; j ++){
             zpred[j+k*rT1+i*rT1*ns] = zpr[j+k*rT1];                                                
         }        
         }

       printR (i, its1); 

       } // end of iteration loop
       PutRNGstate();
            
       free(phi); free(nu); free(sig_e); free(sig_eta); free(sig_l0); free(rho); 
       free(beta); free(mu_l); free(o); free(zpr);
       return;
}       


// Prediction for all sites for all time point "XB"
void z_pr_ar(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *N, double *d, double *d12, double *phip, double *nup,
     double *sig_ep, double *sig_etap, double *sig_l0p, double *rhop, 
     double *betap, double *mu_lp,  double *X, double *valX,
     double *op, int *constant, double *zpred)
{

	 int i, j, l, t, r1, col, rT1, n1, nn, ns, p1, N1;
 
	 col = *constant;
	 r1 = *r;
     rT1 = *rT;
     n1 = *n;
     nn = n1*n1;
     ns = *nsite;
     p1 = *p;
     N1 = *N;

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);   
      
    double *S_eta, *Si_eta, *S_eta12, *S_eta12c, *det; 
    S_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_eta12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_eta12c = (double *) malloc((size_t)((n1)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phip, nup, d, S_eta);
    MInv(S_eta, Si_eta, n, det);    
    covF(cov, n, nsite, phip, nup, d12, S_eta12);

    double *S, *m1, *O11, *O_l0, *XB;
    S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));     
    m1 = (double *) malloc((size_t)((n1)*sizeof(double)));     
    O11 = (double *) malloc((size_t)((n1)*sizeof(double)));     
    O_l0 = (double *) malloc((size_t)((n1*r1*col)*sizeof(double))); 
    XB = (double *) malloc((size_t)((N1*col)*sizeof(double))); 

    for(l=0; l<r1; l++){
       for(i=0; i<nn; i++){
          S[i] = sig_l0p[l]*Si_eta[i];
       }
       for(i=0; i<n1; i++) {        
          m1[i] = mu_lp[l];
       }
       mvrnormal(constant, m1, S, n, O11);
       for(i=0; i<n1; i++) {        
          O_l0[i+l*n1] = O11[i];
       }      
    }
    MProd(betap, constant, p, X, N, XB);

    double *s21, *sig, *m, *sig_0; 
    s21 = (double *) malloc((size_t)((col)*sizeof(double))); 
    sig = (double *) malloc((size_t)((col)*sizeof(double)));     
    m = (double *) malloc((size_t)((col)*sizeof(double))); 
    sig_0 = (double *) malloc((size_t)((col)*sizeof(double)));                                                   


    double *O1, *opp, *XB1, *valX1, *valXB1;
    O1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    opp = (double *) malloc((size_t)((n1)*sizeof(double))); 
    XB1 = (double *) malloc((size_t)((n1)*sizeof(double))); 
    valX1 = (double *) malloc((size_t)((p1*ns)*sizeof(double))); 
    valXB1 = (double *) malloc((size_t)((ns)*sizeof(double))); 

    double *opp1, *oox, *part2, *out, *out1, *mu, *opre; 
    opp1 = (double *) malloc((size_t)((n1)*sizeof(double))); 
    oox = (double *) malloc((size_t)((n1)*sizeof(double))); 
    part2 = (double *) malloc((size_t)((col)*sizeof(double))); 
    out = (double *) malloc((size_t)((col)*sizeof(double))); 
    out1 = (double *) malloc((size_t)((col)*sizeof(double)));                                                   
    mu = (double *) malloc((size_t)((col)*sizeof(double))); 
    opre = (double *) malloc((size_t)((rT1*col)*sizeof(double))); 
   
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
             t=0;
             m[0] = mu_lp[l];
             sig_0[0] = sig_l0p[l];
             mvrnormal(constant, m, sig_0, constant, O1);
             extract_alt_uneqT(l, t, n, r, T, rT, op, opp);
             extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
             extract_X21_uneqT(l, t, nsite, rT, r, T, p, valX, valX1);
//             extract_alt2(l, t, n, rT, T, op, opp); 
//             extract_alt2(l, t, n, rT, T, XB, XB1);
//             extract_X21(l, t, nsite, rT, T, p, valX, valX1);  // X1 = p x n
             MProd(valX1, nsite, p, betap, constant, valXB1);
                   for(i=0; i < n1; i++){
                       opp1[i] = O_l0[i+l*n1];
                   }      
                   for(i=0; i < n1; i++) {
                       oox[i] = opp[i]-rhop[0]*opp1[i]-XB1[i];
                   }         
                   xTay(S_eta12c, Si_eta, oox, n, part2);
             mu[0] = 0.0;
             mvrnormal(constant, mu, sig, constant, out);
             mu[0] = 0.0;
             mvrnormal(constant, mu, sig_ep, constant, out1);     
             opre[t+T2[l]]=rhop[0]*O1[0]+valXB1[j]+part2[0];
             zpred[t+T2[l]+j*(rT1)]=opre[t+T2[l]]+out[0]+out1[0];             
             //rhop[0]*O1[0]+valXB1[j]+part2[0]+out[0]+out1[0]; 

             for(t=1; t < T1[l]; t++) {             
             extract_alt_uneqT(l, t-1, n, r, T, rT, op, opp1);
             extract_alt_uneqT(l, t, n, r, T, rT, op, opp);
             extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
             extract_X21_uneqT(l, t, nsite, rT, r, T, p, valX, valX1);
//             extract_alt2(l, t-1, n, rT, T, op, opp1);
//             extract_alt2(l, t, n, rT, T, op, opp); 
//             extract_alt2(l, t, n, rT, T, XB, XB1);
//             extract_X21(l, t, nsite, rT, T, p, valX, valX1);  // X1 = p x n
             MProd(valX1, nsite, p, betap, constant, valXB1);  // 1 x n 
                   for(i=0; i < n1; i++) {
                       oox[i] = opp[i]-rhop[0]*opp1[i]-XB1[i];
                   }         
                   xTay(S_eta12c, Si_eta, oox, n, part2);
             mu[0] = 0.0;
             mvrnormal(constant, mu, sig, constant, out);
             mu[0] = 0.0;
             mvrnormal(constant, mu, sig_ep, constant, out1);  
    	     opre[t+T2[l]] = rhop[0]*opre[(t-1)+T2[l]]+valXB1[j]+part2[0];		  
             zpred[t+T2[l]+j*(rT1)]=opre[t+T2[l]] + out[0]+ out1[0];
//             zpred[t+T2[l]+j*(rT1)]=rhop[0]*zpred[t+T2[l]+j*(rT1)]+valXB1[j]+part2[0] + out[0]+ out1[0];
             }
	    }	  
    }

    free(opre);
    free(T1); free(T2);
    free(mu); free(out1); free(out); free(part2); free(oox); free(opp1);
    free(valXB1); free(valX1); free(XB1); free(opp); free(O1); free(sig_0); 
    free(m); free(sig); free(s21); free(XB); free(O_l0); free(O11); free(m1); 
    free(S); free(det); free(S_eta12c); free(S_eta12); free(Si_eta); free(S_eta);    

    return;
}




/*

// Prediction for all sites for all time point "XB"
void z_pr_ar(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *N, double *d, double *d12, double *phip, double *nup,
     double *sig_ep, double *sig_etap, double *sig_l0p, double *rhop, 
     double *betap, double *mu_lp,  double *X, double *valX,
     double *op, int *constant, double *zpred)
{
	 int i, j, l, t, r1, col, rT1, n1, nn, ns, p1, N1;
 
	 col = *constant;
	 r1 = *r;
     rT1 = *rT;
     n1 = *n;
     nn = n1*n1;
     ns = *nsite;
     p1 = *p;
     N1 = *N;

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);   
      
    double *S_eta, *Si_eta, *S_eta12, *S_eta12c, *det; 
    S_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_eta12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_eta12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 


    covF(cov, n, n, phip, nup, d, S_eta);
    MInv(S_eta, Si_eta, n, det);    
    covF(cov, n, nsite, phip, nup, d12, S_eta12);

    double *S, *m1, *O11, *O_l0, *XB;
    S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));     
    m1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
    O11 = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
    O_l0 = (double *) malloc((size_t)((n1*r1*col)*sizeof(double))); 
    XB = (double *) malloc((size_t)((N1*col)*sizeof(double))); 

    for(l=0; l<r1; l++){
       for(i=0; i<nn; i++){
          S[i] = sig_l0p[l]*Si_eta[i];
       }
       for(i=0; i<n1; i++) {        
          m1[i] = mu_lp[l];
       }
       mvrnormal(constant, m1, S, n, O11);
       for(i=0; i<n1; i++) {        
          O_l0[i+l*n1] = O11[i];
       }      
    }
    MProd(betap, constant, p, X, N, XB);


    double *s21, *del, *m, *O1, *opp, *XB1, *valX1, *valXB1, *opp1, *oox, *part2;
    s21 = (double *) malloc((size_t)((col)*sizeof(double))); 
    del = (double *) malloc((size_t)((ns*col)*sizeof(double))); 
    m = (double *) malloc((size_t)((col)*sizeof(double))); 
    O1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    opp = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    valX1 = (double *) malloc((size_t)((p1*ns*col)*sizeof(double))); 
    valXB1 = (double *) malloc((size_t)((ns*col)*sizeof(double))); 
    opp1 = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    oox = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    part2 = (double *) malloc((size_t)((col)*sizeof(double))); 

    double *opre, *mu, *sig, *out, *out1, *sig_0;
    opre = (double *) malloc((size_t)((rT1*col)*sizeof(double))); 
    mu = (double *) malloc((size_t)((col)*sizeof(double))); 
    sig = (double *) malloc((size_t)((col)*sizeof(double)));     
    out = (double *) malloc((size_t)((col)*sizeof(double))); 
    out1 = (double *) malloc((size_t)((col)*sizeof(double)));                                                   
    sig_0 = (double *) malloc((size_t)((col)*sizeof(double)));                                                   
    
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
             t=0;
             m[0] = mu_lp[l];
             sig_0[0] = sig_l0p[l];
             mvrnormal(constant, m, sig_0, constant, O1);
             
             extract_alt_uneqT(l, t, n, r, T, rT, op, opp);
             extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
             extract_X21_uneqT(l, t, nsite, rT, r, T, p, valX, valX1);

//             extract_alt2(l, t, n, rT, T, op, opp); 
//             extract_alt2(l, t, n, rT, T, XB, XB1);
//             extract_X21(l, t, nsite, rT, T, p, valX, valX1);  // X1 = p x n

             MProd(valX1, nsite, p, betap, constant, valXB1);

                   for(i=0; i < n1; i++){
                   opp1[i] = O_l0[i+l*n1];
                   }      
                   for(i=0; i < n1; i++) {
                   oox[i] = opp[i]-rhop[0]*opp1[i]-XB1[i];
                   }         
                   xTay(S_eta12c, Si_eta, oox, n, part2);

       	     opre[t+T2[l]] = rhop[0]*O1[0]+valXB1[j]+part2[0];		          	 
             mvrnormal(constant, mu, sig, constant, out);
             mvrnormal(constant, mu, sig_ep, constant, out1);           
             zpred[t+T2[l]+j*(rT1)]=opre[t+T2[l]]+out[0]+out1[0];

             for(t=1; t < T1[l]; t++) {             

             extract_alt_uneqT(l, t-1, n, r, T, rT, op, opp1);
             extract_alt_uneqT(l, t, n, r, T, rT, op, opp);
             extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
             extract_X21_uneqT(l, t, nsite, rT, r, T, p, valX, valX1);

//             extract_alt2(l, t-1, n, rT, T, op, opp1);
//             extract_alt2(l, t, n, rT, T, op, opp); 
//             extract_alt2(l, t, n, rT, T, XB, XB1);
//             extract_X21(l, t, nsite, rT, T, p, valX, valX1);  // X1 = p x n

             MProd(valX1, nsite, p, betap, constant, valXB1);  // 1 x n 

                   for(i=0; i < n1; i++) {
                   oox[i] = opp[i]-rhop[0]*opp1[i]-XB1[i];
                   }         
                   xTay(S_eta12c, Si_eta, oox, n, part2);
    	     opre[t+T2[l]] = rhop[0]*opre[(t-1)+T2[l]]+valXB1[j]+part2[0];		  
             mvrnormal(constant, mu, sig, constant, out);
             mvrnormal(constant, mu, sig_ep, constant, out1);           
             zpred[t+T2[l]+j*(rT1)]=opre[t+T2[l]] + out[0]+ out1[0];
             }
	    }	  
    }

//    double *S_eta, *Si_eta, *S_eta12, *S_eta12c, *det; 
//    double *S, *m1, *O11, *O_l0, *XB;
//    double *s21, *del, *m, *O1, *opp, *XB1, *valX1, *valXB1, *opp1, *oox, *part2;
//    double *opre, *mu, *sig, *out, *out1, *sig_0;
    
    free(sig_0); free(out1);  

    free(T1); free(T1);
    free(S_eta); free(Si_eta); free(S_eta12); free(S_eta12c); free(det);
    free(S); free(m1); free(O11); free(O_l0); free(XB);
    free(s21); free(del); free(m); free(O1); free(opp); free(XB1); 
    free(valX1); free(valXB1); free(opp1); free(oox); free(part2); 
    free(opre); free(mu); free(sig); free(out); 
    return;
}

*/
////////////////////////////// END /////////////////////////////////////////////
