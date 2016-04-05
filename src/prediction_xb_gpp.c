//** File Name 'prediction_xb_gpp.c' **//


#include "main_gpp.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"


////////////////////////////////////////////////////////////

// prediction for spatially varying beta
void z_pr_its_gpp1_sp(int *cov, int *scale, int *its, int *nsite, int *n, 
     int *m, int *r, int *T, int *rT, int *p, int *q, int *nsiterT, 
     double *phi_etap, double *nup, double *dm, double *dnsm, double *wp, 
     double *sig2ep, double *sig2betap, double *betap, double *betasp, 
     double *Xpred, double *Xspred, int *constant, double *betapred,
     double *zpred)
{
     int i, j, its1, nsite1, m1, rT1, p1, q1, col;
     its1 = *its;
     nsite1 = *nsite;
     m1 = *m;
//     r1 =*r;
     rT1 =*rT;
     p1 =*p;
     q1 =*q;
     col =*constant;
     
     double *phi_eta, *nu, *w, *sig2e, *sig2beta, *beta, *betas, *zp, *bp;
     phi_eta = (double *) malloc((size_t)((col)*sizeof(double)));
     nu = (double *) malloc((size_t)((col)*sizeof(double)));     
     w = (double *) malloc((size_t)((rT1*m1)*sizeof(double)));
     sig2e = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2beta = (double *) malloc((size_t)((col)*sizeof(double)));
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     betas = (double *) malloc((size_t)((q1*m1)*sizeof(double)));       
     zp = (double *) malloc((size_t)((nsite1*rT1)*sizeof(double)));
     bp = (double *) malloc((size_t)((nsite1*q1)*sizeof(double)));
          
     GetRNGstate();                                   
     for(i=0; i<its1; i++){

       phi_eta[0] = phi_etap[i];
         if(cov[0]==4){
           nu[0] = nup[i];
         }
         else{
           nu[0] = 0.0;
         }
       for(j=0; j< rT1*m1; j++){
          w[j] = wp[j+i*rT1*m1];
       }    
       sig2e[0] = sig2ep[i];
       sig2beta[0] = sig2betap[i];
       for(j=0; j< p1; j++){
          beta[j] = betap[j+i*p1];
       }    
       for(j=0; j< (q1*m1); j++){
          betas[j] = betasp[j+i*q1*m1];
       }         

       z_pr_gpp1_sp(cov, nsite, n, m, r, T, rT, p, q, nsiterT, phi_eta, nu, 
       dm, dnsm, w, sig2e, sig2beta, beta, betas, Xpred, Xspred, constant, 
       bp, zp);

       for(j=0; j< nsite1*q1; j++){
          betapred[j+i*nsite1*q1] = bp[j];
       }
       for(j=0; j< nsite1*rT1; j++){
           if( *scale == 1){     
           zpred[j+i*nsite1*rT1] = zp[j];
           }
           else if(*scale == 2){
           zpred[j+i*nsite1*rT1] = zp[j]*zp[j];
           }                 
           else if(*scale == 3){
           zpred[j+i*nsite1*rT1] = exp(zp[j]);
           }
           else {
                //;
                //exit(9);
           }                 
//           Rprintf("zpred %2.2f\n", zpred[j]); // ok
       }
     printR(i, its1); 
     } // end of iteration loop
     PutRNGstate();

     free(phi_eta); free(nu); free(w); free(sig2e); 
     free(sig2beta); free(beta); free(betas); free(zp); free(bp);
     return;
}

// approach dnsm x dm     
void z_pr_gpp1_sp(int *cov, int *nsite, int *n, int *m, int *r, int *T, 
     int *rT, int *p, int *q, int *nsiterT, double *phi_etap, double *nup, 
     double *dm, double *dnsm, double *wp, double *sig2ep, double *sig2betap, 
     double *betap, double *betasp, double *Xpred, double *Xspred, 
     int *constant, double *betapred, double *zpred)
{
     int i, j, l, t, m1, nsite1, r1, T1, rT1, q1, col;

     m1 =*m;
     nsite1 =*nsite;
     r1 =*r;
     T1 =*T;
     rT1 =*rT;
     q1 =*q;
     col = *constant;

     double *S_eta, *det, *Snsm, *A, *Aw, *tAw, *XB, *mu, *out;
     S_eta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     Snsm = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     A = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     Aw = (double *) malloc((size_t)((rT1*nsite1)*sizeof(double)));
     tAw = (double *) malloc((size_t)((rT1*nsite1)*sizeof(double)));
     XB = (double *) malloc((size_t)((rT1*nsite1)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));
                        
     covF(cov, m, m, phi_etap, nup, dm, S_eta); // m x m
     covF(cov, m, nsite, phi_etap, nup, dnsm, Snsm); // m x m
     MInv(S_eta, S_eta, m, det); // m x m
     MProd(S_eta, m, m, Snsm, nsite, A); // nsite x m

// wp is a m x rT matrix
     MProd(wp, rT, m, A, nsite, Aw);  // nsite x rT,  
     MTranspose(Aw, rT, nsite, tAw); // rT x nsite
     MProd(betap, constant, p, Xpred, nsiterT, XB); // nsiterT x 1 
// for sp beta
     MProd(betasp, q, m, A, nsite, betapred);  // nsite x q,  

     free(A); free(Aw); 

/*
     double *Snsm_c, *s21, *sigb;
     Snsm_c = (double *) malloc((size_t)((m1)*sizeof(double)));
     s21 = (double *) malloc((size_t)((col)*sizeof(double)));
     sigb = (double *) malloc((size_t)((nsite1)*sizeof(double)));

     for(i=0; i < nsite1; i++){ 
            extn_12(i, m, Snsm, Snsm_c); // m x 1
            xTay(Snsm_c, S_eta, Snsm_c, m, s21);
            if(s21[0] > 1.0){
             s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){
             s21[0] = 1.0-pow(1,-320);
            }
            sigb[i] = sig2betap[0] * (1.0 - s21[0]);
            for(j=0; j < q1; j++){
              out[0] = sigb[i];
              mu[0] = betapred[j*nsite1+i];
              mvrnormal(constant, mu, out, constant, mu);
              betapred[j*nsite1+i] = mu[0];
            }
     }
     free(Snsm_c); free(s21); free(sigb);
*/
     free(S_eta); free(det); free(Snsm); 
     
     double *X1, *Ab1, *XBsp, *XBf;
     X1 = (double *) malloc((size_t)((nsite1*nsite1)*sizeof(double)));
     Ab1 = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     XBsp = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     XBf = (double *) malloc((size_t)((rT1*nsite1)*sizeof(double)));

     mu[0] = 0.0;
     for(l=0; l < r1; l++) {
       for(t=0; t < T1; t++) {             
         for(i=0; i < nsite1; i++){ 
            XBsp[i] = 0.0;
         }
         for(j=0; j<q1; j++){ // for q var
            extract_X_sp2(t, l, j, nsite, r, T, Xspred, X1); // nsite x nsite diagonal mat
            extract_beta_sp(j, nsite, betapred, Ab1); // nsite x 1
            MProd(Ab1, constant, nsite, X1, nsite, Ab1);  // nsite x 1  
            MAdd(XBsp, nsite, constant, Ab1, XBsp);  // nsite x 1
         }
         for(i=0; i < nsite1; i++){ 
            XBf[t+l*T1+i*(rT1)] = XBsp[i];
         }
       } // end of t
     } // end of l

     free(X1); free(Ab1); free(XBsp); 

     mu[0] = 0.0;
     for(i=0; i< nsite1*rT1; i++){
          mvrnormal(constant, mu, sig2ep, constant, out);
          zpred[i] = XB[i] + tAw[i] + XBf[i] + out[0];
     }

     free(tAw); free(XB); free(mu); free(out);  free(XBf);

     return;

}



// the prediction
void z_pr_its_gpp1(int *cov, int *scale, int *its, int *nsite, int *n, int *m, int *r, 
     int *T, int *rT, int *p, int *nsiterT, double *phi_etap, double *nup, double *dm, 
     double *dnsm, double *wp, double *sig2ep, double *betap, double *Xpred, 
     int *constant, double *zpred)
{
     int i, j, its1, nsite1, m1, rT1, p1, col;
     its1 = *its;
     nsite1 = *nsite;
     m1 = *m;
//     r1 =*r;
     rT1 =*rT;
     p1 =*p;
     col =*constant;
     
     double *phi_eta, *nu, *w, *sig2e, *beta, *zp;
     phi_eta = (double *) malloc((size_t)((col)*sizeof(double)));
     nu = (double *) malloc((size_t)((col)*sizeof(double)));     
     w = (double *) malloc((size_t)((rT1*m1)*sizeof(double)));
     sig2e = (double *) malloc((size_t)((col)*sizeof(double)));
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     zp = (double *) malloc((size_t)((nsite1*rT1)*sizeof(double)));
          
     GetRNGstate();                                   
     for(i=0; i<its1; i++){

       phi_eta[0] = phi_etap[i];
         if(cov[0]==4){
           nu[0] = nup[i];
         }
         else{
           nu[0] = 0.0;
         }
       for(j=0; j< rT1*m1; j++){
          w[j] = wp[j+i*rT1*m1];
       }    
       sig2e[0] = sig2ep[i];
       for(j=0; j< p1; j++){
          beta[j] = betap[j+i*p1];
       }    
              
       z_pr_gpp1(cov, nsite, n, m, r, T, rT, p, nsiterT, phi_eta, nu, dm, dnsm, 
       w, sig2e, beta, Xpred, constant, zp);

       for(j=0; j< nsite1*rT1; j++){
           if( *scale == 1){     
           zpred[j+i*nsite1*rT1] = zp[j];
           }
           else if(*scale == 2){
           zpred[j+i*nsite1*rT1] = zp[j]*zp[j];
           }                 
           else if(*scale == 3){
           zpred[j+i*nsite1*rT1] = exp(zp[j]);
           }
           else {
                //;
                //exit(9);
           }                 
//           Rprintf("zpred %2.2f\n", zpred[j]); // ok
       }
     printR(i, its1); 
     } // end of iteration loop
     PutRNGstate();

     free(phi_eta); free(nu); free(w); free(sig2e); free(beta); free(zp);
     return;
}



// approach dnsm x dm     
void z_pr_gpp1(int *cov, int *nsite, int *n, int *m, int *r, int *T, 
     int *rT, int *p, int *nsiterT, double *phi_etap, double *nup, double *dm, 
     double *dnsm, double *wp, double *sig2ep, double *betap, 
     double *Xpred, int *constant, double *zpred)
{
     int i, m1, nsite1, rT1, col;

     m1 =*m;
     nsite1 =*nsite;
//     r1 =*r;
     rT1 =*rT;
     col = *constant;

     double *S_eta, *det, *Snsm, *A, *Aw, *tAw, *XB, *mu, *out;
     S_eta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     Snsm = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));

// slot 1
     A = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     Aw = (double *) malloc((size_t)((rT1*nsite1)*sizeof(double)));
     tAw = (double *) malloc((size_t)((rT1*nsite1)*sizeof(double)));
     XB = (double *) malloc((size_t)((rT1*nsite1)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));

    covF(cov, m, m, phi_etap, nup, dm, S_eta);
    covF(cov, m, nsite, phi_etap, nup, dnsm, Snsm);
    MInv(S_eta, S_eta, m, det);
    MProd(S_eta, m, m, Snsm, nsite, A); // nsite x m

// wp is a m x rT matrix
     MProd(wp, rT, m, A, nsite, Aw);  // nsite x rT,  
     MTranspose(Aw, rT, nsite, tAw); // rT x nsite
     MProd(betap, constant, p, Xpred, nsiterT, XB); // nsiterT x 1 
     
     mu[0] = 0.0;
     for(i=0; i< nsite1*rT1; i++){
          mvrnormal(constant, mu, sig2ep, constant, out);
          zpred[i] = XB[i] + tAw[i] + out[0];
     }
                            
     free(S_eta); free(det); free(Snsm); free(A); 
     free(Aw); free(tAw); free(XB); free(mu); free(out);

     return;
}

   
     

/////////////////////// the end //////////////////////////
