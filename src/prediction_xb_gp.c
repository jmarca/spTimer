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
     int its1, rT1, N1, col, i, j, k, ns, p1;
     its1 = *its;
//     r1 = *r;
//     n1 = *n;
//     rn = r1*n1;
     rT1 = *rT;
     N1 = *N;
     col = *constant;
     ns = *nsite;
     p1 = *p;

//     int *T1; 
//     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
//     for(i=0; i<r1; i++){
//          T1[i] = T[i];
//     }
     
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
       
//       free(T1);     
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
	 int i, j, l, t, r1, col, rT1, n1, ns, p1, N1;
 
	 col = *constant;
	 r1 = *r;
     rT1 = *rT;
     n1 = *n;
//     nn = n1*n1;
     ns = *nsite;
//     nns = n1*ns;
     p1 = *p;
     N1 = *N;
//     valN1 = *valN;

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
    
    double *XB;
    XB = (double *) malloc((size_t)((N1*col)*sizeof(double))); 

    MProd(betap, constant, p, X, N, XB);

    double *s21, *del, *opp, *XB1, *valX1, *valXB1, *oox, *part2;
    s21 = (double *) malloc((size_t)((col)*sizeof(double))); 
    del = (double *) malloc((size_t)((ns*col)*sizeof(double))); 
    opp = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    valX1 = (double *) malloc((size_t)((p1*col)*sizeof(double))); 
    valXB1 = (double *) malloc((size_t)((col)*sizeof(double))); 
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
             for(t=0; t < T1[l]; t++) {             
             extract_alt_uneqT(l, t, n, r, T, rT, op, opp);
             extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
             extract_X41_uneqT(l, t, j, nsite, rT, r, T, p, valX, valX1);
     
//             extract_alt2(l, t, n, rT, T, op, opp); 
//             extract_alt2(l, t, n, rT, T, XB, XB1);
//             extract_X41(l, t, j, nsite, rT, T, p, valX, valX1); // p x 1

             MProd(valX1, constant, p, betap, constant, valXB1);  // 1 x 1 

                   for(i=0; i < n1; i++) {
                   oox[i] = opp[i]-XB1[i];
                   }         
                   xTay(S_eta12c, Si_eta, oox, n, part2);
    	     opre[0] = valXB1[0]+part2[0];		  
             mvrnormal(constant, mu, sig, constant, out);
             mvrnormal(constant, mu, sig_ep, constant, out1);           
             zpred[t+T2[l]+j*(rT1)]=opre[0] + out[0]+ out1[0];
             }
	    }	  
    }

    free(T1); free(T2);
    free(S_eta); free(Si_eta); free(S_eta12); free(S_eta12c); free(det);
    free(XB); free(s21); free(del); free(opp); free(XB1); 
    free(valX1); free(valXB1); free(oox); free(part2); 
    free(opre); free(mu); free(sig); free(out); free(out1);
    return;
}



// Iteration and Prediction of O_lt for all sites and all time points "XB"
// spatial beta
void z_pr_its_gp_sp(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *q, int *N, int *valN, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_betap,
     double *betap, double *betasp, double *X, double *valX, double *Xsp, 
     double *valXsp, double *op, int *constant, double *betapred, double *zpred)
{
     int its1, n1, rT1, N1, col, i, j, k, ns, p1, q1;
     its1 = *its;
//     r1 = *r;
     n1 = *n;
//     T1 = *T;
//     rn = r1*n1;
     rT1 = *rT;
     N1 = *N;
     col = *constant;
     ns = *nsite;
     p1 = *p;
     q1 = *q;
     
     double *phi, *nu, *sig_e, *sig_eta, *sig_beta,  *beta, *betas, *o, *zpr, *bp;

     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_beta = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     betas = (double *) malloc((size_t)((q1*n1)*sizeof(double)));       
     o = (double *) malloc((size_t)((N1*col)*sizeof(double)));       
     zpr = (double *) malloc((size_t)((rT1*col*ns)*sizeof(double)));       
     bp = (double *) malloc((size_t)((ns*q1)*sizeof(double)));

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
         sig_beta[0] = sig_betap[i];
         for(j=0; j < p1; j++) {
         beta[j] = betap[j+i*p1];
         }              
         for(j=0; j< (q1*n1); j++){
         betas[j] = betasp[j+i*q1*n1];
         }         
         for(j=0; j < N1; j++) {
         o[j] = op[j+i*N1];
         }

         z_pr_gp_sp(cov, nsite, n, r, rT, T, p, q, N, valN, d, d12, phi, 
         nu, sig_e, sig_eta, sig_beta, beta, betas, X, valX, Xsp, valXsp, o, 
         constant, bp, zpr);
     
         for(j=0; j< ns*q1; j++){
             betapred[j+i*ns*q1] = bp[j];
         }

         for(k=0; k < ns; k++){       
         for(j=0; j < rT1; j ++){
             zpred[j+k*rT1+i*rT1*ns] = zpr[j+k*rT1];                                                
         }        
         }


       printR (i, its1); 

       } // end of iteration loop
       PutRNGstate();
            
       free(phi); free(nu); free(sig_e); free(sig_eta); free(sig_beta);
       free(beta); free(betas); free(o); free(zpr); free(bp);

       return;
}       



// Prediction for all sites for all time point "XB"
// spatial beta
void z_pr_gp_sp(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *q, int *N, int *valN, double *d, double *d12, double *phip, double *nup, 
     double *sig_ep, double *sig_etap, double *sig_betap, double *betap, double *betasp, 
     double *X, double *valX, double *Xsp, double *valXsp,
     double *op, int *constant, double *betapred, double *zpred)
{
	 int i, j, l, t, r1, T1, col, rT1, n1, ns, p1, q1, N1;
 
	 col = *constant;
	 r1 = *r;
	 T1 = *T;
     rT1 = *rT;
     n1 = *n;
//     nn = n1*n1;
     ns = *nsite;
//     nns = n1*ns;
     p1 = *p;
     q1 = *q;
     N1 = *N;
//     valN1 = *valN;
      
    double *S, *Si, *S_12, *S_12c, *det; 
    S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phip, nup, d, S); // n x n 
    MInv(S, Si, n, det);  // n x n  
    covF(cov, n, nsite, phip, nup, d12, S_12); // ns x n 
    
    double *XB, *XBsp;
    XB = (double *) malloc((size_t)((N1*col)*sizeof(double))); 
    XBsp = (double *) malloc((size_t)((N1*col)*sizeof(double))); 

    MProd(betap, constant, p, X, N, XB); // N x 1
    comb_XB_sp(n, r, T, q, Xsp, betasp, constant, XBsp); // N x 1
    MAdd(XB, N, constant, XBsp, XB);  // N x 1

    free(XBsp);

    double *s21, *del, *opp, *XB1, *valX1, *valXB1, *valXsp1, *valXBsp1, *valXBsp;
    double  *oox, *part2, *betasp0, *opre, *mu, *sig, *sigb, *out, *out1;

    s21 = (double *) malloc((size_t)((col)*sizeof(double))); 
    del = (double *) malloc((size_t)((ns*col)*sizeof(double))); 
    opp = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    valX1 = (double *) malloc((size_t)((p1*col)*sizeof(double))); 
    valXB1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    valXsp1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    valXBsp1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    valXBsp = (double *) malloc((size_t)((col)*sizeof(double))); 
    oox = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    part2 = (double *) malloc((size_t)((col)*sizeof(double))); 

    betasp0 = (double *) malloc((size_t)((col)*sizeof(double))); 
    opre = (double *) malloc((size_t)((col)*sizeof(double))); 
    mu = (double *) malloc((size_t)((col)*sizeof(double))); 
    sig = (double *) malloc((size_t)((col)*sizeof(double)));     
    sigb = (double *) malloc((size_t)((col)*sizeof(double)));         
    out = (double *) malloc((size_t)((col)*sizeof(double))); 
    out1 = (double *) malloc((size_t)((col)*sizeof(double)));                                                   
    
    mu[0] = 0.0;
    for(i=0; i < ns; i++) {
         extn_12(i, n, S_12,S_12c); // n x 1
         xTay(S_12c, Si, S_12c, n, s21); // 1 x 1
         
         if(s21[0] > 1.0){
             s21[0] = 1.0-pow(1,-320);
         }
         if(s21[0] == 1.0){
             s21[0] = 1.0-pow(1,-320);
         }
         sig[0] = sig_etap[0] * (1.0 - s21[0]);
         sigb[0] = sig_betap[0] * (1.0 - s21[0]);

         for(j=0; j<q1; j++){ // for q var
            extract_beta_sp(j, n, betasp, oox); // n x 1
            xTay(S_12c, Si, oox, n, betasp0); // 1 x 1
            mvrnormal(constant, betasp0, sigb, constant, betasp0);
//     Rprintf("   betasp0: %4.4f, \n", betasp0[0]);      
            betapred[j+i*q1] = betasp0[0];
         }

// okay

        for(l=0; l < r1; l++) {
             for(t=0; t < T1; t++) {             
                extract_X41(l, t, i, nsite, rT, T, p, valX, valX1); // p x 1
                MProd(valX1, constant, p, betap, constant, valXB1);  // 1 x 1 
//     Rprintf("   valXB1: %4.4f, valXBsp: %4.4f, \n", valXB1[0], valXBsp[0]);      
                valXBsp[0] = 0.0;

                for(j=0; j<q1; j++){ // for q var
                    betasp0[0] = betapred[j+i*q1]; 
                    extract_X_sp4(t, l, i, j, nsite, r, T, valXsp, valXsp1); // 1 x 1
                    MProd(valXsp1, constant, constant, betasp0, constant, valXBsp1);  // 1 x 1 
                    valXBsp[0] = valXBsp[0] + valXBsp1[0];
                 }

                 extract_alt2(l, t, n, rT, T, op, opp); // n x 1
                 extract_alt2(l, t, n, rT, T, XB, XB1); // n x 1
                 for(j=0; j < n1; j++) {
                   oox[j] = opp[j]-XB1[j];
                 }         
             xTay(S_12c, Si, oox, n, part2); // 1 x 1

//     Rprintf("   valXB1: %4.4f, valXBsp: %4.4f, \n", valXB1[0], valXBsp[0]);      
 
    	     opre[0] = valXB1[0]+valXBsp[0]+part2[0];		  
             mvrnormal(constant, mu, sig, constant, out);
             mvrnormal(constant, mu, sig_ep, constant, out1);           
             zpred[t+l*T1+i*(rT1)]=opre[0] + out[0]+ out1[0];

             }
	    }	  
    }

    free(det);
    free(S); free(Si); free(S_12); free(S_12c); 
    free(XB); free(s21); free(del); free(opp); free(XB1); 
    free(valX1); free(valXB1); free(valXsp1); free(valXBsp1); free(valXBsp);
    free(oox); free(part2); 
    free(betasp0); free(opre); free(mu); free(sig); free(sigb); 
    free(out); free(out1);

    return;
}

// Iteration and Prediction of O_lt for all sites and all time points "XB"
// temporal beta
void z_pr_its_gp_tp(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *u, int *N, int *valN, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_deltap,
     double *sig_op, double *betap, double *rhotp, double *betat0p, double *betatp,  
     double *X, double *valX, double *Xtp, double *valXtp, double *op, 
     int *constant, double *zpred)
{
     int its1, T1, rT1, N1, col, i, j, k, ns, p1, u1;
     its1 = *its;
//     r1 = *r;
//     n1 = *n;
     T1 = *T;
//     rn = r1*n1;
     rT1 = *rT;
     N1 = *N;
     col = *constant;
     ns = *nsite;
     p1 = *p;
     u1 = *u;
     
     double *phi, *nu, *sig_e, *sig_eta, *sig_delta, *sig_o, *beta, *rhot, *betat0, *betat, *o, *zpr;

     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_delta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_o = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     rhot = (double *) malloc((size_t)((u1)*sizeof(double)));       
     betat0 = (double *) malloc((size_t)((u1)*sizeof(double)));       
     betat = (double *) malloc((size_t)((u1*T1)*sizeof(double)));       
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
         sig_delta[0] = sig_deltap[i];
         sig_o[0] = sig_op[i];
         for(j=0; j < p1; j++) {
         beta[j] = betap[j+i*p1];
         }              
         for(j=0; j < u1; j++) {
         rhot[j] = rhotp[j+i*u1];
         }              
         for(j=0; j < u1; j++) {
         betat0[j] = betat0p[j+i*u1];
         }              
         for(j=0; j < u1*T1; j++) {
         betat[j] = betatp[j+i*u1*T1];
         }              
         for(j=0; j < N1; j++) {
         o[j] = op[j+i*N1];
         }

         z_pr_gp_tp(cov, nsite, n, r, rT, T, p, u, N, valN, d, d12, phi, nu, 
         sig_e, sig_eta, sig_delta, sig_o, beta, rhot, betat0, betat, 
         X, valX, Xtp, valXtp, o, constant, zpr);
     
         for(k=0; k < ns; k++){       
         for(j=0; j < rT1; j ++){
             zpred[j+k*rT1+i*rT1*ns] = zpr[j+k*rT1];                                                
         }        
         }

       printR (i, its1); 

       } // end of iteration loop
       PutRNGstate();
            
       free(phi); free(nu); free(sig_e); free(sig_eta); free(sig_delta); 
       free(sig_o); free(beta); free(rhot); free(betat0); free(betat);
       free(o); free(zpr);
       
       return;
}       


// Prediction for all sites for all time point "XB"
// temporal beta, betatp = u x T
void z_pr_gp_tp(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *u, int *N, int *valN, double *d, double *d12, double *phip, double *nup, 
     double *sig_ep, double *sig_etap, double *sig_deltap, double *sig_op,
     double *betap, double *rhotp, double *betat0p, double *betatp, double *X, 
     double *valX, double *Xtp, double *valXtp, double *op, int *constant, 
     double *zpred)
{
	 int i, j, l, t, r1, T1, col, rT1, n1, ns, p1, u1, N1;
 
	 col = *constant;
	 r1 = *r;
	 T1 = *T;
     rT1 = *rT;
     n1 = *n;
//     nn = n1*n1;
     ns = *nsite;
//     nns = n1*ns;
     p1 = *p;
     u1 = *u;
     N1 = *N;
//     valN1 = *valN;
      
    double *S_eta, *Si_eta, *S_eta12, *S_eta12c, *det; 
    S_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_eta12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_eta12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phip, nup, d, S_eta);
    MInv(S_eta, Si_eta, n, det);    
    covF(cov, n, nsite, phip, nup, d12, S_eta12);
    
    double *XB, *XBtp;
    XB = (double *) malloc((size_t)((N1*col)*sizeof(double))); 
    XBtp = (double *) malloc((size_t)((N1*col)*sizeof(double))); 

    MProd(betap, constant, p, X, N, XB); // N x 1
    comb_XB_tp(n, r, T, u, Xtp, betatp, constant, XBtp); // N x 1
    MAdd(XB, N, constant, XBtp, XB);  // N x 1

    free(XBtp);

    double *s21, *del, *opp, *XB1, *valX1, *valXB1, *valXtp1, *betatp1, *valXBtp, *oox, *part2;
    s21 = (double *) malloc((size_t)((col)*sizeof(double))); 
    del = (double *) malloc((size_t)((ns*col)*sizeof(double))); 
    opp = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    valX1 = (double *) malloc((size_t)((p1*col)*sizeof(double))); 
    valXB1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    valXtp1 = (double *) malloc((size_t)((u1)*sizeof(double))); 
    betatp1 = (double *) malloc((size_t)((u1)*sizeof(double))); 
    valXBtp = (double *) malloc((size_t)((col)*sizeof(double))); 
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
             extract_alt2(l, t, n, rT, T, op, opp); // n x 1
             extract_alt2(l, t, n, rT, T, XB, XB1); // n x 1
             extract_X41(l, t, j, nsite, rT, T, p, valX, valX1); // p x 1
             MProd(valX1, constant, p, betap, constant, valXB1);  // 1 x 1 

             extract_X41(l, t, j, nsite, rT, T, u, valXtp, valXtp1); // u x 1
             extract_beta_t(t, T, u, betatp, betatp1); // u x T into u x 1
             MProd(valXtp1, constant, u, betatp1, constant, valXBtp); // 1 x 1

                   for(i=0; i < n1; i++) {
                      oox[i] = opp[i]-XB1[i];
                   }         
                   xTay(S_eta12c, Si_eta, oox, n, part2);
    	     opre[0] = valXB1[0]+ valXBtp[0] + part2[0];		  
             mvrnormal(constant, mu, sig, constant, out);
             mvrnormal(constant, mu, sig_ep, constant, out1);           
             zpred[t+l*T1+j*(rT1)]=opre[0] + out[0]+ out1[0];
             }
	    }	  
    }

    free(S_eta); free(Si_eta); free(S_eta12); free(S_eta12c); free(det);
    free(XB); free(s21); free(del); free(opp); free(XB1); 
    free(valX1); free(valXB1); free(valXtp1); free(betatp1); free(valXBtp);
    free(oox); free(part2); 
    free(opre); free(mu); free(sig); free(out); free(out1);
    
    return;
}


// Iteration and Prediction of O_lt for all sites and all time points "XB"
// spatial and temporal beta
void z_pr_its_gp_sptp(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *q, int *u, int *N, int *valN, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_betap,
     double *sig_deltap, double *sig_op, double *betap, double *betasp,
     double *rhotp, double *betat0p, double *betatp, double *X, double *valX, 
     double *Xsp, double *valXsp, double *Xtp, double *valXtp, double *op, 
     int *constant, double *betapred, double *zpred)
{
     int its1, n1, T1, rT1, N1, col, i, j, k, ns, p1, q1, u1;
     its1 = *its;
//     r1 = *r;
     n1 = *n;
     T1 = *T;
//     rn = r1*n1;
     rT1 = *rT;
     N1 = *N;
     col = *constant;
     ns = *nsite;
     p1 = *p;
     q1 = *q;
     u1 = *u;
     
     double *phi, *nu, *sig_e, *sig_eta, *sig_beta, *sig_delta, *sig_o;
     double *beta, *betas, *rhot, *betat0, *betat, *o, *zpr, *bp;

     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double))); 
     sig_beta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_delta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_o = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     betas = (double *) malloc((size_t)((q1*n1)*sizeof(double)));       
     rhot = (double *) malloc((size_t)((u1)*sizeof(double)));       
     betat0 = (double *) malloc((size_t)((u1)*sizeof(double)));       
     betat = (double *) malloc((size_t)((u1*T1)*sizeof(double)));       
     o = (double *) malloc((size_t)((N1*col)*sizeof(double)));       
     zpr = (double *) malloc((size_t)((rT1*col*ns)*sizeof(double)));       
     bp = (double *) malloc((size_t)((ns*q1)*sizeof(double)));

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
         sig_beta[0] = sig_betap[i];         
         sig_delta[0] = sig_deltap[i];
         sig_o[0] = sig_op[i];
         for(j=0; j < p1; j++) {
         beta[j] = betap[j+i*p1];
         }
         for(j=0; j< (q1*n1); j++){
         betas[j] = betasp[j+i*q1*n1];
         }                       
         for(j=0; j < u1; j++) {
         rhot[j] = rhotp[j+i*u1];
         }              
         for(j=0; j < u1; j++) {
         betat0[j] = betat0p[j+i*u1];
         }              
         for(j=0; j < u1*T1; j++) {
         betat[j] = betatp[j+i*u1*T1];
         }              
         for(j=0; j < N1; j++) {
         o[j] = op[j+i*N1];
         }

         z_pr_gp_sptp(cov, nsite, n, r, rT, T, p, q, u, N, valN, d, d12, phi, nu, 
         sig_e, sig_eta, sig_beta, sig_delta, sig_o, beta, betas, rhot, betat0, betat, 
         X, valX, Xsp, valXsp, Xtp, valXtp, o, constant, bp, zpr);

         for(j=0; j< ns*q1; j++){
             betapred[j+i*ns*q1] = bp[j];
         }

         for(k=0; k < ns; k++){       
         for(j=0; j < rT1; j ++){
             zpred[j+k*rT1+i*rT1*ns] = zpr[j+k*rT1];                                                
         }        
         }

       printR (i, its1); 

       } // end of iteration loop
       PutRNGstate();
            
       free(phi); free(nu); free(sig_e); free(sig_eta); free(sig_beta);
       free(sig_delta); free(sig_o); free(beta); free(betas); free(rhot); 
       free(betat0); free(betat); free(o); free(zpr); free(bp);
       
       return;
}       


// Prediction for all sites for all time point "XB"
// temporal beta, betatp = u x T
// spatial beta, betasp = n x 1
void z_pr_gp_sptp(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *q, int *u, int *N, int *valN, double *d, double *d12, double *phip, 
     double *nup, double *sig_ep, double *sig_etap, double *sig_betap, 
     double *sig_deltap, double *sig_op, double *betap, double *betasp, 
     double *rhotp, double *betat0p, double *betatp, double *X, double *valX, 
     double *Xsp, double *valXsp, double *Xtp, double *valXtp, double *op, 
     int *constant, double *betapred, double *zpred)
{
	 int i, j, l, t, r1, T1, col, rT1, n1, ns, p1, q1, u1, N1;
 
	 col = *constant;
	 r1 = *r;
	 T1 = *T;
     rT1 = *rT;
     n1 = *n;
//     nn = n1*n1;
     ns = *nsite;
//     nns = n1*ns;
     p1 = *p;
     q1 = *q;
     u1 = *u;
     N1 = *N;
//     valN1 = *valN;
      
    double *S, *Si, *S_12, *S_12c, *det; 
    S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phip, nup, d, S);
    MInv(S, Si, n, det);    
    covF(cov, n, nsite, phip, nup, d12, S_12);
    
    double *XB, *XBsp, *XBtp;
    XB = (double *) malloc((size_t)((N1*col)*sizeof(double))); 
    XBtp = (double *) malloc((size_t)((N1*col)*sizeof(double))); 
    XBsp = (double *) malloc((size_t)((N1*col)*sizeof(double))); 

    MProd(betap, constant, p, X, N, XB); // N x 1
    comb_XB_sp(n, r, T, q, Xsp, betasp, constant, XBsp); // N x 1
    MAdd(XB, N, constant, XBsp, XB);  // N x 1
    free(XBsp);
    comb_XB_tp(n, r, T, u, Xtp, betatp, constant, XBtp); // N x 1
    MAdd(XB, N, constant, XBtp, XB);  // N x 1
    free(XBtp);

    double *s21, *del, *opp, *XB1, *valX1, *valXB1, *valXsp1, *valXBsp1, *valXBsp;
    double *valXtp1, *betatp1, *valXBtp, *oox, *part2, *betasp0, *opre, *mu, *sig, *sigb, *out, *out1;

    s21 = (double *) malloc((size_t)((col)*sizeof(double))); 
    del = (double *) malloc((size_t)((ns*col)*sizeof(double))); 
    opp = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    valX1 = (double *) malloc((size_t)((p1*col)*sizeof(double))); 
    valXB1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    valXsp1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    valXBsp1 = (double *) malloc((size_t)((col)*sizeof(double))); 
    valXBsp = (double *) malloc((size_t)((col)*sizeof(double))); 
    valXtp1 = (double *) malloc((size_t)((u1)*sizeof(double))); 
    betatp1 = (double *) malloc((size_t)((u1)*sizeof(double))); 
    valXBtp = (double *) malloc((size_t)((col)*sizeof(double))); 
    oox = (double *) malloc((size_t)((n1*col)*sizeof(double))); 
    part2 = (double *) malloc((size_t)((col)*sizeof(double))); 

    betasp0 = (double *) malloc((size_t)((col)*sizeof(double))); 
    opre = (double *) malloc((size_t)((col)*sizeof(double))); 
    mu = (double *) malloc((size_t)((col)*sizeof(double))); 
    sig = (double *) malloc((size_t)((col)*sizeof(double)));
    sigb = (double *) malloc((size_t)((col)*sizeof(double)));         
    out = (double *) malloc((size_t)((col)*sizeof(double))); 
    out1 = (double *) malloc((size_t)((col)*sizeof(double)));                                                   
    
    mu[0] = 0.0;
    for(j=0; j < ns; j++) {
         extn_12(j, n, S_12,S_12c); // n x 1
         xTay(S_12c, Si, S_12c, n, s21);
         if(s21[0] > 1.0){
             s21[0] = 1.0-pow(1,-320);
         }
         if(s21[0] == 1.0){
             s21[0] = 1.0-pow(1,-320);
         }
         sig[0] = sig_etap[0] * (1.0 - s21[0]);
         sigb[0] = sig_betap[0] * (1.0 - s21[0]);
         for(i=0; i<q1; i++){ // for q var
            extract_beta_sp(i, n, betasp, oox); // n x 1
            xTay(S_12c, Si, oox, n, betasp0); // 1 x 1
            mvrnormal(constant, betasp0, sigb, constant, betasp0);
            betapred[i+j*q1] = betasp0[0];
        }
                  
        for(l=0; l < r1; l++) {
             for(t=0; t < T1; t++) {             
             extract_alt2(l, t, n, rT, T, op, opp); // n x 1
             extract_alt2(l, t, n, rT, T, XB, XB1); // n x 1
             extract_X41(l, t, j, nsite, rT, T, p, valX, valX1); // p x 1
             MProd(valX1, constant, p, betap, constant, valXB1);  // 1 x 1 

             valXBsp[0] = 0.0;
             for(i=0; i<q1; i++){ // for q var
//                    extract_beta_sp(i, n, betasp, oox); // n x 1
//                    xTay(S_12c, Si, oox, n, betasp0); // 1 x 1
//                    mvrnormal(constant, betasp0, sigb, constant, betasp0);
                    betasp0[0] = betapred[i+j*q1]; 
                    extract_X_sp4(t, l, j, i, nsite, r, T, valXsp, valXsp1); // 1 x 1
                    MProd(valXsp1, constant, constant, betasp0, constant, valXBsp1);  // 1 x 1 
                    valXBsp[0] = valXBsp[0] + valXBsp1[0];
             }

             extract_X41(l, t, j, nsite, rT, T, u, valXtp, valXtp1); // u x 1
             extract_beta_t(t, T, u, betatp, betatp1); // u x T into u x 1
             MProd(valXtp1, constant, u, betatp1, constant, valXBtp); // 1 x 1

                   for(i=0; i < n1; i++) {
                   oox[i] = opp[i]-XB1[i];
                   }         
                   xTay(S_12c, Si, oox, n, part2);
    	     opre[0] = valXB1[0]+ valXBsp[0] + valXBtp[0] + part2[0];		  
             mvrnormal(constant, mu, sig, constant, out);
             mvrnormal(constant, mu, sig_ep, constant, out1);           
             zpred[t+l*T1+j*(rT1)]=opre[0] + out[0]+ out1[0];
             }
	    }	  
    }

    free(S); free(Si); free(S_12); free(S_12c); free(det); free(XB); 
    free(s21); free(del); free(opp); free(XB1); free(valX1); free(valXB1); 
    free(valXsp1); free(valXBsp1); free(valXBsp);
    free(valXtp1); free(betatp1); free(valXBtp);
    free(oox); free(part2); free(betasp0);
    free(opre); free(mu); free(sig); free(sigb); free(out); free(out1);
    
    return;
}


////////////////////////////// END /////////////////////////////////////////////
