//** File Name 'forecast_xb_gp.c' **//

#include"main_gp.h"
#include "covariance.h"
#include "common.h"
#include "math.h"
#include "mathematics.h"
#include "randgenerator.h"



// d is the distance for actual locations n x n
// d12 is the distance between pred and observaed locations nsite x n
// nrK = nsite x r x k
void zlt_fore_gp_its(int *cov, int *its, int *K, int *nsite, int *n, int *r, 
     int *p, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *foreX, 
     double *betap, double *wp, int *constant, double *foreZ)
{
     int i, j, its1, n1, ns, r1, K1, p1, col;
     its1 = *its;
     n1 =*n;
     ns =*nsite;
     r1 =*r;
     K1 =*K;
     p1 =*p;
     col =*constant;     

//     unsigned iseed = 44;
//     srand(iseed); 

     double *phi, *nu, *sig_e, *sig_eta, *beta, *w, *fZ;
     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     w = (double *) malloc((size_t)((n1)*sizeof(double)));            
     fZ = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     
     GetRNGstate();     
     for(i=0; i<its1; i++){
     phi[0] = phip[i];
         if(cov[0]==4){
           nu[0] = nup[i];
         }
         else{
           nu[0] = 0.0;
         }
     sig_e[0] = sig_ep[i];
     sig_eta[0] = sig_etap[i];
     for(j=0; j<p1; j++){
        beta[j] = betap[j+i*p1];
     }
     //for(j=0; j<n1; j++){
     //   w[j] = wp[j+i*n1];
     //}
                   
     zlt_fore_gp(cov, K, nsite, n, r, p, rT, T, rK, nrK, d, d12, 
     phi, nu, sig_e, sig_eta, foreX, beta, wp, constant, fZ);
     
     for(j=0; j < ns*r1*K1; j++){       
         foreZ[j+i*ns*r1*K1] = fZ[j];                                                
     }
     printR(i, its1); 
     } // end of iteration loop
     PutRNGstate();
     
     free(phi); free(nu); free(sig_e); free(sig_eta);
     free(beta); free(w); free(fZ);
     
     return;
}

     

// K-step Forecasts without its
void zlt_fore_gp(int *cov, int *K, int *nsite, int *n, int *r, int *p, 
     int *rT, int *T, int *rK, int *nrK, double *d, double *d12, 
     double *phi, double *nu, double *sig_e, double *sig_eta, double *foreX, 
     double *beta, double *w, int *constant, double *foreZ)
{ 
     int l, k, i, K1, r1, n1, ns, col;
     K1 =*K;
     r1 =*r;
     n1 =*n;
     ns =*nsite;
//     nns =n1*ns;
     col =*constant;

    double *S_eta, *Si_eta, *S_eta12, *S_eta12c, *det; 
    S_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_eta12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_eta12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phi, nu, d, S_eta);
    MInv(S_eta, Si_eta, n, det);    
    covF(cov, n, nsite, phi, nu, d12, S_eta12);

     double *mu, *sig, *s21, *XB, *XB1, *eta, *eps, *zfore;
     mu = (double *) malloc((size_t)((col)*sizeof(double))); 
     sig = (double *) malloc((size_t)((col)*sizeof(double)));            
     s21 = (double *) malloc((size_t)((col)*sizeof(double)));      
     XB = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     XB1 = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     eps = (double *) malloc((size_t)((col)*sizeof(double)));       
     zfore = (double *) malloc((size_t)((ns*col)*sizeof(double)));       

    
     MProd(beta, constant, p, foreX, nrK, XB);  // nsiterK x 1     
     for(l=0; l<r1; l++){
       for(k=0; k<1; k++){     
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
//         mvrnormal(constant, mu, sig_e, constant, eps); 
         for(i=0; i<ns; i++){

            extn_12(i, n, S_eta12,S_eta12c); // n x 1
            xTay(S_eta12c, Si_eta, w, n, mu); // 1 x 1 for mean
            xTay(S_eta12c, Si_eta, S_eta12c, n, s21); // 1 x 1
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            mvrnormal(constant, mu, sig, constant, eta); 
            mvrnormal(constant, mu, sig_e, constant, eps); 
            zfore[i] = XB1[i]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
       for(k=1; k<K1; k++){     
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         for(i=0; i<ns; i++){
            extn_12(i, n, S_eta12,S_eta12c);
            xTay(S_eta12c, Si_eta, w, n, mu); // 1 x 1 for mean
            xTay(S_eta12c, Si_eta, S_eta12c, n, s21);
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            mvrnormal(constant, mu, sig, constant, eta); 
            mvrnormal(constant, mu, sig_e, constant, eps);             
            zfore[i] = XB1[i]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
     }

//     free(T1);
     free(S_eta); free(det); 
     free(Si_eta), free(S_eta12), free(S_eta12c); 
     free(mu); free(sig); free(s21); free(XB); free(XB1); 
     free(eta); free(eps); free(zfore);

     return;
}

/*
// for spatially varying beta
// d is the distance for actual locations n x n
// d12 is the distance between pred and observaed locations nsite x n
// nrK = nsite x r x k
void zlt_fore_gp_sp_its(int *cov, int *its, int *K, int *nsite, int *n, int *r, 
     int *p, int *q, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_betap,
     double *foreX, double *foreXsp, double *betap, double *betasp, double *wp, 
     int *constant, double *foreZ)
{
     int i, j, its1, n1, ns, r1, K1, p1, q1, col;
     its1 = *its;
     n1 =*n;
     ns =*nsite;
     r1 =*r;
//     T1 =*T;
     K1 =*K;
     p1 =*p;
     q1 =*q;
     col =*constant;     

     double *phi, *nu, *sig_e, *sig_eta, *sig_beta, *beta, *betas, *fZ;
     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_beta = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     betas = (double *) malloc((size_t)((q1*n1*col)*sizeof(double)));       
     fZ = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     
     GetRNGstate();     
     for(i=0; i<its1; i++){
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
     for(j=0; j<p1; j++){
        beta[j] = betap[j+i*p1];
     }
     for(j=0; j< (q1*n1); j++){
        betas[j] = betasp[j+i*q1*n1];
     }         
                   
     zlt_fore_gp_sp(cov, K, nsite, n, r, p, q, rT, T, rK, nrK, d, d12, 
     phi, nu, sig_e, sig_eta, sig_beta, foreX, foreXsp, beta, betas,
     wp, constant, fZ);

     for(j=0; j < ns*r1*K1; j++){       
         foreZ[j+i*ns*r1*K1] = fZ[j];                                                
     }
     printR(i, its1); 
     } // end of iteration loop
     PutRNGstate();
     
     free(phi); free(nu); free(sig_e); free(sig_eta); free(sig_beta);
     free(beta); free(betas); free(fZ);
     
     return;
}
*/

/*
// K-step Forecasts without its
// for spatial beta
// nrK = nsite*r*K
void zlt_fore_gp_sp(int *cov, int *K, int *nsite, int *n, int *r, int *p, int *q,
     int *rT, int *T, int *rK, int *nrK, double *d, double *d12, 
     double *phi, double *nu, double *sig_e, double *sig_eta, double *sig_beta, 
     double *foreX, double *foreXsp, double *beta, double *betas, double *w, int *constant, 
     double *foreZ)
{ 
     int l, k, i, t, j, T1, K1, r1, n1, ns, col;
     T1 =*T;
     K1 =*K;
     r1 =*r;
     n1 =*n;
     ns =*nsite;
//     nns =n1*ns;
     col =*constant;

    double *S, *Si, *S_12, *S_12c, *det; 
    S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phi, nu, d, S);
    MInv(S, Si, n, det);    
    covF(cov, n, nsite, phi, nu, d12, S_12);

     double *mu, *sig, *sigb, *s21, *XB, *XB1, *eta, *eps, *zfore;
     mu = (double *) malloc((size_t)((col)*sizeof(double))); 
     sig = (double *) malloc((size_t)((col)*sizeof(double)));            
     sigb = (double *) malloc((size_t)((col)*sizeof(double)));            
     s21 = (double *) malloc((size_t)((col)*sizeof(double)));      
     XB = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     XB1 = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     eps = (double *) malloc((size_t)((col)*sizeof(double)));       
     zfore = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     
     double *oox, *betasp0, * foreXsp1, *foreXBsp, *foreXBsp1;
     oox = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     betasp0 = (double *) malloc((size_t)((col)*sizeof(double)));       
     foreXsp1 = (double *) malloc((size_t)((col)*sizeof(double)));       
     foreXBsp = (double *) malloc((size_t)((col)*sizeof(double)));       
     foreXBsp1 = (double *) malloc((size_t)((col)*sizeof(double)));       
         
     MProd(beta, constant, p, foreX, nrK, XB);  // nsiterK x 1     
     
     for(l=0; l<r1; l++){
       for(k=0; k<1; k++){     
         t = (T1-1);
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         mvrnormal(constant, mu, sig_e, constant, eps); 
         for(i=0; i<ns; i++){
            extn_12(i, n, S_12,S_12c); // n x 1
            xTay(S_12c, Si, w, n, mu); // 1 x 1 for mean
            xTay(S_12c, Si, S_12c, n, s21); // 1 x 1
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            sigb[0] = sig_beta[0] * (1.0 - s21[0]);
            foreXBsp[0] = 0.0;

            for(j=0; j< *q; j++){ // for q var
                    extract_beta_sp(j, n, betas, oox); // n x 1
                    xTay(S_12c, Si, oox, n, betasp0); // 1 x 1
//                    mvrnormal(constant, betasp0, sigb, constant, betasp0);
                    extract_X_sp4(k, l, i, j, nsite, r, K, foreXsp, foreXsp1); // 1 x 1
                    MProd(foreXsp1, constant, constant, betasp0, constant, foreXBsp1);  // 1 x 1 
                    foreXBsp[0] = foreXBsp[0] + foreXBsp1[0];
            }

            mvrnormal(constant, mu, sig, constant, eta); 
            mvrnormal(constant, mu, sig_e, constant, eps); 
            zfore[i] = XB1[i]+foreXBsp[0]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
       for(k=1; k<K1; k++){     
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         for(i=0; i<ns; i++){
            extn_12(i, n, S_12,S_12c);
            xTay(S_12c, Si, w, n, mu); // 1 x 1 for mean
            xTay(S_12c, Si, S_12c, n, s21);
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            sigb[0] = sig_beta[0] * (1.0 - s21[0]);
            foreXBsp[0] = 0.0;
          
            for(j=0; j< *q; j++){ // for q var
                    extract_beta_sp(j, n, betas, oox); // n x 1
                    xTay(S_12c, Si, oox, n, betasp0); // 1 x 1
//                    mvrnormal(constant, betasp0, sigb, constant, betasp0);
                    extract_X_sp4(k, l, i, j, nsite, r, K, foreXsp, foreXsp1); // 1 x 1
                    MProd(foreXsp1, constant, constant, betasp0, constant, foreXBsp1);  // 1 x 1 
                    foreXBsp[0] = foreXBsp[0] + foreXBsp1[0];
            }

            mvrnormal(constant, mu, sig, constant, eta); 
            mvrnormal(constant, mu, sig_e, constant, eps);             
            zfore[i] = XB1[i]+foreXBsp[0]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
     }

     free(S); free(Si); free(S_12); free(S_12c); free(det); 
     free(mu); free(sig); free(sigb); free(s21); free(XB);
     free(XB1); free(eta); free(eps); free(zfore); free(oox);
     free(betasp0); free(foreXsp1); free(foreXBsp); free(foreXBsp1);
     
     return;
}
*/

/*
// for temporally varying beta
// d is the distance for actual locations n x n
// d12 is the distance between pred and observaed locations nsite x n
// nrK = nsite x r x k
void zlt_fore_gp_tp_its(int *cov, int *its, int *K, int *nsite, int *n, int *r, 
     int *p, int *u, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_deltap,
     double *sig_op, double *foreX, double *foreXtp, double *betap, double *rhotp,
     double *betat0p, double *betatp, double *wp, int *constant, double *foreZ)
{
     int i, j, its1, ns, r1, T1, K1, p1, u1, col;
     its1 = *its;
//     n1 =*n;
     ns =*nsite;
     r1 =*r;
     T1 =*T;
     K1 =*K;
     p1 =*p;
     u1 =*u;
     col =*constant;     

     double *phi, *nu, *sig_e, *sig_eta, *sig_delta, *sig_o, *beta, *rhot, * betat0, *betat, *fZ;
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
     fZ = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     
     GetRNGstate();     
     for(i=0; i<its1; i++){
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
                   
     zlt_fore_gp_tp(cov, K, nsite, n, r, p, u, rT, T, rK, nrK, d, d12, 
     phi, nu, sig_e, sig_eta, sig_delta, sig_o, foreX, foreXtp, beta, 
     rhot, betat0, betat, wp, constant, fZ);

     for(j=0; j < ns*r1*K1; j++){       
         foreZ[j+i*ns*r1*K1] = fZ[j];                                                
     }
     printR(i, its1); 
     } // end of iteration loop
     PutRNGstate();
     
     free(phi); free(nu); free(sig_e); free(sig_eta); free(sig_delta); free(sig_o);
     free(beta); free(rhot); free(betat0); free(betat); free(fZ);
     
     return;
}
*/

/*
// for temporally varying beta
// K-step Forecasts without its
void zlt_fore_gp_tp(int *cov, int *K, int *nsite, int *n, int *r, int *p, int *u,
     int *rT, int *T, int *rK, int *nrK, double *d, double *d12, double *phi, 
     double *nu, double *sig_e, double *sig_eta, double *sig_delta, double *sig_op,
     double *foreX, double *foreXtp, double *beta, double *rhotp, double *betat0, 
     double *betat, double *w, int *constant, double *foreZ)
{ 
     int l, k, t, i, T1, K1, r1, n1, u1, ns, col;
     T1 =*T;
     K1 =*K;
     r1 =*r;
     n1 =*n;
     u1 = *u;
     ns =*nsite;
//     nns =n1*ns;
     col =*constant;

    double *S_eta, *Si_eta, *S_eta12, *S_eta12c, *det; 
    S_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_eta12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_eta12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phi, nu, d, S_eta);
    MInv(S_eta, Si_eta, n, det);    
    covF(cov, n, nsite, phi, nu, d12, S_eta12);

     double *mu, *sig, *s21, *XB, *XB1, *eta, *eps, *zfore;
     mu = (double *) malloc((size_t)((col)*sizeof(double))); 
     sig = (double *) malloc((size_t)((col)*sizeof(double)));            
     s21 = (double *) malloc((size_t)((col)*sizeof(double)));      
     XB = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     XB1 = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     eps = (double *) malloc((size_t)((col)*sizeof(double)));       
     zfore = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     
     double *betat1, *out, *foreXtp1, *foreXBtp;
     betat1 = (double *) malloc((size_t)((u1)*sizeof(double))); 
     out = (double *) malloc((size_t)((col)*sizeof(double))); 
     foreXtp1 = (double *) malloc((size_t)((u1)*sizeof(double))); 
     foreXBtp = (double *) malloc((size_t)((col)*sizeof(double))); 

     
     MProd(beta, constant, p, foreX, nrK, XB);  // nsiterK x 1   
  
     for(l=0; l<r1; l++){
       for(k=0; k<1; k++){     
         t = (T1-1);
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         extract_beta_t(t, T, u, betat, betat1); // u x T into u x 1
         for(i=0; i< u1; i++){
            out[0] = rhotp[i]*betat1[i];    // u x 1
//            mvrnormal(constant, out, sig_delta, constant, out); 
            betat1[i] = out[0];    // u x 1
         }
         for(i=0; i<ns; i++){
         extract_X41(l, k, i, nsite, rK, K, u, foreXtp, foreXtp1); // nsiterK x u into u x 1
         MProd(foreXtp1, constant, u, betat1, constant, foreXBtp); // 1 x 1

            extn_12(i, n, S_eta12,S_eta12c); // n x 1
            xTay(S_eta12c, Si_eta, w, n, mu); // 1 x 1 for mean
            xTay(S_eta12c, Si_eta, S_eta12c, n, s21); // 1 x 1
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            mvrnormal(constant, mu, sig, constant, eta); 
            mvrnormal(constant, mu, sig_e, constant, eps); 
            zfore[i] = XB1[i]+foreXBtp[0]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
       for(k=1; k<K1; k++){     
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         for(i=0; i< u1; i++){
            out[0] = rhotp[i]*betat1[i];    // u x 1
//            mvrnormal(constant, out, sig_delta, constant, out); 
            betat1[i] = out[0];    // u x 1
         }
         for(i=0; i<ns; i++){
         extract_X41(l, k, i, nsite, rK, K, u, foreXtp, foreXtp1); // u x 1
         MProd(foreXtp1, constant, u, betat1, constant, foreXBtp); // 1 x 1

            extn_12(i, n, S_eta12,S_eta12c);
            xTay(S_eta12c, Si_eta, w, n, mu); // 1 x 1 for mean
            xTay(S_eta12c, Si_eta, S_eta12c, n, s21);
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            mvrnormal(constant, mu, sig, constant, eta); 
            mvrnormal(constant, mu, sig_e, constant, eps);             
            zfore[i] = XB1[i]+foreXBtp[0]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
     }

     free(S_eta); free(det); 
     free(Si_eta), free(S_eta12), free(S_eta12c); 
     free(mu); free(sig); free(s21); free(XB); free(XB1); 
     free(eta); free(eps); free(zfore);
     free(betat1); free(out); free(foreXtp1); free(foreXBtp);
     
     return;
}
*/

/*
// for spatially and temporally varying beta
// d is the distance for actual locations n x n
// d12 is the distance between pred and observaed locations nsite x n
// nrK = nsite x r x k
void zlt_fore_gp_sptp_its(int *cov, int *its, int *K, int *nsite, int *n, int *r, 
     int *p, int *q, int *u, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_betap, 
     double *sig_deltap, double *sig_op, double *foreX, double *foreXsp, double *foreXtp, 
     double *betap, double *betasp, double *rhotp, double *betat0p, double *betatp, 
     double *wp, int *constant, double *foreZ)
{
     int i, j, its1, n1, ns, r1, T1, K1, p1, q1, u1, col;
     its1 = *its;
     n1 =*n;
     ns =*nsite;
     r1 =*r;
     T1 =*T;
     K1 =*K;
     p1 =*p;
     q1 =*q;
     u1 =*u;
     col =*constant;     

     double *phi, *nu, *sig_e, *sig_eta, *sig_beta, *sig_delta, *sig_o;
     double *beta, *betas, *rhot, * betat0, *betat, *fZ;
     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_beta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_delta = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_o = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     betas = (double *) malloc((size_t)((q1*n1*col)*sizeof(double)));       
     rhot = (double *) malloc((size_t)((u1)*sizeof(double)));       
     betat0 = (double *) malloc((size_t)((u1)*sizeof(double)));       
     betat = (double *) malloc((size_t)((u1*T1)*sizeof(double)));       
     fZ = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     
     GetRNGstate();     
     for(i=0; i<its1; i++){
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
                   
     zlt_fore_gp_sptp(cov, K, nsite, n, r, p, q, u, rT, T, rK, nrK, d, d12, 
     phi, nu, sig_e, sig_eta, sig_beta, sig_delta, sig_o, foreX, foreXsp, 
     foreXtp, beta, betas, rhot, betat0, betat, wp, constant, fZ);

     for(j=0; j < ns*r1*K1; j++){       
         foreZ[j+i*ns*r1*K1] = fZ[j];                                                
     }
     printR(i, its1); 
     } // end of iteration loop
     PutRNGstate();
     
     free(phi); free(nu); free(sig_e); free(sig_eta); free(sig_beta);
     free(sig_delta); free(sig_o); free(beta); free(betas); free(rhot); 
     free(betat0); free(betat); free(fZ);
     
     return;
}
*/

/*
// for spatially and temporally varying beta
// K-step Forecasts without its
void zlt_fore_gp_sptp(int *cov, int *K, int *nsite, int *n, int *r, int *p, 
     int *q, int *u, int *rT, int *T, int *rK, int *nrK, double *d, double *d12, 
     double *phi, double *nu, double *sig_e, double *sig_eta, double *sig_beta,
     double *sig_delta, double *sig_op, double *foreX, double *foreXsp, 
     double *foreXtp, double *beta, double *betas, double *rhotp, double *betat0, 
     double *betat, double *w, int *constant, double *foreZ)
{ 
     int l, k, t, i, j, T1, K1, r1, n1, u1, ns, col;
     T1 =*T;
     K1 =*K;
     r1 =*r;
     n1 =*n;
//     q1 =*q;
     u1 =*u;
     ns =*nsite;
//     nns =n1*ns;
     col =*constant;

    double *S, *Si, *S_12, *S_12c, *det; 
    S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

    covF(cov, n, n, phi, nu, d, S);
    MInv(S, Si, n, det);    
    covF(cov, n, nsite, phi, nu, d12, S_12);

     double *mu, *sig, *sigb, *s21, *XB, *XB1, *eta, *eps, *zfore;
     mu = (double *) malloc((size_t)((col)*sizeof(double))); 
     sig = (double *) malloc((size_t)((col)*sizeof(double))); 
     sigb = (double *) malloc((size_t)((col)*sizeof(double)));                 
     s21 = (double *) malloc((size_t)((col)*sizeof(double)));      
     XB = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     XB1 = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     eps = (double *) malloc((size_t)((col)*sizeof(double)));       
     zfore = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     
     double *oox, *betasp0, *betat1, *out, * foreXsp1, *foreXBsp, *foreXBsp1, *foreXtp1, *foreXBtp;
     oox = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     betasp0 = (double *) malloc((size_t)((col)*sizeof(double)));       
     betat1 = (double *) malloc((size_t)((u1)*sizeof(double))); 
     out = (double *) malloc((size_t)((col)*sizeof(double))); 
     foreXsp1 = (double *) malloc((size_t)((col)*sizeof(double)));       
     foreXBsp = (double *) malloc((size_t)((col)*sizeof(double)));       
     foreXBsp1 = (double *) malloc((size_t)((col)*sizeof(double)));       
     foreXtp1 = (double *) malloc((size_t)((u1)*sizeof(double))); 
     foreXBtp = (double *) malloc((size_t)((col)*sizeof(double))); 

     
     MProd(beta, constant, p, foreX, nrK, XB);  // nsiterK x 1   
  
     for(l=0; l<r1; l++){
       for(k=0; k<1; k++){     
         t = (T1-1);
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         extract_beta_t(t, T, u, betat, betat1); // u x T into u x 1
         for(i=0; i< u1; i++){
            out[0] = rhotp[i]*betat1[i];    // u x 1
//            mvrnormal(constant, out, sig_delta, constant, out); 
            betat1[i] = out[0];    // u x 1
         }
         for(i=0; i<ns; i++){
         extract_X41(l, k, i, nsite, rK, K, u, foreXtp, foreXtp1); // nsiterK x u into u x 1
         MProd(foreXtp1, constant, u, betat1, constant, foreXBtp); // 1 x 1

            extn_12(i, n, S_12,S_12c); // n x 1
            xTay(S_12c, Si, w, n, mu); // 1 x 1 for mean
            xTay(S_12c, Si, S_12c, n, s21); // 1 x 1
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            sigb[0] = sig_beta[0] * (1.0 - s21[0]);

            foreXBsp[0] = 0.0;
            for(j=0; j< *q; j++){ // for q var
                    extract_beta_sp(j, n, betas, oox); // n x 1
                    xTay(S_12c, Si, oox, n, betasp0); // 1 x 1
//                    mvrnormal(constant, betasp0, sigb, constant, betasp0);
                    extract_X_sp4(k, l, i, j, nsite, r, K, foreXsp, foreXsp1); // 1 x 1
                    MProd(foreXsp1, constant, constant, betasp0, constant, foreXBsp1);  // 1 x 1 
                    foreXBsp[0] = foreXBsp[0] + foreXBsp1[0];
            }
            mvrnormal(constant, mu, sig, constant, eta); 
            mvrnormal(constant, mu, sig_e, constant, eps); 
            zfore[i] = XB1[i]+foreXBsp[0]+foreXBtp[0]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
       for(k=1; k<K1; k++){     
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         for(i=0; i< u1; i++){
            out[0] = rhotp[i]*betat1[i];    // u x 1
//            mvrnormal(constant, out, sig_delta, constant, out); 
            betat1[i] = out[0];    // u x 1
         }
         for(i=0; i<ns; i++){
         extract_X41(l, k, i, nsite, rK, K, u, foreXtp, foreXtp1); // u x 1
         MProd(foreXtp1, constant, u, betat1, constant, foreXBtp); // 1 x 1

            extn_12(i, n, S_12,S_12c);
            xTay(S_12c, Si, w, n, mu); // 1 x 1 for mean
            xTay(S_12c, Si, S_12c, n, s21);
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            sigb[0] = sig_beta[0] * (1.0 - s21[0]);

            foreXBsp[0] = 0.0;
            for(j=0; j< *q; j++){ // for q var
                    extract_beta_sp(j, n, betas, oox); // n x 1
                    xTay(S_12c, Si, oox, n, betasp0); // 1 x 1
//                    mvrnormal(constant, betasp0, sigb, constant, betasp0);
                    extract_X_sp4(k, l, i, j, nsite, r, K, foreXsp, foreXsp1); // 1 x 1
                    MProd(foreXsp1, constant, constant, betasp0, constant, foreXBsp1);  // 1 x 1 
                    foreXBsp[0] = foreXBsp[0] + foreXBsp1[0];
            }
            mvrnormal(constant, mu, sig, constant, eta); 
            mvrnormal(constant, mu, sig_e, constant, eps);             
            zfore[i] = XB1[i]+foreXBsp[0]+foreXBtp[0]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
     }

     free(S); free(det); free(Si), free(S_12), free(S_12c); free(mu); free(sig); 
     free(sigb); free(s21); free(XB); free(XB1); free(eta); free(eps); free(zfore);
     free(oox); free(betasp0); free(betat1); free(out); free(foreXsp1); free(foreXBsp);
     free(foreXBsp1); free(foreXtp1); free(foreXBtp);

     return;
}
*/

/////////////////////// END ///////////////////////////////////////////////////
