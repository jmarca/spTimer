//** File Name 'forecast_xb_ar.c' **//

#include"main_ar.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"

// K-step Forecasts into the any sites with its
// d is the distance for actual locations n x n
// d12 is the distance between pred/fore and observaed locations nsite x n
// nrK = nsite x r x k

void zlt_fore_ar_its_anysite(int *cov, int *its, int *K, int *nsite, int *n, 
     int *r, int *p, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *rhop, 
     double *foreX, double *betap, double *zpred_exist, double *wp, 
     int *constant, double *foreZ)
{
     int i, j, its1, n1, ns, r1, T1, K1, p1, col;
     its1 = *its;
     n1 =*n;
     ns =*nsite;
     r1 =*r;
     T1 =*T;
     K1 =*K;
     p1 =*p;
     col =*constant;     

//     unsigned iseed = 44;
//     srand(iseed); 

     double *phi, *nu, *sig_e, *sig_eta, *rho, *beta, *z, *fZ, *mu, *w;
     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       

     rho = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     z = (double *) malloc((size_t)((ns*r1*T1*col)*sizeof(double)));       
     fZ = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     mu = (double *) malloc((size_t)((col)*sizeof(double)));       
     w = (double *) malloc((size_t)((n1)*sizeof(double)));            

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
     rho[0] = rhop[i];
     for(j=0; j<p1; j++){
        beta[j] = betap[j+i*p1];
     }
     for(j=0; j<ns*r1*T1; j++){
        z[j] = zpred_exist[j+i*ns*r1*T1];
        mu[0] = z[j];
        mvrnormal(constant, mu, sig_e, constant, mu);
        z[j] = mu[0]; 
     }            
     //for(j=0; j<n1; j++){
     //   w[j] = wp[j+i*n1];
     //}
             
     zlt_fore_ar(cov, K, nsite, n, r, p, rT, T, rK, nrK, d, d12, phi, nu, sig_e, 
     sig_eta, rho, foreX, beta, z, wp, constant, fZ);
     
     for(j=0; j < ns*r1*K1; j++){       
         foreZ[j+i*ns*r1*K1] = fZ[j];                                                
     }
     printR(i, its1); 
     } // end of iteration loop
     PutRNGstate();
     
     free(phi); free(nu); free(sig_e); free(sig_eta); free(rho);
     free(beta); free(z); free(fZ); free(mu); free(w);
     
     return;
}
     

// K-step Forecasts without its
void zlt_fore_ar(int *cov, int *K, int *nsite, int *n, int *r, int *p, 
     int *rT, int *T, int *rK, int *nrK, double *d, double *d12, 
     double *phi, double *nu, double *sig_e, double *sig_eta, double *rho, 
     double *foreX, double *beta, double *z, double *w, 
     int *constant, double *foreZ)
{ 
     int l, k, t, i, T1, K1, r1, n1, ns, nns, col;
     T1 =*T;
     K1 =*K;
     r1 =*r;
     n1 =*n;
     ns =*nsite;
     nns =n1*ns;
     col =*constant;


    double *S_eta, *Si_eta, *S_eta12, *S_eta12c, *det; 
    S_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    Si_eta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));    
    S_eta12 = (double *) malloc((size_t)((n1*ns)*sizeof(double)));     
    S_eta12c = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    det = (double *) malloc((size_t)((col)*sizeof(double))); 

/*      
      // exponential covariance
      if(cov[0] == 1){
        for(i = 0; i < (n1*n1); i++){
          S_eta[i] = exp(-1.0*phi[0]*d[i]);
        }
        MInv(S_eta, Si_eta, n, det);
        for(i=0; i < nns; i++){
          S_eta12[i] = exp(-phi[0]*d12[i]);        
        }
      }
      // gaussian covariance
      if(cov[0] == 2){
        for(i = 0; i < (n1*n1); i++){
          S_eta[i] = exp(-1.0*phi[0]*phi[0]*d[i]*d[i]);
        }
        MInv(S_eta, Si_eta, n, det);
        for(i=0; i < nns; i++){
          S_eta12[i] = exp(-1.0*phi[0]*phi[0]*d12[i]*d12[i]);        
        }
      }
      // spherical covariance
      if(cov[0] == 3){
        for(i = 0; i < (n1*n1); i++){
         if(d[i] > 0 && d[i] <= 1.0/phi[0]){
         S_eta[i] = 1.0-1.5*phi[0]*d[i]+0.5*(phi[0]*d[i])*(phi[0]*d[i])*(phi[0]*d[i]);
         }
         else if(d[i] >= 1.0/phi[0]){
         S_eta[i] = 0.0;
         }
         else{
         S_eta[i] = 1.0;        
         }        
        }
        MInv(S_eta, Si_eta, n, det);
        for(i=0; i < nns; i++){
         if(d12[i] > 0 && d12[i] <= 1.0/phi[0]){
         S_eta12[i] = 1.0-1.5*phi[0]*d12[i]+0.5*(phi[0]*d12[i])*(phi[0]*d12[i])*(phi[0]*d12[i]);
         }
         else if(d12[i] >= 1.0/phi[0]){
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
          S_eta[i] = (1.0+phi[0]*d[i])*exp(-1.0*phi[0]*d[i]);
        }
        MInv(S_eta, Si_eta, n, det);
        for(i=0; i < nns; i++){
          S_eta12[i] = (1.0+phi[0]*d12[i])*exp(-1.0*phi[0]*d12[i]);        
        }
      }
*/

    covF(cov, n, n, phi, nu, d, S_eta);
    MInv(S_eta, Si_eta, n, det);    
    covF(cov, n, nsite, phi, nu, d12, S_eta12);

     double *mu, *sig, *s21, *XB, *zT, *XB1, *eta, *eps, *zfore;
     mu = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig = (double *) malloc((size_t)((col)*sizeof(double)));            
     s21 = (double *) malloc((size_t)((col)*sizeof(double)));      
     XB = (double *) malloc((size_t)((ns*r1*K1*col)*sizeof(double)));       
     zT = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     XB1 = (double *) malloc((size_t)((ns*col)*sizeof(double)));       
     eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     eps = (double *) malloc((size_t)((col)*sizeof(double)));       
     zfore = (double *) malloc((size_t)((ns*col)*sizeof(double)));       

     MProd(beta, constant, p, foreX, nrK, XB);  // nsiterK x 1     
     for(l=0; l<r1; l++){
       for(k=0; k<1; k++){     
         t = (T1-1);
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         extract_alt2(l, t, nsite, rT, T, z, zT); // nsite x 1
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         for(i=0; i<ns; i++){
/*
            if(ns == n1){
            extn_12(i, n, S_eta,S_eta12c); // n x 1
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
            zfore[i] = rho[0]*(zT[i]-eps[0])+XB1[i]+eta[0]+eps[0];                                 
            }
            else{
*/
            extn_12(i, n, S_eta12,S_eta12c); // n x 1
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
            zfore[i] = rho[0]*(zT[i]-eps[0])+XB1[i]+eta[0]+eps[0];                                 
//            }
         }
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
       for(k=1; k<K1; k++){     
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         mvrnormal(constant, mu, sig_e, constant, eps); 
         extract_alt2(l, t, nsite, rT, T, z, zT); // nsite x 1
         extract_alt2(l, k, nsite, rK, K, XB, XB1); // nsite x 1
         for(i=0; i<ns; i++){
/*
            if(ns == n1){
            extn_12(i, n, S_eta,S_eta12c); // n x 1
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
            zfore[i] = rho[0]*(zT[i]-eps[0])+XB1[i]+eta[0]+eps[0];                                 
            }
            else{            
*/
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
            zfore[i] = rho[0]*(zfore[i]-eps[0])+XB1[i]+eta[0]+eps[0];                                 
//            }
         }  
         put_together1(l, k, nsite, r, K, foreZ, zfore);
       }
     }


     free(S_eta); free(Si_eta); free(S_eta12); free(S_eta12c); free(det);  
     free(mu); free(sig); free(s21); free(XB); free(zT); free(XB1); 
     free(eta); free(eps); free(zfore);

     return;
}


/*


// K-step Forecasts into the any sites with its
// 'd' is the distance matrix of the prediction/existing sites
void zlt_fore_ar_its_anysite(int *cov, int *its, int *K, int *n, int *r, int *p, int *rT, 
     int *T, int *rK, int *nrK, double *d, double *phip, double *sig_ep, 
     double *sig_etap, double *rhop, double *foreX, double *betap, 
     double *zpred_exist, int *constant, double *foreZ)
{
     int i, j, its1, n1, r1, T1, K1, p1, col;
     its1 = *its;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     K1 =*K;
     p1 =*p;
     col =*constant;     

     unsigned iseed = 44;
     srand(iseed); 

     double *phi, *sig_e, *sig_eta, *rho, *beta, *z, *fZ, *mu;
     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       

     rho = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     z = (double *) malloc((size_t)((n1*r1*T1*col)*sizeof(double)));       
     fZ = (double *) malloc((size_t)((n1*r1*K1*col)*sizeof(double)));       
     mu = (double *) malloc((size_t)((col)*sizeof(double)));       
     
     for(i=0; i<its1; i++){
     phi[0] = phip[i];         
     sig_e[0] = sig_ep[i];
     sig_eta[0] = sig_etap[i];
     rho[0] = rhop[i];
     for(j=0; j<p1; j++){
        beta[j] = betap[j+i*p1];
     }
     for(j=0; j<n1*r1*T1; j++){
        z[j] = zpred_exist[j+i*n1*r1*T1];
        mu[0] = z[j];
        mvrnormal(constant, mu, sig_e, constant, mu);
        z[j] = mu[0]; 
     }            
              
     zlt_fore_ar(cov, K, n, r, p, rT, T, rK, nrK, d, phi, sig_e, sig_eta, rho, 
     foreX, beta, z, constant, fZ);
     
     for(j=0; j < n1*r1*K1; j++){       
         foreZ[j+i*n1*r1*K1] = fZ[j];                                                
     }
     printR(i, its1); 
     } // end of iteration loop

     free(phi); free(sig_e); free(sig_eta); free(rho);
     free(beta); free(z); free(fZ); free(mu);
     
     return;
}
     

// K-step Forecasts without its
void zlt_fore_ar(int *cov, int *K, int *n, int *r, int *p, int *rT, int *T, int *rK, 
     int *nrK, double *d, double *phi, double *sig_e, double *sig_eta, 
     double *rho, double *foreX, double *beta, double *z, int *constant, 
     double *foreZ)
{ 
     int l, k, t, i, T1, K1, r1, n1, col;
     T1 =*T;
     K1 =*K;
     r1 =*r;
     n1 =*n;
     col =*constant;
     
     double *Sigeta, *det, *Sinv, *Qeta, *I, *mu, *XB, *zT, *XB1, *eta, *eps, *zfore;
     Sigeta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));       
     det = (double *) malloc((size_t)((col)*sizeof(double)));       
     Sinv = (double *) malloc((size_t)((n1*n1)*sizeof(double)));       
     Qeta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));       
     I = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     mu = (double *) malloc((size_t)((col)*sizeof(double)));       
     XB = (double *) malloc((size_t)((n1*r1*K1*col)*sizeof(double)));       
     zT = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     eta = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     eps = (double *) malloc((size_t)((col)*sizeof(double)));       
     zfore = (double *) malloc((size_t)((n1*col)*sizeof(double)));       

     covFormat(cov, n, phi, d, sig_eta, Sigeta, det, Sinv, Qeta);
     
     for(i=0; i<n1; i++){
         I[i] = 0.0;
     }     
    
     mu[0] = 0.0;
     MProd(beta, constant, p, foreX, nrK, XB);  // nrK x 1     
     for(l=0; l<r1; l++){
       for(k=0; k<1; k++){     
//         mvrnormal(constant, I, Sigeta, n, eta); 
         t = (T1-1);
         extract_alt2(l, t, n, rT, T, z, zT); // n x 1
         extract_alt2(l, k, n, rK, K, XB, XB1); // n x 1
            mvrnormal(constant, mu, sig_e, constant, eps); 
         for(i=0; i<n1; i++){
            mvrnormal(constant, mu, sig_eta, constant, eta); 
            zfore[i] = rho[0]*(zT[i]-eps[0])+XB1[i]+eta[0]+eps[0];                                 
//            zfore[i] = rho[0]*zT[i]+XB1[i]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, n, r, K, foreZ, zfore);
       }
       for(k=1; k<K1; k++){     
//         mvrnormal(constant, mu, sig_eta, constant, eta); 
         extract_alt2(l, k, n, rK, K, XB, XB1); // n x 1
            mvrnormal(constant, mu, sig_e, constant, eps); 
         for(i=0; i<n1; i++){
            mvrnormal(constant, mu, sig_eta, constant, eta); 
            zfore[i] = rho[0]*(zfore[i]-eps[0])+XB1[i]+eta[0]+eps[0];                                 
//            zfore[i] = rho[0]*zfore[i]+XB1[i]+eta[0]+eps[0];                                             
         }
         put_together1(l, k, n, r, K, foreZ, zfore);
       }
     }
       
     free(Sigeta); free(det); free(Sinv); free(Qeta); 
     free(I); free(mu); free(XB); free(zT); free(XB1); 
     free(eta); free(eps); free(zfore);

     return;
}




*/

/////////////////////// END ///////////////////////////////////////////////////










