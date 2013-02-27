//** File Name 'forecast_xb_gpp.c' **//

#include "main_gpp.h"
#include "covariance.h"
#include "common.h"
#include <math.h>
#include "mathematics.h"
#include "randgenerator.h"


// K-step Forecast 
// for forecast into the existing sites: 
       // 'n' is the total number of existing sites  
       // 'dnm' is the distance matrix of the existing and knot sites
// for forecast into the prediction sites:
       // 'n' is the total number of prediction sites  
       // 'dnm' is the distance matrix of the prediction and knot sites
// 'dm' is the distance matrix of the knot sites

void zlt_fore_gpp_its(int *cov, int *its, int *K, int *n, int *m, 
     int *r, int *p, int *rT, int *T, int *rK, int *nrK, double *dnm, double *dm, 
     double *phip, double *nup,  
     double *sig_ep, double *sig_etap, double *betap, double *rhop, double *wp, 
     double *foreX, int *constant, double *foreZ)
{
     int i, j, its1, n1, m1, r1, T1, K1, p1, col;
     its1 = *its;
     n1 =*n;
     m1 =*m;
     r1 =*r;
     T1 =*T;
     K1 =*K;
     p1 =*p;
     col =*constant;     

//     unsigned iseed = 44;
//     srand(iseed); 

     double *phi, *nu, *sig_e, *sig_eta, *rho, *beta, *w, *fZ;
     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     nu = (double *) malloc((size_t)((col)*sizeof(double)));            
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     rho = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     w = (double *) malloc((size_t)((m1*r1*T1)*sizeof(double)));       
     fZ = (double *) malloc((size_t)((n1*r1*K1*col)*sizeof(double)));       
     
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
     for(j=0; j<m1*r1*T1; j++){
        w[j] = wp[j+i*m1*r1*T1];
     }            
              
     zlt_fore_gpp(cov, K, n, m, r, p, rT, T, rK, nrK, dnm, dm, phi, nu, sig_e, sig_eta, 
     beta, rho, w, foreX, constant, fZ);
    
     for(j=0; j < n1*r1*K1; j++){       
         foreZ[j+i*n1*r1*K1] = fZ[j];                                                
     }
          printR(i, its1); 
     }// end of iteration loop
     PutRNGstate();
     
     free(phi); free(nu); free(sig_e); free(sig_eta);
     free(rho); free(beta); free(w); free(fZ);
     
     return;
}


// K-step Forecasts without its
void zlt_fore_gpp(int *cov, int *K, int *n, int *m, int *r, int *p, int *rT, int *T, 
     int *rK, int *nrK, double *dnm, double *dm, double *phi, double *nu, double *sig_e, 
     double *sig_eta, double *beta, double *rho, double *wp, double *foreX, 
     int *constant, double *foreZ)
{ 
     int l, k, t, i, T1, K1, r1, n1, m1, col;
     T1 =*T;
     K1 =*K;
     r1 =*r;
     n1 =*n;
     m1 =*m;
     col =*constant;
     
     double *C, *det, *I, *A, *mu, *XB, *XB1;
     double *wp1, *Aw, *eta, *eps, *zfore;
     C = (double *) malloc((size_t)((n1*m1)*sizeof(double)));       
//     Sigeta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));       
//     Sinv = (double *) malloc((size_t)((m1*m1)*sizeof(double)));       
     det = (double *) malloc((size_t)((col)*sizeof(double)));       
     I = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     A = (double *) malloc((size_t)((m1*n1)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));       
     XB = (double *) malloc((size_t)((n1*r1*K1*col)*sizeof(double)));       
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     wp1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));       
     Aw = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     eta = (double *) malloc((size_t)((m1*col)*sizeof(double)));       
     eps = (double *) malloc((size_t)((col)*sizeof(double)));       
     zfore = (double *) malloc((size_t)((n1*col)*sizeof(double)));       

     double *S, *C12c, *s21, *sig;
     S = (double *) malloc((size_t)((m1*m1)*sizeof(double)));            
     C12c = (double *) malloc((size_t)((m1*col)*sizeof(double)));       
     s21 = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig = (double *) malloc((size_t)((col)*sizeof(double)));       

/*     
      // exponential covariance
      if(cov[0] == 1){
        for(i = 0; i < (m1*m1); i++){
          S[i] = exp(-1.0*phi[0]*dm[i]);
//          Sigeta[i] = sig_eta[0]*S[i];
        }
//        MInv(Sigeta, Sinv_eta, m, det);
        for(i=0; i < m1*n1; i++){
          C[i] = exp(-1.0*phi[0]*dnm[i]);        
        }
      }
      // gaussian covariance
      if(cov[0] == 2){
        for(i = 0; i < (m1*m1); i++){
          S[i] = exp(-1.0*phi[0]*phi[0]*dm[i]*dm[i]);
//          Sigeta[i] = sig_eta[0]*S[i];          
        }
//        MInv(Sigeta, Sinv_eta, m, det);
        for(i=0; i < m1*n1; i++){
          C[i] = exp(-1.0*phi[0]*phi[0]*dnm[i]*dnm[i]);        
        }
      }
      // spherical covariance
      if(cov[0] == 3){
        for(i = 0; i < (m1*m1); i++){
         if(dm[i] > 0 && dm[i] <= 1.0/phi[0]){
         S[i] = (1.0-1.5*phi[0]*dm[i]+0.5*(phi[0]*dm[i])*(phi[0]*dm[i])*(phi[0]*dm[i]));
//         Sigeta[i] = sig_eta[0]*S[i];         
         }
         else if(dm[i] >= 1.0/phi[0]){
         S[i] = 0.0;     
//         Sigeta[i] = 0.0;
         }
         else{
         S[i] = 1.0;     
//         Sigeta[i] = 1.0*sig_eta[0];
         }        
        }
//        MInv(Sigeta, Sinv_eta, m, det);
        for(i=0; i < m1*n1; i++){
         if(dnm[i] > 0 && dnm[i] <= 1.0/phi[0]){
         C[i] = 1.0-1.5*phi[0]*dnm[i]+0.5*(phi[0]*dnm[i])*(phi[0]*dnm[i])*(phi[0]*dnm[i]);
         }
         else if(dnm[i] >= 1.0/phi[0]){
         C[i] = 0.0;
         }
         else{
         C[i] = 1.0;        
         }        
        }
      }
      // matern covariance, nu = 3/2
      if(cov[0] == 4){
        for(i = 0; i < (m1*m1); i++){
          S[i] = ((1.0+phi[0]*dm[i])*exp(-1.0*phi[0]*dm[i]));
//          Sigeta[i] = sig_eta[0]*S[i];
        }
//        MInv(Sigeta, Sinv_eta, m, det);
        for(i=0; i < m1*n1; i++){
          C[i] = (1.0+phi[0]*dnm[i])*exp(-1.0*phi[0]*dnm[i]);        
        }
      }
*/
    covF(cov, m, m, phi, nu, dm, S);
    covF(cov, n, m, phi, nu, dnm, C);

   MInv(S, S, m, det); // m x m
   MProd(S, m, m, C, n, A);  // n x m
   
   for(i=0; i<m1; i++){
         I[i] = 0.0;
   }
      
     mu[0] = 0.0;
     MProd(beta, constant, p, foreX, nrK, XB);  // nrK x 1     
     for(l=0; l<r1; l++){
       for(k=0; k<1; k++){     
         t = (T1-1);
         extract_alt2(l, k, n, rK, K, XB, XB1); // n x 1
         for(i=0; i<m1; i++){
             wp1[i] = wp[i+t*m1+l*m1*T1];  
         }     
         MProd(wp1, constant, m, A, n, Aw);  // n x 1
         for(i=0; i<n1; i++){
            extn_12(i, m, C, C12c); // m x 1
            xTay(C12c, S, C12c, m, s21); // 1 x 1
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            mu[0] = 0.0;
            mvrnormal(constant, mu, sig_e, constant, eps); 
            mu[0] = Aw[i];
            mvrnormal(constant, mu, sig, constant, eta);                   
            zfore[i] = XB1[i] + eta[0] + eps[0];                                 
         }
         put_together1(l, k, n, r, K, foreZ, zfore);
       }
       for(k=1; k<K1; k++){     
         for(i=0; i<m1; i++){
             wp1[i] = wp1[i]; // m x 1             
         }     
         MProd(wp1, constant, m, A, n, Aw);  // n x 1
         extract_alt2(l, k, n, rK, K, XB, XB1); // n x 1
         mvrnormal(constant, mu, sig_e, constant, eps); 
         for(i=0; i<n1; i++){
            extn_12(i, m, C, C12c); // m x 1
            xTay(C12c, S, C12c, m, s21); // 1 x 1
            if(s21[0] > 1.0){
                s21[0] = 1.0-pow(1,-320);
            }
            if(s21[0] == 1.0){ 
                s21[0] = 1.0-pow(1,-320);
            }
            sig[0] = sig_eta[0] * (1.0 - s21[0]);
            mu[0] = 0.0;
            mvrnormal(constant, mu, sig_e, constant, eps); 
            mu[0] = Aw[i];
            mvrnormal(constant, mu, sig, constant, eta);                   
            zfore[i] = XB1[i] + eta[0] + eps[0];                                 
         }
         put_together1(l, k, n, r, K, foreZ, zfore);
       }
     }
       
     free(S); free(det); 
     free(C); free(I); 
     free(A); free(mu); free(XB); free(XB1); free(wp1); free(Aw); 
     free(eta); free(eps); free(zfore);
     
     free(C12c); free(s21); free(sig);

     return;
}





/*




// K-step Forecast 
// for forecast into the existing sites: 
       // 'n' is the total number of existing sites  
       // 'dnm' is the distance matrix of the existing and knot sites
// for forecast into the prediction sites:
       // 'n' is the total number of prediction sites  
       // 'dnm' is the distance matrix of the prediction and knot sites
// 'dm' is the distance matrix of the knot sites
void zlt_fore_gpp_its(int *cov, int *its, int *K, int *n, int *m, int *r, int *p, 
     int *rT, int *T, int *rK, int *nrK, double *dnm, double *dm, double *phip, 
     double *sig_ep, double *sig_etap, double *betap, double *rhop, double *wp, 
     double *foreX, int *constant, double *foreZ)
{
     int i, j, its1, n1, m1, r1, T1, K1, p1, col;
     its1 = *its;
     n1 =*n;
     m1 =*m;
     r1 =*r;
     T1 =*T;
     K1 =*K;
     p1 =*p;
     col =*constant;     

     unsigned iseed = 44;
     srand(iseed); 

     double *phi, *sig_e, *sig_eta, *rho, *beta, *w, *fZ;
     phi = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_e = (double *) malloc((size_t)((col)*sizeof(double)));       
     sig_eta = (double *) malloc((size_t)((col)*sizeof(double)));       
     rho = (double *) malloc((size_t)((col)*sizeof(double)));       
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));       
     w = (double *) malloc((size_t)((m1*r1*T1)*sizeof(double)));       
     fZ = (double *) malloc((size_t)((n1*r1*K1*col)*sizeof(double)));       
     
     
     for(i=0; i<its1; i++){

     phi[0] = phip[i];         
     sig_e[0] = sig_ep[i];
     sig_eta[0] = sig_etap[i];
     rho[0] = rhop[i];
     for(j=0; j<p1; j++){
        beta[j] = betap[j+i*p1];
     }
     for(j=0; j<m1*r1*T1; j++){
        w[j] = wp[j+i*m1*r1*T1];
     }            
              
     zlt_fore_gpp(cov, K, n, m, r, p, rT, T, rK, nrK, dnm, dm, phi, sig_e, sig_eta, 
     beta, rho, w, foreX, constant, fZ);
    
     for(j=0; j < n1*r1*K1; j++){       
         foreZ[j+i*n1*r1*K1] = fZ[j];                                                
     }
          printR(i, its1); 
     }// end of iteration loop

     free(phi); free(sig_e); free(sig_eta);
     free(rho); free(beta); free(w); free(fZ);
     
     return;
}


// K-step Forecasts without its
// 'd' is the distance matrix of the existing sites
void zlt_fore_gpp(int *cov, int *K, int *n, int *m, int *r, int *p, int *rT, int *T, 
     int *rK, int *nrK, double *dnm, double *dm, double *phi, double *sig_e, 
     double *sig_eta, double *beta, double *rho, double *wp, double *foreX, 
     int *constant, double *foreZ)
{ 
     int l, k, t, i, T1, K1, r1, n1, m1, col;
     T1 =*T;
     K1 =*K;
     r1 =*r;
     n1 =*n;
     m1 =*m;
     col =*constant;
     
     double *C, *Sigeta, *Sinv_eta, *det, *I, *A, *mu, *XB, *XB1;
     double *wp1, *Aw, *eta, *eps, *zfore;
     C = (double *) malloc((size_t)((n1*m1)*sizeof(double)));       
     Sigeta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));       
     Sinv_eta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));       
     det = (double *) malloc((size_t)((col)*sizeof(double)));       
     I = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     A = (double *) malloc((size_t)((m1*n1)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));       
     XB = (double *) malloc((size_t)((n1*r1*K1*col)*sizeof(double)));       
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     wp1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));       
     Aw = (double *) malloc((size_t)((n1*col)*sizeof(double)));       
     eta = (double *) malloc((size_t)((m1*col)*sizeof(double)));       
     eps = (double *) malloc((size_t)((col)*sizeof(double)));       
     zfore = (double *) malloc((size_t)((n1*col)*sizeof(double)));       

     double *S;
     S = (double *) malloc((size_t)((m1*m1)*sizeof(double)));       

      // exponential covariance
      if(cov[0] == 1){
        for(i = 0; i < (m1*m1); i++){
          Sigeta[i] = sig_eta[0]*exp(-1.0*phi[0]*dm[i]);
          S[i] = Sigeta[i];
        }
        MInv(Sigeta, Sinv_eta, m, det);
        for(i=0; i < m1*n1; i++){
          C[i] = exp(-1.0*phi[0]*dnm[i]);        
        }
      }
      // gaussian covariance
      if(cov[0] == 2){
        for(i = 0; i < (m1*m1); i++){
          Sigeta[i] = sig_eta[0]*exp(-1.0*phi[0]*phi[0]*dm[i]*dm[i]);
          S[i] = Sigeta[i];
        }
        MInv(Sigeta, Sinv_eta, m, det);
        for(i=0; i < m1*n1; i++){
          C[i] = exp(-1.0*phi[0]*phi[0]*dnm[i]*dnm[i]);        
        }
      }
      // spherical covariance
      if(cov[0] == 3){
        for(i = 0; i < (m1*m1); i++){
         if(dm[i] > 0 && dm[i] <= 1.0/phi[0]){
         Sigeta[i] = sig_eta[0]*(1.0-1.5*phi[0]*dm[i]+0.5*(phi[0]*dm[i])*(phi[0]*dm[i])*(phi[0]*dm[i]));
         S[i] = Sigeta[i];
         }
         else if(dm[i] >= 1.0/phi[0]){
         Sigeta[i] = 0.0;
         S[i] = Sigeta[i];
         }
         else{
         Sigeta[i] = 1.0*sig_eta[0];
         S[i] = Sigeta[i];        
         }        
        }
        MInv(Sigeta, Sinv_eta, m, det);
        for(i=0; i < m1*n1; i++){
         if(dnm[i] > 0 && dnm[i] <= 1.0/phi[0]){
         C[i] = 1.0-1.5*phi[0]*dnm[i]+0.5*(phi[0]*dnm[i])*(phi[0]*dnm[i])*(phi[0]*dnm[i]);
         }
         else if(dnm[i] >= 1.0/phi[0]){
         C[i] = 0.0;
         }
         else{
         C[i] = 1.0;        
         }        
        }
      }
      // matern covariance, nu = 3/2
      if(cov[0] == 4){
        for(i = 0; i < (m1*m1); i++){
          Sigeta[i] = sig_eta[0]*((1.0+phi[0]*dm[i])*exp(-1.0*phi[0]*dm[i]));
          S[i] = Sigeta[i];
        }
        MInv(Sigeta, Sinv_eta, m, det);
        for(i=0; i < m1*n1; i++){
          C[i] = (1.0+phi[0]*dnm[i])*exp(-1.0*phi[0]*dnm[i]);        
        }
      }

   MProd(Sinv_eta, m, m, C, n, A);  // n x m

   free(Sigeta);

   for(i=0; i<m1; i++){
         I[i] = 0.0;
   }
      
     mu[0] = 0.0;
     MProd(beta, constant, p, foreX, nrK, XB);  // nrK x 1     
     for(l=0; l<r1; l++){
       for(k=0; k<1; k++){     
         t = (T1-1);
         extract_alt2(l, k, n, rK, K, XB, XB1); // n x 1
//         extract_alt2(l, t, m, rT, T, wp, wp1); // m x 1
//         mvrnormal(constant, I, S, m, eta); 
         for(i=0; i<m1; i++){
             mvrnormal(constant, mu, sig_eta, constant, eta);                   
             wp1[i] =  wp[i+t*m1+l*m1*T1];  
//             wp1[i] = rho[0]*wp1[i] + eta[i]; // m x 1
             wp1[i] = rho[0]*wp1[i] + eta[0]; // m x 1             
         }     

         MProd(wp1, constant, m, A, n, Aw);  // n x 1
            mvrnormal(constant, mu, sig_e, constant, eps); 
         for(i=0; i<n1; i++){
            zfore[i] = XB1[i] + Aw[i] + eps[0];                                 
         }
         put_together1(l, k, n, r, K, foreZ, zfore);
       }
       for(k=1; k<K1; k++){     
//         mvrnormal(constant, I, S, m, eta); 
         for(i=0; i<m1; i++){
             mvrnormal(constant, mu, sig_eta, constant, eta);                   
//             wp1[i] = rho[0]*wp1[i] + eta[i]; // m x 1
             wp1[i] = rho[0]*wp1[i] + eta[0]; // m x 1             
         }     
         MProd(wp1, constant, m, A, n, Aw);  // n x 1
         extract_alt2(l, k, n, rK, K, XB, XB1); // n x 1
            mvrnormal(constant, mu, sig_e, constant, eps); 
         for(i=0; i<n1; i++){
            zfore[i] = XB1[i]+Aw[i]+ eps[0];                                 
         }
         put_together1(l, k, n, r, K, foreZ, zfore);
       }
     }
       
     free(C); free(Sinv_eta); free(det); free(I); 
     free(A); free(mu); free(XB); free(XB1); free(wp1); free(Aw); 
     free(eta); free(eps); free(zfore);
     
     free(S);

     return;
}




*/


////////////////// the end //////////////////////
