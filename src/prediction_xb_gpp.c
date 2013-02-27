//** File Name 'prediction_xb_gpp.c' **//


#include "main_gpp.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"


////////////////////////////////////////////////////////////

// approach w-tilde prediction 
void z_pr_its_gpp_wtilde(int *cov, int *scale, int *its, int *nsite, 
     int *n, int *m, int *r, int *T, int *rT, int *p, int *nsiterT, 
     double *phi_etap, double *dm, double *dnsm, double *w0p, double *wp, 
     double *sig2ep, double *sig2etap, double *sig2lp, double *rhop, 
     double *betap, double *Xpred, int *constant, double *zpred)
{
//     unsigned iseed = 44;
//     srand(iseed); 

     int i, j, its1, nsite1, m1, r1, T1, p1, col;
     its1 = *its;
     nsite1 = *nsite;
     m1 = *m;
     r1 =*r;
     T1 =*T;
     p1 =*p;
     col =*constant;
     
     double *phi_eta, *w0, *w, *sig2e, *sig2eta, *sig2l; 
     double *rho, *beta, *zp;
     phi_eta = (double *) malloc((size_t)((col)*sizeof(double)));
     w0 = (double *) malloc((size_t)((r1*m1)*sizeof(double)));
     w = (double *) malloc((size_t)((r1*T1*m1)*sizeof(double)));
     sig2e = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2eta = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2l = (double *) malloc((size_t)((r1*col)*sizeof(double)));
     rho = (double *) malloc((size_t)((col)*sizeof(double)));
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     zp = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double)));
          
     GetRNGstate();
     for(i=0; i<its1; i++){

       phi_eta[0] = phi_etap[i];
       for(j=0; j<r1*m1; j++){
          w0[j] = w0p[j+i*m1*r1];
       }      
       for(j=0; j< r1*T1*m1; j++){
          w[j] = wp[j+i*r1*T1*m1];
       }    
       sig2e[0] = sig2ep[i];
       sig2eta[0] = sig2etap[i];
       for(j=0; j<r1; j++){
             sig2l[j] = sig2lp[j+i*r1];
       }
       rho[0] = rhop[i];          
       for(j=0; j< p1; j++){
          beta[j] = betap[j+i*p1];
       }    

       z_pr_gpp_wtilde(cov, nsite, n, m, r, T, rT, p, nsiterT, 
       phi_eta, dm, dnsm, w0, w, sig2e, sig2eta, sig2l, rho, 
       beta, Xpred, constant, zp);
     
       for(j=0; j< nsite1*r1*T1; j++){
           if( *scale == 1){     
           zpred[j+i*nsite1*r1*T1] = zp[j];
           }
           else if(*scale == 2){
           zpred[j+i*nsite1*r1*T1] = zp[j]*zp[j];
           }                 
           else if(*scale == 3){
           zpred[j+i*nsite1*r1*T1] = exp(zp[j]);
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
     
     free(phi_eta); free(w0); free(w); free(sig2e);
     free(sig2eta); free(sig2l); free(rho); free(beta);
     free(zp);

     return;
}


// approach dnsm x dm  with w-tilde   
void z_pr_gpp_wtilde(int *cov, int *nsite, 
     int *n, int *m, int *r, int *T, int *rT, int *p, int *nsiterT, 
     double *phi_etap, double *dm, double *dnsm, double *w0p, double *wp, 
     double *sig2ep, double *sig2etap, double *sig2lp, double *rhop, 
     double *betap, double *Xpred, int *constant, double *zpred)
{
     int i, m1, nsite1, r1, T1, col, l, t;

     m1 =*m;
     nsite1 =*nsite;
     r1 =*r;
     T1 =*T;
     col = *constant;

     double *S_eta, *det, *Snsm, *A, *Aw, *Aw0, *XB;
     S_eta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     Snsm = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     A = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     Aw = (double *) malloc((size_t)((r1*T1*nsite1)*sizeof(double)));
     Aw0 = (double *) malloc((size_t)((r1*nsite1)*sizeof(double)));
     XB = (double *) malloc((size_t)((r1*T1*nsite1)*sizeof(double)));


      // exponential covariance
      if(cov[0] == 1){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = exp(-1.0*phi_etap[0]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = exp(-phi_etap[0]*dnsm[i]);        
        }
      }
      // gaussian covariance
      if(cov[0] == 2){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = exp(-1.0*phi_etap[0]*phi_etap[0]*dm[i]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = exp(-phi_etap[0]*phi_etap[0]*dnsm[i]*dnsm[i]);        
        }
      }
      // spherical covariance
      if(cov[0] == 3){
        for(i = 0; i < (m1*m1); i++){
         if(dm[i] > 0 && dm[i] <= 1.0/phi_etap[0]){
         S_eta[i] = 1.0-1.5*phi_etap[0]*dm[i]+0.5*(phi_etap[0]*dm[i])*(phi_etap[0]*dm[i])*(phi_etap[0]*dm[i]);
         }
         else if(dm[i] >= 1.0/phi_etap[0]){
         S_eta[i] = 0.0;
         }
         else{
         S_eta[i] = 1.0;        
         }        
        }
        for(i=0; i < m1*nsite1; i++){
         if(dnsm[i] > 0 && dnsm[i] <= 1.0/phi_etap[0]){
         Snsm[i] = 1.0-1.5*phi_etap[0]*dnsm[i]+0.5*(phi_etap[0]*dnsm[i])*(phi_etap[0]*dnsm[i])*(phi_etap[0]*dnsm[i]);
         }
         else if(dnsm[i] >= 1.0/phi_etap[0]){
         Snsm[i] = 0.0;
         }
         else{
         Snsm[i] = 1.0;        
         }        
        }
      }
      // matern covariance, nu = 3/2
      if(cov[0] == 4){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = (1.0+phi_etap[0]*dm[i])*exp(-1.0*phi_etap[0]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = (1.0+phi_etap[0]*dnsm[i])*exp(-1.0*phi_etap[0]*dnsm[i]);        

        }
      }

     MInv(S_eta, S_eta, m, det);
     MProd(S_eta, m, m, Snsm, nsite, A); // nsite x m
     free(S_eta); free(det); free(Snsm); 
// wp is a m x rT matrix
     MProd(wp, rT, m, A, nsite, Aw);  // nsite x rT,  
// w0p is a m x r matrix
     MProd(w0p, r, m, A, nsite, Aw0);  // nsite x r,   
     free(A); 
     MProd(betap, constant, p, Xpred, nsiterT, XB); // nsiterT x 1 

     double *XB1, *mu, *out, *tmp;
     XB1 = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));
     tmp = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     
     mu[0] = 0.0;
     for(l=0; l < r1; l++) {
       for(t=0; t < T1; t++) {
         extract_alt2(l, t, nsite, rT, T, XB, XB1); // nsite x 1
         for(i=0; i<nsite1; i++){
             if(t==0){
               mvrnormal(constant, mu, sig2lp, constant, out);
               tmp[i] =  Aw[i+t*nsite1+l*nsite1*T1] + out[0];
             }
             else{
               mvrnormal(constant, mu, sig2etap, constant, out);
               tmp[i] =  Aw[i+t*nsite1+l*nsite1*T1] + out[0];
             }                
             mvrnormal(constant, mu, sig2ep, constant, out);
             tmp[i] = XB1[i] + tmp[i] + out[0];
         }         
         put_together1(l, t, nsite, r, T, zpred, tmp);
       }
     }
                            
     free(Aw); free(Aw0); free(XB); 
     free(XB1); free(mu); free(out); free(tmp);

     return;
}




//////////////////////////////////////////////////////////

// approach w prediction      
// dnsm is the nsite x m matrix
// dm is the m x m matrix
// dmns is the m x nsite matrix
void z_pr_its_gpp(int *cov, int *scale, 
     int *its, int *nsite, int *n, int *m, int *r, int *T, 
     int *rT, int *p, int *nsiterT, double *phi_etap, double *phi_0p, 
     double *dm, double *dnsm, double *w0p, double *wp, double *sig2ep,  
     double *sig2etap, double *sig2lp, double *mu_lp, double *rhop, 
     double *betap, double *Xpred, int *constant, double *zpred)
{
     unsigned iseed = 44;
     srand(iseed); 

     int i, j, its1, nsite1, m1, r1, T1, p1, col;
     its1 = *its;
     nsite1 = *nsite;
     m1 = *m;
     r1 =*r;
     T1 =*T;
     p1 =*p;
     col =*constant;
     
     double *phi_eta, *phi_0, *w0, *w, *sig2e, *sig2eta, *sig2l; 
     double *mu_l, *rho, *beta, *zp;
     phi_eta = (double *) malloc((size_t)((col)*sizeof(double)));
     phi_0 = (double *) malloc((size_t)((col)*sizeof(double)));
     w0 = (double *) malloc((size_t)((r1*m1)*sizeof(double)));
     w = (double *) malloc((size_t)((r1*T1*m1)*sizeof(double)));
     sig2e = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2eta = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2l = (double *) malloc((size_t)((r1*col)*sizeof(double)));
     mu_l = (double *) malloc((size_t)((r1*col)*sizeof(double)));
     rho = (double *) malloc((size_t)((col)*sizeof(double)));
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     zp = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double)));
          

     for(i=0; i<its1; i++){

       phi_eta[0] = phi_etap[i];
//           Rprintf("phi_eta %i\n", phi_eta[0]); // ok
       phi_0[0] = phi_0p[i];
//           Rprintf("phi_0 %f\n", phi_0[0]); // ok
       for(j=0; j<r1*m1; j++){
          w0[j] = w0p[j+i*m1*r1];
//           Rprintf("w0 %2.2f\n", w0[j]); // ok
       }      
       for(j=0; j< r1*T1*m1; j++){
          w[j] = wp[j+i*r1*T1*m1];
//           Rprintf("w %2.2f\n", w[j]); // ok
       }    
       sig2e[0] = sig2ep[i];
//           Rprintf("sig2e %2.2f\n", sig2e[0]); // ok
       sig2eta[0] = sig2etap[i];
//           Rprintf("sig2eta %2.2f\n", sig2eta[0]); // ok
       for(j=0; j<r1; j++){
             sig2l[j] = sig2lp[j+i*r1];
//           Rprintf("sig2l %2.2f\n", sig2l[j]); // ok
       }
       for(j=0; j<r1; j++){
             mu_l[j] = mu_lp[j+i*r1];
//           Rprintf("mu_l %2.2f\n", mu_l[j]); // ok
       }
       rho[0] = rhop[i];          
//           Rprintf("rho %2.2f\n", rho[0]); // ok
       for(j=0; j< p1; j++){
          beta[j] = betap[j+i*p1];
//           Rprintf("beta %2.2f\n", beta[j]); // ok
       }    

     z_pr_gpp(cov, nsite, n, m, r, T, rT, p, nsiterT, phi_eta, phi_0,
     dm, dnsm, w0, wp, sig2e, sig2eta, sig2l, mu_l, rho, beta,
     Xpred, constant, zp);
              
       // here, dmns is a m x nsite matrix  
       
       for(j=0; j< nsite1*r1*T1; j++){
           if( *scale == 1){     
           zpred[j+i*nsite1*r1*T1] = zp[j];
           }
           else if(*scale == 2){
           zpred[j+i*nsite1*r1*T1] = zp[j]*zp[j];
           }                 
           else if(*scale == 3){
           zpred[j+i*nsite1*r1*T1] = exp(zp[j]);
           }
           else {
                //;
                //exit(9);
           }                 
//           Rprintf("zpred %2.2f\n", zpred[j]); // ok
       }
     printR(i, its1); 
     } 

     free(phi_eta); free(phi_0); free(w0); free(w); free(sig2e);
     free(sig2eta); free(sig2l); free(mu_l); free(rho); free(beta);
     free(zp);
     
     return;
}

     
void z_pr_gpp(int *cov, int *nsite, int *n, int *m, int *r, int *T, int *rT, 
     int *p, 
     int *nsiterT, double *phi_etap, double *phi_0p, double *dm, double *dnsm, 
     double *w0p, double *wp, double *sig2ep, double *sig2etap,
     double *sig2lp, double *mu_lp, double *rhop, double *betap,
     double *Xpred, int *constant, double *zpred)
{
     int i, m1, nsite1, r1, T1, col;
     m1 =*m;
     nsite1 =*nsite;
     r1 =*r;
     T1 =*T;
     col = *constant;

     double *S_eta, *Snsm, *det, *w0pred, *wpred, *XB, *mu, *out;
     S_eta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     Snsm = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     w0pred = (double *) malloc((size_t)((r1*nsite1)*sizeof(double)));
     wpred = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double)));
     XB = (double *) malloc((size_t)((r1*T1*nsite1)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));

//
      // exponential covariance
      if(cov[0] == 1){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = exp(-1.0*phi_etap[0]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = exp(-phi_etap[0]*dnsm[i]);        
        }
      }
      // gaussian covariance
      if(cov[0] == 2){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = exp(-1.0*phi_etap[0]*phi_etap[0]*dm[i]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = exp(-phi_etap[0]*phi_etap[0]*dnsm[i]*dnsm[i]);        
        }
      }
      // spherical covariance
      if(cov[0] == 3){
        for(i = 0; i < (m1*m1); i++){
         if(dm[i] > 0 && dm[i] <= 1.0/phi_etap[0]){
         S_eta[i] = 1.0-1.5*phi_etap[0]*dm[i]+0.5*(phi_etap[0]*dm[i])*(phi_etap[0]*dm[i])*(phi_etap[0]*dm[i]);
         }
         else if(dm[i] >= 1.0/phi_etap[0]){
         S_eta[i] = 0.0;
         }
         else{
         S_eta[i] = 1.0;        
         }        
        }
        for(i=0; i < m1*nsite1; i++){
         if(dnsm[i] > 0 && dnsm[i] <= 1.0/phi_etap[0]){
         Snsm[i] = 1.0-1.5*phi_etap[0]*dnsm[i]+0.5*(phi_etap[0]*dnsm[i])*(phi_etap[0]*dnsm[i])*(phi_etap[0]*dnsm[i]);
         }
         else if(dnsm[i] >= 1.0/phi_etap[0]){
         Snsm[i] = 0.0;
         }
         else{
         Snsm[i] = 1.0;        
         }        
        }
      }
      // matern covariance, nu = 3/2
      if(cov[0] == 4){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = (1.0+phi_etap[0]*dm[i])*exp(-1.0*phi_etap[0]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = (1.0+phi_etap[0]*dnsm[i])*exp(-1.0*phi_etap[0]*dnsm[i]);        
        }
      }
        
     MInv(S_eta, S_eta, m, det);
     free(det); 
//     for(i=0; i< m1*m1; i++){
//     Rprintf("S_eta %2.2f\n", S_eta[i]); //  ok
//    }
//     for(i=0; i< m1*nsite1; i++){
//     Rprintf("Smns %2.2f\n", Smns[i]); //  ok
//     }
     
     w0_pr_gpp_multi(nsite, m, r, w0p, mu_lp, sig2lp, S_eta, Snsm, constant, 
     w0pred);
//     w0_pr_gpp(nsite, m, r, w0p, mu_lp, sig2lp, S_eta, Snsm, constant, 
//     w0pred);
//     for(i=0; i< nsite1; i++){
//     Rprintf("w0pred %2.2f\n", w0pred[i]); // w0pred not ok
//     w0pred[i]=0.001;
//     }
         
     wlt_pr_gpp_multi(nsite, m, r, T, w0pred, wp, w0p, sig2etap, rhop, S_eta, 
     Snsm, constant, wpred);  // nsitre x rT
//     wlt_pr_gpp(nsite, m, r, T, w0pred, wp, w0p, sig2etap, rhop, S_eta, 
//     Snsm, constant, wpred);  // nsitre x rT
//     for(i=0; i< nsite1*r1*T1; i++){
//     Rprintf("wpred %2.2f\n", wpred[i]); // wpred not ok
//     }
     
     MTranspose(wpred, rT, nsite, wpred); // rT x nsite
    
     MProd(betap, constant, p, Xpred, nsiterT, XB);  // nsiterT x 1
     
     mu[0] = 0.0;
     for(i=0; i< nsite1*r1*T1; i++){
          mvrnormal(constant, mu, sig2ep, constant, out);
          zpred[i] = XB[i] + wpred[i] + out[0];
//    Rprintf("XB %2.2f, wpred %2.2f, er %2.2f\n", XB[i], wpred[i], out[0]); // 
//    Rprintf("zpr %2.2f, XB %2.2f, wpred %2.2f \n", zpred[i], XB[i], wpred[i]); // 
     }

     free(S_eta); free(w0pred); free(Snsm);
     free(wpred); free(XB); free(mu); free(out);
     
     return;
}

  
// multivariate approach                            
void wlt_pr_gpp_multi(int *nsite, int *m, int *r, int *T, 
     double *w0pred, double *wp, double *w0p, 
     double *sig2etap, double *rhop, double *Sinv, double *Snsm,
     int *constant, double *wpred)
{
     int i, k, t, l, nsite1, m1, r1, T1, col;
     nsite1 =*nsite;
     m1 =*m;
     r1 =*r;
     T1 =*T;
     col =*constant;

     double *Smns, *out, *del;
     Smns = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     del = (double *) malloc((size_t)((nsite1*nsite1)*sizeof(double)));

     double *w, *er, *out1, *mu, *pr;
     w = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     er = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     mu = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     pr = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     
// dnsm is nsite x m matrix
// dmns is m x nsite matrix

     MTranspose(Snsm, m, nsite, Smns); // m x nsite  
     MProd(Smns, nsite, m, Sinv, m, out); // m x nsite
     MProd(out, nsite, m, Snsm, nsite, del);  // nsite x nsite
     for(k=0; k<nsite1*nsite1; k++){
       del[k] = sig2etap[0]*(1.0-del[k]);
//     Rprintf("del %2.2f\n", del[k]); // ok
     }
     free(out); free(Smns);
     
     for(l=0; l<r1; l++){
     for(t=0; t<T1; t++){
         if(t == 0){
              for(i=0; i<nsite1; i++){
                w[i] = rhop[0]*w0pred[i+l*nsite1];
              }
              for(i=0; i<m1; i++){
                er[i] = wp[i+t*m1+l*m1*T1] - rhop[0]*w0p[i+l*m1];           
              }
         }
         else{
              for(i=0; i<nsite1; i++){
                w[i] = rhop[0]*w[i];
              }
              for(i=0; i<m1; i++){
                er[i] = wp[i+t*m1+l*m1*T1] - rhop[0]*wp[i+(t-1)*m1+l*m1*T1];           
              }
         }
         MProd(er, constant, m, Sinv, m, out1); // m x 1
         MProd(out1, constant, m, Snsm, nsite, mu); // nsite x 1
         for(i=0; i < nsite1; i++){
            mu[i] = w[i] + mu[i];
//        Rprintf("mu %2.2f\n", mu[i]); // ok
         }
         mvrnormal(constant, mu, del, constant, pr);
         for(i=0; i < nsite1; i++){
            wpred[i+t*nsite1+l*nsite1*T1] = pr[i];
//            wpred[i+t*nsite1+l*nsite1*T1] = mu[i];
//        Rprintf("wpred %2.2f\n", pr[i]); // ok
        }          
     }
     }

     free(del); 
     free(w); free(er); free(out1); free(mu); free(pr);
     return;
}


// multivariate approach
void w0_pr_gpp_multi(int *nsite, int *m, int *r, double *w0p, double *mu_lp, 
     double *sig2lp, double *Sinv, double *Snsm, int *constant, 
     double *w0pred)
{
     int i, k, l, nsite1, m1, r1, col;
     nsite1 =*nsite;
     m1 =*m;
     r1 =*r;
     col =*constant;

     double *Smns, *out, *del;
     Smns = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));     
     out = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     del = (double *) malloc((size_t)((nsite1*nsite1)*sizeof(double)));

     double *er, *out1, *mu, *pr;
     er = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     mu = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     pr = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     

// dnsm is nsite x m matrix
// dmns is m x nsite matrix

     for(l=0; l<r1; l++){
       MTranspose(Snsm, m, nsite, Smns); // m x nsite  
       MProd(Smns, nsite, m, Sinv, m, out); // m x nsite
       MProd(out, nsite, m, Snsm, nsite, del);  // nsite x nsite
       for(k=0; k<nsite1*nsite1; k++){
         del[k] = sig2lp[l]*(1.0-del[k]);
//       Rprintf("del0 %2.2f\n", del[k]); // ok
       }

       for(i=0; i<m1; i++){
            er[i] = w0p[i+l*m1] - mu_lp[l];          
       }
       MProd(er, constant, m, Sinv, m, out1);  // m x 1
       MProd(out1, constant, m, Snsm, nsite, mu); // nsite x 1
       for(i=0; i < nsite1; i++){
            mu[i] = mu_lp[l] + mu[i];
//        Rprintf("mu0 %2.2f\n", mu[i]); // ok
       }
       mvrnormal(constant, mu, del, constant, pr);  
       for(i=0; i < nsite1; i++){
//            w0pred[i+l*nsite1] = pr[i];
            w0pred[i+l*nsite1] = mu[i];            
//        Rprintf("wpred0 %2.2f\n", pr[i]); // ok
       }          
     }

     free(out); free(Smns); free(del);
     free(er); free(out1); free(mu); free(pr);
     return;
}



// univariate approach
void wlt_pr_gpp(int *nsite, int *m, int *r, int *T, 
     double *w0pred, double *wp, double *w0p, 
     double *sig2etap, double *rhop, double *Sinv, double *Snsm,
     int *constant, double *wpred)
{
     int i, k, t, l, nsite1, m1, r1, T1, col;
     nsite1 =*nsite;
     m1 =*m;
     r1 =*r;
     T1 =*T;
     col =*constant;

     double *Sinv_eta, *S12, *S21, *out, *outt, *del, *w, *er;
     double *out1, *tr, *mu, *s2, *pr;
     Sinv_eta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     S12 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     S21 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     outt = (double *) malloc((size_t)((col)*sizeof(double)));
     del = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     
     w = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     er = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     tr = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     mu = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     s2 = (double *) malloc((size_t)((col)*sizeof(double)));
     pr = (double *) malloc((size_t)((col)*sizeof(double)));
     
     for(i=0; i<m1; i++){
      Sinv_eta[i] = Sinv[i]/sig2etap[0];
     }

// dnsm is nsite x m matrix
// dmns is m x nsite matrix

     for(k=0; k<nsite1; k++){
     extn_12(k, m, Snsm, S12); // 1 x m             
     MTranspose(S12, m, constant, S21); // m x 1  
     MProd(S21, constant, m, Sinv_eta, m, out); // m x 1
     MProd(out, constant, m, S12, constant, outt);  // 1 x 1
     del[k] = sig2etap[0]*(1.0 - outt[0]);
         if(del[k] <= 0){
             del[k] = pow(1,-320);
         }
//     Rprintf("del %2.2f\n", del[k]); // ok
     }
     
     for(l=0; l<r1; l++){
     for(t=0; t<T1; t++){
         if(t == 0){
              for(i=0; i<nsite1; i++){
                w[i] = rhop[0]*w0pred[i+l*nsite1];
              }
              for(i=0; i<m1; i++){
                er[i] = wp[i+t*m1+l*m1*T1] - rhop[0]*w0p[i+l*m1];           
              }
         }
         else{
              for(i=0; i<nsite1; i++){
                w[i] = rhop[0]*w[i];
              }
              for(i=0; i<m1; i++){
                er[i] = wp[i+t*m1+l*m1*T1] - rhop[0]*wp[i+(t-1)*m1+l*m1*T1];           
              }
         }
         MProd(er, constant, m, Sinv_eta, m, out1); // m x 1
         MProd(out1, constant, m, Snsm, nsite, tr); // nsite x 1
         for(i=0; i < nsite1; i++){
            mu[0] = w[i] + tr[i];
            s2[0] = del[i];
//            Rprintf("w %2.2f, mu %2.2f, sd %2.2f\n", w[i], mu[0],s2[0]); // ok
            mvrnormal(constant, mu, s2, constant, pr);
            wpred[i+t*nsite1+l*nsite1*T1] = pr[0];
        }          
     }
     }
     free(S12); free(S21); free(out); free(outt); free(del); free(w); 
     free(er); free(out1); free(tr); free(mu); free(s2); free(pr);
     return;
}


// univariate approach
void w0_pr_gpp(int *nsite, int *m, int *r, double *w0p, double *mu_lp, 
     double *sig2lp, double *Sinv, double *Snsm, int *constant, 
     double *w0pred)
{
     int i, k, l, nsite1, m1, r1, col;
     nsite1 =*nsite;
     m1 =*m;
     r1 =*r;
     col =*constant;

     double *Sinv0, *S12, *S21, *out, *del, *er; 
     double *out1, *tr, *mu, *s2, *pr;
     Sinv0 = (double *) malloc((size_t)((m1*m1)*sizeof(double)));     
     S12 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     S21 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     del = (double *) malloc((size_t)((nsite1*nsite1)*sizeof(double)));
     er = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     tr = (double *) malloc((size_t)((nsite1)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     s2 = (double *) malloc((size_t)((col)*sizeof(double)));
     pr = (double *) malloc((size_t)((col)*sizeof(double)));
     

// dnsm is nsite x m matrix
// dmns is m x nsite matrix

     for(l=0; l<r1; l++){
        for(i=0; i<m1; i++){
         Sinv0[i] = Sinv[i]/sig2lp[l];
        }
        for(k=0; k<nsite1; k++){
        extn_12(k, m, Snsm, S12); // 1 x m             
        MTranspose(S12, m, constant, S21); // m x 1  
        MProd(S21, constant, m, Sinv0, m, out);  // m x 1
        MProd(out, constant, m, S12, constant, out);  // 1 x 1
        del[k] = sig2lp[l]*(1.0 - out[0]);
          if(del[k] <= 0){
              del[k] = pow(1,-320);
          }
//     Rprintf("del0 %2.2f\n", del[k]); // ok
        }     
        for(i=0; i<m1; i++){
            er[i] = w0p[i+l*m1] - mu_lp[l];          
        }
        MProd(er, constant, m, Sinv0, m, out1);  // m x 1
        MProd(out1, constant, m, Snsm, nsite, tr); // nsite x 1
        for(i=0; i < nsite1; i++){
            mu[0] = mu_lp[l] + tr[i];
            s2[0] = del[k];
//        Rprintf("mulp %2.2f, mu0 %2.2f, sd0 %22f\n", mu_lp[l], mu[0],s2[0]); // ok
            mvrnormal(constant, mu, s2, constant, pr);  
            w0pred[i+l*nsite1] = pr[0];
        }          
     }

     free(Sinv0); free(S12); free(S21); free(out); free(del);
     free(er); free(out1); free(tr); free(mu); free(s2); free(pr);
     return;
}



///////////////////////////////////////////////////

//
//

// the prediction
void z_pr_its_gpp1(int *cov, int *scale, int *its, int *nsite, int *n, int *m, int *r, 
     int *T, int *rT, int *p, int *nsiterT, double *phi_etap, double *nup, double *dm, 
     double *dnsm, double *wp, double *sig2ep, double *betap, double *Xpred, 
     int *constant, double *zpred)
{
     unsigned iseed = 44;
     srand(iseed); 

     int i, j, its1, nsite1, m1, r1, T1, p1, col;
     its1 = *its;
     nsite1 = *nsite;
     m1 = *m;
     r1 =*r;
     T1 =*T;
     p1 =*p;
     col =*constant;
     
     double *phi_eta, *nu, *w, *sig2e, *beta, *zp;
     phi_eta = (double *) malloc((size_t)((col)*sizeof(double)));
     nu = (double *) malloc((size_t)((col)*sizeof(double)));     
     w = (double *) malloc((size_t)((r1*T1*m1)*sizeof(double)));
     sig2e = (double *) malloc((size_t)((col)*sizeof(double)));
     beta = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     zp = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double)));
          

     for(i=0; i<its1; i++){

       phi_eta[0] = phi_etap[i];
         if(cov[0]==4){
           nu[0] = nup[i];
         }
         else{
           nu[0] = 0.0;
         }
       for(j=0; j< r1*T1*m1; j++){
          w[j] = wp[j+i*r1*T1*m1];
       }    
       sig2e[0] = sig2ep[i];
       for(j=0; j< p1; j++){
          beta[j] = betap[j+i*p1];
       }    
              
       z_pr_gpp1(cov, nsite, n, m, r, T, rT, p, nsiterT, phi_eta, nu, dm, dnsm, 
       w, sig2e, beta, Xpred, constant, zp);

       for(j=0; j< nsite1*r1*T1; j++){
           if( *scale == 1){     
           zpred[j+i*nsite1*r1*T1] = zp[j];
           }
           else if(*scale == 2){
           zpred[j+i*nsite1*r1*T1] = zp[j]*zp[j];
           }                 
           else if(*scale == 3){
           zpred[j+i*nsite1*r1*T1] = exp(zp[j]);
           }
           else {
                //;
                //exit(9);
           }                 
//           Rprintf("zpred %2.2f\n", zpred[j]); // ok
       }
     printR(i, its1); 
     } // end of iteration loop

     free(phi_eta); free(nu); free(w); free(sig2e); free(beta); free(zp);
     return;
}



// approach dnsm x dm     
void z_pr_gpp1(int *cov, int *nsite, int *n, int *m, int *r, int *T, 
     int *rT, int *p, int *nsiterT, double *phi_etap, double *nup, double *dm, 
     double *dnsm, double *wp, double *sig2ep, double *betap, 
     double *Xpred, int *constant, double *zpred)
{
     int i, m1, nsite1, r1, T1, col;

     m1 =*m;
     nsite1 =*nsite;
     r1 =*r;
     T1 =*T;
     col = *constant;

     double *S_eta, *det, *Snsm, *A, *Aw, *tAw, *XB, *mu, *out;
     S_eta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     Snsm = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));

// slot 1
     A = (double *) malloc((size_t)((m1*nsite1)*sizeof(double)));
     Aw = (double *) malloc((size_t)((r1*T1*nsite1)*sizeof(double)));
     tAw = (double *) malloc((size_t)((r1*T1*nsite1)*sizeof(double)));
     XB = (double *) malloc((size_t)((r1*T1*nsite1)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));

/*
      // exponential covariance
      if(cov[0] == 1){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = exp(-1.0*phi_etap[0]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = exp(-phi_etap[0]*dnsm[i]);        
        }
      }
      // gaussian covariance
      if(cov[0] == 2){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = exp(-1.0*phi_etap[0]*phi_etap[0]*dm[i]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = exp(-phi_etap[0]*phi_etap[0]*dnsm[i]*dnsm[i]);        
        }
      }
      // spherical covariance
      if(cov[0] == 3){
        for(i = 0; i < (m1*m1); i++){
         if(dm[i] > 0 && dm[i] <= 1.0/phi_etap[0]){
         S_eta[i] = 1.0-1.5*phi_etap[0]*dm[i]+0.5*(phi_etap[0]*dm[i])*(phi_etap[0]*dm[i])*(phi_etap[0]*dm[i]);
         }
         else if(dm[i] >= 1.0/phi_etap[0]){
         S_eta[i] = 0.0;
         }
         else{
         S_eta[i] = 1.0;        
         }        
        }
        for(i=0; i < m1*nsite1; i++){
         if(dnsm[i] > 0 && dnsm[i] <= 1.0/phi_etap[0]){
         Snsm[i] = 1.0-1.5*phi_etap[0]*dnsm[i]+0.5*(phi_etap[0]*dnsm[i])*(phi_etap[0]*dnsm[i])*(phi_etap[0]*dnsm[i]);
         }
         else if(dnsm[i] >= 1.0/phi_etap[0]){
         Snsm[i] = 0.0;
         }
         else{
         Snsm[i] = 1.0;        
         }        
        }
      }
      // matern covariance, nu = 3/2
      if(cov[0] == 4){
        for(i = 0; i < (m1*m1); i++){
          S_eta[i] = (1.0+phi_etap[0]*dm[i])*exp(-1.0*phi_etap[0]*dm[i]);
        }
        for(i=0; i < m1*nsite1; i++){
          Snsm[i] = (1.0+phi_etap[0]*dnsm[i])*exp(-1.0*phi_etap[0]*dnsm[i]);        

        }
      }
*/      

    covF(cov, m, m, phi_etap, nup, dm, S_eta);
    covF(cov, m, nsite, phi_etap, nup, dnsm, Snsm);
    MInv(S_eta, S_eta, m, det);
    MProd(S_eta, m, m, Snsm, nsite, A); // nsite x m

// wp is a m x rT matrix
     MProd(wp, rT, m, A, nsite, Aw);  // nsite x rT,  
     MTranspose(Aw, rT, nsite, tAw); // rT x nsite
     MProd(betap, constant, p, Xpred, nsiterT, XB); // nsiterT x 1 
     
     mu[0] = 0.0;
     for(i=0; i< nsite1*r1*T1; i++){
          mvrnormal(constant, mu, sig2ep, constant, out);
          zpred[i] = XB[i] + tAw[i] + out[0];
     }
                            
     free(S_eta); free(det); free(Snsm); free(A); 
     free(Aw); free(tAw); free(XB); free(mu); free(out);

     return;
}




/////////////////////// the end //////////////////////////
