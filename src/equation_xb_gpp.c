//** File Name 'equation_xb_gpp.c' **//

#include "main_gpp.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"


// Joint posterior distribution with one phi parameter
void JOINT_onephi_gpp(int *cov, int *spdecay, double *flag, int *n, int *m, 
     int *T, int *r, int *rT, int *p, int *N, double *shape_e, double *shape_eta, 
     double *shape_l, 
     double *phi_a, double *phi_b,
     double *prior_a, double *prior_b, double *mu_beta, 
     double *delta2_beta, double *mu_rho,  double *delta2_rho, double *alpha_l, 
     double *delta2_l, double *phi, double *tau, double *phis, int *phik,
     double *nu, double *dm, double *dnm, int *constant, 
     double *sig2e, double *sig2eta, double *sig2l, double *beta, 
     double *rho, double *mu_l, double *X, double *z, double *w0, double *w,
     double *phip, double *accept, double *nup,
     double *sig2ep, double *sig2etap, double *betap, double *rhop, 
     double *mu_lp, double *sig2lp, double *w0p, double *wp, 
     double *zfit)
{     
     int n1, m1, r1, N1; 
     n1 = *n;
     m1 = *m;
     r1 = *r;
//     p1 = *p;
     N1 = *N;
//     col = *constant;
     
   double *Qeta, *XB, *Sinv, *S, *det, *A, *C, *Aw;

   Qeta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
   XB = (double *) malloc((size_t)((N1)*sizeof(double)));
   Sinv = (double *) malloc((size_t)((m1*m1)*sizeof(double)));   
   S = (double *) malloc((size_t)((m1*m1)*sizeof(double)));   
   det = (double *) malloc((size_t)((1)*sizeof(double)));   
   A = (double *) malloc((size_t)((n1*m1)*sizeof(double)));   
   C = (double *) malloc((size_t)((n1*m1)*sizeof(double)));   
   Aw = (double *) malloc((size_t)((N1)*sizeof(double)));   

   covFormat(cov, m, phi, nu, dm, sig2eta, S, det, Sinv, Qeta);
   covF(cov, n, m, phi, nu, dnm, C);

   MProd(Sinv, m, m, C, n, A);  // n x m
   MProd(beta, constant, p, X, N, XB);

   wlt_gpp(n, m, r, T, rT, p, sig2e, rho, Qeta, A, w0, w, XB, z, 
   constant, wp);     
   w0_gpp(m, r, T, Qeta, sig2l, Sinv, rho, mu_l, wp, constant, w0p);     

// check nu
   if(cov[0]==4){
      nu_gpp_DIS(cov, Qeta, det, phi, nu, m, r, T, rT, dm, rho, sig2eta, 
      mu_l, w0, w, constant, nup); 
   }
   else {
      nup[0] = nu[0];
   }  

// fixed values for phi 
   if(spdecay[0] == 1){
      accept[0] =0.0;
      phip[0] = phi[0];
      covFormat(cov, m, phip, nup, dm, sig2eta, S, det, Sinv, Qeta);
   }
// discrete sampling for phi 
   else if(spdecay[0] == 2){
      phi_gpp_DIS2(cov, Qeta, det, phi, phis, phik, nup, m, r, T, rT, 
      prior_a, prior_b, dm, rho, sig2eta, mu_l, w0p, wp, constant, 
      accept, phip);
      covFormat(cov, m, phip, nup, dm, sig2eta, S, det, Sinv, Qeta);
      MProd(Sinv, m, m, C, n, A);  // n x m
   }
// Random-Walk MH-within-Gibbs sampling for phi
   else if(spdecay[0] == 3){
     double *Qeta2, *det2, *tmp, *phi2;
     Qeta2 = (double *) malloc((size_t)((m1*m1)*sizeof(double))); 
     det2 = (double *) malloc((size_t)((1)*sizeof(double)));   
     tmp = (double *) malloc((size_t)((1)*sizeof(double)));
     phi2 = (double *) malloc((size_t)((1)*sizeof(double)));
     
     if(phi[0] <= 0){
        phi[0] = pow(1,-320);
     }      
     //tmp[0] = -log(phi[0]); 
     //mvrnormal(constant, tmp, tau, constant, phi2);
     //phi2[0]= exp(-phi2[0]);
     tmp[0] = log(phi[0]); 
     mvrnormal(constant, tmp, tau, constant, phi2);
     phi2[0]= exp(phi2[0]);
     covFormat2(cov, m, phi2, nup, dm, sig2eta, det2, Qeta2);
     phi_gpp_MH2(Qeta, Qeta2, det, det2, phi, phi2, m, r, T, rT, 
     prior_a, prior_b, rho, mu_l, w0p, wp, constant, accept, phip);
     free(Qeta2); free(det2);
     free(tmp); free(phi2);
     if(accept[0]==1.0){
       covFormat(cov, m, phip, nup, dm, sig2eta, S, det, Sinv, Qeta);
       MProd(Sinv, m, m, C, n, A);  // n x m
     }
   }
   else {
//     ;
//     exit(9);
   }   

   MProd(wp, rT, m, A, n, Aw);  // n x rT
   rho_gpp(m, r, T, rT, p, mu_rho, delta2_rho, Qeta, w0p, wp, constant, rhop);        
   beta_gpp(n, p, rT, N, mu_beta, delta2_beta, sig2e, X, Aw, z, 
   constant, betap); 
   MProd(betap, constant, p, X, N, XB);
   int j;
   for(j=0; j < r1; j++) {
         mu_lp[j] = mu_l[j];     
   }
//   mu_l_gpp(m, r, sig2l, alpha_l, delta2_l, Sinv, w0p, constant, mu_lp);         
   sig_l_gpp(m, r, shape_l, prior_b, mu_lp, Sinv, w0p, constant, sig2lp);          
   sig_eta_gpp(m, r, T, rT, shape_eta, prior_b, Sinv, rhop, wp, w0p, 
   constant, sig2etap);
   sig_e_gpp(n, rT, N, shape_e, prior_b, XB, Aw, z, constant, sig2ep);
   Z_fit_gpp(flag, n, m, T, r, rT, sig2ep, Aw, XB, z, constant, zfit);
   
   free(Qeta); free(XB); free(Sinv); free(S); free(det); free(A);
   free(C); free(Aw); 
   
   return;
}


// Joint posterior distribution for spatial beta
void JOINT_onephi_sp_gpp(int *intercept, int *cov, int *spdecay, double *flag, 
     int *n, int *m, int *T, int *r, int *rT, int *p, int *q, int *N, 
     double *shape_e, double *shape_eta, double *shape_beta, double *shape_l, 
     double *prior_a, double *prior_b, double *mu_beta, double *delta2_beta, 
     double *mu_rho,  double *delta2_rho, double *alpha_l, double *delta2_l, 
     double *phi, double *tau, double *phis, int *phik, double *nu, double *dm, 
     double *dnm, int *constant, double *sig2e, double *sig2eta, double *sig2beta, 
     double *sig2l, double *beta, double *betas, double *rho, double *mu_l, 
     double *X, double *Xsp, double *z, double *w0, double *w, double *phip, 
     double *accept, double *nup, double *sig2ep, double *sig2etap, 
     double *sig2betap, double *betap, double *betasp, double *rhop, 
     double *mu_lp, double *sig2lp, double *w0p, double *wp, double *zfit)
{     

     int n1, m1, p1, N1; 
     n1 = *n;
     m1 = *m;
//     r1 = *r;
     p1 = *p;
//     q1 = *q;
     N1 = *N;
//     col = *constant;
     
   double *Qeta, *XB, *XBno, *XBsp, *Sinv, *S, *det, *A, *C, *Aw;

   Qeta = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
   XB = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBno = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBsp = (double *) malloc((size_t)((N1)*sizeof(double)));
   Sinv = (double *) malloc((size_t)((m1*m1)*sizeof(double)));   
   S = (double *) malloc((size_t)((m1*m1)*sizeof(double)));   
   det = (double *) malloc((size_t)((1)*sizeof(double)));   
   A = (double *) malloc((size_t)((n1*m1)*sizeof(double)));   
   C = (double *) malloc((size_t)((n1*m1)*sizeof(double)));   
   Aw = (double *) malloc((size_t)((N1)*sizeof(double)));   


   covFormat(cov, m, phi, nu, dm, sig2eta, S, det, Sinv, Qeta);
   covF(cov, n, m, phi, nu, dnm, C);

   if(*intercept == 0){
     int i;
     for(i=0; i<N1; i++){
        XBno[i] = 0.0;
     }
     for(i=0; i<p1; i++){
        beta[i] = 0.0;
     }                  
   }
   else{
     MProd(beta, constant, p, X, N, XBno); // N x 1
   }

   MProd(Sinv, m, m, C, n, A);  // n x m
   comb_XB_sp_gpp(n, m, r, T, q, Xsp, betas, A, constant, XBsp); // N x 1
   MAdd(XBno, N, constant, XBsp, XB);  // N x 1

   wlt_gpp_sp(n, m, r, T, rT, p, sig2e, rho, Qeta, A, w0, w, XB, z, 
   constant, wp);     
   w0_gpp_sp(m, r, T, Qeta, sig2l, Sinv, rho, mu_l, wp, constant, w0p);     

// check nu
   if(cov[0]==4){
      nu_gpp_DIS(cov, Qeta, det, phi, nu, m, r, T, rT, dm, rho, sig2eta, 
      mu_l, w0, w, constant, nup); 
   }
   else {
      nup[0] = nu[0];
   }  

// fixed values for phi 
   if(spdecay[0] == 1){
      accept[0] =0.0;
      phip[0] = phi[0];
      covFormat(cov, m, phip, nup, dm, sig2eta, S, det, Sinv, Qeta);
   }
// discrete sampling for phi 
   else if(spdecay[0] == 2){
      phi_gpp_DIS2(cov, Qeta, det, phi, phis, phik, nup, m, r, T, rT, 
      prior_a, prior_b, dm, rho, sig2eta, mu_l, w0p, wp, constant, 
      accept, phip);
      covFormat(cov, m, phip, nup, dm, sig2eta, S, det, Sinv, Qeta);
      MProd(Sinv, m, m, C, n, A);  // n x m
   }
// Random-Walk MH-within-Gibbs sampling for phi
   else if(spdecay[0] == 3){
     double *Qeta2, *det2, *tmp, *phi2;
     Qeta2 = (double *) malloc((size_t)((m1*m1)*sizeof(double))); 
     det2 = (double *) malloc((size_t)((1)*sizeof(double)));   
     tmp = (double *) malloc((size_t)((1)*sizeof(double)));
     phi2 = (double *) malloc((size_t)((1)*sizeof(double)));
     
     if(phi[0] <= 0){
        phi[0] = pow(1,-320);
     }      
     tmp[0] = log(phi[0]); 
     mvrnormal(constant, tmp, tau, constant, phi2);
     phi2[0]= exp(phi2[0]);
     covFormat2(cov, m, phi2, nup, dm, sig2eta, det2, Qeta2);
     phi_gpp_MH2(Qeta, Qeta2, det, det2, phi, phi2, m, r, T, rT, 
     prior_a, prior_b, rho, mu_l, w0p, wp, constant, accept, phip);
     free(Qeta2); free(det2);
     free(tmp); free(phi2);
     if(accept[0]==1.0){
       covFormat(cov, m, phip, nup, dm, sig2eta, S, det, Sinv, Qeta);
       MProd(Sinv, m, m, C, n, A);  // n x m
     }
   }
   else {
//     ;
//     exit(9);
   }   

   MProd(wp, rT, m, A, n, Aw);  // n x rT
   if(*intercept == 0){
     int i;
     for(i=0; i<N1; i++){
        XBno[i] = 0.0;
     }
     for(i=0; i<p1; i++){
        betap[i] = 0.0;
     }                  
   }
   else{
     beta_gpp_for_sp(n, p, rT, N, mu_beta, delta2_beta, sig2e, X, XBsp, Aw, z, 
     constant, betap); // p x 1
     MProd(betap, constant, p, X, N, XBno); // N x 1
   }

   beta_gpp_sp(n, m, q, r, T, rT, N, sig2beta, Sinv, betas, Xsp, XBno, 
   A, Aw, z, constant, betasp); // m x q
   comb_XB_sp_gpp(n, m, r, T, q, Xsp, betasp, A, constant, XBsp); // N x 1
   MAdd(XBno, N, constant, XBsp, XB);  // N x 1

   free(XBsp); free(XBno);

   rho_gpp_sp(m, r, T, rT, p, mu_rho, delta2_rho, Qeta, w0p, wp, constant, rhop);        
   mu_l_gpp(m, r, sig2l, alpha_l, delta2_l, Sinv, w0p, constant, mu_lp);         
   sig_l_gpp(m, r, shape_l, prior_b, mu_lp, Sinv, w0p, constant, sig2lp);          
   sig_eta_gpp(m, r, T, rT, shape_eta, prior_b, Sinv, rhop, wp, w0p, 
   constant, sig2etap);
   sig_e_gpp(n, rT, N, shape_e, prior_b, XB, Aw, z, constant, sig2ep);
   Z_fit_gpp_sp(flag, n, m, T, r, rT, sig2ep, Aw, XB, z, constant, zfit);

   sig_beta_gpp_sp(m, q, shape_beta, prior_b, betasp, Sinv, constant, sig2betap);

   free(Qeta); free(XB); free(Sinv); free(det); free(S); 
   return;
}




// Posterior distribution for "sig_e"
void sig_e_gpp(int *n, int *rT, int *N, double *shape, double *prior_b, 
     double *XB, double *Aw, double *z, int *constant, double *sig2e)
{
     int i, col, N1;
//     n1 =*n;
     col =*constant;
     N1 =*N;
     
     double *tAw, *er, *out;
     tAw = (double *) malloc((size_t)((N1*col)*sizeof(double)));
     er = (double *) malloc((size_t)((N1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));
     
     double u, v, b;
     u = 0.0;
     b = 0.0;

     MTranspose(Aw, rT, n, tAw); // rt x n
     for(i=0; i<N1; i++){
        er[i] = (z[i] - XB[i] - tAw[i]);
     }
     MProd(er, constant, N, er, constant, out);          
     
     b = *prior_b;
     out[0] = b + 0.5 * out[0];
     u = out[0];

     v = *shape;
     sig2e[0] = rigammaa(v, u);

     free(tAw); free(er); free(out);
     return;                  

}  


// Posterior distribution for "sig_eta"
void sig_eta_gpp(int *m, int *r, int *T, int *rT, double *shape, 
     double *prior_b, double *Sinv_eta, double *rho, double *w, 
     double *w0, int *constant, double *sig2eta)     
{
     int m1, col, t, l, i, r1;
     m1 = *m;
     col = *constant;
     r1 = *r;
//     rT1 = *rT;
     
     double *w1, *w2, *er, *out, u, b, sh, sig[1];

     w1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     w2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     er = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*col)*sizeof(double)));

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);
     
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(t=0; t < T1[l]; t++) {                                
             if(t == 0){
                   for(i=0; i<m1; i++){
                   w2[i] = w0[i+l*m1];
                   }      
             }
             else {       
                    for(i=0; i<m1; i++){
                    w2[i] = w[i+(t-1)*m1+m1*T2[l]];         
//                    w2[i] = w[i+(t-1)*m1+l*m1*T1];         
//             extract_alt2(l, t-1, m, rT, T, w, w2);
                    }
             }

             for(i=0; i<m1; i++){
               w1[i] = w[i+t*m1+m1*T2[l]];         
//             extract_alt2(l, t, m, rT, T, w, w1);
             }
            
             for(i=0; i < m1; i++) {
                 er[i]=w1[i]-rho[0]*w2[i];
             }
             MProd(er, constant, m, Sinv_eta, m, out);
             MProd(out, constant, m, er, constant, out);
             u += out[0];       
//             u += xTay2(er, Sinv_eta, er, m1);
         }
     }
     b = *prior_b;
     u = b + 0.5 * u;
     sh = *shape;
     sig[0] = rigammaa(sh, u);
     *sig2eta = sig[0]; 
    
    free(T1); free(T2);
    free(w1); free(w2); free(er); free(out);
    return;
}


// Posterior distribution for "sig_beta" for spatial varying beta
// betasp = m x q, Sinv = m x m
void sig_beta_gpp_sp(int *m, int *q, double *shape, double *prior_b, 
     double *betasp, double *Sinv, int *constant, double *sig2beta)
{
     int i, j, m1, q1;
     m1 =*m;
     q1 =*q;
//     col =*constant;
     
     double *bt;
     bt = (double *) malloc((size_t)((m1)*sizeof(double)));
          
     double u, v, b;
     u = 0.0;
     v = 0.0;
     b = 0.0;

     for(j=0; j<q1; j++){
         for(i=0; i<m1; i++){
           bt[i] = betasp[i+j*m1];
         }
         u += xTay2(bt, Sinv, bt, m1);
     }      

     b = *prior_b;
     u = b + 0.5 * u;
     v = *shape;
     sig2beta[0] = rigammaa(v, u);

     free(bt); 
     return;                  
}     


// conditional for betasp
// here delta2_beta (pxp diag) = for fixed betas
// sig2beta (1) = for spatial betas 
void beta_gpp_sp(int *n, int *m, int *q, int *r, int *T, int *rT, int *N, 
     double *sig2beta, double *Sinv, double *betas, double *Xsp, 
     double *XBno, double *A, double *Aw, double *z, int *constant, 
     double *betasp) // m x q
{
     int i, j, l, t, q1, n1, m1, r1, T1;
     q1 =*q;
     n1 =*n;
     m1 =*m;
     r1 =*r;
     T1 =*T;

     double *del, *chi, *bt, *xs, *XBrest, *X1, *XA, *tXA, *out;
     double *z1, *XBno1, *XBrest1, *chi1, *det;
     del = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     chi = (double *) malloc((size_t)((m1)*sizeof(double)));     
     bt = (double *) malloc((size_t)((m1*(q1-1))*sizeof(double)));     
     xs = (double *) malloc((size_t)((n1*r1*T1*(q1-1))*sizeof(double)));          
     XBrest = (double *) malloc((size_t)((n1*r1*T1)*sizeof(double)));          
     X1 = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     XA = (double *) malloc((size_t)((n1*m1)*sizeof(double)));
     tXA = (double *) malloc((size_t)((n1*m1)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     z1 = (double *) malloc((size_t)((n1)*sizeof(double)));
     XBno1 = (double *) malloc((size_t)((n1)*sizeof(double)));
     XBrest1 = (double *) malloc((size_t)((n1)*sizeof(double)));
     chi1 = (double *) malloc((size_t)((m1)*sizeof(double)));
     det = (double *) malloc((size_t)((1)*sizeof(double)));
               
     int *qq;
     qq = (int *) malloc((size_t)((1)*sizeof(int)));          
    
     for(j=0; j<q1; j++){
       for(i=0; i<m1*m1; i++){
           del[i] = 0.0;
       }   
       for(i=0; i<m1; i++){
           chi[i] = 0.0;
       }   

       // betas, m x q  === m x (q-1)
       // Xsp, N x q === N x (q-1)
       extract_beta_sp2(j, m, q, betas, bt); // m x (q-1)
       extract_beta_sp2(j, N, q, Xsp, xs); // N x (q-1)
       qq[0] = q1-1;
       comb_XB_sp_gpp(n, m, r, T, qq, xs, bt, A, constant, XBrest); // N x 1

       for(l=0; l<r1; l++){
          for(t=0; t<T1; t++){

           extract_X_sp2(t, l, j, n, r, T, Xsp, X1); // n x n diagonal matrix
           MProd(X1, n, n, A, m, XA);  // n x m 
           MTranspose(XA, m, n, tXA); // m x n
           MProd(XA, m, n, tXA, m, out);   // m x m
           MAdd(del, m, m, out, del);  // m x m
           
           extract_alt2(l, t, n, rT, T, z, z1);  // n x 1
           extract_alt2(l, t, n, rT, T, XBno, XBno1); // n x 1        
           if(q1 > 1){
           extract_alt2(l, t, n, rT, T, XBrest, XBrest1); // n x 1        
             for(i=0; i<n1; i++){
                z1[i] = z1[i]-XBno1[i]-XBrest1[i]-Aw[i+t*n1+l*n1*T1];;
             }
           }
           else {
             for(i=0; i<n1; i++){
                z1[i] = z1[i]-XBno1[i]-Aw[i+t*n1+l*n1*T1];;
             }
           }   
           MProd(z1, constant, n, tXA, m, chi1);  // m x 1
           MAdd(chi, m, constant, chi1, chi);  // m x 1
          }
       }
       for(i=0; i<m1*m1; i++){
           del[i] = del[i] + Sinv[i]*(1.0/sig2beta[0]);
       } 
       MInv(del, del, m, det);
       MProd(chi, constant, m, del, m, chi);  // m x 1    
       mvrnormal(constant, chi, del, m, chi); // m x 1
       for(i=0; i<m1; i++){
           betasp[i+j*m1] = chi[i]; 
       }
     }

     free(qq); 
     free(det); free(chi1); free(XBrest1); free(XBno1);  free(z1); free(out); 
     free(tXA); free(XA); free(X1); free(XBrest); free(xs); free(bt); free(chi); 
     free(del); 
     
     return;
}



// Posterior distribution for "beta"
void beta_gpp(int *n, int *p, int *rT, int *N, double *mu_beta, 
     double *delta2_beta, double *sig2e, double *X, double *Aw, 
     double *z, int *constant, double *beta) 
{
     int i, p1, N1, col;
     p1 =*p;
     N1 =*N;
     col =*constant;
     
     double *delta2, *det, *tX, *XX, *md, *tAw, *zA, *Xz, *mean1;
     delta2 = (double *) malloc((size_t)((p1*p1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     tX = (double *) malloc((size_t)((N1*p1)*sizeof(double)));
     XX = (double *) malloc((size_t)((p1*p1)*sizeof(double)));
     md = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     tAw = (double *) malloc((size_t)((N1*col)*sizeof(double)));
     zA = (double *) malloc((size_t)((N1*col)*sizeof(double)));
     Xz = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     mean1 = (double *) malloc((size_t)((p1*col)*sizeof(double)));
    

     MInv(delta2_beta, delta2, p, det);
     MTranspose(X, p, N, tX);
     MProd(X, p, N, tX, p, XX);
     for(i=0; i<p1*p1; i++){
          XX[i] = XX[i]/sig2e[0] + delta2[i];
     }    
              
     MProd(mu_beta, constant, p, delta2, p, md);
     MTranspose(Aw, rT, n, tAw);  // rT x n
     for(i=0; i<N1; i++){
         zA[i] = z[i] - tAw[i];
     }
     MProd(zA, constant, N, tX, p, Xz);
     for(i=0; i<p1; i++){
          Xz[i] = Xz[i]/sig2e[0] + md[i];
     }    
      
     MInv(XX, XX, p, det);
     MProd(Xz, constant, p, XX, p, mean1);
//     for(i=0; i<p1; i++){
//         beta[i] = mean1[i]; // ok
//     }      
     mvrnormal(constant, mean1, XX, p, beta);
 

     free(delta2); free(det); free(tX); free(XX); free(md);
     free(tAw); free(zA); free(Xz); free(mean1);

     return;
}     



// Posterior distribution for "beta" when spatial beta is used
void beta_gpp_for_sp(int *n, int *p, int *rT, int *N, double *mu_beta, 
     double *delta2_beta, double *sig2e, double *X, double *XBsp, double *Aw, 
     double *z, int *constant, double *beta) 
{
     int i, p1, N1, col;
     p1 =*p;
     N1 =*N;
     col =*constant;
     
     double *delta2, *det, *tX, *XX, *md, *tAw, *zA, *Xz, *mean1;
     delta2 = (double *) malloc((size_t)((p1*p1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     tX = (double *) malloc((size_t)((N1*p1)*sizeof(double)));
     XX = (double *) malloc((size_t)((p1*p1)*sizeof(double)));
     md = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     tAw = (double *) malloc((size_t)((N1*col)*sizeof(double)));
     zA = (double *) malloc((size_t)((N1*col)*sizeof(double)));
     Xz = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     mean1 = (double *) malloc((size_t)((p1*col)*sizeof(double)));
   
     MInv(delta2_beta, delta2, p, det);
     MTranspose(X, p, N, tX);
     MProd(X, p, N, tX, p, XX);
     for(i=0; i<p1*p1; i++){
          XX[i] = XX[i]/sig2e[0] + delta2[i];
     }    
              
     MProd(mu_beta, constant, p, delta2, p, md);
     MTranspose(Aw, rT, n, tAw);  // rT x n
     for(i=0; i<N1; i++){
         zA[i] = z[i] - XBsp[i] - tAw[i]; // N x 1
     }
     MProd(zA, constant, N, tX, p, Xz);
     for(i=0; i<p1; i++){
          Xz[i] = Xz[i]/sig2e[0] + md[i];
     }    
      
     MInv(XX, XX, p, det);
     MProd(Xz, constant, p, XX, p, mean1);
     mvrnormal(constant, mean1, XX, p, beta);

     free(delta2); free(det); free(tX); free(XX); free(md);
     free(tAw); free(zA); free(Xz); free(mean1);
     
     return;
}     



// Posterior distribution for "rho"
void rho_gpp(int *m, int *r, int *T, int *rT, int *p, double *mu_rho, 
     double *delta2, double *Q_eta, double *w0, double *w, 
     int *constant, double *rho)     
{
     int m1, col, t, l, i, r1;
     m1 = *m;
     col = *constant;
     r1 = *r;
//     p1= *p;

     double *w2, *w1, *out, *mu, *s2;
     w2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     w1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     s2 = (double *) malloc((size_t)((col)*sizeof(double)));     

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);
          
     double u, v;
     u = 0.0;         
     v = 0.0;   
      for(l=0; l < r1; l++) {
      for(t=0; t < T1[l]; t++) {
             if(t == 0){
                   for(i=0; i<m1; i++){
                     w2[i] = w0[i+l*m1];
                   }      
             }
             else {       
                   for(i=0; i<m1; i++){
                     w2[i] = w[i+(t-1)*m1+m1*T2[l]];
                   }      
             }
             MProd(w2, constant, m, Q_eta, m, out);
             MProd(out, constant, m, w2, constant, out);             
//             xTay(ot1, Q_eta, ot1, n, out);
             v += out[0];
             for(i=0; i<m1; i++){
                 w1[i] = w[i+t*m1+m1*T2[l]];         
             }
             MProd(w1, constant, m, Q_eta, m, out);
             MProd(out, constant, m, w2, constant, out);             
             u += out[0];
      }
      }
     v = v + (1.0/delta2[0]); // OK
     v = 1.0/v;
     u = u + mu_rho[0]/delta2[0];
     u = u*v;

     *mu = u;
     *s2 = v;

     mvrnormal(constant, mu, s2, constant, out);
/*
     out[0] = fabs(out[0]);
     if(out[0] > 1.0){
         out[0] = 1.0;
     }          
*/
     rho[0] = out[0];
          
     free(T1); free(T2);     
     free(w2); free(w1); free(out); free(mu); free(s2);
     return;
}

// Posterior distribution for "rho"
void rho_gpp_sp(int *m, int *r, int *T, int *rT, int *p, double *mu_rho, 
     double *delta2, double *Q_eta, double *w0, double *w, 
     int *constant, double *rho)     
{
     int m1, col, t, l, i, r1, T1;
     m1 = *m;
     col = *constant;
     r1 = *r;
//     p1= *p;
     T1 =*T;
     
     double *w2, *w1, *out, *mu, *s2;
     w2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     w1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     s2 = (double *) malloc((size_t)((col)*sizeof(double)));     

     double u, v;
     u = 0.0;         
     v = 0.0;   
      for(l=0; l < r1; l++) {
      for(t=0; t < T1; t++) {
             if(t == 0){
                   for(i=0; i<m1; i++){
                     w2[i] = w0[i+l*m1];
                   }      
             }
             else {       
                   for(i=0; i<m1; i++){
                     w2[i] = w[i+(t-1)*m1+l*T1];
                   }      
             }
             MProd(w2, constant, m, Q_eta, m, out);
             MProd(out, constant, m, w2, constant, out);             
//             xTay(ot1, Q_eta, ot1, n, out);
             v += out[0];
             for(i=0; i<m1; i++){
                 w1[i] = w[i+t*m1+l*T1];         
             }
             MProd(w1, constant, m, Q_eta, m, out);
             MProd(out, constant, m, w2, constant, out);             
             u += out[0];
      }
      }
     v = v + (1.0/delta2[0]); // OK
     v = 1.0/v;
     u = u + mu_rho[0]/delta2[0];
     u = u*v;

     *mu = u;
     *s2 = v;

     mvrnormal(constant, mu, s2, constant, out);
     rho[0] = out[0];
          
     free(w2); free(w1); free(out); free(mu); free(s2);
     return;
}



// Posterior distribution for "sig2l"
void sig_l_gpp(int *m, int *r, double *shape, double *prior_b, double *mu_l, 
     double *Sinv_0, double *w0, int *constant, double *sig2l)          
{
     int i, l, m1, r1, col;
     
     r1 = *r;
     m1 = *m;
     col = *constant;

     double *er, *out;
     er =(double *) malloc((size_t)((m1*col)*sizeof(double)));
     out =(double *) malloc((size_t)((m1*col)*sizeof(double)));
//     double n1 = (double) row;     
     double rt, b, u, sh;
     sh = *shape;
     for(l=0; l < r1; l++){
	     rt = 0.0;
	     b = 0.0;
	     u = 0.0;
         for(i=0; i < m1; i++) {
           er[i] = w0[i+l*m1] - mu_l[l];
         }
         MProd(er, constant, m, Sinv_0, m, out);
         MProd(out, constant, m, er, constant, out);
         u = out[0];
         b = *prior_b;
         rt = b + 0.5 * u;            
         sig2l[l] = rigammaa(sh, rt);
     }

     free(er); free(out); 
     return;
}


// Posterior distribution for "mu_l"
void mu_l_gpp(int *m, int *r, double *sig2l, double *alpha_l, double *delta2_l, 
     double *Sinv_0, double *w0, int *constant, double *mu_l)         
{
     int i, l, r1, m1, col;
     r1 =*r;
     m1 =*m;
     col =*constant;
          
     double *I, *out, *del, *w, *chi, *mu;
     I = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     del = (double *) malloc((size_t)((col)*sizeof(double)));
     w = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     chi = (double *) malloc((size_t)((col)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
                              
     for(i=0; i <m1; i++){
           I[i] = 1.0;
     }

     for(l=0; l<r1; l++){
         MProd(I, constant, m, Sinv_0, m, out);
         MProd(out, constant, m, I, constant, del);              
         del[0] = del[0]/sig2l[l] + 1.0/delta2_l[l];
         del[0] = 1.0/del[0];
         for(i=0; i<m1; i++){         
            w[i] = w0[i+l*m1];         
         }
         MProd(w, constant, m, Sinv_0, m, out);
         MProd(out, constant, m, I, constant, chi);              
         chi[0] = chi[0]/sig2l[l] + alpha_l[l]/delta2_l[l];
         
         mu[0] = del[0] * chi[0];
         mvrnormal(constant, mu, del, constant, out);
         mu_l[l] = out[0];
     }

     free(I); free(out); free(del); free(w); free(chi); free(mu);
     return;
}


// Posterior samples for Z_lt obtained using MCMC samples
void Z_fit_gpp(double *flag, int *n, int *m, int *T, int *r, int *rT, 
     double *sig2e, double *Aw, double *XB, double *z, int *constant, 
     double *zfit)
{
     int i, l, t, n1, r1, col;
     n1 = *n;
//     m1 = *m;
     r1 = *r;
     col = *constant;
   
     double *XB1, *er, *z1, *fl, *zz;
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     er = (double *) malloc((size_t)((col)*sizeof(double)));
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     fl = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zz = (double *) malloc((size_t)((n1*col)*sizeof(double)));

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);
       
     for(l=0; l < r1; l++) {
       for(t=0; t < T1[l]; t++) {

         extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
         extract_alt_uneqT(l, t, n, r, T, rT, z, z1);
         extract_alt_uneqT(l, t, n, r, T, rT, flag, fl);         
//         extract_alt2(l, t, n, rT, T, XB, XB1);
//         extract_alt2(l, t, n, rT, T, z, z1);
//         extract_alt2(l, t, n, rT, T, flag, fl);

         er[0] = 0.0;
         mvrnormal(constant, er, sig2e, constant, er);
         for(i=0; i<n1; i++){
             if(fl[i] == 1.0){
             mvrnormal(constant, er, sig2e, constant, er);
             zz[i]=XB1[i]+ Aw[i+t*n1+n1*T2[l]] + er[0];
             } 
             else {
             zz[i] = XB1[i] + Aw[i+t*n1+n1*T2[l]] + er[0]; 
             }
         }         
         put_together1_uneqT(l, t, n, r, T, rT, zfit, zz);
//         put_together1(l, t, n, r, T, zfit, zz);
       }
     }

     free(T1); free(T2);
     free(XB1); free(er); free(z1); free(fl); free(zz);
     return;
}



// Posterior samples for Z_lt obtained using MCMC samples
void Z_fit_gpp_sp(double *flag, int *n, int *m, int *T, int *r, int *rT, 
     double *sig2e, double *Aw, double *XB, double *z, int *constant, 
     double *zfit)
{
     int i, l, t, n1, T1, r1, col;
     n1 = *n;
     T1 =*T;
//     m1 = *m;
     r1 = *r;
     col = *constant;
   
     double *XB1, *er, *z1, *fl, *zz;
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     er = (double *) malloc((size_t)((col)*sizeof(double)));
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     fl = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zz = (double *) malloc((size_t)((n1*col)*sizeof(double)));

      
     for(l=0; l < r1; l++) {
       for(t=0; t < T1; t++) {

         extract_alt2(l, t, n, rT, T, XB, XB1);
         extract_alt2(l, t, n, rT, T, z, z1);
         extract_alt2(l, t, n, rT, T, flag, fl);

         *er = 0.0;
         mvrnormal(constant, er, sig2e, constant, er);
         for(i=0; i<n1; i++){
             if(fl[i] == 1.0){
             mvrnormal(constant, er, sig2e, constant, er);
             zz[i]=XB1[i]+ Aw[i+t*n1+l*T1] + er[0];
             } 
             else {
             zz[i] = XB1[i] + Aw[i+t*n1+l*T1] + er[0]; 
             }
         }         
         put_together1(l, t, n, r, T, zfit, zz);
       }
     }

     free(XB1); free(er); free(z1); free(fl); free(zz);
     return;
}


// Posterior samples for o_lt=XB+Aw obtained using MCMC samples
void o_fit_gpp(double *flag, int *n, int *m, int *T, int *r, int *rT, 
     double *Aw, double *XB, double *z, int *constant, 
     double *zfit)
{
     int i, l, t, n1, r1, col;
     n1 = *n;
//     m1 = *m;
     r1 = *r;
     col = *constant;
   
     double *XB1, *z1, *fl, *zz;
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     fl = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zz = (double *) malloc((size_t)((n1*col)*sizeof(double)));

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);
       
     for(l=0; l < r1; l++) {
       for(t=0; t < T1[l]; t++) {

         extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
         extract_alt_uneqT(l, t, n, r, T, rT, z, z1);
         extract_alt_uneqT(l, t, n, r, T, rT, flag, fl);         
//         extract_alt2(l, t, n, rT, T, XB, XB1);
//         extract_alt2(l, t, n, rT, T, z, z1);
//         extract_alt2(l, t, n, rT, T, flag, fl);

         for(i=0; i<n1; i++){
             if(fl[i] == 1.0){
             zz[i]=XB1[i]+ Aw[i+t*n1+n1*T2[l]];
             } 
             else {
             zz[i] = XB1[i] + Aw[i+t*n1+n1*T2[l]]; 
             }
         }         
         put_together1_uneqT(l, t, n, r, T, rT, zfit, zz);
//         put_together1(l, t, n, r, T, zfit, zz);
       }
     }

     free(T1); free(T2);
     free(XB1); free(z1); free(fl); free(zz);
     return;
}




// approach without simplifications
void w0_gpp(int *m, int *r, int *T, double *Q_eta, double *sig2l, 
     double *Sinv_0, double *rho, double *mu_l, double *w, 
     int *constant, double *w0p)     
{     
     int i, l, m1, r1, col; 
     m1 = *m;
     r1 = *r;
     col = *constant;
     
     double *I, *del, *det, *w1, *pt1, *pt2, *chi, *out;
     I = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     del = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     w1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     pt1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     pt2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     chi = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     
     for(i=0; i<m1; i++){
           I[i] = 1.0;
     }   

     int *T2; 
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     cumsumint(r, T, T2);
                                        
     for(l=0; l<r1; l++){
       for(i=0; i<m1*m1; i++){
           del[i] = Q_eta[i]*rho[0]*rho[0]+Sinv_0[i]/sig2l[l];
       }
       MInv(del, del, m, det);
       
       for(i=0; i<m1; i++){
             w1[i] = w[i+m1*T2[l]];
       }         
       MProd(w1, constant, m, Q_eta, m, pt1);
       MProd(I, constant, m, Sinv_0, m, pt2);
              
       for(i=0; i<m1; i++){
          chi[i]=pt1[i]*rho[0]+pt2[i]*(mu_l[l]/sig2l[l]);
       }
       
       MProd(chi, constant, m, del, m, out);
       mvrnormal(constant, out, del, m, out);

       for(i=0; i<m1; i++){
         w0p[i+l*m1]=out[i];
//         w0p[i+l*m1]=chi[i];
       }
     }

     free(T2);     
     free(I); free(del); free(det); free(w1); free(pt1); 
     free(pt2); free(chi); free(out);
     return;
}     
     


// approach without simplifications
void w0_gpp_sp(int *m, int *r, int *T, double *Q_eta, double *sig2l, 
     double *Sinv_0, double *rho, double *mu_l, double *w, 
     int *constant, double *w0p)     
{     
     int i, l, m1, r1, T1, col; 
     m1 = *m;
     r1 = *r;
     T1 = *T;
     col = *constant;
     
     double *I, *del, *det, *w1, *pt1, *pt2, *chi, *out;
     I = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     del = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     w1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     pt1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     pt2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     chi = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     
     for(i=0; i<m1; i++){
           I[i] = 1.0;
     }   

     for(l=0; l<r1; l++){
       for(i=0; i<m1*m1; i++){
           del[i] = Q_eta[i]*rho[0]*rho[0]+Sinv_0[i]/sig2l[l];
       }
       MInv(del, del, m, det);
       
       for(i=0; i<m1; i++){
             w1[i] = w[i+l*T1];
       }         
       MProd(w1, constant, m, Q_eta, m, pt1);
       MProd(I, constant, m, Sinv_0, m, pt2);
              
       for(i=0; i<m1; i++){
          chi[i]=pt1[i]*rho[0]+pt2[i]*(mu_l[l]/sig2l[l]);
       }
       
       MProd(chi, constant, m, del, m, out);
       mvrnormal(constant, out, del, m, out);

       for(i=0; i<m1; i++){
         w0p[i+l*m1]=out[i];
//         w0p[i+l*m1]=chi[i];
       }
     }

     free(I); free(del); free(det); free(w1); free(pt1); 
     free(pt2); free(chi); free(out);
     return;
}     
     
   
void wlt_gpp(int *n, int *m, int *r, int *T, int *rT, int *p, double *sig2e, 
     double *rho, double *Q_eta, double *A, double *w0, double *w, double *XB, 
     double *z, int *constant, double *wp)     
{
     int i, l, t, r1, n1, m1, mm, col; 
     r1 =*r; n1 =*n; m1 =*m; mm =m1*m1; col =*constant;
     
     double *tA, *AA, *de_tT, *de_T, *det, *w2, *w1;
     double *tr, *tr1, *tr2, *tr3, *z1, *XB1, *chi, *mean1, *ww;

     tA = (double *) malloc((size_t)((m1*n1)*sizeof(double)));
     AA = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     de_tT = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     de_T = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     w2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     w1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));

     tr = (double *) malloc((size_t)((n1*col)*sizeof(double)));

     tr1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     tr2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     tr3 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     chi = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     mean1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     ww = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     
     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);
     
         MTranspose(A, m, n, tA);  // m x n 
         MProd(A, m, n, tA, m, AA); // m x m

// for 1 <= t < T, the delta part
         for(i=0; i < mm; i++) {
            de_tT[i] = ((AA[i]/sig2e[0]) + Q_eta[i] + rho[0]*rho[0]*Q_eta[i]);
         }    
         MInv(de_tT, de_tT, m, det); 

// for t = T, the delta part
         for(i=0; i < mm; i++) {
             de_T[i] = ((AA[i]/sig2e[0]) + Q_eta[i]);       
         }    
         MInv(de_T, de_T, m, det);

// LOOP begains
     for(l=0; l < r1; l++) {
         if(T1[l] == 1){
             // LOOP for t = 1 
             t = 0;
             for(i=0; i<m1; i++){
               w2[i] = w0[i+l*m1];
             }      
             MProd(w2, constant, m, Q_eta, m, tr2);
             extract_alt_uneqT(l, t, n, r, T, rT, z, z1);
             extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
//             extract_alt2(l, t, n, rT, T, z, z1);
//             extract_alt2(l, t, n, rT, T, XB, XB1);
             for(i=0; i<n1; i++){
                  tr[i] = z1[i] - XB1[i];
             }         
             MProd(tr, constant, n, tA, m, tr1);
             for(i=0; i<m1; i++){
                 chi[i] = tr1[i]/sig2e[0]+tr2[i]*rho[0]; // fixed
             }
             MProd(chi, constant, m, de_tT, m, mean1);                     
             mvrnormal(constant, mean1, de_tT, m, ww);     // random generator
             for(i=0; i<m1; i++) {
                wp[i+t*m1+m1*T2[l]] = ww[i];     
             }
         } // end of first if
         else {
             // LOOP for 1 <= t < T
             for(t=0; t < (T1[l]-1); t++) {          
                if(t == 0){
                    for(i=0; i<m1; i++){
                    w2[i] = w0[i+l*m1];
                    }      
                }
                else {       
                    for(i=0; i<m1; i++){
                    w2[i] = w[i+(t-1)*m1+m1*T2[l]];         
                    }
                }
                for(i=0; i<m1; i++){
                    w1[i] = w[i+(t+1)*m1+m1*T2[l]];         
                }
                MProd(w2, constant, m, Q_eta, m, tr2);
                MProd(w1, constant, m, Q_eta, m, tr3);
                extract_alt_uneqT(l, t, n, r, T, rT, z, z1);
                extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
//               extract_alt2(l, t, n, rT, T, z, z1);
//               extract_alt2(l, t, n, rT, T, XB, XB1);
                for(i=0; i<n1; i++){
                    tr[i] = z1[i] - XB1[i];
                }         
                MProd(tr, constant, n, tA, m, tr1);
                for(i=0; i<m1; i++){
                    chi[i] = tr1[i]/sig2e[0]+tr2[i]*rho[0]+tr3[i]*rho[0];
                }
                MProd(chi, constant, m, de_tT, m, mean1);                     
                mvrnormal(constant, mean1, de_tT, m, ww);     // random generator
                for(i=0; i<m1; i++) {
                   wp[i+t*m1+m1*T2[l]] = ww[i];     
                }         
             } // end of for loop: 1 <= t < T
             // LOOP for t = T, 
             t = (T1[l]-1);
             for(i=0; i<m1; i++){
                w2[i] = w[i+(t-1)*m1+m1*T2[l]];         
             }
             MProd(w2, constant, m, Q_eta, m, tr2);
             extract_alt_uneqT(l, t, n, r, T, rT, z, z1);
             extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
//             extract_alt2(l, t, n, rT, T, z, z1);
//             extract_alt2(l, t, n, rT, T, XB, XB1);
             for(i=0; i<n1; i++){
                  tr[i] = z1[i] - XB1[i];
             }         
             MProd(tr, constant, n, tA, m, tr1);
             for(i=0; i<m1; i++){
                 chi[i] = tr1[i]/sig2e[0]+tr2[i]*rho[0];
             }
             MProd(chi, constant, m, de_T, m, mean1);                     
             mvrnormal(constant, mean1, de_T, m, ww);     // random generator
             for(i=0; i<m1; i++) {
                wp[i+t*m1+m1*T2[l]] = ww[i];     
             }
         } // End of else loop        
     } // End of loop year "l"

     free(T1); free(T2);
     free(tA); free(AA); free(de_tT); free(de_T); free(det);
     free(w2); free(w1); free(tr); free(tr1); free(tr2); free(tr3); 
     free(z1); free(XB1); free(chi); free(mean1); free(ww);
     return;
} 


   
void wlt_gpp_sp(int *n, int *m, int *r, int *T, int *rT, int *p, double *sig2e, 
     double *rho, double *Q_eta, double *A, double *w0, double *w, double *XB, 
     double *z, int *constant, double *wp)     
{
     int i, l, t, r1, n1, m1, mm, col, T1; 
     r1 =*r; n1 =*n; m1 =*m; mm =m1*m1; T1 =*T; col =*constant;
     
     double *tA, *AA, *de_tT, *de_T, *det, *w2, *w1;
     double *tr, *tr1, *tr2, *tr3, *z1, *XB1, *chi, *mean1, *ww;

     tA = (double *) malloc((size_t)((m1*n1)*sizeof(double)));
     AA = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     de_tT = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     de_T = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     w2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     w1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));

     tr = (double *) malloc((size_t)((n1*col)*sizeof(double)));

     tr1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     tr2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     tr3 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     chi = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     mean1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     ww = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     
         MTranspose(A, m, n, tA);  // m x n 
         MProd(A, m, n, tA, m, AA); // m x m

// for 1 <= t < T, the delta part
         for(i=0; i < mm; i++) {
            de_tT[i] = ((AA[i]/sig2e[0]) + Q_eta[i] + rho[0]*rho[0]*Q_eta[i]);
         }    
         MInv(de_tT, de_tT, m, det); 

// for t = T, the delta part
         for(i=0; i < mm; i++) {
             de_T[i] = ((AA[i]/sig2e[0]) + Q_eta[i]);       
         }    
         MInv(de_T, de_T, m, det);

// LOOP begains
     for(l=0; l < r1; l++) {
             // LOOP for 1 <= t < T
             for(t=0; t < (T1-1); t++) {          
                if(t == 0){
                    for(i=0; i<m1; i++){
                    w2[i] = w0[i+l*m1];
                    }      
                }
                else {       
                    for(i=0; i<m1; i++){
                    w2[i] = w[i+(t-1)*m1+l*T1];         
                    }
                }
                for(i=0; i<m1; i++){
                    w1[i] = w[i+(t+1)*m1+l*T1];         
                }
                MProd(w2, constant, m, Q_eta, m, tr2);
                MProd(w1, constant, m, Q_eta, m, tr3);
               extract_alt2(l, t, n, rT, T, z, z1);
               extract_alt2(l, t, n, rT, T, XB, XB1);
                for(i=0; i<n1; i++){
                    tr[i] = z1[i] - XB1[i];
                }         
                MProd(tr, constant, n, tA, m, tr1);
                for(i=0; i<m1; i++){
                    chi[i] = tr1[i]/sig2e[0]+tr2[i]+tr3[i];
                }
                MProd(chi, constant, m, de_tT, m, mean1);                     
                mvrnormal(constant, mean1, de_tT, m, ww);     // random generator
                for(i=0; i<m1; i++) {
                   wp[i+t*m1+l*T1] = ww[i];     
                }         
             } // end of for loop: 1 <= t < T
             // LOOP for t = T, 
             t = (T1-1);
             for(i=0; i<m1; i++){
                w2[i] = w[i+(t-1)*m1+l*T1];         
             }
             MProd(w2, constant, m, Q_eta, m, tr2);
             extract_alt2(l, t, n, rT, T, z, z1);
             extract_alt2(l, t, n, rT, T, XB, XB1);
             for(i=0; i<n1; i++){
                  tr[i] = z1[i] - XB1[i];
             }         
             MProd(tr, constant, n, tA, m, tr1);
             for(i=0; i<m1; i++){
                 chi[i] = tr1[i]/sig2e[0]+tr2[i];
             }
             MProd(chi, constant, m, de_T, m, mean1);                     
             mvrnormal(constant, mean1, de_T, m, ww);     // random generator
             for(i=0; i<m1; i++) {
                wp[i+t*m1+l*T1] = ww[i];     
             }
     } // End of loop year "l"

     free(tA); free(AA); free(de_tT); free(de_T); free(det);
     free(w2); free(w1); free(tr); free(tr1); free(tr2); free(tr3); 
     free(z1); free(XB1); free(chi); free(mean1); free(ww);
     return;
} 



// phi sample random walk 
void phi_gpp_MH2(double *Qeta1, double *Qeta2, double *det1, double *det2,
     double *phi1, double *phi2,   
     int *m, int *r, int *T, int *rT, double *prior_a, double *prior_b, 
     double *rho, double *mu_l, double *w0, double *w, int *constant, 
     double *accept, double *phip) 
{
     int m1, col, t, l, i, r1, rT1;
     m1 = *m;
     col = *constant;
     r1 = *r;
     rT1 = *rT;

     double *w2, *out, *er;
     w2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     er = (double *) malloc((size_t)((m1*col)*sizeof(double)));

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);
     
     double u, v;
     u = 0.0;         
     v = 0.0;   

     for(l=0; l < r1; l++) {
     for(t=0; t < T1[l]; t++) {
             if(t == 0){
                   for(i=0; i<m1; i++){
                     w2[i] = w0[i+l*m1];
    		         er[i] = w[i+t*m1+m1*T2[l]]-rho[0]*w2[i];	
                   }      
             }
             else {       
                   for(i=0; i<m1; i++){
                     w2[i] = w[i+(t-1)*m1+m1*T2[l]];
  		             er[i] = w[i+t*m1+m1*T2[l]]-rho[0]*w2[i];	
                   }      
             }
             MProd(er, constant, m, Qeta1, m, out);
             MProd(out, constant, m, er, constant, out);             
             u += out[0];
             MProd(er, constant, m, Qeta2, m, out);
             MProd(out, constant, m, er, constant, out);             
             v += out[0];
     }
     }

     free(T1); free(T2); 
     free(w2); free(out);

     u =  0.5 * u;
     v =  0.5 * v;

     double *ratio, *U;
     ratio = (double *) malloc((size_t)((col)*sizeof(double)));             
     U = (double *) malloc((size_t)((col)*sizeof(double)));         

     double a, b;
     a = *prior_a;
     b = *prior_b;

     if(det1[0] <= 0){
        det1[0] = pow(1,-320);
     }
     if(det2[0] <= 0){
        det2[0] = pow(1,-320);
     }                   
     if(phi1[0] <= 0){
        phi1[0] = pow(1,-320);
     }        
     if(phi2[0] <= 0){
        phi2[0] = pow(1,-320);
     }      


     if(phi2[0] < 0.0001){
          phip[0] = phi1[0];
          accept[0] = 0.0;
     }
     //else if(phi2[0] > 0.9999){
     //     phip[0] = phi1[0];
     //     accept[0] = 0.0;
     //}
     else{

// with Gamma prior    
     double tr1, tr2; //, u1;
   
     tr1 = (a-1.0)*log(phi1[0])-b*phi1[0]-0.5*rT1*log(det1[0])- u; 
     tr2 = (a-1.0)*log(phi2[0])-b*phi2[0]-0.5*rT1*log(det2[0])- v; 

//     double tr;
//     tr = tr2 + exp(tr2) - tr1 - exp(tr1);
//     Rprintf("for phi1: %f for phi2: %f\n", tr1, tr2);
//     *ratio = exp(tr2 + exp(tr2) - tr1 - exp(tr1));
//       u1 = tr2 - tr1 ; /* log density ratio */
//       u 1+= log(phi2[0]) - log(phi1[0]);
//       *ratio = exp(u1); 
     *ratio = exp(tr2 - tr1 + log(phi2[0]) - log(phi1[0]));
//     *ratio = exp(tr);
     ratio_fnc(ratio, constant, U);
     if(U[0] < ratio[0]){
          phip[0] = phi2[0];
          accept[0] = 1.0;
     }             
     else {
          phip[0] = phi1[0];
          accept[0] = 0.0;
     }     
     }

     free(ratio); free(U);
     return;   
}


// phi sample random walk 
void phi_gpp_DIS2(int *cov, double *Qeta1, double *det1, double *phi1,    
     double *phis, int *phik, double *nu, int *m, int *r, int *T, int *rT, 
     double *prior_a, double *prior_b, double *dm, double *rho, 
     double *sig2eta, double *mu_l, double *w0, double *w, int *constant, 
     double *accept, double *phip) 
{
    
     int row, col, i, pk;
     row = *m;
     col = *constant;
//     r1 = *r;
//     rT1 = *rT;
//     N1 = row*rT1;
     pk = *phik;

     double *phitmp, *pden, *Qeta, *det, *out;
     phitmp = (double *) malloc((size_t)((col)*sizeof(double)));             
     pden = (double *) malloc((size_t)((pk)*sizeof(double)));             
     Qeta = (double *) malloc((size_t)((row*row)*sizeof(double)));             
     det = (double *) malloc((size_t)((col)*sizeof(double)));             
     out = (double *) malloc((size_t)((col)*sizeof(double))); 

     double u;
     u =0.0;     
          
     for(i=0; i< *phik; i++){
        phitmp[0] = phis[i];
        covFormat2(cov, m, phitmp, nu, dm, sig2eta, det, Qeta);
     	phiden_gpp(phitmp, Qeta, det, m, r, T, rT, prior_a, prior_b,
        rho, w0, w, constant, out);
        pden[i] = out[0];
        u += out[0];
     }     
     free(phitmp); free(Qeta); free(det); free(out);

     double *pprob, *U, *tr2;
     pprob = (double *) malloc((size_t)((pk)*sizeof(double)));             
     U = (double *) malloc((size_t)((col)*sizeof(double)));
     tr2 = (double *) malloc((size_t)((col)*sizeof(double)));             

     pprob[0] = pden[0]/u;         

     for(i=0; i< (pk-1); i++){
        pprob[i+1] = pprob[i] + pden[i+1]/u;
     }
     runif_val(constant, constant, U);
     if ( U[0] >  pprob[0]){
     i = 0 ;
     do{
       i = i + 1;
       } while ( ( U[0] > pprob[i] ) & ( i< pk - 1 ) ) ;   
     }
     else i=0;
     tr2[0] = pden[i];  
     
     free(pprob);

     double *ratio, *tr1;
     ratio = (double *) malloc((size_t)((col)*sizeof(double)));             
     tr1 = (double *) malloc((size_t)((col)*sizeof(double)));         

     phiden_gpp(phi1, Qeta1, det1, m, r, T, rT, prior_a, prior_b,
     rho, w0, w, constant, tr1);
     ratio[0] = exp(tr2[0] + exp(tr2[0]) - tr1[0] - exp(tr1[0]));
     ratio_fnc(ratio, constant, U);
     if(U[0] < ratio[0]){
          phip[0] = phis[i];
          accept[0] = 0.0;
     }             
     else {
        phip[0] = phi1[0];
        accept[0] = 0.0;
     }     
 
     free(ratio); free(tr2); free(tr1); free(pden); free(U);
     return;
}



// phi density for the gpp
void phiden_gpp(double *phi, double *Qeta, double *det, int *m, int *r, 
     int *T, int *rT, double *prior_a, double *prior_b, double *rho, 
     double *w0, double *w, int *constant, double *out)
{
     int m1, col, t, l, i, r1, rT1;
     m1 = *m;
     col = *constant;
     r1 = *r;
     rT1 = *rT;

     double *w2, *er;
     w2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     er = (double *) malloc((size_t)((m1*col)*sizeof(double)));

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);

     double u, a, b;
     u = 0.0;         
     for(l=0; l < r1; l++) {
     for(t=0; t < T1[l]; t++) {
             if(t == 0){
                   for(i=0; i<m1; i++){
                     w2[i] = w0[i+l*m1];
                     er[i] = w[i+t*m1+m1*T2[l]]-rho[0]*w2[i];	
                   }      
             }
             else {       
                   for(i=0; i<m1; i++){
                     w2[i] = w[i+(t-1)*m1+m1*T2[l]];
        		     er[i] = w[i+t*m1+m1*T2[l]]-rho[0]*w2[i];	
                   }      
             }
             u += xTay2(er, Qeta, er, m1);
     }
     }
     free(T1); free(T2);
     free(w2); free(er);

     u =  0.5 * u;
     a = *prior_a;
     b = *prior_b;
     if(det[0] <= 0){
        det[0] = pow(1,-320);
     }
     if(phi[0] <= 0){
        phi[0] = pow(1,-320);
     }
     double tr;
     tr = 0.0;        
     tr = (a-1.0)*log(phi[0])-b*phi[0]-0.5*rT1*log(det[0])-u; 
     out[0] = tr;

     return;
}



// nu sample random walk 
void nu_gpp_DIS(int *cov, double *Qeta1, double *det1, double *phi, double *nu1, 
     int *m, int *r, int *T, int *rT, double *dm, double *rho, double *sig2eta, 
     double *mu_l, double *w0, double *w, int *constant, double *nup) 
{
     int row, col, i;
     row = *m;
     col = *constant;
//     r1 = *r;
//     rT1 = *rT;
//     N1 = row*rT1;

     int nuk;

     nuk=30;
     double *nus;
     nus = (double *) malloc((size_t)((nuk)*sizeof(double)));             
     nus[0]=0.05; nus[1]=0.10; nus[2]=0.15; nus[3]=0.20; nus[4]=0.25;  
     nus[5]=0.30; nus[6]=0.35; nus[7]=0.40; nus[8]=0.45; nus[9]=0.50;  
     nus[10]=0.55; nus[11]=0.60; nus[12]=0.65; nus[13]=0.70; nus[14]=0.75; 
     nus[15]=0.80; nus[16]=0.85; nus[17]=0.90; nus[18]=0.95; nus[19]=1.0; 
     nus[20]=1.05; nus[21]=1.10; nus[22]=1.15; nus[23]=1.20; nus[24]=1.25; 
     nus[25]=1.30; nus[26]=1.35; nus[27]=1.40; nus[28]=1.45; nus[29]=1.50;      


//     nuk=20;
//     double *nus;
//     nus = (double *) malloc((size_t)((nuk)*sizeof(double)));             
//     nus[0]=0.05; nus[1]=0.10; nus[2]=0.15; nus[3]=0.20; nus[4]=0.25;  
//     nus[5]=0.30; nus[6]=0.35; nus[7]=0.40; nus[8]=0.45; nus[9]=0.50;  
//     nus[10]=0.55; nus[11]=0.60; nus[12]=0.65; nus[13]=0.70; nus[14]=0.75; 
//     nus[15]=0.80; nus[16]=0.85; nus[17]=0.90; nus[18]=0.95; nus[19]=1.0; 


     double *nutmp, *pden, *Qeta, *det, *out;
     nutmp = (double *) malloc((size_t)((col)*sizeof(double)));             
     pden = (double *) malloc((size_t)((nuk)*sizeof(double)));             
     Qeta = (double *) malloc((size_t)((row*row)*sizeof(double)));             
     det = (double *) malloc((size_t)((col)*sizeof(double)));             
     out = (double *) malloc((size_t)((col)*sizeof(double))); 
     double u;
     u =0.0;     
          
     for(i=0; i< nuk; i++){
        nutmp[0] = nus[i];
        covFormat2(cov, m, phi, nutmp, dm, sig2eta, det, Qeta);
        nuden_gpp(Qeta, det, m, r, T, rT, rho, w0, w, constant, out);
        pden[i] = out[0];
        u += out[0];
     }     
     free(nutmp); free(Qeta); free(det); free(out);

     double *pprob, *U, *tr2;
     pprob = (double *) malloc((size_t)((nuk)*sizeof(double)));             
     U = (double *) malloc((size_t)((col)*sizeof(double)));
     tr2 = (double *) malloc((size_t)((col)*sizeof(double)));             

     pprob[0] = pden[0]/u;         

     for(i=0; i< (nuk-1); i++){
        pprob[i+1] = pprob[i] + pden[i+1]/u;
     }
     runif_val(constant, constant, U);
     if ( U[0] >  pprob[0]){
     i = 0 ;
     do{
       i = i + 1;
       } while ( ( U[0] > pprob[i] ) & ( i< nuk - 1 ) ) ;   
     }
     else i=0;
     tr2[0] = pden[i];  
     
     free(pprob);

     double *ratio, *tr1;
     ratio = (double *) malloc((size_t)((col)*sizeof(double)));             
     tr1 = (double *) malloc((size_t)((col)*sizeof(double)));         

     nuden_gpp(Qeta1, det1, m, r, T, rT, rho, w0, w, constant, tr1);
     ratio[0] = exp(tr2[0] + exp(tr2[0]) - tr1[0] - exp(tr1[0]));
     ratio_fnc(ratio, constant, U);
     if(U[0] < ratio[0]){
          nup[0] = nus[i];
     }             
     else {
        nup[0] = nu1[0];
     }     
 
     free(ratio); free(tr2); free(tr1); free(pden); free(U);
     return;
}



// nu density for the gpp
void nuden_gpp(double *Qeta, double *det, int *m, int *r, int *T, int *rT, 
     double *rho, double *w0, double *w, int *constant, double *out)
{
     int m1, col, t, l, i, r1, rT1;
     m1 = *m;
     col = *constant;
     r1 = *r;
     rT1 = *rT;

     double *w2, *er;
     w2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));
     er = (double *) malloc((size_t)((m1*col)*sizeof(double)));

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     cumsumint(r, T, T2);
          
     double u;
     u = 0.0;         
     for(l=0; l < r1; l++) {
     for(t=0; t < T1[l]; t++) {
             if(t == 0){
                   for(i=0; i<m1; i++){
                     w2[i] = w0[i+l*m1];
                     er[i] = w[i+t*m1+m1*T2[l]]-rho[0]*w2[i];	
                   }      
             }
             else {       
                   for(i=0; i<m1; i++){
                     w2[i] = w[i+(t-1)*m1+m1*T2[l]];
                     er[i] = w[i+t*m1+m1*T2[l]]-rho[0]*w2[i];	
                   }      
             }
             u += xTay2(er, Qeta, er, m1);
     }
     }
     free(T1); free(T2);
     free(w2); free(er);

     u =  0.5 * u;
     if(det[0] <= 0){
        det[0] = pow(1,-320);
     }
     double tr;
     tr = 0.0;        
     tr = -0.5*rT1*log(det[0])-u; // uniform prior log(1) = 0
     out[0] = tr;

     return;
}



// phi sampling random walk another approach
void phi_gpp_MH(int *cov, double *phi2, double *nu, double *dm, double *dnm,
     double *Sinv1, double *det1, double *phi1, double *A1,   
     int *n, int *m, int *r, int *T, int *rT, double *prior_a, double *prior_b, 
     double *rho, double *mu_l, double *w0, double *w, double *z, double *XB, 
     int *constant, double *accept, double *phip) 
{
     int n1, m1, r1, rT1, col, i, l, t; 
     n1 = *n;
     m1 = *m;
     r1 = *r;
     rT1 =*rT;
     col = *constant;
     
   double *Sinv2, *det2, *C2, *A2; 
   Sinv2 = (double *) malloc((size_t)((m1*m1)*sizeof(double)));
   det2 = (double *) malloc((size_t)((1)*sizeof(double)));
   C2 = (double *) malloc((size_t)((m1*n1)*sizeof(double)));
   A2 = (double *) malloc((size_t)((m1*n1)*sizeof(double)));

    covF(cov, m, m, phi2, nu, dm, Sinv2);
    MInv(Sinv2, Sinv2, n, det2);    
    covF(cov, n, m, phi2, nu, dnm, C2);
    MProd(Sinv2, m, m, C2, n, A2);  // n x m

     double *tA1, *tA2, *XB1, *z1, *zx, *er1, *er2; 
     tA1 = (double *) malloc((size_t)((m1*n1)*sizeof(double)));   
     tA2 = (double *) malloc((size_t)((m1*n1)*sizeof(double)));   
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));   
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));   
     zx = (double *) malloc((size_t)((n1*col)*sizeof(double)));   
     er1 = (double *) malloc((size_t)((m1*col)*sizeof(double)));   
     er2 = (double *) malloc((size_t)((m1*col)*sizeof(double)));   

     double *AA1,*AA2;
     AA1 = (double *) malloc((size_t)((m1*m1)*sizeof(double)));   
     AA2 = (double *) malloc((size_t)((m1*m1)*sizeof(double)));   
     
     int *T1; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     
     double a, b, u, v;        
     u = 0.0;
     v = 0.0;
     a = *prior_a;
     b = *prior_b;
     MTranspose(A1, m, n, tA1);  // m x n 
     MProd(A1, m, n, tA1, m, AA1); // m x m 

     MTranspose(A2, m, n, tA2);  // m x n 
     MProd(A2, m, n, tA2, m, AA2); // m x m 

     for(l=0; l < r1; l++) {
         for(t=0; t < T1[l]; t++) {                                
           extract_alt_uneqT(l, t, n, r, T, rT, z, z1);
           extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
//           extract_alt2(l, t, n, rT, T, XB, XB1); // n x 1
//           extract_alt2(l, t, n, rT, T, z, z1); // n x 1

           for(i=0; i<n1; i++){
             zx[i] = z1[i] - XB1[i]; //  n x 1                 
           }      
           MProd(zx, constant, n, tA1, m, er1); // m x 1
//          MProd(er1, constant, m, AA1, m, er1); // m x 1
           MProd(zx, constant, n, tA2, m, er2); // m x 1
//           MProd(er2, constant, m, AA2, m, er2); // m x 1

           u += xTay2(er1, Sinv1, er2, m1);
           v += xTay2(er2, Sinv2, er2, m1);
         }
     }
     u =  0.5 * u;
     v =  0.5 * v;

    free(T1); 
    free(Sinv2); free(C2); free(A2); free(tA1); free(tA2);
    free(XB1); free(z1); free(zx); free(er1); free(er2);
    free(AA1); free(AA2);

     double *ratio, *U;
     ratio = (double *) malloc((size_t)((col)*sizeof(double)));             
     U = (double *) malloc((size_t)((col)*sizeof(double)));         

     double tr1, tr2;
     tr1 = 0.0;
     tr2 = 0.0;

     if(det1[0] <= 0){
        det1[0] = pow(1,-320);
     }
     if(det2[0] <= 0){
        det2[0] = pow(1,-320);
     }                   
     if(phi1[0] <= 0){
        phi1[0] = pow(1,-320);
     }        
     if(phi2[0] <= 0){
        phi2[0] = pow(1,-320);
     }      
          

     if(phi2[0] < 0.0001){
          phip[0] = phi1[0];
          accept[0] = 0.0;
     }
     //else if(phi2[0] > 0.9999){
     //     phip[0] = phi1[0];
     //     accept[0] = 0.0;
     //}
     else{
// with Gamma prior    
     tr1 = (a-1.0)*log(phi1[0])-b*phi1[0]-0.5*rT1*log(det1[0])- u;      
     tr2 = (a-1.0)*log(phi2[0])-b*phi2[0]-0.5*rT1*log(det2[0])- v;        
//     Rprintf("for phi1: %f for phi2: %f\n", tr1, tr2);
     ratio[0] = exp(tr2 + exp(tr2) - tr1 - exp(tr1));     
     ratio_fnc(ratio, constant, U);
       if(U[0] < ratio[0]){
          phip[0] = phi2[0];
          accept[0] = 1.0;
       }             
       else {
          phip[0] = phi1[0];
          accept[0] = 0.0;
       }     
     }

    free(ratio); free(U); free(det2); 
    return;   

}

///////////////////// THE END ///////////////////////////
