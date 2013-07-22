//** File Name 'equation_xb_ar.c' **//

#include "main_ar.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"


// Joint posterior distribution
void JOINT_ar(int *n, int *T, int *r, int *rT, int *p, int *N, 
     int *cov, int *spdecay,
     double *shape_e, double *shape_eta, double *shape_0,  
     double *prior_a, double *prior_b, double *prior_sig, double *phi, 
     double *tau, double *phis, int *phik, double *nu, double *d, int *constant, 
     double *sig_e, double *sig_eta, 
     double *sig_l0, double *mu_l,  
     double *rho, double *beta, double *X, double *z, double *o, 
     double *phip, double *accept, double *nup, 
     double *sig_ep, double *sig_etap, double *rhop, double *betap, 
     double *mu_lp, double *sig_l0p, double *op, double *w)
{     
     int n1, nn, r1, p1, rn, N1, col; 
     n1 = *n;
     nn = n1*n1;
     r1 = *r;
     p1 = *p;
     rn = r1 *n1;
     N1 = *N;
     col = *constant;
     
    int i;
     
   double *Qeta, *XB, *Sinv, *Qeta2, *det, *det2, *S, *O_l0, *thetap, *tmp, *phi2;
   Qeta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
   XB = (double *) malloc((size_t)((N1)*sizeof(double)));
   Sinv = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   Qeta2 = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   det = (double *) malloc((size_t)((1)*sizeof(double)));   
   det2 = (double *) malloc((size_t)((1)*sizeof(double)));   
   S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   O_l0 = (double *) malloc((size_t)((r1*n1)*sizeof(double)));
   thetap = (double *) malloc((size_t)((p1+1)*sizeof(double)));
   tmp = (double *) malloc((size_t)((1)*sizeof(double)));   
   phi2 = (double *) malloc((size_t)((1)*sizeof(double)));   

   
   covFormat(cov, n, phi, nu, d, sig_eta, S, det, Sinv, Qeta);   
   MProd(beta, constant, p, X, N, XB);
   o0_ar(n, r, T, rT, p, sig_eta, sig_l0, rho, mu_l, Sinv, XB, o, 
   constant, O_l0);

// check nu 
   if(cov[0]==4){
      nu_ar_DIS(cov, Qeta, det, phi, nu, n, r, T, rT, N, d, sig_eta, rho, 
      mu_l, O_l0, XB, o, constant, nup);
   }
   else {
      nup[0] = nu[0];
   }  

// fixed values for phi 
   if(spdecay[0] == 1){
    accept[0] =0.0;
    phip[0] = phi[0];
    covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
   }
// discrete sampling for phi 
   else if(spdecay[0] == 2){

     phi_ar_DIS(cov, Qeta, det, phi, phis, phik, nup,  n, r, T, rT, N, 
     prior_a, prior_b, d, sig_eta, rho, mu_l, O_l0, XB, o, constant, 
     accept, phip);
     covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
//     if(accept[0] == 1.0){    
//     covFormat(cov, n, phip, nu, d, sig_eta, S, det, Sinv, Qeta);   
//     }
   }
// Random-Walk MH-within-Gibbs sampling for phi
   else if(spdecay[0] == 3){
   double *Qeta2, *det2, *tmp, *phi2;
   Qeta2 = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   det2 = (double *) malloc((size_t)((1)*sizeof(double)));   
   tmp = (double *) malloc((size_t)((1)*sizeof(double)));   
   phi2 = (double *) malloc((size_t)((1)*sizeof(double)));   

      if(phi[0] <= 0){
        phi[0] = pow(1,-320);
      }      
      tmp[0] = -log(phi[0]); 
      mvrnormal(constant, tmp, tau, constant, phi2);
      phi2[0]= exp(-phi2[0]);
      //phi2[0]= exp(-phi2[0]-exp(-phi2[0]));      

      covFormat2(cov, n, phi2, nup, d, sig_eta, det2, Qeta2);
//     Rprintf("det2 %f phi2 %f\n", det2[0], phi2[0]);

    phi_ar_MH(Qeta, Qeta2, det, det2, phi, phi2, n, r, T, rT, p, N, 
    prior_a, prior_b, rho, mu_l, O_l0, XB, o, constant, accept, phip);
    if(accept[0] == 1.0){
    covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
    }
     free(Qeta2); free(det2); free(tmp); free(phi2);
   }
   else {
     //abort();
     //exit(9);
   }   

   theta_ar(n, r, T, rT, p, prior_sig, Qeta, O_l0, X, o, constant, 
   thetap); 
   rhop[0] = thetap[0];
   for(i=0; i<p1; i++){
       betap[i] = thetap[i+1];
   }      
   MProd(betap, constant, p, X, N, XB);

   mu_l_ar(n, r, sig_l0, Sinv, O_l0, constant, mu_lp);        
   sig_0l_ar(n, r, shape_0, prior_b, phip, mu_lp, O_l0, Sinv, constant, sig_l0p);          
   sig2_ar(n, r, T, rT, p, shape_e, shape_eta, prior_b, S, Sinv,
   rhop, O_l0, XB, o, z, constant, sig_ep, sig_etap);     
   o0_ar(n, r, T, rT, p, sig_etap, sig_l0p, rhop, mu_lp, Sinv, XB, o, 
   constant, O_l0);
   o_ar(n, r, T, rT, p, sig_ep, sig_etap, rhop, S, Qeta, O_l0, XB, 
   z, o, constant, op);
   w_ar(n, r, T, rT, p, O_l0, X, o, thetap, constant, w);
     
   free(Qeta); free(XB); free(Sinv); free(det); 
   free(S); free(O_l0); free(thetap);

   return;
}

// w
void w_ar(int *n, int *r, int *T, int *rT, int *p, double *O_l0, 
     double *X, double *o, double *thetap, int *constant, double *w) 
{
     int t, l, i, n1, p1, r1, T1, col;
     n1 =*n;
     p1 =*p;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     int *p1n;
     p1n = (int *) malloc((size_t)((col)*sizeof(int)));
     
     *p1n = (p1+1);
     
     double *ot, *ot1, *out, *w1, *X1, *y, *mu;
     ot = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
     ot1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
     w1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
     X1 = (double *) malloc((size_t)((n1*p1)*sizeof(double)));
     y = (double *) malloc((size_t)((n1*(p1+1))*sizeof(double)));
     mu = (double *) malloc((size_t)(((p1+1)*col)*sizeof(double)));     
                    
     for(l=0; l<r1; l++){
     for(t=0; t<T1; t++){
          if(t == 0){
                 for(i=0; i<n1; i++){
                  ot1[i] = O_l0[i+l*n1];
                 }      
          }
          else {       
            extract_alt2(l, t-1, n, rT, T, o, ot1);  // n x 1
          }
          extract_X(t, l, n, r, T, p, X, X1);        // n x p

          for(i=0; i<n1; i++){
              y[i] = ot1[i];     
          }              
          for(i=0; i<n1*p1; i++){
              y[i+n1] = X1[i];
          }     
          MProd(thetap, constant, p1n, y, n, out);   // n x 1
          extract_alt2(l, t, n, rT, T, o, ot);  // n x 1
          for(i=0; i<n1; i++){
              w1[i] = ot[i]-out[i];     
          }              
     put_together1(l, t, n, r, T, w, w1);
     }
     }

     free(p1n);
     free(ot); free(ot1); free(out); free(w1); 
     free(X1); free(y); free(mu);
     
     return;
}     



// Posterior distribution for 2 "sig"s
void sig2_ar(int *n, int *r,  int *T, int *rT, int *p,
     double *shape_e, double *shape_eta, double *prior_b, 
     double *S, double *Sinv, double *rho, double *O_l0, 
     double *XB, double *o, double *z, int *constant, 
     double *sig2ep, double *sig2etap)
{
     int t, l, i, n1, r1, T1, p1, col;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     p1 =*p;
     col =*constant;
     
    double *A, *SinvA, *det, *o2, *z1, *o1, *XB1, *ov, *oz, *out1, *out2, *out3, *tmp;
    A = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
    SinvA = (double *) malloc((size_t)((n1*n1)*sizeof(double)));     
    det = (double *) malloc((size_t)((col)*sizeof(double)));
    o2 = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
    z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    o1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));          
    XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
    ov = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
    oz = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    out1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    out2 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
    out3 = (double *) malloc((size_t)((col)*sizeof(double)));
    tmp = (double *) malloc((size_t)((col)*sizeof(double)));
                        
// Create an indentity matrix
     for(i=0; i <n1*n1; i++){
         if(i == (i+n1*i)){
             A[i] = 1.0;
         }
         else{
             A[i] = 0.0;
         }       
//         Sinv[i] = exp(-(d[i]*phi[0]));
         SinvA[i] = S[i] + A[i];
     }
//     MInv(Sinv, Sinv, n, det);
     MInv(SinvA, SinvA, n, det); 
     
     double u, v, sh1, sh2, b, sige[1], siget[1];
     v = 0.0;      
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(t=0; t < T1; t++) {                                
             if(t == 0){
                   for(i=0; i<n1; i++){
                   o2[i] = O_l0[i+l*n1];
                   }      
             }
             else {       
             extract_alt2(l, t-1, n, rT, T, o, o2);
             }
             extract_alt2(l, t, n, rT, T, z, z1);             
             extract_alt2(l, t, n, rT, T, o, o1);
             extract_alt2(l, t, n, rT, T, XB, XB1);

             for(i=0; i < n1; i++) {
                 ov[i] = o1[i]-rho[0]*o2[i]-XB1[i];
                 out3[0] = z1[i]-o1[i];  
                 tmp[0] = 0.005;
                 mvrnormal(constant, out3, tmp, constant, out3);                              
                 oz[i]=out3[0]; 
             }
             MProd(ov, constant, n, Sinv, n, out1);
             MProd(out1, constant, n, ov, constant, out1);
             u += out1[0];       

             MProd(oz, constant, n, oz, constant, out2);
             v += out2[0];          
         }
     }
     b = *prior_b;
     v = b + 0.5 * v;
     u = b + 0.5 * u;
    
     sh1 = *shape_e;
     sh2 = *shape_eta;

     sige[0] = rigammaa(sh1, v);
     *sig2ep = sige[0];    
     siget[0] = rigammaa(sh2, u);
     *sig2etap = siget[0];         

    
     free(A); free(SinvA); free(det); free(o2); free(z1); 
     free(o1); free(XB1); free(ov); free(oz); free(out1); free(out2);
     free(out3); free(tmp);

     return;
}



// Posterior distribution for "sig_e"
void sig_e_ar(int *n, int *r, int *T, int *rT, double *shape, double *prior_b, 
     double *z, double *o, int *constant, double *sig2e)
{
     int i, t, l, n1, r1, T1, col;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     double *z1, *o1, *zo, *zzoo;
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     o1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zo = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zzoo = (double *) malloc((size_t)((col)*sizeof(double)));
     
     double u, v, b, sige[1];
     u = 0.0;
     v = 0.0;

     b = 0.0;
     for(l=0; l<r1; l++){
        for(t=0; t<T1; t++){
             extract_alt2(l, t, n, rT, T, o, o1);
             extract_alt2(l, t, n, rT, T, z, z1);
             for(i=0; i<n1; i++){
                 zo[i] = z1[i]-o1[i];
             }
             MProd(zo, constant, n, zo, constant, zzoo);
             u += zzoo[0];          
        }
     }      
     b = *prior_b;
     u = b + 0.5 * u;
     v = *shape;
     sige[0] = rigammaa(v, u);
     *sig2e = sige[0];

     free(z1); free(o1); free(zo); free(zzoo);
     return;                  

}     


// Posterior distribution for "sig_eta"
void sig_eta_ar(int *n, int *r,  int *T, int *rT, int *p, double *phi,
     double *shape, double *prior_b, double *Sinv, double *rho, 
     double *O_l0, double *XB, double *o, int *constant, 
     double *sig2eta)     
{
     double *ov, *o1, *o2, *out, u, b, sh, sig[1];
     int row, col, l, i, j, r1, T1, p1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     p1= *p;
     
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     o2 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     out = (double *) malloc((size_t)((row*col)*sizeof(double)));
          
     double *XB1;
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1; i++) {                                
             if(i == 0){
                   for(j=0; j<row; j++){
                   o2[j] = O_l0[j+l*row];
                   }      
             }
             else {       
             extract_alt2(l, i-1, n, rT, T, o, o2);
             }

             extract_alt2(l, i, n, rT, T, o, o1);
             extract_alt2(l, i, n, rT, T, XB, XB1);
            
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-rho[0]*o2[j]-XB1[j];
             }
             MProd(ov, constant, n, Sinv, n, out);
             MProd(out, constant, n, ov, constant, out);
             u += out[0];       
//             u += xTay2(ov, Sinv, ov, row);
         }
     }
     b = *prior_b;
     u = b + 0.5 * u;
     sh = *shape;
     sig[0] = rigammaa(sh, u);
     *sig2eta = sig[0]; 
    
     free(ov); free(o1); free(o2); free(XB1); free(out);
     return;
}



// Posterior distribution for "beta"
void beta_ar(int *n, int *r, int *T, int *rT, int *p, 
     double *prior_sig, double *Q_eta, double *rho, 
     double *O_l0, double *X, double *o, int *constant, double *beta) 
{
     double *ot, *ot1, *ot2, *ox, u, v;
     int row, col, l, i, j, k, r1, T1, p1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     p1 = *p;

     ot = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ot1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ot2 = (double *) malloc((size_t)((row*col)*sizeof(double)));     
     ox = (double *) malloc((size_t)((row*col)*sizeof(double)));

        
     double *out, *X1;
     out = (double *) malloc((size_t)((row*col)*sizeof(double)));
     X1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     double *mu, *s2, *bt;
     mu = (double *) malloc((size_t)((col*col)*sizeof(double)));
     s2 = (double *) malloc((size_t)((col*col)*sizeof(double)));     
     bt = (double *) malloc((size_t)((col*col)*sizeof(double)));     

     for(k=0; k < p1; k++){
     u = 0.0;               
     v = 0.0;                        
            for(l=0; l < r1; l++) {
            for(i=0; i < T1; i++) {
            extract_X3(l, i, k, n, rT, T, p, X, X1);
            MProd(X1, constant, n, Q_eta, n, out);
            MProd(out, constant, n, X1, constant, out);
//           xTay(X1, Q_eta, X1, n, out);
            u += out[0];
            }
            }
            
            
            for(l=0; l < r1; l++) {
            for(i=0; i < T1; i++) {

              if(i == 0){
                   for(j=0; j<row; j++){
                   ot2[j] = O_l0[j+l*row];
                   }      
              }
              else {       
              extract_alt2(l, i-1, n, rT, T, o, ot2);
              }

              extract_alt2(l, i, n, rT, T, o, ot);
              extract_X3(l, i, k, n, rT, T, p, X, X1);
            
              for(j=0; j< row; j++) {
                ox[j] = ot[j] - rho[0]*ot2[j];
              }    
              MProd(ox, constant, n, Q_eta, n, out);
              MProd(out, constant, n, X1, constant, out);
//              xTay(X1, Q_eta, ox, n, out);
              v += out[0];
            }
            }
        u = 1.0/(u + (1.0/prior_sig[0])); // OK      
        v = (u * v);
     *mu = v;
     *s2 = u;
      mvrnormal(constant, mu, s2, constant, bt);
      beta[k] = bt[0];
     } 

     free(ot); free(ot1);  free(ot2); free(ox); free(out); free(X1); 
     free(mu); free(s2); free(bt); 
     return;
}     



// Posterior distribution for "theta"
void theta_ar(int *n, int *r, int *T, int *rT, int *p, double *prior_sig, 
     double *Q_eta, double *O_l0, double *X, double *o, int *constant, 
     double *thetap) 
{
     int t, l, i, n1, p1, r1, T1, col;
     n1 =*n;
     p1 =*p;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     int *p1n;
     p1n = (int *) malloc((size_t)((col)*sizeof(int)));
     
     *p1n = (p1+1);
     
     double *del, *chi, *ot1, *X1, *y, *ty, *out, *tyQy, *ot2, *tyQo, *det, *mu, *I;
     del = (double *) malloc((size_t)(((p1+1)*(p1+1))*sizeof(double)));
     chi = (double *) malloc((size_t)(((p1+1)*col)*sizeof(double)));     
     ot1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     X1 = (double *) malloc((size_t)((n1*p1)*sizeof(double)));
     y = (double *) malloc((size_t)((n1*(p1+1))*sizeof(double)));
     ty = (double *) malloc((size_t)((n1*(p1+1))*sizeof(double)));
     out = (double *) malloc((size_t)((n1*(p1+1))*sizeof(double)));
     tyQy = (double *) malloc((size_t)(((p1+1)*(p1+1))*sizeof(double)));
     ot2 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     tyQo = (double *) malloc((size_t)(((p1+1)*col)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     mu = (double *) malloc((size_t)(((p1+1)*col)*sizeof(double)));     
     I = (double *) malloc((size_t)(((p1+1)*(p1+1))*sizeof(double)));                         
     
     for(i=0; i<(p1+1)*(p1+1); i++){
           del[i] = 0.0;
     }   
     for(i=0; i<(p1+1); i++){
           chi[i] = 0.0;
     }   
     
     for(l=0; l<r1; l++){
     for(t=0; t<T1; t++){
          if(t == 0){
                 for(i=0; i<n1; i++){
                  ot1[i] = O_l0[i+l*n1];
                 }      
          }
          else {       
            extract_alt2(l, t-1, n, rT, T, o, ot1);  // n x 1
          }
          extract_X(t, l, n, r, T, p, X, X1);        // n x p

          for(i=0; i<n1; i++){
              y[i] = ot1[i];     
          }              
          for(i=0; i<n1*p1; i++){
              y[i+n1] = X1[i];
          }     
          //for(i=0; i<n1*(p1+1); i++){
          //    y[i] = y[i];
          //}     
          MTranspose(y, p1n, n, ty);         // (p+1) x n
          MProd(y, p1n, n, Q_eta, n, out);   // n x (p+1)
          MProd(out, p1n, n, ty, p1n, tyQy); // (p+1) x (p+1)
          MAdd(del, p1n, p1n, tyQy, del);  // (p+1) x (p+1)
          extract_alt2(l, t, n, rT, T, o, ot2);   // n x 1
          MProd(ot2, constant, n, Q_eta, n, out); // n x 1
          MProd(out, constant, n, ty, p1n, tyQo);  // (p+1) x 1
          MAdd(chi, p1n, constant, tyQo, chi);  // (p+1) x 1

     }
     }

     IdentityM(p1n, I);
     for(i=0; i<(p1+1)*(p1+1); i++){
     del[i] = del[i] + I[i]*(1.0/prior_sig[0]);
     }
     //for(i=0; i<(p1+1); i++){
     //chi[i] = chi[i] + prior_mu[0]/prior_sig[0];
     //}

     free(I);
     MInv(del, del, p1n, det);
     MProd(chi, constant, p1n, del, p1n, mu);  // (p+1) x 1      
     mvrnormal(constant, mu, del, p1n, thetap);
//     for(i=0; i<(p1+1); i++){
//      thetap[i] = mu[i];
//     }


     free(p1n);
     free(del); free(chi); free(ot1); free(X1); free(y); free(ty);
     free(out); free(tyQy); free(ot2); free(tyQo); free(det); free(mu);
     
     return;
}     



// Posterior distribution for "rho"
void rho_ar(int *n, int *r, int *T, int *rT, int *p, 
     double *prior_sig, double *Q_eta, double *O_l0, 
     double *XB, double *o, int *constant, double *rho)     
{
     double *ot, *ot1, u, v;
     int row, col, l, i, j, r1, T1, nn, p1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     nn = row*row;
     p1= *p;

     ot = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ot1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     double *term1, *out;
     term1 = (double *) malloc((size_t)((row*col)*sizeof(double)));     
     out = (double *) malloc((size_t)((row*col)*sizeof(double)));


     double *XB1;
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     
     v = 0.0;   
      for(l=0; l < r1; l++) {
      for(i=0; i < T1; i++) {
             if(i == 0){
                   for(j=0; j<row; j++){
                   ot1[j] = O_l0[j+l*row];
                   }      
             }
             else {       
             extract_alt2(l, i-1, n, rT, T, o, ot1);
             }
             MProd(ot1, constant, n, Q_eta, n, out);
             MProd(out, constant, n, ot1, constant, out);             
//             xTay(ot1, Q_eta, ot1, n, out);
             v += out[0];
      }
      }
     u = 0.0;         
      for(l=0; l < r1; l++) {
      for(i=0; i < T1; i++) {
               
             if(i == 0){
                   for(j=0; j<row; j++){
                   ot1[j] = O_l0[j+l*row];
                   }      
             }
             else {       
             extract_alt2(l, i-1, n, rT, T, o, ot1);
             }

             extract_alt2(l, i, n, rT, T, o, ot);
             extract_alt2(l, i, n, rT, T, XB, XB1);
             
             for(j=0; j < row; j++) {
                   term1[j] = (ot[j] - XB1[j]);      
             }
             MProd(term1, constant, n, Q_eta, n, out);
             MProd(out, constant, n, ot1, constant, out);              
//             xTay(ot1, Q_eta, term1, n, out1);
             u += out[0];
         }
     }
     v = 1.0/(v + (1.0/prior_sig[0])); // OK
     u = u*v;

     double *mu, *s2, *out2;
     mu = (double *) malloc((size_t)((col*col)*sizeof(double)));
     s2 = (double *) malloc((size_t)((col*col)*sizeof(double)));     
     out2 = (double *) malloc((size_t)((col*col)*sizeof(double)));     

     *mu = u;
     *s2 = v;

     mvrnormal(constant, mu, s2, constant, out2);
     rho[0] = out2[0];
          
     free(ot); free(ot1); free(out); free(mu); free(s2);
     free(term1); free(XB1); free(out2); 
     return;
}


// Posterior distribution for "sig_0l"
void sig_0l_ar(int *n, int *r, double *shape, double *prior_b, double *phi,
     double *mu_l, double *O_l0, double *Sinv, int *constant, double *sig0)          
{
     double sh, rt, *Ql, *I, *out, b, u;
     int i, l, r1, row, col, nn;
     
     sh = *shape;
     r1 = *r;
     row = *n;
     nn = row * row;
     col = *constant;
     
     Ql = (double *) malloc((size_t)((nn*col)*sizeof(double)));
     out = (double *) malloc((size_t)((row*col)*sizeof(double)));
     I =(double *) malloc((size_t)((row*col)*sizeof(double)));
//     double n1 = (double) row;     

     for(l=0; l < r1; l++){
	 rt = 0.0;
	 b = 0.0;
	 u = 0.0;
     
     for(i=0; i < row; i++) {
           I[i] = O_l0[i+l*row] - mu_l[l];
      }
     MProd(I, constant, n, Sinv, n, out);
     MProd(out, constant, n, I, constant, out);
     u = out[0];
     b = *prior_b;
     rt = b + 0.5 * u;            
     sig0[l] = rigammaa(sh, rt);
     }

     free(Ql); free(out); free(I);
     return;
}


// Posterior distribution for "mu_l"
void mu_l_ar (int *n, int *r, double *sig_l0, double *Sinv, double *O_l0, 
     int *constant, double *mu_l)         
{
     int i, l, r1, n1, nn, col;
     r1 =*r;
     n1 =*n;
     nn =n1*n1;
     col =*constant;
          
     double *I, *S, *s2, *m, *mu, *out;
     I = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     s2 = (double *) malloc((size_t)((col)*sizeof(double)));
     m = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((n1*col)*sizeof(double)));
                              
     for(i=0; i <n1; i++){
           I[i] = 1.0;
     }

     for(l=0; l<r1; l++){
         for(i=0; i<nn; i++){
            S[i] = Sinv[i]/sig_l0[l];
         }
         MProd(I, constant, n, S, n, out);
         MProd(out, constant, n, I, constant, s2);              
         s2[0] = 1.0/s2[0];
         for(i=0; i<n1; i++){         
            m[i] = O_l0[i+l*n1];         
         }
         mean(n, m, mu);
         mvrnormal(constant, mu, s2, constant, out);
         mu_l[l] = out[0];
     }

     free(I); free(S); free(s2); free(m); free(mu); free(out);
     
     return;
}



// Posterior samples for Z_lt obtained using MCMC samples of O_lt
void Z_fitfnc(int *its, int *N, double *sig_ep, double *o_p, 
                int *constant, double *z_p)
{
     unsigned iseed = 44;
     srand(iseed); 
     
     int i, j, N1, its1, col;
     N1 = *N;
     its1 = *its;
     col = *constant;
     
     double *mu, *sig, *out;
     mu = (double *) malloc((size_t)((col)*sizeof(double)));
     sig = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));
     
     for(i=0; i<its1; i++){
         for(j=0; j<N1; j++){
             mu[0]=o_p[j+i*N1];
             sig[0]=sig_ep[i];     
             mvrnormal(constant, mu, sig, constant, out);
             z_p[j+i*N1]=out[0];
         }
     }
     free(mu); free(sig); free(out);            
     return;
}



// sampling o_lt

// Posterior distribution for "o", the true ozone values
// you need to put the sqrt of values in "z" and "o"
void o0_ar(int *n, int *r, int *T, int *rT, int *p, double *sig_eta, 
     double *sig_l0, double *rho, double *mu_l, double *Sinv, 
     double *XB, double *o, int *constant, double *o0post)     
{     
     int i, l, n1, r1, nn, col; 
     n1 = *n;
     r1 = *r;
     nn = n1*n1;
     col = *constant;
     
     double *del, *det, *o1, *XB1, *m, *out, *I;

     del = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     o1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     m = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     out = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     I = (double *) malloc((size_t)((n1*n1)*sizeof(double)));

     IdentityM(n, I);
                                        
     for(l=0; l<r1; l++){
       for(i=0; i<nn; i++){
           del[i] = Sinv[i]*(rho[0]*rho[0]/sig_eta[0]+1.0/sig_l0[l]);
           del[i] = del[i] + I[i]*(1/10000);
       }
       MInv(del, del, n, det);
       
           extract_alt2(l, 0, n, rT, T, o, o1);
           extract_alt2(l, 0, n, rT, T, XB, XB1);

       for(i=0; i<n1; i++){
       m[i]=(rho[0]*(o1[i]-XB1[i])*(1.0/sig_eta[0]) + 
           mu_l[l]/sig_l0[l])/(rho[0]*rho[0]/sig_eta[0]+1.0/sig_l0[l]);
       }
//       mvrnormal(constant, m, del, n, out);

       for(i=0; i<n1; i++){
         o0post[i+l*n1]=m[i]; //out[i];
       }
     }
     
     free(del); free(det); free(o1); free(XB1); 
     free(m); free(out); free(I);
     return;
}     
     

void o_ar(int *n, int *r, int *T, int *rT, int *p, double *sig_e, 
     double *sig_eta, double *rho, double *S, double *Q_eta, 
     double *O_l0, double *XB, double *z, double *o, 
     int *constant, double *opost)     
{
     int i, l, t, r1, nn, row, T1, col, p1; 
     r1 = *r;     row = *n;     T1 = *T;     nn = row * row;  col = *constant;
     p1 = *p;
     
     double *term1, *term2, *o_1, *de_tT, *d_tT, *de_T, *delT; 
     double *det1, *det2, *mean1, *XB1, *XB2;
     
     o_1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     term1 = (double *) malloc((size_t)((nn)*sizeof(double)));
     term2 = (double *) malloc((size_t)((row)*sizeof(double)));
     
     de_tT = (double *) malloc((size_t)((nn*r1*(T1-2))*sizeof(double)));
     d_tT = (double *) malloc((size_t)((nn*r1*(T1-2))*sizeof(double)));
     de_T = (double *) malloc((size_t)((nn*r1)*sizeof(double)));
     delT = (double *) malloc((size_t)((nn*r1)*sizeof(double)));
     det1 = (double *) malloc((size_t)((col)*sizeof(double)));
     det2 = (double *) malloc((size_t)((col)*sizeof(double)));     
     mean1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB2 = (double *) malloc((size_t)((row*col)*sizeof(double)));

    
     double *zT, *zt, *o_2, *I, *o1, *II;
     zT = (double *) malloc((size_t)((row*col)*sizeof(double)));          
     zt = (double *) malloc((size_t)((row*col)*sizeof(double)));               
     o_2 = (double *) malloc((size_t)((row*col)*sizeof(double)));          
     I = (double *) malloc((size_t)((row*col)*sizeof(double)));          
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     II = (double *) malloc((size_t)((nn)*sizeof(double)));


     IdentityM(n, II);
                    
// for 1 < t < T, the delta part
         for(i=0; i < nn; i++) {
            de_tT[i] = ((1.0/sig_e[0]) + Q_eta[i] + rho[0]*rho[0]*Q_eta[i]);
            de_tT[i] = de_tT[i] + II[i]*(1/10000);
         }    
         MInv(de_tT, d_tT, n, det1); 

// for t = T, the delta part
         for(i=0; i < nn; i++) {
             de_T[i] = ((1.0/sig_e[0]) + Q_eta[i]);       
             de_T[i] = de_T[i] + II[i]*(1/10000);
         }    
         MInv(de_T, delT, n, det2);

         free(II);
         
// term 1
         for(i=0; i < nn; i++) {
//             term1[i] = (sig_eta[0]/sig_e[0])* exp(-(d[i]*phi[0]));       
             term1[i] = (sig_eta[0]/sig_e[0])* S[i];       
//             opost[i]=term1[i];  // OK
         }    
         for(i=0; i < row; i++) {
              I[i]= 1.0;
         }    
         MProd(I, constant, n, term1, n, term2);
//         MProd(I, constant, n, term1, n, opost);    // OK     

        free(I);   
// LOOP
     for(l=0; l < r1; l++) {

// for 1 <= t < T, 
         for(t=0; t < (T1-1); t++) {          
            if(t == 0){
                   for(i=0; i<row; i++){
                   o_1[i] = O_l0[i+l*row];
                   }      
             }
             else {       
             extract_alt2(l, t-1, n, rT, T, o, o_1);
             }
              
              extract_alt2(l, (t+1), n, rT, T, o, o_2);
              extract_alt2(l, t, n, rT, T, XB, XB1);
              extract_alt2(l, (t+1), n, rT, T, XB, XB2);

          extract_alt2(l, t, n, rT, T, z, zT);
          MProd(zT, constant, n, term1, n, zt);
                            
           for(i=0; i < row; i++) {
              mean1[i] = (XB1[i]+rho[0]*o_1[i]+rho[0]*o_2[i]-
                       rho[0]*XB2[i]+zt[i])/(1+rho[0]*rho[0]+term2[i]);     
           }                  
              mvrnormal(constant, mean1, d_tT, n, o_1);     // random generator
         put_together1(l, t, n, r, T, opost, o_1);
//       put_together1(l, t, n, r, T, opost, mean1); // to check
         }

// for t = T, 
         t = (T1-1);
          extract_alt2(l, t, n, rT, T, z, zT);
          MProd(zT, constant, n, term1, n, zt);

              extract_alt2(l, (t-1), n, rT, T, o, o_1);
              extract_alt2(l, (t-1), n, rT, T, XB, XB1);

          for(i=0; i < row; i++) {
            mean1[i] = (rho[0]*o_1[i] + XB1[i] + zt[i])/(1+term2[i]);     

           }                  
              mvrnormal(constant, mean1, delT, n, o_1);     // random generator
         put_together1(l, t, n, r, T, opost, o_1);
//       put_together1(l, t, n, r, T, opost, mean1); // check

         } // End of loop year
         
       

       free(o_1); free(term1); free(de_tT); free(d_tT); 
       free(de_T); free(delT); free(det1); free(det2); free(mean1); free(zT); 
       free(zt); free(o_2); free(I); free(o1); free(XB1); free(XB2);

       return;
} 




// Phi sample random walk
void phi_ar_MH(double *Sinv1, double *Sinv2, double *det1, double *det2,
     double *phi1, double *phi2,
     int *n, int *r, int *T, int *rT, int *p, int *N, 
     double *prior_a, double *prior_b, double *rho,  
     double *mu_l, double *O_l0, double *XB, double *o, 
     int *constant, double *accept, double *phip)
{
     
     double *ov, *o1, *o2, u, v, w, x, a, b;
     int row, col, l, i, j, r1, T1, p1, N1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     p1= *p;
     N1 = row*r1*T1;
     rT1 = *rT;
     
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     o2 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     
     double *XB1; 
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     double *ratio, *U; 
     ratio = (double *) malloc((size_t)((col)*sizeof(double)));             
     U = (double *) malloc((size_t)((col)*sizeof(double)));         

     u = 0.0;
     v = 0.0;
     w = 0.0;
     x = 0.0;     

     for(l=0; l < r1; l++) {
         for(i=0; i < T1; i++) {                                
             if(i == 0){
                   for(j=0; j<row; j++){
                   o2[j] = O_l0[j+l*row];
                   }      
             }
             else {       
             extract_alt2(l, i-1, n, rT, T, o, o2);
             }

             extract_alt2(l, i, n, rT, T, o, o1);
             extract_alt2(l, i, n, rT, T, XB, XB1);
            
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-rho[0]*o2[j]-XB1[j];
             }       
             u += xTay2(ov, Sinv1, ov, row);
             v += xTay2(ov, Sinv2, ov, row);
         }
//         for(j=0; j < row; j++) {
//             ov[j]=O_l0[j+l*row] - mu_l[l];
//         }       
//         w += xTay2(ov, Sinv1, ov, row);
//         x += xTay2(ov, Sinv2, ov, row);
     }
     
     a = *prior_a;
     b = *prior_b;
     u =  0.5 * u;
     v =  0.5 * v;
//     u = b + 0.5 * u;
//     v = b + 0.5 * v;
//     w = b + 0.5 * w;
//     x = b + 0.5 * x;

     free(o1); free(o2); free(ov); free(XB1); 

     double tr1, tr2;
     tr1 = 0.0;
     tr2 = 0.0;

//     double gN, ga;
//     gN = sp_gamma(N1/2+a);
//     ga = sp_gamma(a);

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


     if(phi2[0] < 0.0010){
          phip[0] = phi1[0];
          accept[0] = 0.0;
     }
     else if(phi2[0] > 0.9999){
          phip[0] = phi1[0];
          accept[0] = 0.0;
     }
     else{    
// with Gamma prior    
     tr1 = (a-1.0)*log(phi1[0])-b*phi1[0]-0.5*rT1*log(det1[0])-u; 
     tr2 = (a-1.0)*log(phi2[0])-b*phi2[0]-0.5*rT1*log(det2[0])-v; 
//     Rprintf("for u: %f for v: %f\n", u, v);
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
     
     free(ratio); free(U);
     
     return;
}     



// Discrete Phi sample
void phi_ar_DIS(int *cov, double *Qeta1, double *det1, double *phi1, 
     double *phis, int *phik, double *nu, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *d, double *sig2eta, 
     double *rho, double *mu_l, double *O_l0, double *XB, double *o, 
     int *constant, double *accept, double *phip)
{
     int row, col, i, r1, T1, N1, rT1, pk;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     N1 = row*r1*T1;
     rT1 = *rT;
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
        covFormat2(cov, n, phitmp, nu, d, sig2eta, det, Qeta);
        phidens_ar(phitmp, Qeta, det, n, r, T, rT, N, prior_a, prior_b, XB, 
        rho, O_l0, o, constant, out);
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
     phidens_ar(phi1, Qeta1, det1, n, r, T, rT, N, prior_a, prior_b, XB, 
     rho, O_l0, o, constant, tr1);

     ratio[0] = exp(tr2[0] - tr1[0]);
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


void phidens_ar(double *phi, double *Qeta, double *det, int *n, int *r, 
     int *T, int *rT, int *N, double *prior_a, double *prior_b, double *XB, 
     double *rho, double *O_l0, double *o, int *constant, double *out)
{

     int row, col, l, i, j, r1, T1, N1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     N1 = row*r1*T1;
     rT1 = *rT;
     
     double *ov, *o1, *o2, *XB1;     
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     o2 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     double u, a, b;
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1; i++) {                                
             if(i == 0){
                   for(j=0; j<row; j++){
                   o2[j] = O_l0[j+l*row];
                   }      
             }
             else {       
             extract_alt2(l, i-1, n, rT, T, o, o2);
             }
             extract_alt2(l, i, n, rT, T, o, o1);
             extract_alt2(l, i, n, rT, T, XB, XB1);
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-rho[0]*o2[j]-XB1[j];
             }       
             u += xTay2(ov, Qeta, ov, row);
         }
     }
     free(ov); free(o1); free(o2); free(XB1);

     a = *prior_a;
     b = *prior_b;
     u =  0.5 * u;
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


// Discrete nu sample
void nu_ar_DIS(int *cov, double *Qeta1, double *det1, double *phi, double *nu1, 
     int *n, int *r, int *T, int *rT, int *N, double *d, double *sig2eta, 
     double *rho, double *mu_l, double *O_l0, double *XB, double *o, 
     int *constant, double *nup)
{
     int row, col, i, r1, T1, N1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     N1 = row*r1*T1;
     rT1 = *rT;
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
//     nus[20]=1.25; nus[21]=1.50; nus[22]=1.75; nus[23]=2.00; nus[24]=2.5; 
//     nus[25]=3.00; nus[26]=4.00; nus[27]=5.00; nus[28]=10.00; nus[29]=20.00;      

/*
     nuk=20;
     double *nus;
     nus = (double *) malloc((size_t)((nuk)*sizeof(double)));             
     nus[0]=0.05; nus[1]=0.10; nus[2]=0.15; nus[3]=0.20; nus[4]=0.25;  
     nus[5]=0.30; nus[6]=0.35; nus[7]=0.40; nus[8]=0.45; nus[9]=0.50;  
     nus[10]=0.55; nus[11]=0.60; nus[12]=0.65; nus[13]=0.70; nus[14]=0.75; 
     nus[15]=0.80; nus[16]=0.85; nus[17]=0.90; nus[18]=0.95; nus[19]=1.0; 
*/

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
        covFormat2(cov, n, phi, nutmp, d, sig2eta, det, Qeta);
        nudens_ar(Qeta, det, n, r, T, rT, N, XB, rho, O_l0, o, constant, out);
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
     nudens_ar(Qeta1, det1, n, r, T, rT, N, XB, rho, O_l0, o, constant, tr1);
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


void nudens_ar(double *Qeta, double *det, int *n, int *r, int *T, int *rT, 
     int *N, double *XB, double *rho, double *O_l0, double *o, int *constant, 
     double *out)
{
     int row, col, l, i, j, r1, T1, N1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     N1 = row*r1*T1;
     rT1 = *rT;
     
     double *ov, *o1, *o2, *XB1;     
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     o2 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     double u;
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1; i++) {                                
             if(i == 0){
                   for(j=0; j<row; j++){
                   o2[j] = O_l0[j+l*row];
                   }      
             }
             else {       
             extract_alt2(l, i-1, n, rT, T, o, o2);
             }
             extract_alt2(l, i, n, rT, T, o, o1);
             extract_alt2(l, i, n, rT, T, XB, XB1);
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-rho[0]*o2[j]-XB1[j];
             }       
             u += xTay2(ov, Qeta, ov, row);
         }
     }
     free(ov); free(o1); free(o2); free(XB1);

     u =  0.5 * u;
     if(det[0] <= 0){
        det[0] = pow(1,-320);
     }
     double tr;
     tr = 0.0;        
     tr = -0.5*rT1*log(det[0])-u; 
     out[0] = tr;

     return;
}




///////////////////////////// THE END /////////////////////////////////////////
