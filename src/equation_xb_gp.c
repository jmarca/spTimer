//** File Name 'equation_xb_gp.c' **//
//** Gibbs codes are in this file **//

#include "main_gp.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"


// Joint posterior distribution
void JOINT_gp(int *n, int *T, int *r, int *rT, int *p, int *N, 
     int *cov, int *spdecay, double *shape_e, double *shape_eta,
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik, double *nu,
     double *d, double *sig_e, double *sig_eta, 
     double *beta, double *X, double *z, double *o, int *constant,
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *betap, double *op)
{     
     int n1, nn, r1, p1, rn, N1, col; 
     n1 = *n;
     nn = n1*n1;
     r1 = *r;
     p1 = *p;
     rn = r1 *n1;
     N1 = *N;
     col = *constant;
     
   double *Qeta, *XB, *Sinv, *det, *S;
   Qeta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
   XB = (double *) malloc((size_t)((N1)*sizeof(double)));
   Sinv = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   det = (double *) malloc((size_t)((1)*sizeof(double)));   
   S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
// Rprintf("   phi: %4.4f, cov: %i\n", phi[0], cov[0]);         
   covFormat(cov, n, phi, nu, d, sig_eta, S, det, Sinv, Qeta);   
   MProd(beta, constant, p, X, N, XB);
// check nu
   if(cov[0]==4){
      nu_gp_DIS(cov, Qeta, det, phi, nu, n, r, T, rT, N, d, sig_eta, XB, o, 
      constant, nup);
//      Rprintf("   nu: %4.4f, nup: %4.4f \n", nu[0], nup[0]);
//      covFormat(cov, n, phi, nup, d, sig_eta, S, det, Sinv, Qeta);   
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
     phi_gp_DIS(cov, Qeta, det, phi, phis, phik, nup, n, r, T, rT, N, 
     prior_a, prior_b, d, sig_eta, XB, o, constant, accept, phip);
     covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
//       if(accept[0] == 1.0){    
//        covFormat(cov, n, phip, nu, d, sig_eta, S, det, Sinv, Qeta);   
//       }
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
// Rprintf("   phi: %4.4f, tmp: %4.4f, cov: %i\n", phi[0], tmp[0], cov[0]);      
      mvrnormal(constant, tmp, tau, constant, phi2);
      phi2[0]= exp(-phi2[0]);      
// Rprintf("   phi: %4.4f, tmp: %4.4f, cov: %i\n", phi[0], tmp[0], cov[0]);      
      covFormat(cov, n, phi2, nup, d, sig_eta, S, det2, Sinv, Qeta2);   
// Rprintf("   phi: %4.4f, phi2: %4.4f, cov: %i\n", phi[0], phi2[0], cov[0]);
//     randow-walk M  
     phi_gp_MH(Qeta, Qeta2, det, det2, phi, phi2, n, r, T, rT, N, 
     prior_a, prior_b, XB, o, constant, accept, phip);
// Rprintf("   phi: %4.4f, phi2: %4.4f, cov: %i\n", phi[0], phi2[0], cov[0]);
     if(accept[0] == 1.0){    
         covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
     }
     free(Qeta2); free(det2); free(tmp); free(phi2);
   }
   else {
     //;
//     exit(9);
   }   
   beta_gp(n, r, T, rT, p, prior_mubeta, prior_sigbeta, Qeta, X, o,
   constant, betap);
   MProd(betap, constant, p, X, N, XB);
   sig_e_gp(n, r, T, rT, N, shape_e, prior_b, XB, o, z, sig_e, constant, sig_ep);
   sig_eta_gp(n, r, T, rT, shape_eta, prior_b, Sinv, XB, o, constant, sig_etap);
   o_gp(n, r, T, rT, p, prior_omu, prior_osig, sig_e, sig_etap, S, Qeta, 
   XB, z, constant, op);     
           
   free(Qeta); free(XB); free(Sinv); free(det); free(S); 
   return;
}


// Posterior distribution for "sig_e"
void sig_e_gp(int *n, int *r, int *T, int *rT, int *N, double *shape, 
     double *prior_b, double *XB, double *o, double *z, double *sig_e,
     int *constant, double *sig2e)
{
     int i, t, l, n1, r1, T1, col;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     double *z1, *o1, *zo, *zzoo, *XB1, *tmp;
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     o1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zo = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zzoo = (double *) malloc((size_t)((col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
     tmp = (double *) malloc((size_t)((col)*sizeof(double)));     
          
     double u, v, b, sige[1];
     u = 0.0;
     v = 0.0;
     b = 0.0;
     for(l=0; l<r1; l++){
        for(t=0; t<T1; t++){
             extract_alt2(l, t, n, rT, T, o, o1);
             extract_alt2(l, t, n, rT, T, z, z1);
             extract_alt2(l, t, n, rT, T, XB, XB1);
             for(i=0; i<n1; i++){
//                 zzoo[0] = o1[i]-XB1[i];
                 zzoo[0] = z1[i]-o1[i];
                 tmp[0] = 0.005;
                 mvrnormal(constant, zzoo, tmp, constant, zzoo);                              
//                 zo[i] = z1[i]-(o1[i]+zzoo[0]);
                 zo[i] = zzoo[0];
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

     free(z1); free(o1); free(zo); free(zzoo); free(XB1);free(tmp);
     return;                  
}     



// Posterior distribution for "sig_eta"
void sig_eta_gp(int *n, int *r,  int *T, int *rT, double *shape, 
     double *prior_b, double *Sinv, double *XB, double *o, int *constant, 
     double *sig2eta)     
{
     double *ov, *o1, *out, u, b, sh, sig[1];
     int row, col, l, i, j, r1, T1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     out = (double *) malloc((size_t)((row*col)*sizeof(double)));
          
     double *XB1;
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1; i++) {                                
             extract_alt2(l, i, n, rT, T, o, o1);
             extract_alt2(l, i, n, rT, T, XB, XB1);
            
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-XB1[j];
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
    
     free(ov); free(o1); free(XB1); free(out);
     return;
}


// Posterior distribution for "theta"
void beta_gp(int *n, int *r, int *T, int *rT, int *p, double *prior_mu,
     double *prior_sig, double *Qeta, double *X, double *o, int *constant, 
     double *betap) 
{
     int t, l, i, n1, p1, r1, T1, col;
     n1 =*n;
     p1 =*p;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     
     double *del, *chi, *ot1, *X1, *tX1, *out, *tX1QX1, *tX1Qo, *det, *mu, *I;
     del = (double *) malloc((size_t)((p1*p1)*sizeof(double)));
     chi = (double *) malloc((size_t)((p1*col)*sizeof(double)));     
     ot1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     X1 = (double *) malloc((size_t)((n1*p1)*sizeof(double)));
     tX1 = (double *) malloc((size_t)((n1*p1)*sizeof(double)));
     out = (double *) malloc((size_t)((n1*p1)*sizeof(double)));
     tX1QX1 = (double *) malloc((size_t)((p1*p1)*sizeof(double)));
     tX1Qo = (double *) malloc((size_t)((p1*col)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     mu = (double *) malloc((size_t)((p1*col)*sizeof(double)));     
     I = (double *) malloc((size_t)((p1*p1)*sizeof(double)));                         
     
     for(i=0; i<p1*p1; i++){
           del[i] = 0.0;
     }   
     for(i=0; i<p1; i++){
           chi[i] = 0.0;
     }   
 
     for(l=0; l<r1; l++){
     for(t=0; t<T1; t++){
          extract_X(t, l, n, r, T, p, X, X1);    // n x p
          MTranspose(X1, p, n, tX1);         // p x n
          MProd(X1, p, n, Qeta, n, out);   // n x p
          MProd(out, p, n, tX1, p, tX1QX1); // pxp
          MAdd(del, p, p, tX1QX1, del);  // pxp

          extract_alt2(l, t, n, rT, T, o, ot1);  // n x 1
          MProd(ot1, constant, n, Qeta, n, out); // n x 1
          MProd(out, constant, n, tX1, p, tX1Qo);  // p x 1
          MAdd(chi, p, constant, tX1Qo, chi);  // p x 1

     }
     }

     IdentityM(p, I);
     for(i=0; i<p1*p1; i++){
     del[i] = del[i] + I[i]*(1.0/prior_sig[0]);     
     }
     free(I); 
     for(i=0; i<p1; i++){
     chi[i] = chi[i] + prior_mu[0]/prior_sig[0];
     }
              
     MInv(del, del, p, det);
     MProd(chi, constant, p, del, p, mu);  // p x 1      
     mvrnormal(constant, mu, del, p, betap);


     free(del); free(chi); free(ot1); free(X1); free(tX1);
     free(out); free(tX1QX1); free(tX1Qo); free(det); free(mu);
     
     return;
}     


// conditional posterior for o_lt
void o_gp(int *n, int *r, int *T, int *rT, int *p, double *prior_omu,
     double *prior_osig, double *sig_e, double *sig_eta, double *S, 
     double *Qeta, double *XB, double *z, int *constant, double *opost)     
{
     int i, l, t, r1, nn, row, T1, col, p1; 
     r1 = *r;     row = *n;     T1 = *T;     nn = row * row;  col = *constant;
     p1 = *p;
     
     double *o_1, *de_tT, *det1, *chi_tT, *mean1, *XB1, *QXB1, *zT;
     
     o_1 = (double *) malloc((size_t)((row)*sizeof(double)));
     de_tT = (double *) malloc((size_t)((nn)*sizeof(double)));
     det1 = (double *) malloc((size_t)((col)*sizeof(double)));
     chi_tT = (double *) malloc((size_t)((row)*sizeof(double)));
     mean1 = (double *) malloc((size_t)((row)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row)*sizeof(double)));
     QXB1 = (double *) malloc((size_t)((row)*sizeof(double)));    
     zT = (double *) malloc((size_t)((row)*sizeof(double)));          
               

// for 1 <= t <= T, the delta part
         for(i=0; i < nn; i++) {
            de_tT[i] = (1.0/sig_e[0]) + Qeta[i] + 1.0/prior_osig[0];
         }    
         MInv(de_tT, de_tT, n, det1); // n x n

     double *term1, *I, *term2, *zt;
     term1 = (double *) malloc((size_t)((nn)*sizeof(double)));
     I = (double *) malloc((size_t)((row)*sizeof(double)));
     term2 = (double *) malloc((size_t)((row)*sizeof(double)));
     zt = (double *) malloc((size_t)((row)*sizeof(double)));               
         
// term1 and term2
         for(i=0; i < nn; i++) {
             term1[i] = (sig_eta[0]/sig_e[0])* S[i];       
         }    
         for(i=0; i < row; i++) {
              I[i]= 1.0;
         }    
         MProd(I, constant, n, term1, n, term2);
 
     for(l=0; l < r1; l++) {
     for(t=0; t < T1; t++) {          
           extract_alt2(l, t, n, rT, T, XB, XB1);
           extract_alt2(l, t, n, rT, T, z, zT);
           MProd(zT, constant, n, term1, n, zt);
           for(i=0; i < row; i++) {
              mean1[i] = (XB1[i]+zt[i])/(1+term2[i]) + prior_omu[0];     
           }                  
           mvrnormal(constant, mean1, de_tT, n, o_1);     // random generator
           put_together1(l, t, n, r, T, opost, o_1);
     }
     } // End of loop year

         
     free(o_1); free(de_tT); free(det1); free(chi_tT);
     free(mean1); free(XB1); free(QXB1); free(zT);
     free(term1); free(I); free(term2); free(zt);

     return;
} 


// Rnadom-walk metropolis for phi
void phi_gp_MH(double *Qeta1, double *Qeta2, double *det1, double *det2,
     double *phi1, double *phi2, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *XB, double *o, int *constant, 
     double *accept, double *phip)
{
     
     int row, col, l, i, j, r1, T1, N1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     N1 = row*r1*T1;
     rT1 = *rT;
     
     double *ov, *o1, *XB1, *ratio, *U; 
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ratio = (double *) malloc((size_t)((col)*sizeof(double)));             
     U = (double *) malloc((size_t)((col)*sizeof(double)));         
    
     double u, v, a, b;
     u = 0.0;
     v = 0.0;

     for(l=0; l < r1; l++) {
         for(i=0; i < T1; i++) {                                
             extract_alt2(l, i, n, rT, T, o, o1);
             extract_alt2(l, i, n, rT, T, XB, XB1);
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-XB1[j];
             }       
             u += xTay2(ov, Qeta1, ov, row);
             v += xTay2(ov, Qeta2, ov, row);
         }
     }
     
     a = *prior_a;
     b = *prior_b;
     u =  0.5 * u;
     v =  0.5 * v;

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


     if(phi2[0] < 0.0010){
          phip[0] = phi1[0];
          accept[0] = 0.0;
     }
     else if(phi2[0] > 0.9999){
          phip[0] = phi1[0];
          accept[0] = 0.0;
     }
     else{    
     tr1 = (a-1.0)*log(phi1[0])-b*phi1[0]-0.5*rT1*log(det1[0])- u; 
     tr2 = (a-1.0)*log(phi2[0])-b*phi2[0]-0.5*rT1*log(det2[0])- v; 
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

     free(o1); free(ov); free(XB1); free(ratio); free(U);
     
     return;
}     


// Discrete sampling for phi
void phi_gp_DIS(int *cov, double *Qeta1, double *det1, double *phi1,  
     double *phis, int *phik, double *nu, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *d,
     double *sig2eta, double *XB, double *o, int *constant, double *accept, 
     double *phip)
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
        phidens_gp(phitmp, Qeta, det, n, r, T, rT, N, prior_a, prior_b, XB, 
        o, constant, out);
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
     phidens_gp(phi1, Qeta1, det1, n, r, T, rT, N, prior_a, prior_b, XB, 
     o, constant, tr1);
     ratio[0] = exp(tr2[0] + exp(tr2[0]) - tr1[0] - exp(tr1[0]));
     ratio_fnc(ratio, constant, U);
     if(U[0] < ratio[0]){
        phip[0] = phis[i];
     }             
     else {
        phip[0] = phi1[0];
     }     
     accept[0] = 0.0;

     free(ratio); free(tr2); free(tr1); free(pden); free(U);

     
     return;
}     


void phidens_gp(double *phi, double *Qeta, double *det, int *n, int *r, 
     int *T, int *rT, int *N, double *prior_a, double *prior_b, double *XB, 
     double *o, int *constant, double *out)
{
     int row, col, l, i, j, r1, T1, N1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     N1 = *N;
     rT1 = *rT;

     double *ov, *o1, *XB1; 
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     double u, a, b;
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1; i++) {                                
             extract_alt2(l, i, n, rT, T, o, o1);
             extract_alt2(l, i, n, rT, T, XB, XB1);
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-XB1[j];
             }       
             u += xTay2(ov, Qeta, ov, row);
         }
     }
     free(o1); free(ov); free(XB1); 
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

// Discrete sampling for nu
void nu_gp_DIS(int *cov, double *Qeta1, double *det1, double *phi,  
     double *nu1, int *n, int *r, int *T, int *rT, int *N, 
     double *d,  double *sig2eta, double *XB, double *o, int *constant, 
     double *nup)
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
/*
     nuk=10;
     double *nus;
     nus = (double *) malloc((size_t)((nuk)*sizeof(double)));             
     nus[0]=0.10; nus[1]=0.20; nus[2]=0.30; nus[3]=0.40; nus[4]=0.50;  
     nus[5]=0.60; nus[6]=0.70; nus[7]=0.80; nus[8]=0.90; nus[9]=1.0;  
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
        nudens_gp(Qeta, det, n, r, T, rT, N, XB, o, constant, out);
        pden[i] = out[0];
        u += out[0];
     }     
     free(nutmp); free(Qeta); free(det); free(out);
     
     double *pprob, *U, *tr2; 
     pprob = (double *) malloc((size_t)((nuk)*sizeof(double)));             
     U = (double *) malloc((size_t)((col)*sizeof(double)));
     tr2 = (double *) malloc((size_t)((col)*sizeof(double)));             

// cumulative prob
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
//     nup[0] = nus[i];
   
     tr2[0] = pden[i];  
     
     free(pprob);
     
     double *ratio, *tr1;
     ratio = (double *) malloc((size_t)((col)*sizeof(double)));             
     tr1 = (double *) malloc((size_t)((col)*sizeof(double)));             
     nudens_gp(Qeta1, det1, n, r, T, rT, N, XB, o, constant, tr1);
     ratio[0] = exp(tr2[0] + exp(tr2[0]) - tr1[0] - exp(tr1[0]));
     ratio_fnc(ratio, constant, U);
     if(U[0] < ratio[0]){
        nup[0] = nus[i];
     }             
     else {
        nup[0] = nu1[0];
     }     
//     Rprintf("   i: %i, ratio: %4.4f, U: %4.4f, nup: %4.4f \n", i, ratio[i],U[0],nup[0]);
      
     free(ratio); free(tr2); free(tr1); free(pden); free(U);
     free(nus);
     
     return;
}     

void nudens_gp(double *Qeta, double *det, int *n, int *r, int *T, int *rT, 
     int *N, double *XB, double *o, int *constant, double *out)
{
     int row, col, l, i, j, r1, T1, N1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     N1 = *N;
     rT1 = *rT;

     double *ov, *o1, *XB1; 
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     double u;
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1; i++) {                                
             extract_alt2(l, i, n, rT, T, o, o1);
             extract_alt2(l, i, n, rT, T, XB, XB1);
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-XB1[j];
             }       
             u += xTay2(ov, Qeta, ov, row);
         }
     }
     free(o1); free(ov); free(XB1); 

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

////////////////////// THE END //////////////////////////
