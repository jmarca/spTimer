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
     double *phi_a, double *phi_b,
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik, double *nu,
     double *d, double *sig_e, double *sig_eta, 
     double *beta, double *X, double *z, double *o, int *constant,
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *betap, double *op)
{     
     int n1, N1; 
     n1 = *n;
//     nn = n1*n1;
//     r1 = *r;
//     p1 = *p;
//     rn = r1 *n1;
     N1 = *N;
//     col = *constant;
     
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
     phi_a, phi_b, d, sig_eta, XB, o, constant, accept, phip);
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
      tmp[0] = log(phi[0]); 
// Rprintf("   phi: %4.4f, tmp: %4.4f, cov: %i\n", phi[0], tmp[0], cov[0]);      
      mvrnormal(constant, tmp, tau, constant, phi2);
      phi2[0]= exp(phi2[0]);      
// Rprintf("   phi: %4.4f, tmp: %4.4f, cov: %i\n", phi[0], tmp[0], cov[0]);      
      covFormat(cov, n, phi2, nup, d, sig_eta, S, det2, Sinv, Qeta2);   
// Rprintf("   phi: %4.4f, phi2: %4.4f, cov: %i\n", phi[0], phi2[0], cov[0]);
//     randow-walk M  
     phi_gp_MH(Qeta, Qeta2, det, det2, phi, phi2, n, r, T, rT, N, 
     phi_a, phi_b, XB, o, constant, accept, phip);
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
   sig_e_gp(n, r, T, rT, N, shape_e, prior_b, o, z, constant, sig_ep);
   sig_eta_gp(n, r, T, rT, shape_eta, prior_b, Sinv, XB, o, constant, sig_etap);
   o_gp(n, r, T, rT, prior_omu, prior_osig, sig_e, sig_etap, S, Qeta, 
   XB, z, constant, op);     
           
   free(Qeta); free(XB); free(Sinv); free(det); free(S); 
   return;
}


// Joint posterior distribution for temporal and spatial beta (dyn)
void JOINTsptp_gp(int *intercept, int *n, int *T, int *r, int *rT, int *p, int *q, int *u, int *N, 
     int *cov, int *spdecay, int *rhocheck,
     double *shape_e, double *shape_eta, double *shape_beta, double *shape_del, double *shape_0,
     double *prior_a, double *prior_b, double *prior_mubeta, double *prior_sigbeta, 
     double *prior_omu, double *prior_osig, 
     double *phi, double *tau, double *phis, int *phik, double *nu,
     double *d, double *sig_e, double *sig_eta, double *sig_beta, double *sigdelta, double *sig0,
     double *beta, double *betas, double *betat, double *rho, double *X, double *Xsp, double *Xtp, 
     double *z, double *o, int *constant,
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *sigbetap, double *sigdeltap, double *sig0p, double *rhop, 
     double *betap, double *betasp, double *betat0p, double *betatp, double *op)
{     
     int n1, p1, u1, N1; 
     n1 = *n;
//     nn = n1*n1;
//     r1 = *r;
     p1 = *p;
//     q1 = *q;
     u1 = *u;
//     rn = r1 *n1;
     N1 = *N;
//     col = *constant;
     
   double *Qeta, *XB, *XBtmp, *XBno, *XBsp, *XBtp, *Sinv, *det, *S, *G;
   Qeta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
   XB = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBtmp = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBno = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBsp = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBtp = (double *) malloc((size_t)((N1)*sizeof(double)));
   Sinv = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   det = (double *) malloc((size_t)((1)*sizeof(double)));   
   S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   G = (double *) malloc((size_t)((u1*u1)*sizeof(double)));   

   covFormat(cov, n, phi, nu, d, sig_eta, S, det, Sinv, Qeta);   
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
   comb_XB_sp(n, r, T, q, Xsp, betas, constant, XBsp); // N x 1
   MAdd(XBno, N, constant, XBsp, XBtmp);  // N x 1
   comb_XB_tp(n, r, T, u, Xtp, betat, constant, XBtp); // N x 1
   MAdd(XBtmp, N, constant, XBtp, XB);  // N x 1

//  int i; 
//  for(i=0; i< *N; i++){
//     Rprintf("  o: %4.4f, XB: %4.4f, XBno: %4.4f, XBsp: %4.4f, XBtp: %4.4f,\n", o[i], XB[i], XBno[i], XBsp[i], XBtp[i]);      
//  }
   
// check nu
   if(cov[0]==4){

      nu_gp_DIS_sptp(cov, Qeta, det, phi, nu, n, r, T, rT, N, d, sig_eta, XB, o, 
      constant, nup);

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

     phi_gp_DIS_sptp(cov, Qeta, det, phi, phis, phik, nup, n, r, T, rT, N, 
     prior_a, prior_b, d, sig_eta, XB, o, constant, accept, phip);

     covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
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
      tmp[0] = log(phi[0]); 
      mvrnormal(constant, tmp, tau, constant, phi2);
      phi2[0]= exp(phi2[0]);      
      covFormat(cov, n, phi2, nup, d, sig_eta, S, det2, Sinv, Qeta2);   
//     randow-walk M  

     phi_gp_MH_sptp(Qeta, Qeta2, det, det2, phi, phi2, n, r, T, rT, N, 
     prior_a, prior_b, XB, o, constant, accept, phip);

     if(accept[0] == 1.0){    
         covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
     }
     free(Qeta2); free(det2); free(tmp); free(phi2);
   }
   else {

   }   

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
     beta_gp_for_sptp(n, r, T, rT, p, prior_mubeta, prior_sigbeta, Qeta, X, 
     XBsp, XBtp, o, constant, betap); 
     MProd(betap, constant, p, X, N, XBno); // N x 1
   }

   beta_gp_sp(n, r, T, rT, q, N, prior_mubeta, sig_beta, betas, Qeta, Sinv, 
   Xsp, XBno, o, constant, betasp); 
   comb_XB_sp(n, r, T, q, Xsp, betasp, constant, XBsp); // N x 1
   MAdd(XBno, N, constant, XBsp, XBtmp);  // N x 1

   sig_beta_gp_sp(n, q, shape_beta, prior_b, betasp, Sinv, constant, sigbetap);

   free(XBsp); free(XBno);

   if(rhocheck[0] == 0){
     IdentityMX(u, rho, G);
     beta_gp_tp(n, r, T, rT, u, sig0, sigdelta, Qeta, G, betat, XBtmp, Xtp, o, 
     constant, betat0p, betatp);
     rho_gp_tp(u, T, prior_mubeta, prior_sigbeta, sigdelta, betat0p, betatp, 
     constant, rhop);
   }
   else {
     int i;
     for(i=0; i<u1; i++){   
     rhop[i] = rho[i];   
     }
     IdentityMX(u, rhop, G);
     beta_gp_tp(n, r, T, rT, u, sig0, sigdelta, Qeta, G, betat, XBtmp, Xtp, o, 
     constant, betat0p, betatp);
   }

//  int i, T1;
//  T1 = *T; 
//  for(i=0; i< u1*T1; i++){
//     Rprintf("   betatp: %4.4f, \n", betatp[i]);      
//  }
   comb_XB_tp(n, r, T, u, Xtp, betatp, constant, XBtp); // N x 1
   MAdd(XBtmp, N, constant, XBtp, XB);  // N x 1

   sig_0_gp_tp(u, shape_0, prior_b, betat0p, constant, sig0p);
   sig_del_gp_tp(u, T, shape_del, prior_b, betat0p, betatp, G, 
   constant, sigdeltap);

//  int i; 
//  for(i=0; i< *N; i++){
//     Rprintf("  o: %4.4f, XB: %4.4f, XBno: %4.4f, XBsp: %4.4f, XBtmp: %4.4f, XBtp: %4.4f,\n", o[i], XB[i], XBno[i], XBsp[i], XBtmp[i], XBtp[i]);      
//  }

   free(XBtp); free(XBtmp);

//   sig_beta_gp_sp(n, q, shape_beta, prior_b, betasp, Sinv, constant, sigbetap);
//   sig_0_gp_tp(u, shape_0, prior_b, betat0p, constant, sig0p);
//   sig_del_gp_tp(u, T, shape_del, prior_b, betat0p, betatp, G, 
//   constant, sigdeltap);

   sig_e_gp_sptp(n, r, T, rT, N, shape_e, prior_b, o, z, constant, sig_ep);

   sig_eta_gp_sptp(n, r, T, rT, shape_eta, prior_b, Sinv, XB, o, constant, sig_etap);

   o_gp_sptp(n, r, T, rT, prior_omu, prior_osig, sig_e, sig_eta, S, Qeta, 
   XB, z, constant, op);     
      
   free(Qeta); free(XB); free(Sinv); free(det); free(S); free(G);
   
   return;
}



// Joint posterior distribution for temporal beta (dyn)
void JOINTtp_gp(int *intercept, int *n, int *T, int *r, int *rT, int *p, int *u, int *N, 
     int *cov, int *spdecay, int *rhocheck,
     double *shape_e, double *shape_eta, double *shape_del, double *shape_0,
     double *prior_a, double *prior_b, double *prior_mubeta, double *prior_sigbeta, 
     double *prior_omu, double *prior_osig, 
     double *phi, double *tau, double *phis, int *phik, double *nu,
     double *d, double *sig_e, double *sig_eta, double *sigdelta, double *sig0,
     double *beta, double *betat, double *rho, double *X, double *Xtp, 
     double *z, double *o, int *constant,
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *sigdeltap, double *sig0p, double *rhop, double *betap, 
     double *betat0p, double *betatp, double *op)
{     
     int n1, p1, u1, N1; 
     n1 = *n;
//     nn = n1*n1;
//     r1 = *r;
     p1 = *p;
     u1 = *u;
//     rn = r1 *n1;
     N1 = *N;
//     col = *constant;
     
   double *Qeta, *XB, *XBno, *XBtp, *Sinv, *det, *S, *G;
   Qeta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
   XB = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBno = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBtp = (double *) malloc((size_t)((N1)*sizeof(double)));
   Sinv = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   det = (double *) malloc((size_t)((1)*sizeof(double)));   
   S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   G = (double *) malloc((size_t)((u1*u1)*sizeof(double)));   
   
   covFormat(cov, n, phi, nu, d, sig_eta, S, det, Sinv, Qeta);   
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
   IdentityMX(u, rho, G);
   beta_gp_tp(n, r, T, rT, u, sig0, sigdelta, Qeta, G, betat, XBno, Xtp, o, 
   constant, betat0p, betatp);
   comb_XB_tp(n, r, T, u, Xtp, betatp, constant, XBtp); // N x 1
   MAdd(XBno, N, constant, XBtp, XB);  // N x 1

//  int i; 
//  for(i=0; i< *N; i++){
//     Rprintf("  o: %4.4f, XB: %4.4f, XBno: %4.4f, XBtp: %4.4f, \n", o[i], XB[i], XBno[i], XBtp[i]);      
//  }

// check nu
   if(cov[0]==4){
      nu_gp_DIS_sptp(cov, Qeta, det, phi, nu, n, r, T, rT, N, d, sig_eta, XB, o, 
      constant, nup);
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
     phi_gp_DIS_sptp(cov, Qeta, det, phi, phis, phik, nup, n, r, T, rT, N, 
     prior_a, prior_b, d, sig_eta, XB, o, constant, accept, phip);
     covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
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
      tmp[0] = log(phi[0]); 
      mvrnormal(constant, tmp, tau, constant, phi2);
      phi2[0]= exp(phi2[0]);      
      covFormat(cov, n, phi2, nup, d, sig_eta, S, det2, Sinv, Qeta2);   
//     randow-walk M  
     phi_gp_MH_sptp(Qeta, Qeta2, det, det2, phi, phi2, n, r, T, rT, N, 
     prior_a, prior_b, XB, o, constant, accept, phip);
     if(accept[0] == 1.0){    
         covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
     }
     free(Qeta2); free(det2); free(tmp); free(phi2);
   }
   else {

   }   

   if(rhocheck[0] == 0){
     rho_gp_tp(u, T, prior_mubeta, prior_sigbeta, sigdelta, betat0p, betatp, 
     constant, rhop);
     IdentityMX(u, rhop, G);
   }
   else {
     int i;
     for(i=0; i<u1; i++){   
     rhop[i] = rho[i];   
     }
     IdentityMX(u, rhop, G);
   }        

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
     beta_gp(n, r, T, rT, p, prior_mubeta, prior_sigbeta, Qeta, X, o,
     constant, betap);
     MProd(betap, constant, p, X, N, XBno); // N x 1
   }

   beta_gp_tp(n, r, T, rT, u, sig0, sigdelta, Qeta, G, betatp, XBno, Xtp, o, 
   constant, betat0p, betatp); // betatp = u x T
   comb_XB_tp(n, r, T, u, Xtp, betatp, constant, XBtp); // N x 1
   MAdd(XBno, N, constant, XBtp, XB);  // N x 1

   free(XBtp); free(XBno);           

   sig_0_gp_tp(u, shape_0, prior_b, betat0p, constant, sig0p);
   sig_del_gp_tp(u, T, shape_del, prior_b, betat0p, betatp, G, 
   constant, sigdeltap);

   sig_e_gp_sptp(n, r, T, rT, N, shape_e, prior_b, o, z, constant, sig_ep);
   sig_eta_gp_sptp(n, r, T, rT, shape_eta, prior_b, Sinv, XB, o, constant, sig_etap);
   o_gp_sptp(n, r, T, rT, prior_omu, prior_osig, sig_e, sig_etap, S, Qeta, 
   XB, z, constant, op);     

   free(Qeta); free(XB); free(Sinv); free(det); free(S); free(G);
   
   return;
}



// Joint posterior distribution for spatial beta
void JOINTsp_gp(int *intercept, int *n, int *T, int *r, int *rT, int *p, int *q, int *N, 
     int *cov, int *spdecay, double *shape_e, double *shape_eta, double *shape_beta,
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik, double *nu,
     double *d, double *sig_e, double *sig_eta, double *sig_beta,
     double *beta, double *betas, double *X, double *Xsp, double *z, double *o, 
     int *constant,
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *sig_betasp, double *betap, double *betasp, double *op)
{     
     int n1, p1, N1; 
     n1 = *n;
//     nn = n1*n1;
//     r1 = *r;
     p1 = *p;
//     rn = r1 *n1;
     N1 = *N;
//     col = *constant;
     
   double *Qeta, *XB, *XBno, *XBsp, *Sinv, *det, *S;
   Qeta = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
   XB = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBno = (double *) malloc((size_t)((N1)*sizeof(double)));
   XBsp = (double *) malloc((size_t)((N1)*sizeof(double)));
   Sinv = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   
   det = (double *) malloc((size_t)((1)*sizeof(double)));   
   S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));   

   covFormat(cov, n, phi, nu, d, sig_eta, S, det, Sinv, Qeta);   
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
   comb_XB_sp(n, r, T, q, Xsp, betas, constant, XBsp); // N x 1
   MAdd(XBno, N, constant, XBsp, XB);  // N x 1

//  int i; 
//  for(i=0; i< *N; i++){
//     Rprintf("  o: %4.4f, XB: %4.4f, XBno: %4.4f, XBsp: %4.4f, \n", o[i], XB[i], XBno[i], XBsp[i]);      
//  }

// check nu
   if(cov[0]==4){
      nu_gp_DIS_sptp(cov, Qeta, det, phi, nu, n, r, T, rT, N, d, sig_eta, XB, o, 
      constant, nup);
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
     phi_gp_DIS_sptp(cov, Qeta, det, phi, phis, phik, nup, n, r, T, rT, N, 
     prior_a, prior_b, d, sig_eta, XB, o, constant, accept, phip);
     covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
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
      covFormat(cov, n, phi2, nup, d, sig_eta, S, det2, Sinv, Qeta2);   
//     randow-walk M  
     phi_gp_MH_sptp(Qeta, Qeta2, det, det2, phi, phi2, n, r, T, rT, N, 
     prior_a, prior_b, XB, o, constant, accept, phip);
     if(accept[0] == 1.0){    
         covFormat(cov, n, phip, nup, d, sig_eta, S, det, Sinv, Qeta);   
     }
     free(Qeta2); free(det2); free(tmp); free(phi2);
   }
   else {
   }

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
     beta_gp_for_sp(n, r, T, rT, p, prior_mubeta, prior_sigbeta, Qeta, X, 
     XBsp, o, constant, betap); 
     MProd(betap, constant, p, X, N, XBno); // N x 1
   }
   
   beta_gp_sp(n, r, T, rT, q, N, prior_mubeta, sig_beta, betas, Qeta, Sinv, 
   Xsp, XBno, o, constant, betasp); 
   
   comb_XB_sp(n, r, T, q, Xsp, betasp, constant, XBsp); // N x 1
   MAdd(XBno, N, constant, XBsp, XB);  // N x 1

   free(XBsp); free(XBno);

//   double *shape_beta;
//   shape_beta = (double *) malloc((size_t)((1)*sizeof(double)));   
//   free(shape_beta);  
   
   sig_e_gp_sptp(n, r, T, rT, N, shape_e, prior_b, o, z, constant, sig_ep);
   sig_eta_gp_sptp(n, r, T, rT, shape_eta, prior_b, Sinv, XB, o, constant, sig_etap);
   o_gp_sptp(n, r, T, rT, prior_omu, prior_osig, sig_e, sig_eta, S, Qeta, 
   XB, z, constant, op); 
   
//   beta_gp_for_sp(n, r, T, rT, p, prior_mubeta, prior_sigbeta, Qeta, X, 
//   XBsp, op, constant, betap);        
//   shape_beta[0] = (n[0]*q[0])/2.0 + prior_a[0];
   sig_beta_gp_sp(n, q, shape_beta, prior_b, betasp, Sinv, constant, sig_betasp);

/*
  int i; 
  for(i=0; i< *N; i++){
     Rprintf("   z: %4.4f, o: %4.4f, \n", z[i], o[i]);      
  }
   
  int i; 
  for(i=0; i< *N; i++){
     Rprintf("   XB: %4.4f, \n", XB[i]);      
  }
*/
   free(Qeta); free(XB); free(Sinv); free(det); free(S); 
   return;
}


// Posterior distribution for "sig_e"
void sig_e_gp(int *n, int *r, int *T, int *rT, int *N, double *shape, 
     double *prior_b, double *o, double *z, int *constant, double *sig2e)
{
     int i, t, l, n1, r1, col;
     n1 =*n;
     r1 =*r;
//     T1 =*T;
     col =*constant;
     
     double *z1, *o1, *zo, *zzoo, *tmp;
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     o1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zo = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zzoo = (double *) malloc((size_t)((col)*sizeof(double)));
     tmp = (double *) malloc((size_t)((col)*sizeof(double)));     
          
     double u, v, b;
     u = 0.0;
     v = 0.0;
     b = 0.0;

     int *T1; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }

     for(l=0; l<r1; l++){
        for(t=0; t< T1[l]; t++){
             extract_alt_uneqT(l, t, n, r, T, rT, o, o1);
             extract_alt_uneqT(l, t, n, r, T, rT, z, z1);
//             extract_alt2(l, t, n, rT, T, o, o1);
//             extract_alt2(l, t, n, rT, T, z, z1);
             for(i=0; i<n1; i++){
                 zzoo[0] = z1[i]-o1[i];
                 tmp[0] = 0.005;
                 mvrnormal(constant, zzoo, tmp, constant, zzoo);                              
                 zo[i] = zzoo[0];
             }
             MProd(zo, constant, n, zo, constant, zzoo);
             u += zzoo[0];          
        }
     }      

     b = *prior_b;
     u = b + 0.5 * u;
     v = *shape;
     sig2e[0] = rigammaa(v, u);

     free(T1);
     free(z1); free(o1); free(zo); free(zzoo); free(tmp);
     return;                  
}     

// Posterior distribution for "sig_e"
// for sp tp model
void sig_e_gp_sptp(int *n, int *r, int *T, int *rT, int *N, double *shape, 
     double *prior_b, double *o, double *z, int *constant, double *sig2e)
{
     int i, t, l, n1, r1, T1, col;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     double *z1, *o1, *zo, *zzoo, *tmp;
     z1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     o1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zo = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     zzoo = (double *) malloc((size_t)((col)*sizeof(double)));
     tmp = (double *) malloc((size_t)((col)*sizeof(double)));     
          
     double u, v, b;
     u = 0.0;
     v = 0.0;
     b = 0.0;

     for(l=0; l<r1; l++){
        for(t=0; t< T1; t++){
             extract_alt2(l, t, n, rT, T, o, o1);
             extract_alt2(l, t, n, rT, T, z, z1);
             for(i=0; i<n1; i++){
                 zzoo[0] = z1[i]-o1[i];
                 tmp[0] = 0.005;
                 mvrnormal(constant, zzoo, tmp, constant, zzoo);                              
                 zo[i] = zzoo[0];
             }
             MProd(zo, constant, n, zo, constant, zzoo);
             u += zzoo[0];          
        }
     }      

     b = *prior_b;
     u = b + 0.5 * u;
     v = *shape;
     sig2e[0] = rigammaa(v, u);

     free(z1); free(o1); free(zo); free(zzoo); free(tmp);
     return;                  
}     


// Posterior distribution for "sig_eta"
void sig_eta_gp(int *n, int *r,  int *T, int *rT, double *shape, 
     double *prior_b, double *Sinv, double *XB, double *o, int *constant, 
     double *sig2eta)     
{
     double *ov, *o1, *out, u, b, sh;
     int row, col, l, i, j, r1;
     row = *n;
     col = *constant;
     r1 = *r;
//     T1 = *T;
     
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     out = (double *) malloc((size_t)((row*col)*sizeof(double)));
          
     double *XB1;
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     int *T1; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
               
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1[l]; i++) {                                
         extract_alt_uneqT(l, i, n, r, T, rT, o, o1);
         extract_alt_uneqT(l, i, n, r, T, rT, XB, XB1);
//             extract_alt2(l, i, n, rT, T, o, o1);
//             extract_alt2(l, i, n, rT, T, XB, XB1);
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
     sig2eta[0] = rigammaa(sh, u);
     
     free(T1);
     free(ov); free(o1); free(XB1); free(out);
     return;
}


// Posterior distribution for "sig_eta"
// for sp tp model
void sig_eta_gp_sptp(int *n, int *r,  int *T, int *rT, double *shape, 
     double *prior_b, double *Sinv, double *XB, double *o, int *constant, 
     double *sig2eta)     
{
     double *ov, *o1, *out, u, b, sh;
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
     sig2eta[0] = rigammaa(sh, u);
     
     free(ov); free(o1); free(XB1); free(out);
     return;
}

// Posterior distribution for "sig_beta" for spatial varying beta
// betasp = n x q
void sig_beta_gp_sp(int *n, int *q, double *shape, double *prior_b, double *betasp, 
     double *Sinv, int *constant, double *sig2beta)
{
     int i, j, n1, q1;
     n1 =*n;
     q1 =*q;
//     col =*constant;
     
     double *bt;
     bt = (double *) malloc((size_t)((n1)*sizeof(double)));
          
     double u, v, b;
     u = 0.0;
     v = 0.0;
     b = 0.0;

     for(j=0; j<q1; j++){
         for(i=0; i<n1; i++){
           bt[i] = betasp[i+j*n1];
         }
         u += xTay2(bt, Sinv, bt, n1);
     }      

     b = *prior_b;
     u = b + 0.5 * u;
     v = *shape;
     sig2beta[0] = rigammaa(v, u);

     free(bt); 
     return;                  
}     


// Posterior distribution for "sig_0": dyna
// initial variance for temporal parameter gamma, gam_0 = u x 1
void sig_0_gp_tp(int *u, double *shape, double *prior_b, double *gam_0, 
     int *constant, double *sig0)          
{
     int col;
//     u1 = *u;
     col = *constant;

     double sh, rt, *out, b, v;
     out = (double *) malloc((size_t)((col)*sizeof(double)));
     sh = 0.0;
     rt = 0.0;
     b = 0.0;
     v = 0.0;
     MProd(gam_0, constant, u, gam_0, constant, out);
     v = out[0];
     b = *prior_b;
     rt = b + 0.5 * v;            
     sh = *shape;
     sig0[0] = rigammaa(sh, rt);

     free(out);
     return;
}

// Posterior distribution for "sig_del": dyna
// initial variance for temporal parameter gamma, gam = u x T
void sig_del_gp_tp(int *u, int *T, double *shape, double *prior_b, double *gam_0, 
     double *gam, double *G, int *constant, double *sigdelta)          
{
     int i, t, T1, u1, col;
     T1 = *T;
     u1 = *u;
     col = *constant;

     double v, *out, *out1;
     out = (double *) malloc((size_t)((u1*col)*sizeof(double)));
     out1 = (double *) malloc((size_t)((col)*sizeof(double)));
     v = 0.0;
     for(t=0; t < T1; t++) {          
            if(t == 0){
               MProd(gam_0, constant, u, G, u, out); // u x 1
               for(i=0; i<u1; i++){
                   out[i] = gam[i+t*u1] - out[i];     
               }
            }
            else {
               for(i=0; i<u1; i++){
                   out[i] = gam[i+(t-1)*u1];
               }                          
               MProd(out, constant, u, G, u, out); // u x 1
               for(i=0; i<u1; i++){
                   out[i] = gam[i+t*u1] - out[i];     
               }
            }
            MProd(out,constant,u,out,constant, out1); // 1 x 1
            v += out1[0];
     } // End of loop year
         
     double sh, rt, b;
     sh = 0.0;
     rt = 0.0;
     b = 0.0;
     sh = *shape;
     b = *prior_b;
     rt = b + 0.5 * v;            
     sigdelta[0] = rigammaa(sh, rt);

     free(out); free(out1);
     return;
}


// conditional posterior for betatp
// temporal dynamics only
// Xtp = nrT x u; betat = u x T; 
void beta_gp_tp(int *n, int *r, int *T, int *rT, int *u, double *sig0, 
     double *sigdelta, double *Qeta, double *G, double *betat, 
     double *XB, double *Xtp, double *o, int *constant, double *betat0p, 
     double *betatp)     
{
     int i, t, l, n1, r1, T1, u1, uu, col;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     u1 =*u;          
     uu = u1*u1;
     col =*constant;

     double *GG, *I, *del, *chi, *mn, *W1, *tW1, *out, *tW1QW1, *XB1, *o1, *det, *check;
     GG = (double *) malloc((size_t)((u1*u1)*sizeof(double)));
     I = (double *) malloc((size_t)((u1*u1)*sizeof(double)));
     del = (double *) malloc((size_t)((u1*u1)*sizeof(double)));
     chi = (double *) malloc((size_t)((u1*col)*sizeof(double)));
     mn = (double *) malloc((size_t)((u1*col)*sizeof(double)));
     W1 = (double *) malloc((size_t)((u1*n1)*sizeof(double)));
     tW1 = (double *) malloc((size_t)((u1*n1)*sizeof(double)));
     out = (double *) malloc((size_t)((u1*n1)*sizeof(double)));
     tW1QW1 = (double *) malloc((size_t)((u1*u1)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((n1)*sizeof(double)));
     o1 = (double *) malloc((size_t)((n1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     check = (double *) malloc((size_t)((n1)*sizeof(double)));
               
     MProd(G, u, u, G, u, GG); // u x u
     IdentityM(u, I); // u x u diagonal matrix

// for initial value of gamma0     
     for(i=0; i<uu; i++){
         del[i] = GG[i]/sigdelta[0] + (1.0/sig0[0])*I[i]; // u x u
     }
     MInv(del, del, u, det); // u x u
     for(i=0; i<u1; i++){
         chi[i] =  betat[i]; // u  x 1
     }
     MProd(chi, constant, u, G, u, chi); // u x 1    
     for(i=0; i<u1; i++){
         chi[i] = chi[i]+ (1.0/sigdelta[0]);
     }               
     MProd(chi, constant, u, del, u, mn);  // u x 1    
     mvrnormal(constant, mn, del, u, mn); // u x 1
     for(i=0; i<u1; i++){
         betat0p[i] = mn[i]; 
     }

// for other gamma     
     for(t=0; t < T1; t++) {  
       for(i=0; i<uu; i++){
         del[i] = 0.0;
       }
       for(i=0; i<u1; i++){
         chi[i]=0.0;
       }                    
// for delta part
       for(l=0; l<r1; l++){
         extract_X(t, l, n, r, T, u, Xtp, W1); // extract nrT x u into n x u matrix
         MTranspose(W1, u, n, tW1);         // u x n
         MProd(W1, u, n, Qeta, n, out);   // n x u
         MProd(out, u, n, tW1, u, tW1QW1); // u x u
         MAdd(del, u, u, tW1QW1, del);  // u x u
       }
       for(i=0; i < uu; i++) {
         if(t == (T1-1)){       
            del[i] =  del[i]+(1.0/sigdelta[0])*(I[i]);      
         }
         else {        
            del[i] =  del[i]+(1.0/sigdelta[0])*(I[i]+GG[i]);
         }
       }    
       MInv(del, del, u, det); // u x u
// for chi part
       for(l=0; l<r1; l++){
         extract_alt2(l, t, n, rT, T, XB, XB1); // n x 1
         extract_alt2(l, t, n, rT, T, o, o1); // n x 1
         for(i=0; i < n1; i++) {
             check[i] =  o1[i] - XB1[i]; // n x 1
         }    
         extract_X(t, l, n, r, T, u, Xtp, W1); // extract nrT x u into n x u matrix
         MProd(W1, u, n, Qeta, n, out);   // n x u
         MProd(check, constant, n, out, u, mn);   // u x 1
         MAdd(chi, u, constant, mn, chi);  // u x 1
       }
       for(i=0; i<u1; i++){
           if(t == 0){
              mn[i] = (betat[i+t*u1] + betat0p[i])/sigdelta[0]; // u x 1
           }
           else if(t == (T1-1)){
              mn[i] = (betat[i+(t-1)*u1])/sigdelta[0]; // u x 1
           }
           else {
              mn[i] = (betat[i+(t-1)*u1] + betat[i+t*u1])/sigdelta[0]; // u x 1
           }                
       }
       MProd(mn, constant, u, G, u, mn); // u x 1
       for(i=0; i<u1; i++){
         chi[i] = chi[i]+ mn[i];
       }               
       MProd(chi, constant, u, del, u, mn);  // u x 1    
       mvrnormal(constant, mn, del, u, mn); // n x 1
       for(i=0; i<u1; i++){
         betatp[i+t*u1] = mn[i]; 
       }
     }
// end 

     free(GG); free(I); free(del); free(chi); free(mn); free(W1); free(tW1);
     free(out); free(tW1QW1); free(XB1); free(o1); free(det); free(check);

     return;
} 

// rho_k for dynamic .. tp
void rho_gp_tp(int *u, int *T, double *prior_mu, double *prior_sig, 
     double *sigdelta, double *gam0, double *gamma, int *constant, double *rhop)
{
     int i, t, u1, T1;
     u1 = *u;
     T1 = *T;
     
     double *del, *chi, *vv, *out1, *out2;
     del = (double *) malloc((size_t)((1)*sizeof(double)));
     chi = (double *) malloc((size_t)((1)*sizeof(double)));     
     vv = (double *) malloc((size_t)((1)*sizeof(double)));     
     out1 = (double *) malloc((size_t)((1)*sizeof(double)));     
     out2 = (double *) malloc((size_t)((1)*sizeof(double)));     
          
     for(i=0; i<u1; i++){
         out1[0] = 0.0;
         out2[0] = 0.0;
         for(t=0; t<T1; t++){
             if(t==0){
               vv[0] = gam0[i]*gam0[i];
               MAdd(out1, constant, constant, vv, out1);  // 1 x 1             
               vv[0] = gam0[i]*gamma[i+t*u1];
               MAdd(out2, constant, constant, vv, out2);  // 1 x 1             
             }
             else {
               vv[0] = gamma[i+t*u1]*gamma[i+t*u1];
               MAdd(out1, constant, constant, vv, out1);  // 1 x 1             
               vv[0] = gamma[i+(t-1)*u1]*gamma[i+t*u1];
               MAdd(out2, constant, constant, vv, out2);  // 1 x 1 
             }            
         } // end of t     
                         
         del[0] = 1.0/((1.0/sigdelta[0])*out1[0] + (1.0/prior_sig[0]));
         chi[0] = (1.0/sigdelta[0])*out2[0] + (prior_mu[0]/prior_sig[0]);

//     Rprintf("   out2: %4.4f, out1: %4.4f, chi: %4.4f, del: %4.4f, \n", out2[0], out1[0], chi[0], del[0]);      

         chi[0] = chi[0]*del[0];
         mvrnormal(constant, chi, del, constant, chi); // 1 x 1
         if(chi[0] > 1.0){
              chi[0] = 1.0;
         }
         else if(chi[0] < -1.0){
              chi[0] =  -1.0;
         }
         else{
              chi[0] = chi[0];     
         }   
         rhop[i] = chi[0];
     } // end of i

     free(del); free(chi); free(vv); free(out1); free(out2);
     return;
}      



// Posterior distribution for spatial "beta"
void beta_gp_sp(int *n, int *r, int *T, int *rT, int *q, int *N, double *prior_mu,
     double *prior_sig2betasp, double *betas, double *Qeta, double *Sinv, 
     double *Xsp, double *XB, double *o, int *constant, double *betap) 
{
     int t, l, i, j, n1, q1, r1, T1, col;
     n1 =*n;
     q1 =*q;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     double *del, *chi, *ot1, *XB1, *X1, *out, *det, *mu, *I, *tmp;
     del = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     chi = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
     ot1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));
     X1 = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     out = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     det = (double *) malloc((size_t)((col)*sizeof(double)));
     mu = (double *) malloc((size_t)((n1*col)*sizeof(double)));     
     I = (double *) malloc((size_t)((n1*n1)*sizeof(double)));          
     tmp = (double *) malloc((size_t)((n1*n1)*sizeof(double)));          

     double *bt, *xs, *XBrest, *XB1rest; 
     bt = (double *) malloc((size_t)((n1*(q1-1))*sizeof(double)));          
     xs = (double *) malloc((size_t)((n1*r1*T1*(q1-1))*sizeof(double)));          
     XBrest = (double *) malloc((size_t)((n1*r1*T1)*sizeof(double)));          
     XB1rest = (double *) malloc((size_t)((n1)*sizeof(double)));          

     int *qq;
     qq = (int *) malloc((size_t)((col)*sizeof(int)));          

     IdentityM(n, I); // n x n
     
     for(j=0; j<q1; j++){
       for(i=0; i<n1*n1; i++){
           del[i] = 0.0;
           tmp[i] = 0.0;
       }   
       for(i=0; i<n1; i++){
           chi[i] = 0.0;
       }   

       // betas, n x q  === n x (q-1)
       // Xsp, N x q === N x (q-1)
       extract_beta_sp2(j, n, q, betas, bt); // n x (q-1)
       extract_beta_sp2(j, N, q, Xsp, xs); // N x (q-1)
       qq[0] = q1-1;
       comb_XB_sp(n, r, T, qq, xs, bt, constant, XBrest); // N x 1

       for(l=0; l<r1; l++){
          for(t=0; t<T1; t++){
           extract_X_sp2(t, l, j, n, r, T, Xsp, X1); // n x n diagonal matrix
           MProd(X1, n, n, Qeta, n, out);   // n x n
           MProd(out, n, n, X1, n, out);   // n x n
           MAdd(del, n, n, out, del);  // n x n

           extract_X_sp2(t, l, j, n, r, T, Xsp, X1); // n x n diagonal matrix
           MProd(X1, n, n, X1, n, out);   // n x n
           MAdd(tmp, n, n, out, tmp);  // n x n

           extract_alt2(l, t, n, rT, T, o, ot1);  // n x 1
           extract_alt2(l, t, n, rT, T, XB, XB1); // n x 1        
           if(q1 > 1){
           extract_alt2(l, t, n, rT, T, XBrest, XB1rest); // n x 1        
             for(i=0; i<n1; i++){
                ot1[i] = ot1[i]-XB1[i]-XB1rest[i];
             }
           }
           else {
             for(i=0; i<n1; i++){
                ot1[i] = ot1[i]-XB1[i];
             }
           }   
           MProd(ot1, constant, n, X1, n, ot1);  // n x 1
           MAdd(chi, n, constant, ot1, chi);  // n x 1
          }
       }
       for(i=0; i<n1*n1; i++){
           del[i] = del[i] + Sinv[i]*(1.0/prior_sig2betasp[0]);
           tmp[i] = tmp[i] + I[i]*(1.0/prior_sig2betasp[0]);
       } 
       for(i=0; i<n1; i++){
           chi[i] = chi[i]  + prior_mu[0]/prior_sig2betasp[0];
       }
       MInv(del, del, n, det);
       MInv(tmp, tmp, n, det);     
       MProd(chi, constant, n, tmp, n, mu);  // n x 1    
//       mvrnormal(constant, mu, del, n, mu); // n x 1
       for(i=0; i<n1; i++){
           betap[i+j*n1] = mu[i]; 
       }
     }
     
     free(del); free(chi); free(ot1); free(XB1); free(X1); 
     free(out); free(det); free(mu); free(I); free(tmp); 
     free(bt); free(xs); free(qq); free(XBrest); free(XB1rest); 
     
     return;
// for(i=0; i<n1*n1; i++){
//     Rprintf("   X: %4.4f, del: %4.4f,  delinv: %4.4f, \n", X1[i], del[i], out[i]);      
//   }
//  for(i=0; i<n1*n1; i++){
//     Rprintf("   X: %4.4f, \n", X1[i]);      
//  }


}     


// Posterior distribution for "beta" when spatial and temporal beta is used
void beta_gp_for_sptp(int *n, int *r, int *T, int *rT, int *p, double *prior_mu,
     double *prior_sig, double *Qeta, double *X, double *XBsp, double *XBtp, 
     double *o, int *constant, double *betap) 
{
     int t, l, i, n1, p1, r1, T1, col;
     n1 =*n;
     p1 =*p;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     
     double *del, *chi, *ot1, *X1, *tX1, *out, *tX1QX1, *tX1Qo, *det, *mu, *I, *XBsp1, *XBtp1;
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
     XBsp1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));                         
     XBtp1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));                         
          
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
          extract_alt2(l, t, n, rT, T, XBsp, XBsp1);  // n x 1
          extract_alt2(l, t, n, rT, T, XBtp, XBtp1);  // n x 1
          for(i=0; i<n1; i++){
              ot1[i] = ot1[i] - XBsp1[i] - XBtp1[i];
          }     
          MProd(ot1, constant, n, Qeta, n, out); // n x 1
          MProd(out, constant, n, tX1, p, tX1Qo);  // p x 1
          MAdd(chi, p, constant, tX1Qo, chi);  // p x 1

     }
     }
 
     free(ot1); free(X1); free(tX1); free(out); free(tX1QX1); 
     free(tX1Qo); free(XBsp1); free(XBtp1);

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

     free(del); free(chi); free(mu); free(det); 
     
     return;
}     

// Posterior distribution for "beta" when spatial beta is used
void beta_gp_for_sp(int *n, int *r, int *T, int *rT, int *p, double *prior_mu,
     double *prior_sig, double *Qeta, double *X, double *XBsp, double *o, 
     int *constant, double *betap) 
{
     int t, l, i, n1, p1, r1, T1, col;
     n1 =*n;
     p1 =*p;
     r1 =*r;
     T1 =*T;
     col =*constant;
     
     
     double *del, *chi, *ot1, *X1, *tX1, *out, *tX1QX1, *tX1Qo, *det, *mu, *I, *XBsp1;
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
     XBsp1 = (double *) malloc((size_t)((n1*col)*sizeof(double)));                         
          
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
          extract_alt2(l, t, n, rT, T, XBsp, XBsp1);  // n x 1
          for(i=0; i<n1; i++){
              ot1[i] = ot1[i] - XBsp1[i];
          }     
          MProd(ot1, constant, n, Qeta, n, out); // n x 1
          MProd(out, constant, n, tX1, p, tX1Qo);  // p x 1
          MAdd(chi, p, constant, tX1Qo, chi);  // p x 1

     }
     }
 
     free(ot1); free(X1); free(tX1); free(out); free(tX1QX1); 
     free(tX1Qo); free(XBsp1);

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

     free(del); free(chi); free(mu); free(det); 
     
     return;
}     

// Posterior distribution for "beta" not spatial
void beta_gp(int *n, int *r, int *T, int *rT, int *p, double *prior_mu,
     double *prior_sig, double *Qeta, double *X, double *o, int *constant, 
     double *betap) 
{
     int t, l, i, n1, p1, r1, col;
     n1 =*n;
     p1 =*p;
     r1 =*r;
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
 
     int *T1; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }

     for(l=0; l<r1; l++){
     for(t=0; t<T1[l]; t++){
              //// upto this poit

          extract_X_uneqT(t, l, n, r, T, rT, p, X, X1);
//          extract_X(t, l, n, r, T, p, X, X1);    // n x p
          MTranspose(X1, p, n, tX1);         // p x n
          MProd(X1, p, n, Qeta, n, out);   // n x p
          MProd(out, p, n, tX1, p, tX1QX1); // pxp
          MAdd(del, p, p, tX1QX1, del);  // pxp

          extract_alt_uneqT(l, t, n, r, T, rT, o, ot1);
//          extract_alt2(l, t, n, rT, T, o, ot1);  // n x 1
          MProd(ot1, constant, n, Qeta, n, out); // n x 1
          MProd(out, constant, n, tX1, p, tX1Qo);  // p x 1
          MAdd(chi, p, constant, tX1Qo, chi);  // p x 1

     }
     }

     free(T1);
     free(ot1); free(X1); free(tX1);
     free(out); free(tX1QX1); free(tX1Qo);
     
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
     mvrnormal(constant, mu, del, p, betap); // p x 1


     free(del); free(chi); free(det); free(mu);
     
     return;
}     



// conditional posterior for o_lt
void o_gp(int *n, int *r, int *T, int *rT, double *prior_omu,
     double *prior_osig, double *sig_e, double *sig_eta, double *S, 
     double *Qeta, double *XB, double *z, int *constant, double *opost)     
{
     int i, l, t, r1, nn, row, col; 
     r1 = *r;     row = *n;    nn = row * row;  col = *constant;
     
     double *o_1, *de_tT, *det1, *chi_tT, *mean1, *XB1, *QXB1, *zT;
     
     o_1 = (double *) malloc((size_t)((row)*sizeof(double)));
     de_tT = (double *) malloc((size_t)((nn)*sizeof(double)));
     det1 = (double *) malloc((size_t)((col)*sizeof(double)));
     chi_tT = (double *) malloc((size_t)((row)*sizeof(double)));
     mean1 = (double *) malloc((size_t)((row)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row)*sizeof(double)));
     QXB1 = (double *) malloc((size_t)((row)*sizeof(double)));    
     zT = (double *) malloc((size_t)((row)*sizeof(double)));          
               
     int *T1; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
     
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
     for(t=0; t < T1[l]; t++) {          

           extract_alt_uneqT(l, t, n, r, T, rT, XB, XB1);
           extract_alt_uneqT(l, t, n, r, T, rT, z, zT);
//           extract_alt2(l, t, n, rT, T, XB, XB1);
//           extract_alt2(l, t, n, rT, T, z, zT);
           MProd(zT, constant, n, term1, n, zt);
           for(i=0; i < row; i++) {
              mean1[i] = (XB1[i]+zt[i])/(1+term2[i]) + prior_omu[0];     
           }                  
           mvrnormal(constant, mean1, de_tT, n, o_1);     // random generator
           put_together1_uneqT(l, t, n, r, T, rT, opost, o_1);
//           put_together1(l, t, n, r, T, opost, o_1);
     }
     } // End of loop year

     free(T1);    
     free(o_1); free(de_tT); free(det1); free(chi_tT);
     free(mean1); free(XB1); free(QXB1); free(zT);
     free(term1); free(I); free(term2); free(zt);

     return;
} 



// conditional posterior for o_lt
// for sp tp model
void o_gp_sptp(int *n, int *r, int *T, int *rT, double *prior_omu,
     double *prior_osig, double *sig_e, double *sig_eta, double *S, 
     double *Qeta, double *XB, double *z, int *constant, double *opost)     
{
     int i, l, t, r1, T1, nn, row, col; 
     r1 = *r;   T1= *T;  row = *n;    nn = row * row;  col = *constant;
     
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
     
     int row, col, l, i, j, r1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     rT1 = *rT;
//     N1 = row*rT1;
     
     double *ov, *o1, *XB1, *ratio, *U; 
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ratio = (double *) malloc((size_t)((col)*sizeof(double)));             
     U = (double *) malloc((size_t)((col)*sizeof(double)));         

     int *T1; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }
         
     double u, v, a, b;
     u = 0.0;
     v = 0.0;

     for(l=0; l < r1; l++) {
         for(i=0; i < T1[l]; i++) {                                
             extract_alt_uneqT(l, i, n, r, T, rT, o, o1);
             extract_alt_uneqT(l, i, n, r, T, rT, XB, XB1);           
//             extract_alt2(l, i, n, rT, T, o, o1);
//             extract_alt2(l, i, n, rT, T, XB, XB1);
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
     //else if(phi2[0] > 0.9999){
     //     phip[0] = phi1[0];
     //     accept[0] = 0.0;
     //}
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

     free(T1);
     free(o1); free(ov); free(XB1); free(ratio); free(U);
     
     return;
}     


// Rnadom-walk metropolis for phi
// for sp tp model
void phi_gp_MH_sptp(double *Qeta1, double *Qeta2, double *det1, double *det2,
     double *phi1, double *phi2, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *XB, double *o, int *constant, 
     double *accept, double *phip)
{
     
     int row, col, l, i, j, r1, N1, rT1, T1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 =*T;
     rT1 = *rT;
     N1 = row*rT1;
     
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
     //else if(phi2[0] > 0.9999){
     //     phip[0] = phi1[0];
     //     accept[0] = 0.0;
     //}
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
    
     int row, col, i, pk;
     row = *n;
     col = *constant;
//     r1 = *r;
//     rT1 = *rT;
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
     int row, col, l, i, j, r1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
//     N1 = *N;
     rT1 = *rT;

     double *ov, *o1, *XB1; 
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     int *T1; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }

     double u, a, b;
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1[l]; i++) {                                
             extract_alt_uneqT(l, i, n, r, T, rT, o, o1);
             extract_alt_uneqT(l, i, n, r, T, rT, XB, XB1);  
//             extract_alt2(l, i, n, rT, T, o, o1);
//             extract_alt2(l, i, n, rT, T, XB, XB1);
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-XB1[j];
             }       
             u += xTay2(ov, Qeta, ov, row);
         }
     }
     
     free(T1);
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

// Discrete sampling for phi
// for sptp
void phi_gp_DIS_sptp(int *cov, double *Qeta1, double *det1, double *phi1,  
     double *phis, int *phik, double *nu, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *d,
     double *sig2eta, double *XB, double *o, int *constant, double *accept, 
     double *phip)
{
    
     int row, col, i, pk;
     row = *n;
     col = *constant;
//     r1 = *r;
//     rT1 = *rT;
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

        phidens_gp_sptp(phitmp, Qeta, det, n, r, T, rT, N, prior_a, prior_b, XB, 
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
     phidens_gp_sptp(phi1, Qeta1, det1, n, r, T, rT, N, prior_a, prior_b, XB, 
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

// for sp tp models
void phidens_gp_sptp(double *phi, double *Qeta, double *det, int *n, int *r, 
     int *T, int *rT, int *N, double *prior_a, double *prior_b, double *XB, 
     double *o, int *constant, double *out)
{
     int row, col, l, i, j, r1, T1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
//     N1 = *N;
     T1 = *T;
     rT1 =*rT;

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
    
     int row, col, i;
     row = *n;
     col = *constant;
//     r1 = *r;
//     rT1 = *rT;
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

/*
     nuk=20;
     double *nus;
     nus = (double *) malloc((size_t)((nuk)*sizeof(double)));             
     nus[0]=0.05; nus[1]=0.10; nus[2]=0.15; nus[3]=0.20; nus[4]=0.25;  
     nus[5]=0.30; nus[6]=0.35; nus[7]=0.40; nus[8]=0.45; nus[9]=0.50;  
     nus[10]=0.55; nus[11]=0.60; nus[12]=0.65; nus[13]=0.70; nus[14]=0.75; 
     nus[15]=0.80; nus[16]=0.85; nus[17]=0.90; nus[18]=0.95; nus[19]=1.0; 

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
     int row, col, l, i, j, r1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     rT1 = *rT;

     double *ov, *o1, *XB1; 
     o1 = (double *) malloc((size_t)((row*col)*sizeof(double)));
     ov = (double *) malloc((size_t)((row*col)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((row*col)*sizeof(double)));

     int *T1; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     for(i=0; i<r1; i++){
          T1[i] = T[i];
     }

     double u;
     u = 0.0;
     for(l=0; l < r1; l++) {
         for(i=0; i < T1[l]; i++) {                                
             extract_alt_uneqT(l, i, n, r, T, rT, o, o1);
             extract_alt_uneqT(l, i, n, r, T, rT, XB, XB1);  
//             extract_alt2(l, i, n, rT, T, o, o1);
//             extract_alt2(l, i, n, rT, T, XB, XB1);
             for(j=0; j < row; j++) {
                 ov[j]=o1[j]-XB1[j];
             }       
             u += xTay2(ov, Qeta, ov, row);
         }
     }
     
     free(T1);
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


// Discrete sampling for nu
// for sp tp models
void nu_gp_DIS_sptp(int *cov, double *Qeta1, double *det1, double *phi,  
     double *nu1, int *n, int *r, int *T, int *rT, int *N, 
     double *d,  double *sig2eta, double *XB, double *o, int *constant, 
     double *nup)
{
    
     int row, col, i;
     row = *n;
     col = *constant;
//     r1 = *r;
//     T1 = *T;
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

/*
     nuk=20;
     double *nus;
     nus = (double *) malloc((size_t)((nuk)*sizeof(double)));             
     nus[0]=0.05; nus[1]=0.10; nus[2]=0.15; nus[3]=0.20; nus[4]=0.25;  
     nus[5]=0.30; nus[6]=0.35; nus[7]=0.40; nus[8]=0.45; nus[9]=0.50;  
     nus[10]=0.55; nus[11]=0.60; nus[12]=0.65; nus[13]=0.70; nus[14]=0.75; 
     nus[15]=0.80; nus[16]=0.85; nus[17]=0.90; nus[18]=0.95; nus[19]=1.0; 

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
        
        nudens_gp_sptp(Qeta, det, n, r, T, rT, N, XB, o, constant, out);
        
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

     nudens_gp_sptp(Qeta1, det1, n, r, T, rT, N, XB, o, constant, tr1);
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
// for sp and tp model
void nudens_gp_sptp(double *Qeta, double *det, int *n, int *r, int *T, int *rT, 
     int *N, double *XB, double *o, int *constant, double *out)
{
     int row, col, l, i, j, r1, T1, rT1;
     row = *n;
     col = *constant;
     r1 = *r;
     T1 = *T;
     rT1 =*rT;

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

/*



void sig_e_gp2 (int *N, double *shape, double *prior_b, double *o, double *z,
     int *constant, double *sig2eps)
{
     int i, N1, col;
     N1 = *N;
     col = *constant;
     
     double *oz, *out, *out1;
     oz = (double *) malloc((size_t)((N1)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));
     out1 = (double *) malloc((size_t)((col)*sizeof(double)));
               
     for(i=0; i< *N; i++){
          out1[0] = 0.005;
          out[0] = o[i];
          mvrnormal(constant, out, out1, constant, out);                              
          oz[i] = z[i] - out[0];
     }    
     MProd(oz, constant, N, oz, constant, out);     
     
     double u, v, b;
     u = 0.0;
     v = 0.0;
     b = 0.0;
     
     b = *prior_b;
     u = b + 0.5 * out[0];
     v = *shape;
     sig2eps[0] = rigammaa(v, u);

     free(oz); free(out); free(out1);
     return;
}


*/
