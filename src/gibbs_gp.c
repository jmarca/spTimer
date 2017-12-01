//** File Name 'gibbs_gp.c' **//

#include "main_gp.h"
#include "covariance.h"
#include "common.h"
#include "math.h"
#include "mathematics.h"
#include "randgenerator.h"
//#include "Print.h"


  
// The programme for GIBBS SAMPLING with XB and missing values
// with all summary values (mean, variance/sd, low2.5, up97.5)
// also the predictions into another sites
// output into the txt files
void GIBBS_sumpred_txt_gp(int *aggtype, double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *N, int *report,
     int *cov, int *spdecay, double *shape_e, double *shape_eta,  
     double *phi_a, double *phi_b,     
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik,
     double *d, double *sig_e, double *sig_eta, 
     double *beta, double *X, double *z, double *o, int *constant, 
     int *nsite, int *valN, double *d12, double *valX, int *transform, 
     double *accept_f, double *gof, double *penalty)    
{
//     unsigned iseed = 44;
//     srand(iseed); 
     
     int its1, col, i, j, r1, rT1, p1, N1, rep1, nsite1, brin, trans1;
     its1 = *its;
     col = *constant;
//     n1 = *n;
     r1 = *r;
     rT1 = *rT;
     p1 = *p;
     N1 = *N;
//     nr = n1 * r1;
     rep1 = *report;
     nsite1 = *nsite;
     brin = *burnin;
     trans1 = *transform;

     double *phip, *sig_ep, *sig_etap, *betap;
     double *op;
     double *phi1, *sig_e1, *sig_eta1, *beta1;
     double *o1;
     double *oo, *ot, *acc;

     double accept1, *mn_rep, *var_rep;
     accept1 = 0.0;
     mn_rep = (double *) malloc((size_t)((N1)*sizeof(double)));
     var_rep = (double *) malloc((size_t)((N1)*sizeof(double)));
     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      
         
     double *pr_mn, *pr_var;
     pr_mn = (double *) malloc((size_t)((nsite1*rT1)*sizeof(double))); 
     pr_var = (double *) malloc((size_t)((nsite1*rT1)*sizeof(double))); 
     for(j=0; j<nsite1*rT1; j++){
          pr_mn[j] = 0.0;
          pr_var[j] = 0.0;
     }         

     phip = (double *) malloc((size_t)((col)*sizeof(double)));          
     sig_ep = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_etap = (double *) malloc((size_t)((col)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1)*sizeof(double)));          
     op = (double *) malloc((size_t)((N1)*sizeof(double)));

     phi1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1)*sizeof(double)));     
     o1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     
     oo = (double *) malloc((size_t)((col)*sizeof(double)));
     ot = (double *) malloc((size_t)((col)*sizeof(double)));
     acc = (double *) malloc((size_t)((col)*sizeof(double)));

     double *zp, *anf;
     zp = (double *) malloc((size_t)((nsite1*rT1)*sizeof(double)));      
     anf = (double *) malloc((size_t)((nsite1*r1)*sizeof(double)));      

     double *nu, *nup;
     nu = (double *) malloc((size_t)((col)*sizeof(double)));
     nup = (double *) malloc((size_t)((col)*sizeof(double)));     
     nu[0] = 0.5;
              
       ext_sige(phi, phi1);
       ext_sige(sig_e, sig_e1);
       ext_sigeta(sig_eta, sig_eta1);
       ext_beta(p, beta, beta1);
       ext_o(N, o, o1);

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=o1[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z[j] =  ot[0];
          }
          else {
          z[j] = z[j];
          }     
     }

     FILE *parafile;
     parafile = fopen("OutGP_Values_Parameter.txt", "w");
     FILE *zpfile;
     zpfile = fopen("OutGP_Stats_FittedValue.txt", "w");
     FILE *predfile;
     predfile = fopen("OutGP_Values_Prediction.txt", "w");
     FILE *predfilestat;
     predfilestat = fopen("OutGP_Stats_PredValue.txt", "w");

     int type1;
     type1= *aggtype;

     FILE *textan;
     // none
     if(type1==0){
       textan = fopen("OutGP_NONE.txt", "w");
     }
     // annual average value
     if(type1==1){
       textan = fopen("OutGP_Annual_Average_Prediction.txt", "w");
     }
     // annual 4th highest value
     if(type1==2){
       textan = fopen("OutGP_Annual_4th_Highest_Prediction.txt", "w");
     }
     // annual w126 option
     if(type1==3){
       textan = fopen("OutGP_Annual_w126_Prediction.txt", "w");
     }

     GetRNGstate();                 
     for(i=0; i < its1; i++) {

     JOINT_gp(n, T, r, rT, p, N, cov, spdecay, shape_e, shape_eta, phi_a, phi_b,
     prior_a, prior_b, prior_mubeta, prior_sigbeta, prior_omu, prior_osig,
     phi1, tau, phis, phik, nu, d, sig_e1, sig_eta1, beta1, X, z, o1, constant,
     phip, acc, nup, sig_ep, sig_etap, betap, op);

     z_pr_gp(cov, nsite, n, r, rT, T, p, N, valN, d, d12, phip, nup,
     sig_ep, sig_etap, betap, X, valX, op, constant, zp);
 
     accept1 += acc[0];

     for(j=0; j < p1; j++){
     fprintf(parafile, "%f ", betap[j]);
     }         
     fprintf(parafile, "%f ", sig_ep[0]);
     fprintf(parafile, "%f ", sig_etap[0]);
     fprintf(parafile, "%f ", phip[0]);
     if(cov[0]==4){
      fprintf(parafile, "%f ", nup[0]);
     }     
     fprintf(parafile, "\n");


// for pmcc, fitted    
     for(j=0; j < N1; j++){
         if(i >= brin){  
          oo[0] = op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          mn_rep[j] += ot[0];
          var_rep[j] += ot[0]*ot[0];
         } 
     }

// prediction samples
     for(j=0; j<(nsite1*rT1); j++){
         if(trans1 == 0){
         if(i >= brin){  
           zp[j] = zp[j];
           fprintf(predfile, "%f ", zp[j]);
          pr_mn[j] += zp[j];
          pr_var[j] += zp[j]*zp[j];
         } 
         }
         if(trans1 == 1){
         if(i >= brin){  
           zp[j] = zp[j]*zp[j];
           fprintf(predfile, "%f ", zp[j]);
          pr_mn[j] += zp[j];
          pr_var[j] += zp[j]*zp[j];
         } 
         }
         if(trans1 == 2){
         if(i >= brin){  
           zp[j] = exp(zp[j]);
           fprintf(predfile, "%f ", zp[j]);
          pr_mn[j] += zp[j];
          pr_var[j] += zp[j]*zp[j];
         } 
         }                                
     }
     fprintf(predfile, "\n");

     if(cov[0]==4){
     GP_para_printRnu(i, its1, rep1, p1, accept1, phip, nup, sig_ep, sig_etap, betap);  
     }                   
     else {              
     GP_para_printR (i, its1, rep1, p1, accept1, phip, sig_ep, sig_etap, betap);  
     }

     if(i >= brin){
         annual_aggregate_uneqT(aggtype, nsite, r, T, rT, zp, anf);
//         annual_aggregate(aggtype, nsite, r, T, zp, anf);
         for(j=0; j<(nsite1*r1); j++){
           fprintf(textan, "%f ", anf[j]);
         }
         fprintf(textan, "\n");
     } // end of loop i >= brin

       ext_sige(phip, phi1);
       ext_sige(nup, nu);        
       ext_sige(sig_ep, sig_e1);
       ext_sige(sig_etap, sig_eta1);
       ext_beta(p, betap, beta1);              
       
// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z[j] =  ot[0];
          }
          else {
          z[j] = z[j];
          }     
     }

     } // end of iteration loop
     PutRNGstate();

     fclose(parafile);
     fclose(predfile);
     fclose(textan);

     free(phip); free(nu); free(nup); free(sig_ep); free(sig_etap); free(betap); 
     free(op); free(phi1); free(sig_e1); 
     free(sig_eta1); free(beta1);  
     free(o1); free(oo); free(ot); free(acc);    

     free(zp); free(anf); 

     accept_f[0] = accept1;

     int iit;
     iit = 0;
     iit = its1 - brin;

     double pen, go;
     pen = 0.0;
     go =0.0;

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit;
          var_rep[j] = var_rep[j]/iit;
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          fprintf(zpfile, "%f , %f \n", mn_rep[j], sqrt(var_rep[j]));
     }
     fclose(zpfile);

// pmcc          
     for(j=0; j < N1; j++){
         if (flag[j] == 1.0){
          mn_rep[j] = 0.0;
          var_rep[j] = 0.0;
         }
         else{
          mn_rep[j] = mn_rep[j];
          var_rep[j] = var_rep[j];
          mn_rep[j] = (mn_rep[j] - z[j])*(mn_rep[j] - z[j]);          
         }
         pen += var_rep[j]; 
         go += mn_rep[j];
     }
     free(mn_rep); free(var_rep);

     penalty[0] = pen;
     gof[0] = go;       


// predicted mean and sd
     for(j=0; j < nsite1*rT1; j++){
          pr_mn[j] = pr_mn[j]/iit;
          pr_var[j] = pr_var[j]/iit;
          pr_var[j] = pr_var[j] - pr_mn[j]*pr_mn[j];
          fprintf(predfilestat, "%f , %f \n", pr_mn[j], sqrt(pr_var[j]));
     }
     fclose(predfilestat);
     free(pr_mn); free(pr_var);


//     Rprintf("\n---------------------------------------------------------\n");
     
     return;
}




// The programme for GIBBS SAMPLING with XB and missing values
void GIBBS_gp(double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *N, int *report,
     int *cov, int *spdecay, int *ft, double *shape_e, double *shape_eta,   
     double *phi_a, double *phi_b,
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik,
     double *d, double *sig_e, double *sig_eta, 
     double *beta, double *X, double *z, double *o, int *constant, 
     double *phipf, double *accept, double *nupf, double *sig_epf, 
     double *sig_etapf, double *betapf, double *opf, double *zlt_mean_sd, 
     double *gof, double *penalty)
{
//     unsigned iseed = 44;
//     srand(iseed); 
     
     int its1, brin, col, i, j, p1, N1, rep1;
     double *phip, *sig_ep, *sig_etap, *betap, *op;
     double *phi1, *sig_e1, *sig_eta1, *beta1, *o1;
     double *z1, *oo, *ot, *acc;
    
     its1 = *its;
     brin = *burnin;
     col = *constant;
//     n1 = *n;
//     r1 = *r;
     p1 = *p;
     N1 = *N;
//     nr = n1 * r1;
     rep1 = *report;
          
     double accept1, mn_rep[N1], var_rep[N1];
     accept1 = 0.0;
     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      


     phip = (double *) malloc((size_t)((col)*sizeof(double)));          
     sig_ep = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_etap = (double *) malloc((size_t)((col)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1)*sizeof(double)));          
     op = (double *) malloc((size_t)((N1)*sizeof(double)));

     phi1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1)*sizeof(double)));     
     o1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     
     z1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     oo = (double *) malloc((size_t)((col)*sizeof(double)));
     ot = (double *) malloc((size_t)((col)*sizeof(double)));
     acc = (double *) malloc((size_t)((col)*sizeof(double)));

     double *nu, *nup;
     nu = (double *) malloc((size_t)((col)*sizeof(double)));
     nup = (double *) malloc((size_t)((col)*sizeof(double)));     
     nu[0] = 0.5;
                    
       ext_sige(phi, phi1);
       ext_sige(sig_e, sig_e1);
       ext_sigeta(sig_eta, sig_eta1);
       ext_beta(p, beta, beta1);
       ext_o(N, o, o1);
       ext_o(N, z, z1);

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=o1[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z1[j] =  ot[0];
          }
          else {
          z1[j] = z1[j];
          }     
     }


          
     GetRNGstate();       
     for(i=0; i < its1; i++) {

     JOINT_gp(n, T, r, rT, p, N, cov, spdecay, shape_e, shape_eta, phi_a, phi_b,
     prior_a, prior_b, prior_mubeta, prior_sigbeta, prior_omu, prior_osig,
     phi1, tau, phis, phik, nu, d, sig_e1, sig_eta1, beta1, X, z1, o1, constant,
     phip, acc, nup, sig_ep, sig_etap, betap, op);

     accept1 += acc[0];

     phipf[i] = phip[0];
     nupf[i] = nup[0];
     sig_epf[i] = sig_ep[0];
     sig_etapf[i] = sig_etap[0];
     for(j=0; j < p1; j++){
         betapf[j+i*p1] = betap[j];
     }         
     for(j=0; j < N1; j++) {
         opf[j+i*N1] = op[j];     
     }

       ext_sige(phip, phi1);
       ext_sige(nup, nu);       
       ext_sige(sig_ep, sig_e1);
       ext_sige(sig_etap, sig_eta1);
       ext_beta(p, betap, beta1);              
//       ext_o(N, op, o1);
       
// for pmcc   
     for(j=0; j < N1; j++){
          if(i >= brin){    
          oo[0] = op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          // Three options: ft: 0=NONE, 1=SQRT, 2=LOG
		  if(ft[0]==0){
              mn_rep[j] += ot[0];
              var_rep[j] += ot[0]*ot[0];
		  }
		  else{
			  if(ft[0]==1){
              mn_rep[j] += ot[0]*ot[0];
              var_rep[j] += ot[0]*ot[0]*ot[0]*ot[0];
			  }
			  else{
              mn_rep[j] += exp(ot[0]);
              var_rep[j] += exp(ot[0])*exp(ot[0]);
			  }
		  }
          }
     }

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z1[j] =  ot[0];
          }
          else {
          z1[j] = z1[j];
          }     
     }

     if(cov[0]==4){
     GP_para_printRnu(i, its1, rep1, p1, accept1, phip, nup, sig_ep, sig_etap, betap);  
     }                   
     else {              
     GP_para_printR (i, its1, rep1, p1, accept1, phip, sig_ep, sig_etap, betap);  
     }

     } // end of iteration loop
     PutRNGstate();     
     
     accept[0] = accept1;

     double pen, go;
     pen = 0;
     go =0;
     int iit;
     iit = 0;
     iit = its1-brin;

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit;
          var_rep[j] = var_rep[j]/iit;
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          zlt_mean_sd[j] = mn_rep[j];
          zlt_mean_sd[j+N1] = sqrt(var_rep[j]);
     }
// pmcc          
     for(j=0; j < N1; j++){
         if (flag[j] == 1.0){
          mn_rep[j] = 0.0;
          var_rep[j] = 0.0;
         }
         else{
          mn_rep[j] = mn_rep[j];
          var_rep[j] = var_rep[j];
          mn_rep[j] = (mn_rep[j] - z1[j])*(mn_rep[j] - z1[j]);          
         }
         pen += var_rep[j]; 
         go += mn_rep[j];
     }

     gof[0] = go;
     penalty[0] = pen;
          
     free(phip); free(nu); free(nup); free(sig_ep); free(sig_etap); free(betap); 
     free(op); free(phi1); free(sig_e1); free(sig_eta1); free(beta1); free(o1); 
     free(z1); free(oo); free(ot); free(acc);
         
     return;
}


/*
// The programme for GIBBS SAMPLING with XB and missing values
// spatially varying covariates
void GIBBSsp_gp(int *intercept, double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *q, int *N, int *report,
     int *cov, int *spdecay, double *shape_e, double *shape_eta, double *shape_beta,  
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik,
     double *d, double *sig_e, double *sig_eta, double *sig_beta,
     double *beta, double *betas, double *X, double *Xsp, double *z, double *o, 
     int *constant, 
     double *phipf, double *accept, double *nupf, double *sig_epf, 
     double *sig_etapf, double *sig_betaspf, double *betapf, double *betaspf, 
     double *opf, double *zlt_mean_sd, double *gof, double *penalty)
{
     
     int its1, brin, col, i, j, n1, r1, p1, q1, N1, rep1;
     double *phip, *sig_ep, *sig_etap, *sig_betap, *betap, *betasp, *op;
     double *phi1, *sig_e1, *sig_eta1, * sig_beta1, *beta1, *o1;
     double *z1, *oo, *ot, *acc;
    
     its1 = *its;
     brin = *burnin;
     col = *constant;
     n1 = *n;
     r1 = *r;
//     T1 = *T;
     p1 = *p;
     q1 = *q;
     N1 = *N;
//     nr = n1 * r1;
     rep1 = *report;
          
     double accept1, mn_rep[N1], var_rep[N1];
     accept1 = 0.0;
     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      


     phip = (double *) malloc((size_t)((col)*sizeof(double)));          
     sig_ep = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_etap = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_betap = (double *) malloc((size_t)((col)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1)*sizeof(double))); 
     betasp = (double *) malloc((size_t)((n1*q1)*sizeof(double)));               
     op = (double *) malloc((size_t)((N1)*sizeof(double)));

     phi1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_beta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1)*sizeof(double)));  
     o1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     
     z1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     oo = (double *) malloc((size_t)((col)*sizeof(double)));
     ot = (double *) malloc((size_t)((col)*sizeof(double)));
     acc = (double *) malloc((size_t)((col)*sizeof(double)));

     double *nu, *nup;
     nu = (double *) malloc((size_t)((col)*sizeof(double)));
     nup = (double *) malloc((size_t)((col)*sizeof(double)));     
     nu[0] = 0.5;
                    
       ext_sige(phi, phi1);
       ext_sige(sig_e, sig_e1);
       ext_sigeta(sig_eta, sig_eta1);
       ext_sigeta(sig_beta, sig_beta1);
       ext_beta(p, beta, beta1);
       ext_o(N, o, o1);
       ext_o(N, z, z1);

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=o1[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z1[j] =  ot[0];
          }
          else {
          z1[j] = z1[j];
          }     
     }
          
     GetRNGstate();       
     for(i=0; i < its1; i++) {

     JOINTsp_gp(intercept, n, T, r, rT, p, q, N, cov, spdecay, shape_e, 
     shape_eta, shape_beta, prior_a, prior_b, prior_mubeta, prior_sigbeta, 
     prior_omu, prior_osig, phi1, tau, phis, phik, nu, d, sig_e1, sig_eta1, 
     sig_beta1, beta1, betas, X, Xsp, z1, o1, constant, phip, acc, nup, 
     sig_ep, sig_etap, sig_betap, betap, betasp, op);
     
     accept1 += acc[0];

     phipf[i] = phip[0];
     nupf[i] = nup[0];
     sig_epf[i] = sig_ep[0];
     sig_etapf[i] = sig_etap[0];
     sig_betaspf[i] = sig_betap[0];
     for(j=0; j < p1; j++){
         betapf[j+i*p1] = betap[j];
     }         
     for(j=0; j < n1*q1; j++){
         betaspf[j+i*n1*q1] = betasp[j];
     }         
     for(j=0; j < N1; j++) {
         opf[j+i*N1] = op[j];     
     }

       ext_sige(phip, phi1);
       ext_sige(nup, nu);       
       ext_sige(sig_ep, sig_e1);
       ext_sige(sig_etap, sig_eta1);
       ext_sige(sig_betap, sig_beta1);
       ext_beta(p, betap, beta1);   
//       for(j=0; j < n1; j++){
//         betas[j] = betasp[j];
//       }         
           
              
// for pmcc   
     for(j=0; j < N1; j++){
          if(i >= brin){    
          oo[0] = op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          mn_rep[j] += ot[0];
          var_rep[j] += ot[0]*ot[0];
          }
     }

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z1[j] =  ot[0];
          }
          else {
          z1[j] = z1[j];
          }     
     }

     if(cov[0]==4){
     GPsp_para_printRnu(i, its1, rep1, p1, accept1, phip, nup, sig_ep, sig_etap, sig_betap, betap); 
     }                   
     else {              
     GPsp_para_printR (i, its1, rep1, p1, accept1, phip, sig_ep, sig_etap, sig_betap, betap);  
     }

     } // end of iteration loop
     PutRNGstate();     

//      Rprintf("---------------------------------------------------------------\n");
//      Rprintf(" ## Model used spatially varying parameters \n");
//      Rprintf(" ## Spatially varying parameters are omitted in the display ");
//      Rprintf("\n---------------------------------------------------------------\n");
     
     accept[0] = accept1;

     double pen, go;
     pen = 0;
     go =0;
     int iit;
     iit = 0;
     iit = its1-brin;

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit;
          var_rep[j] = var_rep[j]/iit;
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          zlt_mean_sd[j] = mn_rep[j];
          zlt_mean_sd[j+N1] = sqrt(var_rep[j]);
     }

// pmcc          
     for(j=0; j < N1; j++){
         if (flag[j] == 1.0){
          mn_rep[j] = 0.0;
          var_rep[j] = 0.0;
         }
         else{
          mn_rep[j] = mn_rep[j];
          var_rep[j] = var_rep[j];
          mn_rep[j] = (mn_rep[j] - z1[j])*(mn_rep[j] - z1[j]);          
         }
         pen += var_rep[j]; 
         go += mn_rep[j];
     }

     gof[0] = go;
     penalty[0] = pen;
          
     free(phip); free(nu); free(nup); free(sig_ep); free(sig_etap); free(sig_betap); free(betap); 
     free(betasp); free(op); free(phi1); free(sig_e1); free(sig_eta1); free(sig_beta1); 
     free(beta1); free(o1); free(z1); free(oo); free(ot); free(acc);

     return;
}
*/

/*
// The programme for GIBBS SAMPLING with XB and missing values
// temporally varying covariates
void GIBBStp_gp(int *intercept, double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *u, int *N, int *report,
     int *cov, int *spdecay, int *rhocheck,
     double *shape_e, double *shape_eta, double *shape_del, double *shape_0,  
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik,
     double *d, double *sig_e, double *sig_eta, double *sig_delta, double *sig_0,
     double *beta, double *betat, double *rho, double *X, double *Xtp, 
     double *z, double *o, int *constant, 
     double *phipf, double *accept, double *nupf, double *sig_epf, 
     double *sig_etapf, double *sig_deltapf, double *sig_0pf, double *rhopf, 
     double *betapf, double *betat0pf, double *betatpf, 
     double *opf, double *zlt_mean_sd, double *gof, double *penalty)
{
     
     int its1, brin, col, i, j, n1, r1, T1, p1, u1, N1, rep1;
     double *phip, *sig_ep, *sig_etap, *sig_deltap, *sig_0p, *rhop, *betap, * betat0p, *betatp, *op;
     double *phi1, *sig_e1, *sig_eta1, *sig_delta1, *sig_01, *rho1, *beta1, *o1;
     double *z1, *oo, *ot, *acc;
    
     its1 = *its;
     brin = *burnin;
     col = *constant;
     n1 = *n;
     r1 = *r;
     T1 = *T;
     p1 = *p;
     u1 = *u;
     N1 = *N;
//     nr = n1 * r1;
     rep1 = *report;
          
     double accept1, mn_rep[N1], var_rep[N1];
     accept1 = 0.0;
     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      


     phip = (double *) malloc((size_t)((col)*sizeof(double)));          
     sig_ep = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_etap = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_deltap = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_0p = (double *) malloc((size_t)((col)*sizeof(double)));
     rhop = (double *) malloc((size_t)((u1)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1)*sizeof(double))); 
     betat0p = (double *) malloc((size_t)((T1)*sizeof(double)));               
     betatp = (double *) malloc((size_t)((T1*u1)*sizeof(double)));               
     op = (double *) malloc((size_t)((N1)*sizeof(double)));

     phi1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_delta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_01 = (double *) malloc((size_t)((col)*sizeof(double)));
     rho1 = (double *) malloc((size_t)((u1)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1)*sizeof(double)));  
     o1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     
     z1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     oo = (double *) malloc((size_t)((col)*sizeof(double)));
     ot = (double *) malloc((size_t)((col)*sizeof(double)));
     acc = (double *) malloc((size_t)((col)*sizeof(double)));

     double *nu, *nup;
     nu = (double *) malloc((size_t)((col)*sizeof(double)));
     nup = (double *) malloc((size_t)((col)*sizeof(double)));     
     nu[0] = 0.5;
                    
       ext_sige(phi, phi1);
       ext_sige(sig_e, sig_e1);
       ext_sigeta(sig_eta, sig_eta1);
       ext_sigeta(sig_delta, sig_delta1);
       ext_sigeta(sig_0, sig_01);       
       ext_beta(p, beta, beta1);
       ext_beta(u, rho, rho1);
       ext_o(N, o, o1);
       ext_o(N, z, z1);

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=o1[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z1[j] =  ot[0];
          }
          else {
          z1[j] = z1[j];
          }     
     }

         
     GetRNGstate();       
     for(i=0; i < its1; i++) {

     JOINTtp_gp(intercept, n, T, r, rT, p, u, N, cov, spdecay, rhocheck,
     shape_e, shape_eta, shape_del, shape_0,
     prior_a, prior_b, prior_mubeta, prior_sigbeta, prior_omu, prior_osig, 
     phi1, tau, phis, phik, nu, d, sig_e1, sig_eta1, sig_delta1, sig_01,
     beta1, betat, rho1, X, Xtp, z1, o1, constant, phip, acc, nup, sig_ep, 
     sig_etap, sig_deltap, sig_0p, rhop, betap, betat0p, betatp, op);

     accept1 += acc[0];

     phipf[i] = phip[0];
     nupf[i] = nup[0];
     sig_epf[i] = sig_ep[0];
     sig_etapf[i] = sig_etap[0];
     sig_deltapf[i] = sig_deltap[0];
     sig_0pf[i] = sig_0p[0];
     for(j=0; j < u1; j++){
         rhopf[j+i*u1] = rhop[j];
     }
     for(j=0; j < p1; j++){
         betapf[j+i*p1] = betap[j];
     }
     for(j=0; j < u1; j++){
         betat0pf[j+i*u1] = betat0p[j];
     }
     for(j=0; j < u1*T1; j++){
         betatpf[j+i*u1*T1] = betatp[j];
     }         
     for(j=0; j < N1; j++) {
         opf[j+i*N1] = op[j];     
     }

       ext_sige(phip, phi1);
       ext_sige(nup, nu);       
       ext_sige(sig_ep, sig_e1);
       ext_sige(sig_etap, sig_eta1);
       ext_sige(sig_deltap, sig_delta1);
       ext_beta(u, rhop, rho1);              
//       for(j=0; j < u1; j++){
//         rho1[j] = rhop[j];
//         Rprintf("  rhop: %4.4f,", rhop[j]);      
//       }         
//         Rprintf("\n");      
       ext_beta(p, betap, beta1);  
       for(j=0; j < u1*T1; j++){
         betat[j] = betatp[j];
       }         
                   

              
// for pmcc   
     for(j=0; j < N1; j++){
          if(i >= brin){    
          oo[0] = op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          mn_rep[j] += ot[0];
          var_rep[j] += ot[0]*ot[0];
          }
     }

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z1[j] =  ot[0];
          }
          else {
          z1[j] = z1[j];
          }     
     }

     if(cov[0]==4){
     GPtp_para_printRnu(i, its1, rep1, p1, u1, accept1, phip, nup, sig_ep, sig_etap, sig_deltap, sig_0p, rhop, betap); 
     }                   
     else {              
     GPtp_para_printR (i, its1, rep1, p1, u1, accept1, phip, sig_ep, sig_etap, sig_deltap, sig_0p, rhop, betap);  
     }

     } // end of iteration loop
     PutRNGstate();     

     accept[0] = accept1;

     double pen, go;
     pen = 0;
     go =0;
     int iit;
     iit = 0;
     iit = its1-brin;

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit;
          var_rep[j] = var_rep[j]/iit;
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          zlt_mean_sd[j] = mn_rep[j];
          zlt_mean_sd[j+N1] = sqrt(var_rep[j]);
     }

// pmcc          
     for(j=0; j < N1; j++){
         if (flag[j] == 1.0){
          mn_rep[j] = 0.0;
          var_rep[j] = 0.0;
         }
         else{
          mn_rep[j] = mn_rep[j];
          var_rep[j] = var_rep[j];
          mn_rep[j] = (mn_rep[j] - z1[j])*(mn_rep[j] - z1[j]);          
         }
         pen += var_rep[j]; 
         go += mn_rep[j];
     }

     gof[0] = go;
     penalty[0] = pen;
          
     free(nu); free(nup); 
     free(phip); free(sig_ep); free(sig_etap); free(sig_deltap); free(sig_0p);
     free(rhop); free(betap); free(betat0p); free(betatp); free(op);
     free(phi1); free(sig_e1); free(sig_eta1); free(sig_delta1); free(sig_01);
     free(rho1); free(beta1); free(o1); free(z1); free(oo); free(ot); free(acc);


     return;
}
*/


/*
// The programme for GIBBS SAMPLING with XB and missing values
// temporally varying covariates
void GIBBSsptp_gp(int *intercept, double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *q, int *u, int *N, int *report,
     int *cov, int *spdecay, int *rhocheck,
     double *shape_e, double *shape_eta, double *shape_beta, double *shape_del, double *shape_0,  
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik,
     double *d, double *sig_e, double *sig_eta, double *sig_beta, double *sig_delta, double *sig_0,
     double *beta, double *betas, double *betat, double *rho, double *X, double *Xsp, double *Xtp, 
     double *z, double *o, int *constant, 
     double *phipf, double *accept, double *nupf, double *sig_epf, 
     double *sig_etapf, double *sig_betapf, double *sig_deltapf, double *sig_0pf, double *rhopf, 
     double *betapf, double *betaspf, double *betat0pf, double *betatpf, 
     double *opf, double *zlt_mean_sd, double *gof, double *penalty)
{
     
     int its1, brin, col, i, j, n1, r1, T1, p1, q1, u1, N1, rep1;
     double *phip, *sig_ep, *sig_etap, *sig_betap, *sig_deltap, *sig_0p, *betap, *betasp, * betat0p, *betatp, *op, *rhop; 
     double *phi1, *sig_e1, *sig_eta1, *sig_beta1, *sig_delta1, *sig_01, *beta1, *o1, *rho1;
     double *z1, *oo, *ot, *acc;
    
     its1 = *its;
     brin = *burnin;
     col = *constant;
     n1 = *n;
     r1 = *r;
     T1 = *T;
     p1 = *p;
     q1 =*q;
     u1 = *u;
     N1 = *N;
//     nr = n1 * r1;
     rep1 = *report;
          
     double accept1, mn_rep[N1], var_rep[N1];
     accept1 = 0.0;
     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      


     phip = (double *) malloc((size_t)((col)*sizeof(double)));          
     sig_ep = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_etap = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_betap = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_deltap = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_0p = (double *) malloc((size_t)((col)*sizeof(double)));
     rhop = (double *) malloc((size_t)((u1)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1)*sizeof(double))); 
     betasp = (double *) malloc((size_t)((n1*q1)*sizeof(double)));  
     betat0p = (double *) malloc((size_t)((u1)*sizeof(double)));               
     betatp = (double *) malloc((size_t)((T1*u1)*sizeof(double)));               
     op = (double *) malloc((size_t)((N1)*sizeof(double)));

     phi1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_beta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_delta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_01 = (double *) malloc((size_t)((col)*sizeof(double)));
     rho1 = (double *) malloc((size_t)((u1)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1)*sizeof(double)));  
     o1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     
     z1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     oo = (double *) malloc((size_t)((col)*sizeof(double)));
     ot = (double *) malloc((size_t)((col)*sizeof(double)));
     acc = (double *) malloc((size_t)((col)*sizeof(double)));

     double *nu, *nup;
     nu = (double *) malloc((size_t)((col)*sizeof(double)));
     nup = (double *) malloc((size_t)((col)*sizeof(double)));     
     nu[0] = 0.5;
                    
       ext_sige(phi, phi1);
       ext_sige(sig_e, sig_e1);
       ext_sige(sig_eta, sig_eta1);
       ext_sige(sig_beta, sig_beta1);
       ext_sige(sig_delta, sig_delta1);
       ext_sige(sig_0, sig_01);       
       ext_beta(p, beta, beta1);
       ext_beta(u, rho, rho1);
       ext_o(N, o, o1);
       ext_o(N, z, z1);

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=o1[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z1[j] =  ot[0];
          }
          else {
          z1[j] = z1[j];
          }     
     }
          
     GetRNGstate();       
     for(i=0; i < its1; i++) {

//     Rprintf("  sig2eps1: %4.4f, sig2eta1: %4.4f, sig2beta1: %4.4f, sig2delta1: %4.4f, sig2op1: %4.4f,\n", sig_e1[0], sig_eta1[0], sig_beta1[0], sig_delta1[0], sig_01[0]);      

     JOINTsptp_gp(intercept, n, T, r, rT, p, q, u, N, cov, spdecay, rhocheck,
     shape_e, shape_eta, shape_beta, shape_del, shape_0,
     prior_a, prior_b, prior_mubeta, prior_sigbeta, prior_omu, prior_osig, 
     phi1, tau, phis, phik, nu, d, sig_e1, sig_eta1, sig_beta1, sig_delta1, sig_01,
     beta1, betas, betat, rho, X, Xsp, Xtp, z1, o1, constant,
     phip, acc, nup, sig_ep, sig_etap, sig_betap, sig_deltap, sig_0p, rhop, 
     betap, betasp, betat0p, betatp, op);

//    Rprintf("  sig2eps: %4.4f, sig2eta: %4.4f, sig2beta: %4.4f, sig2delta: %4.4f, sig2op: %4.4f,\n", sig_ep[0], sig_etap[0], sig_betap[0], sig_deltap[0], sig_0p[0]);      

     accept1 += acc[0];

     phipf[i] = phip[0];
     nupf[i] = nup[0];
     sig_epf[i] = sig_ep[0];
     sig_etapf[i] = sig_etap[0];
     sig_betapf[i] = sig_betap[0];     
     sig_deltapf[i] = sig_deltap[0];
     sig_0pf[i] = sig_0p[0];
     for(j=0; j < u1; j++){
         rhopf[j+i*u1] = rhop[j];
     }
     for(j=0; j < p1; j++){
         betapf[j+i*p1] = betap[j];
     }
     for(j=0; j < n1*q1; j++){
         betaspf[j+i*n1*q1] = betasp[j];
     }         
     for(j=0; j < u1; j++){
         betat0pf[j+i*u1] = betat0p[j];
     }
     for(j=0; j < u1*T1; j++){
         betatpf[j+i*u1*T1] = betatp[j];
     }         
     for(j=0; j < N1; j++) {
         opf[j+i*N1] = op[j];     
     }

//     Rprintf("  sig2eps: %4.4f, sig2eta: %4.4f, sig2beta: %4.4f, sig2delta: %4.4f, sig2op: %4.4f,\n", sig_ep[0], sig_etap[0], sig_betap[0], sig_deltap[0], sig_0p[0]);      

       ext_sige(phip, phi1);
       ext_sige(nup, nu);       
       ext_sige(sig_ep, sig_e1);
       ext_sige(sig_etap, sig_eta1);
       ext_sige(sig_betap, sig_beta1);
       ext_sige(sig_deltap, sig_delta1);
       ext_beta(u, rhop, rho1);              
       ext_beta(p, betap, beta1);  
       for(j=0; j < u1*T1; j++){
         betat[j] = betatp[j];
       }         

              
// for pmcc   
     for(j=0; j < N1; j++){
          if(i >= brin){    
          oo[0] = op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          mn_rep[j] += ot[0];
          var_rep[j] += ot[0]*ot[0];
          }
     }

// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          oo[0]=op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          z1[j] =  ot[0];
          }
          else {
          z1[j] = z1[j];
          }     
     }

     if(cov[0]==4){
     GPsptp_para_printRnu(i, its1, rep1, p1, u1, accept1, phip, nup, sig_ep, sig_etap, 
     sig_betap, sig_deltap, sig_0p, rhop, betap); 
     }                   
     else {              
     GPsptp_para_printR (i, its1, rep1, p1, u1, accept1, phip, sig_ep, sig_etap, 
     sig_betap, sig_deltap, sig_0p, rhop, betap);  
     }

     } // end of iteration loop
     PutRNGstate();     

     accept[0] = accept1;

     double pen, go;
     pen = 0;
     go =0;
     int iit;
     iit = 0;
     iit = its1-brin;

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit;
          var_rep[j] = var_rep[j]/iit;
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          zlt_mean_sd[j] = mn_rep[j];
          zlt_mean_sd[j+N1] = sqrt(var_rep[j]);
     }

// pmcc          
     for(j=0; j < N1; j++){
         if (flag[j] == 1.0){
          mn_rep[j] = 0.0;
          var_rep[j] = 0.0;
         }
         else{
          mn_rep[j] = mn_rep[j];
          var_rep[j] = var_rep[j];
          mn_rep[j] = (mn_rep[j] - z1[j])*(mn_rep[j] - z1[j]);          
         }
         pen += var_rep[j]; 
         go += mn_rep[j];
     }

     gof[0] = go;
     penalty[0] = pen;
          
     free(nu); free(nup); 
     free(phip); free(sig_ep); free(sig_etap); free(sig_betap); free(sig_deltap); free(sig_0p);
     free(betap); free(betasp); free(betat0p); free(betatp); free(op);
     free(phi1); free(sig_e1); free(sig_eta1); free(sig_beta1); free(sig_delta1); free(sig_01);
     free(beta1); free(o1); free(z1); free(oo); free(ot); free(acc);

     free(rho1); free(rhop);

     return;
}
*/

/////////////////////////// THE END ///////////////////////////////////////////
