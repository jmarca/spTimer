//** File Name 'gibbs_gp.c' **//

#include "main_gp.h"
#include "covariance.h"
#include "common.h"
#include "math.h"
#include "mathematics.h"
#include "randgenerator.h"
#include "Print.h"


  
// The programme for GIBBS SAMPLING with XB and missing values
// with all summary values (mean, variance/sd, low2.5, up97.5)
// also the predictions into another sites
// output into the txt files
void GIBBS_sumpred_txt_gp(int *aggtype, double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *N, int *report,
     int *cov, int *spdecay, double *shape_e, double *shape_eta,  
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
     
     int its1, col, i, j, n1, r1, T1, p1, N1, nr, rep1, nsite1, brin, trans1;
     its1 = *its;
     col = *constant;
     n1 = *n;
     r1 = *r;
     T1 = *T;
     p1 = *p;
     N1 = *N;
     nr = n1 * r1;
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
     pr_mn = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double))); 
     pr_var = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double))); 
     for(j=0; j<nsite1*r1*T1; j++){
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
     zp = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double)));      
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

     JOINT_gp(n, T, r, rT, p, N, cov, spdecay, shape_e, shape_eta,
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
     for(j=0; j<(nsite1*r1*T1); j++){
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
         annual_aggregate(aggtype, nsite, r, T, zp, anf);
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
     for(j=0; j < nsite1*r1*T1; j++){
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
     int *cov, int *spdecay, double *shape_e, double *shape_eta,   
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
     
     int its1, brin, col, i, j, n1, r1, T1, p1, N1, nr, rep1;
     double *phip, *sig_ep, *sig_etap, *betap, *op;
     double *phi1, *sig_e1, *sig_eta1, *beta1, *o1;
     double *z1, *oo, *ot, *acc;
    
     its1 = *its;
     brin = *burnin;
     col = *constant;
     n1 = *n;
     r1 = *r;
     T1 = *T;
     p1 = *p;
     N1 = *N;
     nr = n1 * r1;
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

     JOINT_gp(n, T, r, rT, p, N, cov, spdecay, shape_e, shape_eta,
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


/////////////////////////// THE END ///////////////////////////////////////////
