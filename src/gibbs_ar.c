//** File Name 'gibbs_ar.c' **//

#include "main_ar.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"
//#include "Print.h"

   
// The programme for GIBBS SAMPLING with XB and missing values
// with all summary values (mean, variance/sd, low2.5, up97.5)
// also the predictions into another sites
// output into the txt files
void GIBBS_sumpred_txt_ar(int *aggtype, double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *N, int *report,
     int *cov, int *spdecay, double *shape_e, double *shape_eta, double *shape_0,
     double *phi_a, double *phi_b,  
     double *prior_a, double *prior_b, double *prior_sig, double *phi, 
     double *tau, double *phis, int *phik, double *d, int *constant, 
     double *sig_e, double *sig_eta, double *sig_0, double *mu_l,  
     double *rho, double *beta, double *X, double *z, double *o, 
     int *nsite, int *predN, double *d12, double *predX, int *transform,
     double *accept_f, double *gof, double *penalty)    
{
//     unsigned iseed = 44;
//     srand(iseed); 
     
     int its1, col, i, j, r1, rT1, p1, N1, rep1, nsite1, brin, trans1;
     its1 = *its;
     col = *constant;
//     n1 = *n;
     r1 = *r;
     rT1 =*rT;
     p1 = *p;
     N1 = *N;
//     nr = n1 * r1;
     rep1 = *report;
     nsite1 = *nsite;
     brin = *burnin;
     trans1 = *transform;

     double *phip, *sig_ep, *sig_etap, *rhop, *betap;
     double *mu_lp, *sig_0p, *op;
     double *phi1, *sig_e1, *sig_eta1, *rho1, *beta1;
     double *mu_l1, *sig_01, *o1, *w;
     double *oo, *ot, *acc;

     double accept1, *mn_rep, *var_rep; //*tr_mn, *tr_var, 
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
     rhop = (double *) malloc((size_t)((col)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1)*sizeof(double)));          
     mu_lp = (double *) malloc((size_t)((r1)*sizeof(double)));    
     sig_0p = (double *) malloc((size_t)((r1)*sizeof(double))); 
     op = (double *) malloc((size_t)((N1)*sizeof(double)));

     phi1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     rho1 = (double *) malloc((size_t)((col)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1)*sizeof(double)));     
     mu_l1 = (double *) malloc((size_t)((r1)*sizeof(double)));    
     sig_01 = (double *) malloc((size_t)((r1)*sizeof(double))); 
     o1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     w = (double *) malloc((size_t)((N1)*sizeof(double)));
          
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
       ext_rho(rho, rho1);
       ext_beta(p, beta, beta1);
       ext_mul(r, mu_l, mu_l1);
       ext_sigl(r, sig_0, sig_01);
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
     parafile = fopen("OutAR_Values_Parameter.txt", "w");
     FILE *zpfile;
     zpfile = fopen("OutAR_Stats_FittedValue.txt", "w");
     FILE *predfile;
     predfile = fopen("OutAR_Values_Prediction.txt", "w");
     FILE *predfilestat;
     predfilestat = fopen("OutAR_Stats_PredValue.txt", "w");

     int type1;
     type1= *aggtype;
   
     FILE *textan;
     // none
     if(type1==0){
       textan = fopen("OutAR_NONE.txt", "w");
     }
     // annual average value
     if(type1==1){
       textan = fopen("OutAR_Annual_Average_Prediction.txt", "w");
     }
     // annual 4th highest value
     if(type1==2){
       textan = fopen("OutAR_Annual_4th_Highest_Prediction.txt", "w");
     }
     // annual w126 option
     if(type1==3){
       textan = fopen("OutAR_Annual_w126_Prediction.txt", "w");
     }

     GetRNGstate();            
     for(i=0; i < its1; i++) {

     JOINT_ar(n, T, r, rT, p, N, cov, spdecay, shape_e, shape_eta, shape_0,
     phi_a, phi_b,  
     prior_a, prior_b, prior_sig, phi1, tau, phis, phik, nu, d, constant, 
     sig_e1, sig_eta1, sig_01, mu_l1, rho1, beta1, X, z, o1, 
     phip, acc, nup, sig_ep, sig_etap, rhop, betap, mu_lp, sig_0p, op, w);

     z_pr_ar(cov, nsite, n, r, rT, T, p, N, d, d12, phip, nup, sig_ep, 
     sig_etap, sig_0p, rhop, betap, mu_lp, X, predX, op, constant, zp);
     
     accept1 += acc[0];

     for(j=0; j < p1; j++){
     fprintf(parafile, "%f ", betap[j]);
     }         
     fprintf(parafile, "%f ", rhop[0]);
     fprintf(parafile, "%f ", sig_ep[0]);
     fprintf(parafile, "%f ", sig_etap[0]);
     fprintf(parafile, "%f ", phip[0]);
     if(cov[0]==4){
      fprintf(parafile, "%f ", nup[0]);
     }     
     for(j=0; j < r1; j++) {
     fprintf(parafile, "%f ", sig_0p[j]);
     }
     for(j=0; j < r1; j++) {  
     fprintf(parafile, "%f ", mu_lp[j]);
     }

     fprintf(parafile, "\n");


// for pmcc, true fitted    
     for(j=0; j < N1; j++){
         if(i >= brin){  
          oo[0] = op[j];
          mvrnormal(constant, oo, sig_e1, constant, ot);  
          mn_rep[j] += ot[0];
          var_rep[j] += ot[0]*ot[0];
//     Rprintf("itt: %i mn %f\n", i+1, mn_rep[j]);
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
     para_printRnu (i, its1, rep1, p1, accept1, phip, nup, rhop, sig_ep, sig_etap, betap); 
     }                   
     else {              
     para_printR (i, its1, rep1, p1, accept1, phip, rhop, sig_ep, sig_etap, betap); 
     }

      if(i >= brin){
        annual_aggregate_uneqT(aggtype, nsite, r, T, rT, zp, anf);
//        annual_aggregate(aggtype, nsite, r, T, zp, anf);
  	     for(j=0; j<(nsite1*r1); j++){
           fprintf(textan, "%f ", anf[j]);
         }
         fprintf(textan, "\n");

      } // end of loop i >= brin

       ext_sige(phip, phi1);
       ext_sige(nup, nu);       
       ext_sige(sig_ep, sig_e1);
       ext_sige(sig_etap, sig_eta1);
       ext_rho(rhop, rho1); 
       ext_beta(p, betap, beta1);              
       ext_mul(r, mu_lp, mu_l1);
       ext_sigl(r, sig_0p, sig_01);
       
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

     free(phip); free(nu); free(nup); free(sig_ep); free(sig_etap); free(rhop); 
     free(betap); free(mu_lp); free(sig_0p); free(op); free(phi1); free(sig_e1); 
     free(sig_eta1); free(rho1); free(beta1); free(mu_l1); free(sig_01); 
     free(o1); free(w); free(oo); free(ot); free(acc);    

     free(zp); free(anf); 

     accept_f[0] = accept1;

     int iit;
     iit = 0;
     iit = its1 - brin;

//     Rprintf("itt: %i \n", itt[0]);

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
//     Rprintf("penalty: %4.4f, gof: %4.4f \n", var_rep[j], mn_rep[j]);
         pen += var_rep[j]; 
         go += mn_rep[j];
     }
     free(mn_rep); free(var_rep);

//     Rprintf("penalty: %4.4f, gof: %4.4f \n", pen, go);

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
void GIBBS_ar(double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *N, int *report,
     int *cov, int *spdecay, int *ft, 
	 double *shape_e, double *shape_eta, double *shape_0, double *phi_a, double *phi_b,
     double *prior_a, double *prior_b, double *prior_sig, double *phi, 
     double *tau, double *phis, int *phik, double *d, int *constant, 
     double *sig_e, double *sig_eta, double *sig_0, double *mu_l,  
     double *rho, double *beta, double *X, double *z, double *o, 
     double *phipf, double *accept, double *nupf, double *sig_epf, double *sig_etapf, 
     double *rhopf, double *betapf, double *mu_lpf, double *sig_0pf, 
     double *opf, double *wf, double *zlt_mean_sd, double *gof, double *penalty)
{
//     unsigned iseed = 44;
//     srand(iseed); 
     
     int its1, brin, col, i, j, r1, p1, N1, rep1;
     double *phip, *sig_ep, *sig_etap, *rhop, *betap;
     double *mu_lp, *sig_0p, *op;
     double *phi1, *sig_e1, *sig_eta1, *rho1, *beta1;
     double *mu_l1, *sig_01, *o1, *w;
     double *z1, *oo, *ot, *acc;

    
     its1 = *its;
     brin = *burnin;
     col = *constant;
//     n1 = *n;
     r1 = *r;
//     T1 = *T;
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
     rhop = (double *) malloc((size_t)((col)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1)*sizeof(double)));          
     mu_lp = (double *) malloc((size_t)((r1)*sizeof(double)));    
     sig_0p = (double *) malloc((size_t)((r1)*sizeof(double))); 
     op = (double *) malloc((size_t)((N1)*sizeof(double)));

     phi1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     rho1 = (double *) malloc((size_t)((col)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1)*sizeof(double)));     
     mu_l1 = (double *) malloc((size_t)((r1)*sizeof(double)));    
     sig_01 = (double *) malloc((size_t)((r1)*sizeof(double))); 
     o1 = (double *) malloc((size_t)((N1)*sizeof(double)));
     w = (double *) malloc((size_t)((N1)*sizeof(double)));
          
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
       ext_rho(rho, rho1);
       ext_beta(p, beta, beta1);
       ext_mul(r, mu_l, mu_l1);
       ext_sigl(r, sig_0, sig_01);
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

     JOINT_ar(n, T, r, rT, p, N, cov, spdecay, shape_e, shape_eta, shape_0,  
     phi_a, phi_b,
     prior_a, prior_b, prior_sig, phi1, tau, phis, phik, nu, d, constant, 
     sig_e1, sig_eta1, sig_01, mu_l1, rho1, beta1, X, z1, o1, 
     phip, acc, nup, sig_ep, sig_etap, rhop, betap, mu_lp, sig_0p, op, w);
     
     accept1 += acc[0];
     

     phipf[i] = phip[0];
     nupf[i] = nup[0];     
     sig_epf[i] = sig_ep[0];
     sig_etapf[i] = sig_etap[0];
     rhopf[i] = rhop[0];
     for(j=0; j < r1; j++) {
         sig_0pf[j+i*r1] = sig_0p[j];     
     }
     for(j=0; j < p1; j++){
         betapf[j+i*p1] = betap[j];
     }         
     for(j=0; j < r1; j++) {
         mu_lpf[j+i*r1] = mu_lp[j];     
     }
     for(j=0; j < N1; j++) {
         opf[j+i*N1] = op[j];     
         wf[j+i*N1] = w[j];     
     }

       ext_sige(phip, phi1);
       ext_sige(nup, nu);       
       ext_sige(sig_ep, sig_e1);
       ext_sige(sig_etap, sig_eta1);
       ext_rho(rhop, rho1); 
       ext_beta(p, betap, beta1);              
       ext_mul(r, mu_lp, mu_l1);
       ext_sigl(r, sig_0p, sig_01);

       
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
     para_printRnu (i, its1, rep1, p1, accept1, phip, nup, rhop, sig_ep, sig_etap, betap); 
     }                   
     else {              
     para_printR (i, its1, rep1, p1, accept1, phip, rhop, sig_ep, sig_etap, betap); 
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
          
     free(phip); free(nu); free(nup); free(sig_ep); free(sig_etap); free(rhop); 
     free(betap); free(mu_lp); free(sig_0p); free(op);
     free(phi1); free(sig_e1); free(sig_eta1); free(rho1); free(beta1); 
     free(mu_l1); free(sig_01); free(o1); free(w);
     free(z1); free(oo); free(ot); free(acc);
         
     return;
}


/////////////////////////// THE END ///////////////////////////////////////////
