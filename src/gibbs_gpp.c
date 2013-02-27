
#include "main_gpp.h"
#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"
#include "Print.h"


// with all summary values (mean, variance/sd, low2.5, up97.5)
// also the predictions into another sites
// one phi value
// by defalt text output: as big data files
void GIBBS_sumpred_gpp(int *aggtype, int *cov, int *spdecay, double *flag, 
     int *its, int *burnin,
     int *n, int *m, int *T, int *r, int *rT, int *p, int *N, int *report,
     double *shape_e, double *shape_eta, double *shape_l,  
     double *prior_a, double *prior_b, double *mu_beta, double *delta2_beta,
     double *mu_rho,  double *delta2_rho, double *alpha_l, double *delta2_l,
     double *phi, double *tau, double *phis, int *phik,
     double *dm, double *dnm, int *constant, 
     double *sig2e, double *sig2eta, double *sig2l, double *beta, 
     double *rho, double *mu_l, double *X, double *z, double *w0, double *w,
     int *nsite, int *nsiterT, double *dnsm, double *Xpred, int *transform, 
     double *accept_f, double *gof, double *penalty)    
{
//     unsigned iseed = 44;
//     srand(iseed); 
     
     int its1, col, i, j, n1, m1, r1, T1, p1, N1, rep1, nsite1, brin, trans1;
     its1 = *its;
     col = *constant;
     n1 = *n;
     m1 = *m;
     r1 = *r;
     T1 = *T;
     p1 = *p;
     N1 = *N;
     rep1 = *report;
     nsite1 = *nsite;
     brin = *burnin;
     trans1 = *transform;

     double accept1, *mn_rep, *var_rep, *mn_pred, *var_pred;
     mn_rep = (double *) malloc((size_t)((N1)*sizeof(double)));          
     var_rep = (double *) malloc((size_t)((N1)*sizeof(double)));          
     mn_pred = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double)));          
     var_pred = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double)));          

     accept1 = 0.0;

     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      
     for(j=0; j<nsite1*r1*T1; j++){
        mn_pred[j] = 0.0;
        var_pred[j] = 0.0;
     }      

     double *phip, *sig2ep, *sig2etap, *rhop, *betap;
     double *mu_lp, *sig2lp, *wp, *w0p, *zfit, *ofit, *acc;
     phip = (double *) malloc((size_t)((col)*sizeof(double)));          
     sig2ep = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2etap = (double *) malloc((size_t)((col)*sizeof(double)));
     rhop = (double *) malloc((size_t)((col)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1*col)*sizeof(double)));          
     mu_lp = (double *) malloc((size_t)((r1*col)*sizeof(double)));    
     sig2lp = (double *) malloc((size_t)((r1*col)*sizeof(double))); 
     wp = (double *) malloc((size_t)((m1*r1*T1)*sizeof(double))); 
     w0p = (double *) malloc((size_t)((m1*r1)*sizeof(double))); 
     zfit = (double *) malloc((size_t)((N1*col)*sizeof(double)));
     ofit = (double *) malloc((size_t)((N1*col)*sizeof(double)));     
     acc = (double *) malloc((size_t)((col)*sizeof(double)));

     double *phi1, *sig2e1, *sig2eta1, *rho1, *beta1;
     phi1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     rho1 = (double *) malloc((size_t)((col)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1*col)*sizeof(double)));     

     double *zp, *anf; //, tmp[T1];
     zp = (double *) malloc((size_t)((nsite1*r1*T1)*sizeof(double)));      
     anf = (double *) malloc((size_t)((nsite1*r1)*sizeof(double)));      

     double *nu, *nup;
     nu = (double *) malloc((size_t)((col)*sizeof(double)));
     nup = (double *) malloc((size_t)((col)*sizeof(double)));     
     nu[0] = 0.5;
          
       ext_sige(phi, phi1);
       ext_sige(sig2e, sig2e1);
       ext_sigeta(sig2eta, sig2eta1);
       ext_rho(rho, rho1);
       ext_beta(p, beta, beta1);


     FILE *parafile;
     parafile = fopen("OutGPP_Values_Parameter.txt", "w");
     FILE *fitfile;
     fitfile = fopen("OutGPP_Stats_FittedValue.txt", "w");
     FILE *prdfile;
     prdfile = fopen("OutGPP_Stats_PredValue.txt", "w");
     FILE *predfile;
     predfile = fopen("OutGPP_Values_Prediction.txt", "w");

     int type1;
     type1= *aggtype;
   
     FILE *textan;
     // none
     if(type1==0){
       textan = fopen("OutGPP_NONE.txt", "w");
     }
     // annual average value
     if(type1==1){
       textan = fopen("OutGPP_Annual_Average_Prediction.txt", "w");
     }
     // annual 4th highest value
     if(type1==2){
       textan = fopen("OutGPP_Annual_4th_Highest_Prediction.txt", "w");
     }
     // annual w126 option
     if(type1==3){
       textan = fopen("OutGPP_Annual_w126_Prediction.txt", "w");
     }

     GetRNGstate();            
     for(i=0; i < its1; i++) {

     JOINT_onephi_gpp(cov, spdecay, flag, n, m, T, r, rT, p, N, 
     shape_e, shape_eta, shape_l, prior_a, prior_b, mu_beta, delta2_beta, 
     mu_rho, delta2_rho, alpha_l, delta2_l, phi1, tau, phis, phik, nu, dm, dnm, 
     constant, sig2e1, sig2eta1, sig2l, beta1, rho1, mu_l, X, z, w0, w, 
     phip, acc, nup, sig2ep, sig2etap, betap, rhop, mu_lp, sig2lp, w0p, wp, zfit);

     z_pr_gpp1(cov, nsite, n, m, r, T, rT, p, nsiterT, phip, nup, dm, dnsm, 
     wp, sig2ep, betap, Xpred, constant, zp);

     accept1 += acc[0];

     for(j=0; j < p1; j++){
     fprintf(parafile, "%f ", betap[j]);
     }         
     fprintf(parafile, "%f ", rhop[0]);
     fprintf(parafile, "%f ", sig2ep[0]);
     fprintf(parafile, "%f ", sig2etap[0]);
     fprintf(parafile, "%f ", phip[0]);
     if(cov[0]==4){
      fprintf(parafile, "%f ", nup[0]);
     }     
     fprintf(parafile, "\n");


// for pmcc, fitted sum values    
     for(j=0; j < N1; j++){
         if(i >= brin){  
          mn_rep[j] += zfit[j];
          var_rep[j] += zfit[j]*zfit[j];
         } 
     }


// prediction samples
     if(i >= brin){
     for(j=0; j<(nsite1*r1*T1); j++){
         if(trans1 == 0){
           zp[j] = zp[j];
           fprintf(predfile, "%f ", zp[j]);
         }
         if(trans1 == 1){
           zp[j] = zp[j]*zp[j];
           fprintf(predfile, "%f ", zp[j]);
         }
         if(trans1 == 2){
           zp[j] = exp(zp[j]);
           fprintf(predfile, "%f ", zp[j]);
         }
     }                            
     fprintf(predfile, "\n");
     }	

// prediction samples summary
     for(j=0; j<(nsite1*r1*T1); j++){
         if(i >= brin){  
          mn_pred[j] += zp[j];
          var_pred[j] += zp[j]*zp[j];
         }
     }

       ext_sige(phip, phi1);
       ext_sige(sig2ep, sig2e1);
       ext_sigeta(sig2etap, sig2eta1); 
       ext_rho(rhop, rho1); 
       ext_beta(p, betap, beta1);              

     if(cov[0]==4){
     para_printRnu (i, its1, rep1, p1, accept1, phip, nup, rhop, sig2ep, sig2etap, betap); 
     }                   
     else {              
     para_printR (i, its1, rep1, p1, accept1, phip, rhop, sig2ep, sig2etap, betap); 
     }

      if(i >= brin){
         annual_aggregate(aggtype, nsite, r, T, zp, anf);
  	     for(j=0; j<(nsite1*r1); j++){
           fprintf(textan, "%f ", anf[j]);
         }
         fprintf(textan, "\n");
/*
       // annual average value
       if(type1==1){
         annual_aggregate(aggtype, nsite, r, T, zp, anf);
  	     for(j=0; j<(nsite1*r1); j++){
           fprintf(textan, "%f ", anf[j]);
         }
         fprintf(textan, "\n");
       }
       // annual 4th highest value
       if(type1==2){
         annual_aggregate(aggtype, nsite, r, T, zp, anf);
  	     for(j=0; j<(nsite1*r1); j++){
           fprintf(textan, "%f ", anf[j]);
         }
         fprintf(textan, "\n");
       }     
       // annual w126 
       if(type1==3){
         annual_aggregate(aggtype, nsite, r, T, zp, anf);
  	     for(j=0; j<(nsite1*r1); j++){
           fprintf(textan, "%f ", anf[j]);
         }
         fprintf(textan, "\n");
       }     
*/
      } // end of loop i >= brin
     } // end of iteration loop
     PutRNGstate();

     fclose(parafile);
     fclose(predfile);
     fclose(textan);

    free(phip); free(sig2ep); free(sig2etap); free(rhop); free(betap);
    free(mu_lp); free(sig2lp); free(wp); free(w0p); free(zfit); free(ofit); 
    free(acc); 
    free(phi1); free(sig2e1); free(sig2eta1); free(rho1); free(beta1);
    free(zp); free(anf); 

     Rprintf(" Calculating Summary Statistics \n ... ");

     accept_f[0] = accept1;

//     int iit;
//     iit = 0;
//     iit = its1 - brin;

	int *iit;
    iit = (int *) malloc((size_t)((col)*sizeof(int)));
	iit[0] = its[0] - burnin[0];

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit[0];
          var_rep[j] = var_rep[j]/iit[0];
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          fprintf(fitfile, "%f , %f \n", mn_rep[j], sqrt(var_rep[j]));
     }
     fclose(fitfile);

// pred mean and sd
     for(j=0; j < nsite1*r1*T1; j++){
          mn_pred[j] = mn_pred[j]/iit[0];
          var_pred[j] = var_pred[j]/iit[0];
          var_pred[j] = var_pred[j] - mn_pred[j]*mn_pred[j];
          fprintf(prdfile, "%f , %f \n", mn_pred[j], sqrt(var_pred[j]));
     }
     fclose(prdfile);
     free(mn_pred); free(var_pred); free(iit); 	
   

     double pen, go;
     pen =0;
     go =0;
     
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

     penalty[0] = pen;
     gof[0] = go;       
     free(mn_rep); free(var_rep); 

     Rprintf(" ... ");
     Rprintf("\n---------------------------------------------------------\n");
       
     return;
}




// with zfit summary values (mean and variance/sd)
// with one phi parameter
void GIBBS_zfitsum_onephi_gpp(int *cov, int *spdecay, double *flag, int *its, int *burnin,
     int *n, int *m, int *T, int *r, int *rT, int *p, int *N, int *report,
     double *shape_e, double *shape_eta, double *shape_l,  
     double *prior_a, double *prior_b, double *mu_beta, double *delta2_beta,
     double *mu_rho,  double *delta2_rho, double *alpha_l, double *delta2_l,
     double *phi_eta, double *tau_eta, double *phis, int *phik, 
     double *dm, double *dnm, int *constant, 
     double *sig2e, double *sig2eta, double *sig2l, double *beta, 
     double *rho, double *mu_l, double *X, double *z, double *w0, double *w,
     int *transform, double *phi_etaf, double *accept_etaf, double *nupf,
     double *sig2ef, double *sig2etaf, double *betaf, double *rhof, 
     double *mu_lf, double *sig2lf, double *w0f, double *wf, 
     double *gof, double *penalty, double *z_mean_sd)
{
//     unsigned iseed = 44;
//     srand(iseed); 
     
     int its1, col, i, j, n1, m1, r1, T1, p1, N1, rep1, brin, trans1;
     its1 = *its;
     col = *constant;
     n1 = *n;
     m1 = *m;
     r1 = *r;
     T1 = *T;
     p1 = *p;
     N1 = *N;
     rep1 = *report;
     brin = *burnin;
     trans1 = *transform;

     double accept1, mn_rep[N1], var_rep[N1];
     accept1 = 0.0;
     for(j=0; j<N1; j++){
        mn_rep[j] = 0.0;
        var_rep[j] = 0.0;
     }      

     double *phi_etap, *sig2ep, *sig2etap, *rhop, *betap;
     double *mu_lp, *sig2lp, *wp, *w0p, *zfit;
     double *phi_eta1, *sig2e1, *sig2eta1, *rho1, *beta1;
     double *mu_l1, *sig2l1; //, *w01, *w1;
     double *acc_eta, *out; // *z1;
          
     phi_etap = (double *) malloc((size_t)((col)*sizeof(double)));          
     sig2ep = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2etap = (double *) malloc((size_t)((col)*sizeof(double)));
     rhop = (double *) malloc((size_t)((col)*sizeof(double)));
     betap = (double *) malloc((size_t)((p1*col)*sizeof(double)));          
     mu_lp = (double *) malloc((size_t)((r1*col)*sizeof(double)));    
     sig2lp = (double *) malloc((size_t)((r1*col)*sizeof(double))); 
     wp = (double *) malloc((size_t)((m1*r1*T1)*sizeof(double))); 
     w0p = (double *) malloc((size_t)((m1*r1)*sizeof(double))); 
     zfit = (double *) malloc((size_t)((N1*col)*sizeof(double)));

     phi_eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2e1 = (double *) malloc((size_t)((col)*sizeof(double)));
     sig2eta1 = (double *) malloc((size_t)((col)*sizeof(double)));
     rho1 = (double *) malloc((size_t)((col)*sizeof(double)));
     beta1 = (double *) malloc((size_t)((p1*col)*sizeof(double)));     
     mu_l1 = (double *) malloc((size_t)((r1*col)*sizeof(double)));    
     sig2l1 = (double *) malloc((size_t)((r1*col)*sizeof(double))); 
//     w01 = (double *) malloc((size_t)((m1*r1)*sizeof(double))); 
//     w1 = (double *) malloc((size_t)((m1*r1*T1)*sizeof(double))); 
     
     acc_eta = (double *) malloc((size_t)((col)*sizeof(double)));
     out = (double *) malloc((size_t)((col)*sizeof(double)));     
//     z1 = (double *) malloc((size_t)((N1)*sizeof(double)));

     double *nu, *nup;
     nu = (double *) malloc((size_t)((col)*sizeof(double)));
     nup = (double *) malloc((size_t)((col)*sizeof(double)));     
     nu[0] = 2.0/3.0;
          
       ext_sige(phi_eta, phi_eta1);
       ext_sige(sig2e, sig2e1);
       ext_sigeta(sig2eta, sig2eta1);
       ext_rho(rho, rho1);
       ext_beta(p, beta, beta1);
       ext_mul(r, mu_l, mu_l1);
       ext_sigl(r, sig2l, sig2l1);
//       ext_o(N, z, z1);
            
/*
// for missing
     for(j=0; j < N1; j++){
          if (flag[j] == 1.0){
          out[0]=z[j];
          mvrnormal(constant, out, sig2e1, constant, out);  
          z1[j] =  out[0];
          }
          else {
          z1[j] = z[j];
          }     
     }
*/

     GetRNGstate();
     for(i=0; i < its1; i++) {

     JOINT_onephi_gpp(cov, spdecay, flag, n, m, T, r, rT, p, N, 
     shape_e, shape_eta, shape_l, prior_a, prior_b, mu_beta, delta2_beta, 
     mu_rho, delta2_rho, alpha_l, delta2_l, phi_eta1, tau_eta, phis, phik, nu, 
     dm, dnm, constant, sig2e1, sig2eta1, sig2l1, beta1, rho1, mu_l1, X, z, 
     w0, w, phi_etap, acc_eta, nup, sig2ep, sig2etap, betap, rhop, mu_lp, sig2lp, 
     w0p, wp, zfit);

     accept1 += acc_eta[0];
        
     phi_etaf[i] = phi_etap[0];
     nupf[i] = nup[0];
     
     sig2ef[i] = sig2ep[0];
     sig2etaf[i] = sig2etap[0];
     rhof[i] = rhop[0];
     for(j=0; j < r1; j++) {
         sig2lf[j+i*r1] = sig2lp[j];     
     }
     for(j=0; j < p1; j++){
         betaf[j+i*p1] = betap[j];
     }         
     for(j=0; j < r1; j++) {
         mu_lf[j+i*r1] = mu_lp[j];     
     }
     for(j=0; j<m1*r1; j++){
         w0f[j+i*m1*r1] = w0p[j];
     }     
     for(j=0; j<m1*r1*T1; j++){
         wf[j+i*m1*r1*T1] = wp[j];
     }     


// for pmcc, fitted sum values    
     for(j=0; j < N1; j++){
         if(i >= brin){  
          mn_rep[j] += zfit[j];
          var_rep[j] += zfit[j]*zfit[j];
         } 
     }

       ext_sige(phi_etap, phi_eta1);
       ext_sige(nup, nu);       
       ext_sige(sig2ep, sig2e1);
       ext_sigeta(sig2etap, sig2eta1); 
       ext_rho(rhop, rho1); 
       ext_beta(p, betap, beta1);              

     if(cov[0]==4){
     para_printRnu (i, its1, rep1, p1, accept1, phi_etap, nup, rhop, sig2ep, sig2etap, betap); 
     }                   
     else {              
     para_printR (i, its1, rep1, p1, accept1, phi_etap, rhop, sig2ep, sig2etap, betap); 
     }
//     printR(i, its1); 
     
     } // end of iteration loop
     PutRNGstate();

     free(phi_etap); free(nu); free(nup); free(sig2ep); free(sig2etap); 
     free(rhop); free(betap); free(mu_lp); free(sig2lp); free(wp); free(w0p); 
     free(zfit); free(phi_eta1); free(sig2e1); free(sig2eta1); free(rho1); 
     free(beta1); free(mu_l1); free(sig2l1); free(acc_eta); free(out); 
     //free(z1); // free(w01); free(w1); 

     accept_etaf[0] = accept1;

     int iit;
     iit = 0;
     iit = its1 - brin;

// fitted zlt, mean and sd
     for(j=0; j < N1; j++){
          mn_rep[j] = mn_rep[j]/iit;
          var_rep[j] = var_rep[j]/iit;
          var_rep[j] = var_rep[j] - mn_rep[j]*mn_rep[j];
          z_mean_sd[j] = mn_rep[j];
          z_mean_sd[N1+j] = sqrt(var_rep[j]);
     }

     double pen, go;
     pen = 0;
     go =0;

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
     penalty[0] = pen;
     gof[0] = go;       

     
     return;
}



/////////////////////////// THE END ////////////////////////////
