//** File Name 'main_gpp.h' **//

#include "header.h"


/**************************** From "gibbs_ar.c" *********************************/
/*****************************************************************************/

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
     double *accept_f, double *gof, double *penalty);   


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
     double *gof, double *penalty, double *z_mean_sd);

  
/******************** From "equation_xb_gpp.c" file ***************************/
/*****************************************************************************/

void JOINT_onephi_gpp(int *cov, int *spdecay, double *flag, int *n, int *m, 
     int *T, int *r, int *rT, int *p, int *N, double *shape_e, double *shape_eta, 
     double *shape_l, double *prior_a, double *prior_b, double *mu_beta, 
     double *delta2_beta, double *mu_rho,  double *delta2_rho, double *alpha_l, 
     double *delta2_l, double *phi, double *tau, double *phis, int *phik,
     double *nu, double *dm, double *dnm, int *constant, 
     double *sig2e, double *sig2eta, double *sig2l, double *beta, 
     double *rho, double *mu_l, double *X, double *z, double *w0, double *w,
     double *phip, double *accept, double *nup,
     double *sig2ep, double *sig2etap, double *betap, double *rhop, 
     double *mu_lp, double *sig2lp, double *w0p, double *wp, 
     double *zfit);


void sig_e_gpp(int *n, int *rT, int *N, double *shape, double *prior_b, 
     double *XB, double *Aw, double *z, int *constant, double *sig2e);
     
void sig_eta_gpp(int *m, int *r,  int *T, int *rT, double *shape, 
     double *prior_b, double *Sinv_eta, double *rho, double *w, 
     double *w0, int *constant, double *sig2eta); 

void beta_gpp(int *n, int *p, int *rT, int *N, double *mu_beta, 
     double *delta2_beta, double *sig2e, double *X, double *Aw, 
     double *z, int *constant, double *beta);

     
void rho_gpp(int *m, int *r, int *T, int *rT, int *p, double *mu_rho, 
     double *delta2, double *Q_eta, double *w0, double *w, 
     int *constant, double *rho);
     
void sig_l_gpp(int *m, int *r, double *shape, double *prior_b, double *mu_l, 
     double *Sinv_0, double *w0, int *constant, double *sig2l);

void mu_l_gpp(int *m, int *r, double *sig2l, double *alpha_l, double *delta2_l, 
     double *Sinv_0, double *w0, int *constant, double *mu_l);
          
void Z_fit_gpp(double *flag, int *n, int *m, int *T, int *r, int *rT, 
     double *sig2e, double *Aw, double *XB, double *z, int *constant, 
     double *zfit);
     
void o_fit_gpp(double *flag, int *n, int *m, int *T, int *r, int *rT, 
     double *Aw, double *XB, double *z, int *constant, 
     double *zfit);

/*****************************************************************************/

void wlt_gpp_s(int *n, int *m, int *r, int *T, int *rT, int *p, double *dm, 
     double *phi_eta, double *sig2e, double *sig2eta, double *rho, 
     double *Q_eta, double *A, double *w0, double *w, double *XB, double *z, 
     int *constant, double *wp);

void wlt_gpp(int *n, int *m, int *r, int *T, int *rT, int *p, double *sig2e, 
     double *rho, double *Q_eta, double *A, double *w0, double *w, double *XB, 
     double *z, int *constant, double *wp);

void w0_gpp(int *m, int *r, int *T, double *Q_eta, double *sig2l, 
     double *Sinv_0, double *rho, double *mu_l, double *w, 
     int *constant, double *w0p);
 
     
void ww_gpp(int *n, int *m, int *r, int *T, int *rT, int *p, int *N, 
     double *phi_eta, double *phi_0, double *dm, double *dnm, double *sig2e, 
     double *sig2eta, double *sig2l, double *mu_l, double *rho, double *beta, 
     double *w, double *X, double *z, int *constant, double *wp);

/*****************************************************************************/

void phi_gpp_DIS2(int *cov, double *Qeta1, double *det1, double *phi1,    
     double *phis, int *phik, double *nu, int *m, int *r, int *T, int *rT, 
     double *prior_a, double *prior_b, double *dm, double *rho, 
     double *sig2eta, double *mu_l, double *w0, double *w, int *constant, 
     double *accept, double *phip);

void phiden_gpp(double *phi, double *Qeta, double *det, int *m, int *r, 
     int *T, int *rT, double *prior_a, double *prior_b, double *rho, 
     double *w0, double *w, int *constant, double *out);

void nu_gpp_DIS(int *cov, double *Qeta1, double *det1, double *phi, double *nu1, 
     int *m, int *r, int *T, int *rT, double *dm, double *rho, double *sig2eta, 
     double *mu_l, double *w0, double *w, int *constant, double *nup);

void nuden_gpp(double *Qeta, double *det, int *m, int *r, int *T, int *rT, 
     double *rho, double *w0, double *w, int *constant, double *out);
     
void phi_gpp_MH2(double *Qeta1, double *Qeta2, double *det1, double *det2,
     double *phi1, double *phi2,   
     int *m, int *r, int *T, int *rT, double *prior_a, double *prior_b, 
     double *rho, double *mu_l, double *w0, double *w, int *constant, 
     double *accept, double *phip);

void phi_gpp_MH(int *cov, double *phi2, double *nu, double *dm, double *dnm,
     double *Sinv1, double *det1, double *phi1, double *A1,   
     int *n, int *m, int *r, int *T, int *rT, double *prior_a, double *prior_b, 
     double *rho, double *mu_l, double *w0, double *w, double *z, double *XB, 
     int *constant, double *accept, double *phip);
     

/******************* From "prediction_xb_gpp.c" file **************************/
/*****************************************************************************/
void z_pr_its_gpp1(int *cov, int *scale, int *its, int *nsite, int *n, int *m, int *r, 
     int *T, int *rT, int *p, int *nsiterT, double *phi_etap, double *nup, double *dm, 
     double *dnsm, double *wp, double *sig2ep, double *betap, double *Xpred, 
     int *constant, double *zpred);

void z_pr_gpp1(int *cov, int *nsite, int *n, int *m, int *r, int *T, 
     int *rT, int *p, int *nsiterT, double *phi_etap, double *nup, double *dm, 
     double *dnsm, double *wp, double *sig2ep, double *betap, 
     double *Xpred, int *constant, double *zpred);

void z_pr_its_gpp_wtilde(int *cov, int *scale, int *its, int *nsite, 
     int *n, int *m, int *r, int *T, int *rT, int *p, int *nsiterT, 
     double *phi_etap, double *dm, double *dnsm, double *w0p, double *wp, 
     double *sig2ep, double *sig2eta, double *sig2lp, double *rhop, 
     double *betap, double *Xpred, int *constant, double *zpred);
     
void z_pr_gpp_wtilde(int *cov, int *nsite, 
     int *n, int *m, int *r, int *T, int *rT, int *p, int *nsiterT, 
     double *phi_etap, double *dm, double *dnsm, double *w0p, double *wp, 
     double *sig2ep, double *sig2etap, double *sig2lp, double *rhop, 
     double *betap, double *Xpred, int *constant, double *zpred);
          
void z_pr_its_gpp(int *cov, int *scale, 
     int *its, int *nsite, int *n, int *m, int *r, int *T, 
     int *rT, int *p, int *nsiterT, double *phi_etap, double *phi_0p, 
     double *dm, double *dnsm, double *w0p, double *wp, double *sig2ep,  
     double *sig2etap, double *sig2lp, double *mu_lp, double *rhop, 
     double *betap, double *Xpred, int *constant, double *zpred);
     
void z_pr_gpp(int *cov, int *nsite, int *n, int *m, int *r, int *T, int *rT, 
     int *p, 
     int *nsiterT, double *phi_etap, double *phi_0p, double *dm, double *dmns, 
     double *w0p, double *wp, double *sig2ep, double *sig2etap,
     double *sig2lp, double *mu_lp, double *rhop, double *betap,
     double *Xpred, int *constant, double *zpred);

// multivariate approach                            
void wlt_pr_gpp_multi(int *nsite, int *m, int *r, int *T, 
     double *w0pred, double *wp, double *w0p, 
     double *sig2etap, double *rhop, double *Sinv, double *Snsm,
     int *constant, double *wpred);

// multivariate approach
void w0_pr_gpp_multi(int *nsite, int *m, int *r, double *w0p, double *mu_lp, 
     double *sig2lp, double *Sinv, double *Snsm, int *constant, 
     double *w0pred);
 
void wlt_pr_gpp(int *nsite, int *m, int *r, int *T, 
     double *w0pred, double *wp, double *w0p, 
     double *sig2eta, double *rhop, double *Sinv, double *S12,
     int *constant, double *wpred);

void w0_pr_gpp(int *nsite, int *m, int *r, double *w0p, double *mu_lp, 
     double *sig2lp, double *Sinv, double *S12, int *constant, 
     double *w0pred);


/************************* From forecast_xb_ar.c ***************************/
/*****************************************************************************/

void zlt_fore_gpp_its(int *cov, int *its, int *K, int *n, int *m, 
     int *r, int *p, int *rT, int *T, int *rK, int *nrK, double *dnm, double *dm, 
     double *phip, double *nup,  
     double *sig_ep, double *sig_etap, double *betap, double *rhop, double *wp, 
     double *foreX, int *constant, double *foreZ);
            
void zlt_fore_gpp(int *cov, int *K, int *n, int *m, int *r, int *p, int *rT, int *T, 
     int *rK, int *nrK, double *dnm, double *dm, double *phi, double *nu, double *sig_e, 
     double *sig_eta, double *beta, double *rho, double *wp, double *foreX, 
     int *constant, double *foreZ);


/////////////////////////////////////////////////////////////////////////////////




