//** File Name 'main_gp.h' **//
// Gaussian processes (GP) models

#include "header.h"

/******************** From "gibbs_gp.c" file ***************************/
/*****************************************************************************/

void GIBBS_sumpred_txt_gp(int *aggtype, double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *N, int *report,
     int *cov, int *spdecay, double *shape_e, double *shape_eta,  
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik,
     double *d, double *sig_e, double *sig_eta, 
     double *beta, double *X, double *z, double *o, int *constant, 
     int *nsite, int *valN, double *d12, double *valX, int *transform, 
     double *accept_f, double *gof, double *penalty);

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
     double *gof, double *penalty);
     

/******************** From "equation_xb_gp.c" file ***************************/
/*****************************************************************************/

void JOINT_gp(int *n, int *T, int *r, int *rT, int *p, int *N, 
     int *cov, int *spdecay, double *shape_e, double *shape_eta,
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik, double *nu,
     double *d, double *sig_e, double *sig_eta, 
     double *beta, double *X, double *z, double *o, int *constant,
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *betap, double *op);


void sig_e_gp(int *n, int *r, int *T, int *rT, int *N, double *shape, 
     double *prior_b, double *XB, double *o, double *z, double *sig_e,
     int *constant, double *sig2e);
     
void sig_eta_gp(int *n, int *r,  int *T, int *rT, double *shape, 
     double *prior_b, double *Sinv, double *XB, double *o, int *constant, 
     double *sig2eta);

void beta_gp(int *n, int *r, int *T, int *rT, int *p, double *prior_mu,
     double *prior_sig, double *Qeta, double *X, double *o, int *constant, 
     double *betap);

void o_gp(int *n, int *r, int *T, int *rT, int *p, double *prior_omu,
     double *prior_osig, double *sig_e, double *sig_eta, double *S, 
     double *Qeta, double *XB, double *z, int *constant, double *opost);
     
void phi_gp_MH(double *Qeta1, double *Qeta2, double *det1, double *det2,
     double *phi1, double *phi2, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *XB, double *o, int *constant, 
     double *accept, double *phip);
     
void phi_gp_DIS(int *cov, double *Qeta1, double *det1, double *phi1,  
     double *phis, int *phik, double *nu, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *d,
     double *sig2eta, double *XB, double *o, int *constant, double *accept, 
     double *phip);
         
void phidens_gp(double *phi, double *Qeta, double *det, int *n, int *r, 
     int *T, int *rT, int *N, double *prior_a, double *prior_b, double *XB, 
     double *o, int *constant, double *out);

void nu_gp_DIS(int *cov, double *Qeta1, double *det1, double *phi,  
     double *nu1, int *n, int *r, int *T, int *rT, int *N, 
     double *d,  double *sig2eta, double *XB, double *o, int *constant, 
     double *nup);

void nudens_gp(double *Qeta, double *det, int *n, int *r, int *T, int *rT, 
     int *N, double *XB, double *o, int *constant, double *out);

/******************** From "prediction_xb_gp.c" file ***************************/
/*****************************************************************************/


void z_pr_its_gp(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *N, int *valN, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *betap,   
     double *X, double *valX, double *op, int *constant, double *zpred);
     
void z_pr_gp(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *N, int *valN, double *d, double *d12, double *phip, double *nup, 
     double *sig_ep, double *sig_etap, double *betap, double *X, double *valX,
     double *op, int *constant, double *zpred);
     

/******************** From "forecast_xb_gp.c" file ***************************/
/*****************************************************************************/

void zlt_fore_gp_its(int *cov, int *its, int *K, int *nsite, int *n, int *r, 
     int *p, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *foreX, 
     double *betap, double *wp, int *constant, double *foreZ);
     
    
void zlt_fore_gp(int *cov, int *K, int *nsite, int *n, int *r, int *p, 
     int *rT, int *T, int *rK, int *nrK, double *d, double *d12, 
     double *phi, double *nu, double *sig_e, double *sig_eta, double *foreX, 
     double *beta, double *w, int *constant, double *foreZ);
     
