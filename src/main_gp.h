//** File Name 'main_gp.h' **//
// Gaussian processes (GP) models

#include "header.h"

/******************** From "gibbs_gp.c" file ***************************/
/*****************************************************************************/

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
     double *accept_f, double *gof, double *penalty);

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
     double *gof, double *penalty);

/*
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
     double *opf, double *zlt_mean_sd, double *gof, double *penalty);

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
     double *opf, double *zlt_mean_sd, double *gof, double *penalty);
     
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
     double *opf, double *zlt_mean_sd, double *gof, double *penalty);
*/     
      
/******************** From "equation_xb_gp.c" file ***************************/
/*****************************************************************************/

void JOINT_gp(int *n, int *T, int *r, int *rT, int *p, int *N, 
     int *cov, int *spdecay, double *shape_e, double *shape_eta,
     double *phi_a, double *phi_b,
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik, double *nu,
     double *d, double *sig_e, double *sig_eta, 
     double *beta, double *X, double *z, double *o, int *constant,
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *betap, double *op);

/*
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
     double *betap, double *betasp, double *betat0p, double *betatp, double *op);
     
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
     double *betat0p, double *betatp, double *op);
          
void JOINTsp_gp(int *intercept, int *n, int *T, int *r, int *rT, int *p, int *q, int *N, 
     int *cov, int *spdecay, double *shape_e, double *shape_eta, double *shape_beta,
     double *prior_a, double *prior_b, double *prior_mubeta, 
     double *prior_sigbeta, double *prior_omu, double *prior_osig,
     double *phi, double *tau, double *phis, int *phik, double *nu,
     double *d, double *sig_e, double *sig_eta, double *sig_beta,
     double *beta, double *betas, double *X, double *Xsp, double *z, double *o, 
     int *constant,
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *sig_betasp, double *betap, double *betasp, double *op);
*/
          
void sig_e_gp(int *n, int *r, int *T, int *rT, int *N, double *shape, 
     double *prior_b, double *o, double *z, int *constant, double *sig2e);

//void sig_e_gp_sptp(int *n, int *r, int *T, int *rT, int *N, double *shape, 
//     double *prior_b, double *o, double *z, int *constant, double *sig2e);
          
void sig_eta_gp(int *n, int *r,  int *T, int *rT, double *shape, 
     double *prior_b, double *Sinv, double *XB, double *o, int *constant, 
     double *sig2eta);

/*
void sig_eta_gp_sptp(int *n, int *r,  int *T, int *rT, double *shape, 
     double *prior_b, double *Sinv, double *XB, double *o, int *constant, 
     double *sig2eta);
     
void sig_beta_gp_sp(int *n, int *q, double *shape, double *prior_b, double *betasp, 
     double *Sinv, int *constant, double *sig2beta);

void sig_0_gp_tp(int *u, double *prior_a, double *prior_b, double *gam_0, 
     int *constant, double *sig0);

void sig_del_gp_tp(int *u, int *T, double *shape, double *prior_b, double *gam_0, 
     double *gam, double *G, int *constant, double *sigdelta);

void beta_gp_tp(int *n, int *r, int *T, int *rT, int *u, double *sig0, 
     double *sigdelta, double *Qeta, double *G, double *betat, 
     double *XB, double *Xtp, double *o, int *constant, double *betat0p, 
     double *betatp);
     
void rho_gp_tp(int *u, int *T, double *prior_mu, double *prior_sig, 
     double *sigdelta, double *gam0, double *gamma, int *constant, double *rhop);
*/
                         
void beta_gp(int *n, int *r, int *T, int *rT, int *p, double *prior_mu,
     double *prior_sig, double *Qeta, double *X, double *o, int *constant, 
     double *betap);

/*     
void beta_gp_for_sp(int *n, int *r, int *T, int *rT, int *p, double *prior_mu,
     double *prior_sig, double *Qeta, double *X, double *XBsp, double *o, 
     int *constant, double *betap); 
void beta_gp_for_sptp(int *n, int *r, int *T, int *rT, int *p, double *prior_mu,
     double *prior_sig, double *Qeta, double *X, double *XBsp, double *XBtp, 
     double *o, int *constant, double *betap);
void beta_gp_sp(int *n, int *r, int *T, int *rT, int *q, int *N, double *prior_mu,
     double *prior_sig2betasp, double *betas, double *Qeta, double *Sinv, 
     double *Xsp, double *XB, double *o, int *constant, double *betap); 
*/
     
void o_gp(int *n, int *r, int *T, int *rT, double *prior_omu,
     double *prior_osig, double *sig_e, double *sig_eta, double *S, 
     double *Qeta, double *XB, double *z, int *constant, double *opost);

//void o_gp_sptp(int *n, int *r, int *T, int *rT, double *prior_omu,
//     double *prior_osig, double *sig_e, double *sig_eta, double *S, 
//     double *Qeta, double *XB, double *z, int *constant, double *opost);
          
void phi_gp_MH(double *Qeta1, double *Qeta2, double *det1, double *det2,
     double *phi1, double *phi2, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *XB, double *o, int *constant, 
     double *accept, double *phip);

//void phi_gp_MH_sptp(double *Qeta1, double *Qeta2, double *det1, double *det2,
//     double *phi1, double *phi2, int *n, int *r, int *T, int *rT, int *N, 
//     double *prior_a, double *prior_b, double *XB, double *o, int *constant, 
//     double *accept, double *phip);
          
void phi_gp_DIS(int *cov, double *Qeta1, double *det1, double *phi1,  
     double *phis, int *phik, double *nu, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *d,
     double *sig2eta, double *XB, double *o, int *constant, double *accept, 
     double *phip);
         
void phidens_gp(double *phi, double *Qeta, double *det, int *n, int *r, 
     int *T, int *rT, int *N, double *prior_a, double *prior_b, double *XB, 
     double *o, int *constant, double *out);

/*
void phi_gp_DIS_sptp(int *cov, double *Qeta1, double *det1, double *phi1,  
     double *phis, int *phik, double *nu, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *d,
     double *sig2eta, double *XB, double *o, int *constant, double *accept, 
     double *phip);
     
void phidens_gp_sptp(double *phi, double *Qeta, double *det, int *n, int *r, 
     int *T, int *rT, int *N, double *prior_a, double *prior_b, double *XB, 
     double *o, int *constant, double *out);
*/
     
void nu_gp_DIS(int *cov, double *Qeta1, double *det1, double *phi,  
     double *nu1, int *n, int *r, int *T, int *rT, int *N, 
     double *d,  double *sig2eta, double *XB, double *o, int *constant, 
     double *nup);

void nudens_gp(double *Qeta, double *det, int *n, int *r, int *T, int *rT, 
     int *N, double *XB, double *o, int *constant, double *out);

/*	 
void nu_gp_DIS_sptp(int *cov, double *Qeta1, double *det1, double *phi,  
     double *nu1, int *n, int *r, int *T, int *rT, int *N, 
     double *d,  double *sig2eta, double *XB, double *o, int *constant, 
     double *nup);
void nudens_gp_sptp(double *Qeta, double *det, int *n, int *r, int *T, int *rT, 
     int *N, double *XB, double *o, int *constant, double *out);
*/
     
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
	 
/*     
void z_pr_its_gp_sp(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *q, int *N, int *valN, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_betap,
     double *betap, double *betasp, double *X, double *valX, double *Xsp, 
     double *valXsp, double *op, int *constant, double *betapred, double *zpred);
void z_pr_gp_sp(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *q, int *N, int *valN, double *d, double *d12, double *phip, double *nup, 
     double *sig_ep, double *sig_etap, double *sig_betap, double *betap, double *betasp, 
     double *X, double *valX, double *Xsp, double *valXsp,
     double *op, int *constant, double *betapred, double *zpred);
void z_pr_its_gp_tp(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *u, int *N, int *valN, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_deltap,
     double *sig_op, double *betap, double *rhotp, double *betat0p, double *betatp,  
     double *X, double *valX, double *Xtp, double *valXtp, double *op, 
     int *constant, double *zpred);
void z_pr_gp_tp(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *u, int *N, int *valN, double *d, double *d12, double *phip, double *nup, 
     double *sig_ep, double *sig_etap, double *sig_deltap, double *sig_op,
     double *betap, double *rhotp, double *betat0p, double *betatp, double *X, 
     double *valX, double *Xtp, double *valXtp, double *op, int *constant, 
     double *zpred);
void z_pr_its_gp_sptp(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *q, int *u, int *N, int *valN, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_betap,
     double *sig_deltap, double *sig_op, double *betap, double *betasp,
     double *rhotp, double *betat0p, double *betatp, double *X, double *valX, 
     double *Xsp, double *valXsp, double *Xtp, double *valXtp, double *op, 
     int *constant, double *betapred, double *zpred);
void z_pr_gp_sptp(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *q, int *u, int *N, int *valN, double *d, double *d12, double *phip, 
     double *nup, double *sig_ep, double *sig_etap, double *sig_betap, 
     double *sig_deltap, double *sig_op, double *betap, double *betasp, 
     double *rhotp, double *betat0p, double *betatp, double *X, double *valX, 
     double *Xsp, double *valXsp, double *Xtp, double *valXtp, double *op, 
     int *constant, double *betapred, double *zpred);
*/
     
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
	 
/*
void zlt_fore_gp_sp_its(int *cov, int *its, int *K, int *nsite, int *n, int *r, 
     int *p, int *q, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_betap,
     double *foreX, double *foreXsp, double *betap, double *betasp, double *wp, 
     int *constant, double *foreZ);
          
void zlt_fore_gp_sp(int *cov, int *K, int *nsite, int *n, int *r, int *p, int *q,
     int *rT, int *T, int *rK, int *nrK, double *d, double *d12, 
     double *phi, double *nu, double *sig_e, double *sig_eta, double *sig_beta, 
     double *foreX, double *foreXsp, double *beta, double *betas, double *w, int *constant, 
     double *foreZ);

void zlt_fore_gp_tp_its(int *cov, int *its, int *K, int *nsite, int *n, int *r, 
     int *p, int *u, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_deltap,
     double *sig_op, double *foreX, double *foreXtp, double *betap, double *rhotp,
     double *betat0p, double *betatp, double *wp, int *constant, double *foreZ);
     
void zlt_fore_gp_tp(int *cov, int *K, int *nsite, int *n, int *r, int *p, int *u,
     int *rT, int *T, int *rK, int *nrK, double *d, double *d12, double *phi, 
     double *nu, double *sig_e, double *sig_eta, double *sig_delta, double *sig_op,
     double *foreX, double *foreXtp, double *beta, double *rhotp, double *betat0, 
     double *betat, double *w, int *constant, double *foreZ);

void zlt_fore_gp_sptp_its(int *cov, int *its, int *K, int *nsite, int *n, int *r, 
     int *p, int *q, int *u, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_betap, 
     double *sig_deltap, double *sig_op, double *foreX, double *foreXsp, double *foreXtp, 
     double *betap, double *betasp, double *rhotp, double *betat0p, double *betatp, 
     double *wp, int *constant, double *foreZ);
     
void zlt_fore_gp_sptp(int *cov, int *K, int *nsite, int *n, int *r, int *p, 
     int *q, int *u, int *rT, int *T, int *rK, int *nrK, double *d, double *d12, 
     double *phi, double *nu, double *sig_e, double *sig_eta, double *sig_beta,
     double *sig_delta, double *sig_op, double *foreX, double *foreXsp, 
     double *foreXtp, double *beta, double *betas, double *rhotp, double *betat0, 
     double *betat, double *w, int *constant, double *foreZ);
*/
          
/*****************************************************************************/
     
