//** File Name 'main_ar.h' **//

#include "header.h"

/**************************** From "gibbs_ar.c" *********************************/
/*****************************************************************************/

void GIBBS_sumpred_txt_ar(int *aggtype, double *flag, int *its, int *burnin,
     int *n, int *T, int *r, int *rT, int *p, int *N, int *report,
     int *cov, int *spdecay, double *shape_e, double *shape_eta, double *shape_0,  
     double *phi_a, double *phi_b,
     double *prior_a, double *prior_b, double *prior_sig, double *phi, 
     double *tau, double *phis, int *phik, double *d, int *constant, 
     double *sig_e, double *sig_eta, double *sig_0, double *mu_l,  
     double *rho, double *beta, double *X, double *z, double *o, 
     int *nsite, int *predN, double *d12, double *predX, int *transform,
     double *accept_f, double *gof, double *penalty);
     
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
     double *opf, double *wf, double *zlt_mean_sd, double *gof, double *penalty);


  
/******************** From "equation_xb_ar.c" file ***************************/
/*****************************************************************************/


void JOINT_ar(int *n, int *T, int *r, int *rT, int *p, int *N, 
     int *cov, int *spdecay,
     double *shape_e, double *shape_eta, double *shape_0,  
     double *phi_a, double *phi_b,
     double *prior_a, double *prior_b, double *prior_sig, double *phi, 
     double *tau, double *phis, int *phik, double *nu, double *d, int *constant, 
     double *sig_e, double *sig_eta, 
     double *sig_l0, double *mu_l,  
     double *rho, double *beta, double *X, double *z, double *o, 
     double *phip, double *nup, double *accept,
     double *sig_ep, double *sig_etap, double *rhop, double *betap, 
     double *mu_lp, double *sig_l0p, double *op, double *w);

void JOINTsp_ar(int *n, int *T, int *r, int *rT, int *p, int *q, int *N, 
     int *cov, int *spdecay, double *shape_e, double *shape_eta, double *shape_0,  
     double *prior_a, double *prior_b, double *prior_sig, double *phi, 
     double *tau, double *phis, int *phik, double *nu, double *d, int *constant, 
     double *sig_e, double *sig_eta, double *sig_l0, double *mu_l, double *rho, 
     double *beta, double *betas, double *X, double *Xsp, double *z, double *o, 
     double *phip, double *accept, double *nup, double *sig_ep, double *sig_etap, 
     double *rhop, double *betap, double *betasp, double *mu_lp, double *sig_l0p, 
     double *op, double *w);
     
void w_ar(int *n, int *r, int *T, int *rT, int *p, double *O_l0, 
     double *X, double *o, double *thetap, int *constant, double *w);

void sig2_ar(int *n, int *r,  int *T, int *rT, int *p,
     double *shape_e, double *shape_eta, double *prior_b, 
     double *S, double *Sinv, double *rho, double *O_l0, 
     double *XB, double *o, double *z, int *constant, 
     double *sig2ep, double *sig2etap);

void sig_e_ar(int *n, int *r, int *T, int *rT, double *shape, double *prior_b, 
     double *z, double *o, int *constant, double *sig2e);

void sig_eta_ar(int *n, int *r,  int *T, int *rT, int *p, double *phi,
     double *shape, double *prior_b, double *Sinv, double *rho, 
     double *O_l0, double *XB, double *o, int *constant, 
     double *sig2eta);

void theta_ar(int *n, int *r, int *T, int *rT, int *p, double *prior_sig, 
     double *Q_eta, double *O_l0, double *X, double *o, int *constant, 
     double *thetap);

void theta_ar_for_sp(int *n, int *r, int *T, int *rT, int *p, double *prior_sig, 
     double *Q_eta, double *O_l0, double *X, double *XBsp, double *o, 
     int *constant, double *thetap);

void beta_ar_sp(int *n, int *r, int *T, int *rT, int *q, //double *prior_mu,
     double *prior_sig, double *Q_eta, double *rho, double *O_l0, double *Xsp, 
     double *XB, double *o, int *constant, double *betap);
          
void beta_ar(int *n, int *r, int *T, int *rT, int *p, 
     double *prior_sig, double *Q_eta, double *rho, 
     double *O_l0, double *X, double *o, int *constant, double *beta);
     
void rho_ar(int *n, int *r, int *T, int *rT, int *p, 
     double *prior_sig, double *Q_eta, double *O_l0, 
     double *XB, double *o, int *constant, double *rho);
     
void sig_0l_ar(int *n, int *r, double *shape, double *prior_b, double *phi,
     double *mu_l, double *O_l0, double *Sinv, int *constant, double *sig0);

void mu_l_ar (int *n, int *r, double *sig_0, double *Sinv, double *O_l0, 
     int *constant, double *mu_l);
     
void Z_fitfnc(int *its, int *N, double *sig_ep, double *o_p, 
                int *constant, double *z_p);


/*****************************************************************************/

void o0_ar(int *n, int *r, int *T, int *rT, int *p, double *sig_eta, 
     double *sig_l0, double *rho, double *mu_l, double *Sinv, 
     double *XB, double *o, int *constant, double *o0post);

void o_ar(int *n, int *r, int *T, int *rT, int *p, double *sig_e, 
     double *sig_eta, double *rho, double *S, double *Q_eta, 
     double *O_l0, double *XB, double *z, double *o, 
     int *constant, double *opost);
          

/*****************************************************************************/

void phi_ar_MH(double *Sinv1, double *Sinv2, double *det1, double *det2,
     double *phi1, double *phi2,
     int *n, int *r, int *T, int *rT, int *p, int *N, 
     double *prior_a, double *prior_b, double *rho,  
     double *mu_l, double *O_l0, double *XB, double *o, 
     int *constant, double *accept, double *phip);

void phi_ar_DIS(int *cov, double *Qeta1, double *det1, double *phi1, 
     double *phis, int *phik, double *nu, int *n, int *r, int *T, int *rT, int *N, 
     double *prior_a, double *prior_b, double *d, double *sig2eta, 
     double *rho, double *mu_l, double *O_l0, double *XB, double *o, 
     int *constant, double *accept, double *phip);

void phidens_ar(double *phi, double *Qeta, double *det, int *n, int *r, 
     int *T, int *rT, int *N, double *prior_a, double *prior_b, double *XB, 
     double *rho, double *O_l0, double *o, int *constant, double *out);

void nu_ar_DIS(int *cov, double *Qeta1, double *det1, double *phi1, double *nu, 
     int *n, int *r, int *T, int *rT, int *N, double *d, double *sig2eta, 
     double *rho, double *mu_l, double *O_l0, double *XB, double *o, 
     int *constant, double *nup);

void nudens_ar(double *Qeta, double *det, int *n, int *r, int *T, int *rT, 
     int *N, double *XB, double *rho, double *O_l0, double *o, int *constant, 
     double *out);
     
/******************* From "dic_ar.c" file **************************/
/*****************************************************************************/

void deviance_ar(int *n, int *T, int *r, int *rT, int *p, int *N, 
     double *Sinv, double *sig_e, double *sig_eta, double *rho, double *XB, 
     double *z, double *o, double *O_l0, int *constant, double *deviance);



/******************* From "prediction_xb_ar.c" file **************************/
/*****************************************************************************/

void z_pr_its_ar(int *cov, int *its, int *nsite, int *n, int *r, int *rT, 
     int *T, int *p, int *N, double *d, double *d12, 
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *sig_l0p, 
     double *rhop, double *betap, double *mu_lp,  
     double *X, double *valX, double *op, int *constant, double *zpred);

void z_pr_ar(int *cov, int *nsite, int *n, int *r, int *rT, int *T, int *p, 
     int *N, double *d, double *d12, double *phip, double *nup,
     double *sig_ep, double *sig_etap, double *sig_l0p, double *rhop, 
     double *betap, double *mu_lp,  double *X, double *valX,
     double *op, int *constant, double *zpred);
     

/************************* From forecast_xb_ar.c ***************************/
/*****************************************************************************/

void zlt_fore_ar_its_anysite(int *cov, int *its, int *K, int *nsite, int *n, 
     int *r, int *p, int *rT, int *T, int *rK, int *nrK, double *d, double *d12,
     double *phip, double *nup, double *sig_ep, double *sig_etap, double *rhop, 
     double *foreX, double *betap, double *zpred_exist, double *wp, 
     int *constant, double *foreZ);

void zlt_fore_ar(int *cov, int *K, int *nsite, int *n, int *r, int *p, 
     int *rT, int *T, int *rK, int *nrK, double *d, double *d12, 
     double *phi, double *nu, double *sig_e, double *sig_eta, double *rho, 
     double *foreX, double *beta, double *z, double *w, 
     int *constant, double *foreZ);
     
/////////////////////////////////////////////////////////////////////////////////




