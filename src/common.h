
#include "header.h"

/*****************************************************************************/
/*****************************************************************************/

/**************************** From "common.c" file ***************************/
/*****************************************************************************/

void ratio_fnc(double *ratio, int *constant, double *U);

void sort_4th(double *sample, int *n, int *r, int *T, double *an4th);

int sort_fnc(const void *x, const void *y);

void test_C_read(double *input, int *T, double *out);

void test_RW(int *its, int *constant);

void GP_para_printR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *beta);
void GP_para_printRnu (int i, int iteration, int report, int p, double accept, 
     double *phi, double *nu, double *sig2e, double *sig2eta, double *beta);

void GPsp_para_printR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *sig2beta, double *beta); 
void GPsp_para_printRnu (int i, int iteration, int report, int p, double accept, 
     double *phi, double *nu, double *sig2e, double *sig2eta, double *sig2beta, 
     double *beta);

void GPtp_para_printR (int i, int iteration, int report, int p, int u, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *sig2delta, double *sig20,
     double *rho, double *beta);
void GPtp_para_printRnu (int i, int iteration, int report, int p, int u, double accept, 
     double *phi, double *nu, double *sig2e, double *sig2eta, double *sig2delta, 
     double *sig20, double *rho, double *beta);

void GPsptp_para_printR (int i, int iteration, int report, int p, int u, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *sig2beta, double *sig2delta, 
     double *sig20, double *rho, double *beta);
void GPsptp_para_printRnu (int i, int iteration, int report, int p, int u, double accept, 
     double *phi, double *nu, double *sig2e, double *sig2eta, double *sig2beta, 
     double *sig2delta, double *sig20, double *rho, double *beta);

void GPPsp_para_printRnu (int i, int iteration, int report, int p, double accept, 
     double *phi, double *nu, double *rho, double *sig2e, double *sig2eta, 
     double *sig2beta, double *beta);
void GPPsp_para_printR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *rho, double *sig2e, double *sig2eta, 
     double *sig2beta, double *beta);
                    
void para_printR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *rho, double *sig2e, double *sig2eta, double *beta);
void para_printRnu (int i, int iteration, int report, int p, double accept, 
     double *phi, double *nu, double *rho, double *sig2e, double *sig2eta, 
     double *beta);

void para_print_spTR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *beta);
     
void printR (int i, int iteration);

void ext_sumstat(int i, int *its, double *x, double *alt);
     
void ext_sumstat_burnin(int i, int *its, int *burnin, double *x, double *alt);

void extn_12(int j, int *n, double *S_12, double *S_12c);

void comb_XB_tp(int *n, int *r, int *T, int *p, double *Xtp, double *betatp, 
     int *constant, double *XB);
void comb_XB_sp(int *n, int *r, int *T, int *q, double *Xsp, double *betasp, 
     int *constant, double *XB);
void comb_XB_sp_gpp(int *n, int *m, int *r, int *T, int *q, double *Xsp, 
     double *betasp, double *A, int *constant, double *XB);

void extract_X(int t, int l, int *n, int *r, int *T, int *p, 
     double *x, double *alt);
void extract_X_uneqT(int t, int l, int *n, int *r, int *T, int *rT, int *p, 
     double *x, double *alt);     
void extract_X_sp2(int t, int l, int j, int *n, int *r, int *T,  
     double *x, double *alt);
void extract_X_sp3(int t, int l, int j, int *n, int *r, int *T, 
     double *x, double *alt);
void extract_X_sp4(int t, int l, int i, int j, int *n, int *r, int *T, 
     double *x, double *alt);     
void extract_X_sp(int t, int l, int i, int j, int *n, int *r, int *T, 
     double *x, double *alt);

void extract_beta_sp(int j, int *n, double *betasp, double *alt);
void extract_beta_sp2(int j, int *n, int *q, double *betasp, double *alt);

void extract_beta_t(int t, int *T, int *p, double *beta, double *alt);
void extract_beta_l(int l, int *r, int *p, double *beta, double *alt);

void extract_X5(int t, int *n, int *r, int *T, int *p, double *x, double *alt);
void extract_X4(int i, int l, int *n, int *r, int *p, double *x, double *alt);
void extract_X41(int l, int t, int i, int *n, int *rT, int *T, int *p, 
     double *x, double *alt);
void extract_X41_uneqT(int l, int t, int i, int *n, int *rT, int *r, int *T, 
     int *p, double *x, double *alt);

     
void ext(double *IN, double *OUT);
void ext_wlt(int *m, int *r, int *T, double *wp, double *w);

void ext_sige(double *sig_ep, double *sig_e);
void ext_sigeta(double *sig_etap, double *sig_eta);
void ext_rho(double *rhop, double *rho);
void ext_beta(int *p, double *betap, double *beta);
void ext_xil(int *r, double *xi_lp, double *xi_l);
void ext_gaml(int *n, int *r, double *gamma_lp, double *gamma_l);
void ext_mul(int *r, double *mu_lp, double *mu_l);
void ext_sigl(int *r, double *sig_lp, double *sig_l);
void ext_o(int *N, double *op, double *o);
void ext_v(int *M, double *vp, double *v); 


void extract_X21(int l, int t, int *n, int *rT, int *T, int *p, 
     double *x, double *alt);
void extract_X21_uneqT(int l, int t, int *n, int *rT, int *r, int *T, int *p, 
     double *x, double *alt);
     
void extract_X2(int l, int t, int *n, int *rT, int *T, int *p, 
     double *x, double *alt);

void extract_X3(int l, int t, int k, int *n, int *rT, int *T, int *p, 
     double *x, double *alt);

void extract_X3_uneqT(int l, int t, int k, int *n, int *r, int *rT, int *T, 
     int *p, double *x, double *alt);
     
void extract_alt(int *l, int *t, int *n, int *rT, int *T, double *x, 
     double *alt);
     
void extract_alt2(int l, int t, int *n, int *rT, int *T, double *x, 
     double *alt);

void extract_alt_uneqT(int l, int t, int *n, int *r, int *T, int *rT, 
     double *x, double *alt);

void put_together(int *n, int *r, int *T, double *x, double *out);

void put_together1(int l, int t, int *n, int *r, int *T, 
     double *x, double *alt);

void put_together1_uneqT(int l, int t, int *n, int *r, int *T, int *rT, 
     double *x, double *alt);

void put(int t, int *n, int *T, double *x, double *alt);
void extract(int t, int *n, int *T, double *x, double *alt);

/////////////////////////////////////////////////////////////////////
