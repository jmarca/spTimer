
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

void extract_X(int t, int l, int *n, int *r, int *T, int *p, 
     double *x, double *alt);

void extract_beta_t(int t, int *T, int *p, double *beta, double *alt);
void extract_beta_l(int l, int *r, int *p, double *beta, double *alt);

void extract_X5(int t, int *n, int *r, int *T, int *p, double *x, double *alt);
void extract_X4(int i, int l, int *n, int *r, int *p, double *x, double *alt);

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
     
void extract_X2(int l, int t, int *n, int *rT, int *T, int *p, 
     double *x, double *alt);

void extract_X3(int l, int t, int k, int *n, int *rT, int *T, int *p, 
     double *x, double *alt);
     
void extract_alt(int *l, int *t, int *n, int *rT, int *T, double *x, 
     double *alt);
     
void extract_alt2(int l, int t, int *n, int *rT, int *T, double *x, 
     double *alt);

void put_together(int *n, int *r, int *T, double *x, double *out);

void put_together1(int l, int t, int *n, int *r, int *T, 
     double *x, double *alt);

void put(int t, int *n, int *T, double *x, double *alt);
void extract(int t, int *n, int *T, double *x, double *alt);

/////////////////////////////////////////////////////////////////////
