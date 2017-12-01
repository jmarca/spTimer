#include "header.h"

/*****************************************************************************/
/*****************************************************************************/


/*********************** From "randgenerator.c" file *************************/
/*****************************************************************************/

double drand48 (void);
double myrunif(int constant);
double rexpon (double scale);
double rnormal(double mu, double sigma);
double rgammaa (double shape, double rate);
double rigammaa (double shape, double rate);

void runif_val(int *n, int *constant, double *out);
void rexp_val(int *n, double *scale, double *result);
void rnorm_val(int *n, double *m, double *s, double *result);
void rgamma_val(int *n, double *shape, double *rate, double *result);
void rigamma_val(int *n, double *shape, double *rate, double *result);

double rnorm_for_multivariate(void);
void chol_for_multivariate(double *s, int *n, double *ltrg);
void mvrnormal(int *n, double *mu, double *s2, int *p, double *y);

/*****************************************************************************/
/*****************************************************************************/
