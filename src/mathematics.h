
#include "header.h"

/*****************************************************************************/
/*****************************************************************************/

/****************************** From "math.c" file ***************************/
/*****************************************************************************/

double w126_from_daily(); 
double max_pos(); 

void annual_aggregate(int *aggtype, int *n, int *r, int *T, double *z, 
     double *out);
void sum_stat(int *its, int *N, double *tX, int *constant, 
     double *Mean, double *Median, double *Var, double *Low, double *Up);
void sum_stat2(int *its, int *burnin, int *N, double *tX, int *constant, 
         double *Mean, double *SD, double *Low, double *Up);
void stats(int *n, double *x, double *ave, double *sd, double *low, double *up);
void variance(int *n, double *x, double *var);
void stdeviation(int *n, double *x, double *sd);
void mean( int *n, double *x, double *ave);
void median(int *n, double *x, double *med);
void range(int *n, double *x, double *low, double *up);
void maximum(int *n, double *x, double *maxi);
void minimum(int *n, double *x, double *mini);
double fabs(double x);
void absol(double *a, double *b, double *out);

void GeoDist_miles(int *n, double *Lat, double *Long, double *out);
void GeoDist_km(int *n, double *Lat, double *Long, double *out);

double geodeticdistance (double *point1, double *point2);
double sqr(double x);

//void matprint2(double *x, int m, int n);
void IdentityM(int *n, double *I);
void MAdd(double *x, int *xrow, int *xcol, double *y, double *out);
//void MAdd2(double *x, int *xrow, int *xcol, double *y, double *out);
void MSub(double *x, int *xrow, int *xcol, double *y, double *out);
void MProd(double *y, int *nycol, int *nyrow, double *x, 
     int *nxrow, double *out);
//void MProd3(double *y, int *nycol, int *nyrow, double *x, 
//     int *nxrow, double *out);
          
void MTranspose(double *x, int *ncol, int *nrow, double *tx);

void QuadMat(double *a, int *n, double *x, int *p, double *out);
void QuadMat2(double *a, int *n, double *x, int *p, double *out);
void xTay(double *x, double *A, double *y, int *p, double *xTAy);
double xTay2(double *x, double *A, double *y, int p);

void determinant_val(double mat[], int *order, double *out);
double determinant(double mat[], int order);
double* submatrix(double mat[], int m, int n, int order);

void MInv(double *S, double *inv, int *pp, double *det);
void sq_rt(double *a, double *t, int p, double *det);
void trans_pose (double *a, double *b, int m, int n);
void tinv_mat (double *t, double *c, int p);
void mat_mul(double *ma1, double *ma2, double *mul, double no, double po, 
     double qo);

void MProd2(double *ma1, double *ma2, double *mul, double *no,
     double *po, double *qo);
     
/*****************************************************************************/
/*****************************************************************************/
