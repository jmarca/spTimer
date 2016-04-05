//** File Name 'math.c' **//

#include "header.h"
#include "covariance.h"
#include "common.h"
#include "math.h"
#include "mathematics.h"
#include "randgenerator.h"

#ifdef USE_MKL
#include <mkl_blas.h>
#include <mkl_lapack.h>
#include <mkl_vml_functions.h>
#include <mkl_vml_defines.h>
#endif


//double w126_from_daily(); 
//double max_pos(); 
double max_pos(double *x, int n)
{ 
  int i; 
  double v=-9999999.0;

  for (i=0; i<n; i++) { 
    if (v < x[i]) v = x[i]; 
  }
  return(v); 
}


double w126_from_daily(double *day) 
{ 
  // day should have 214 values 
  // may not contain missing values as well

  int l;
  double v, apr, may, jun, jul, aug, sep, oct, tot5[5];

     apr=0.0; may=0.0; jun=0.0; jul=0.0; 
     aug = 0.0; sep=0.0; oct=0.0; 

     for (l=0; l<30; l++) 
       apr += day[l]; 
     for (l=30; l<61; l++) 
       may += day[l]; 
     for (l=61; l<91; l++) 
       jun += day[l]; 
     for (l=91; l<122; l++) 
       jul += day[l]; 
     for (l=122; l<153; l++)
       aug += day[l]; 
     for (l=153; l<183; l++) 
       sep += day[l];
     for (l=183; l<214; l++) 
       oct += day[l];
     
     tot5[0] = apr + may + jun; 
     tot5[1] = may + jun + jul; 
     tot5[2] = jun + jul + aug;
     tot5[3] = jul + aug + sep;
     tot5[4] = aug + sep + oct;
     v = max_pos(tot5, 5); 
     return (v); 
}


/***************************************************/
// annual summary stat

void annual_aggregate(int *aggtype, int *n, int *r, int *T, double *z, 
     double *out)
{
      int j, l, t, n1, r1, T1, type1;
      n1= *n;
      r1= *r;
      T1= *T;
      type1= *aggtype;

      double *tmp, *ave;
      tmp = (double *) malloc((size_t)((T1)*sizeof(double)));      
      ave = (double *) malloc((size_t)((1)*sizeof(double)));      
            
	  for(j=0; j<n1; j++){
	    for(l=0; l<r1; l++){
	        for(t=0; t<T1; t++){
                tmp[t] = z[t+l*T1+j*T1*r1];
            }
            // NONE
            if(type1==0){
 		       out[l+j*r1] = 0.0; 
            }         
            // annual average value
            if(type1==1){
               mean(T, tmp, ave);                 
 		       out[l+j*r1] = ave[0]; 
            }         
            // annual 4th highest value
            if(type1==2){
 		        qsort(tmp, T1, sizeof(double), sort_fnc);
		        out[l+j*r1] = tmp[T1-4];
	        }
	        // annual w126 option
            if(type1==3){
   	            out[l+j*r1] = w126_from_daily(tmp);
            }         
	    }
      }
      free(tmp); free(ave);
      return;
}

// with unequal T
void annual_aggregate_uneqT(int *aggtype, int *n, int *r, int *T, int *rT, 
     double *z, double *out)
{
      int j, l, t, n1, r1, rT1, type1;
      n1= *n;
      r1= *r;
      rT1= *rT;
      type1= *aggtype;

      double *tmp, *ave;
      tmp = (double *) malloc((size_t)((rT1)*sizeof(double)));      
      ave = (double *) malloc((size_t)((1)*sizeof(double)));      

     int *T1, *T2; 
     T1 = (int *) malloc((size_t)((r1)*sizeof(int)));
     T2 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     for(j=0; j<r1; j++){
          T1[j] = T[j];
     }
     cumsumint(r, T, T2);  
     
	  for(j=0; j<n1; j++){
	    for(l=0; l<r1; l++){
	        for(t=0; t<T1[l]; t++){
                tmp[t] = z[t+T2[l]+j*rT1];
            }
            // NONE
            if(type1==0){
 		       out[l+j*r1] = 0.0; 
            }         
            // annual average value
            if(type1==1){
               mean(T, tmp, ave);                 
 		       out[l+j*r1] = ave[0]; 
            }         
            // annual 4th highest value
            if(type1==2){
 		        qsort(tmp, T1[l], sizeof(double), sort_fnc);
		        out[l+j*r1] = tmp[T1[l]-4];
	        }
	        // annual w126 option
            if(type1==3){
   	            out[l+j*r1] = w126_from_daily(tmp);
            }         
	    }
      }
      
      free(T1); free(T2);
      free(tmp); free(ave);
      return;
}



/**************************************************/
// Summary Statistics for the MCMC samples
// X is : N x its matrix
// tX is : its x N matrix
void sum_stat(int *its, int *N, double *tX, int *constant, 
     double *Mean, double *Median, double *Var, double *Low, double *Up)
{
     int i, col, its1;
     col =*constant;
     its1 =*its;
     
     double *x, *ave, *med, *var, *low, *up;
     x = (double *) malloc((size_t)((its1)*sizeof(double)));
     ave = (double *) malloc((size_t)((col)*sizeof(double)));
     med = (double *) malloc((size_t)((col)*sizeof(double)));
     var = (double *) malloc((size_t)((col)*sizeof(double)));
     low = (double *) malloc((size_t)((col)*sizeof(double)));
     up = (double *) malloc((size_t)((col)*sizeof(double)));
                                        
     for(i=0; i<*N; i++){
       ext_sumstat(i, its, tX, x); // its x 1
       mean(its, x, ave);
       Mean[i] = ave[0];
       median(its, x, med);
       Median[i] = med[0];
       variance(its, x, var);
       Var[i] = var[0];
       range(its, x, low, up);
       Low[i] = low[0];
       Up[i] = up[0];
     }

     free(x); free(ave); free(med); free(var); free(low); free(up);
     
     return;
}       

// Summary Statistics for the MCMC samples
// X is : N x its matrix
// tX is : its x N matrix
void sum_stat2(int *its, int *burnin, int *N, double *tX, int *constant, 
         double *Mean, double *SD, double *Low, double *Up)
{
     int i, col, its1, burn1;
     col =*constant;
     its1 =*its;
     burn1 =*burnin;
     
     double *x, *ave, *sd, *low, *up;
     x = (double *) malloc((size_t)((its1-burn1)*sizeof(double)));
     ave = (double *) malloc((size_t)((col)*sizeof(double)));
     sd = (double *) malloc((size_t)((col)*sizeof(double)));
     low = (double *) malloc((size_t)((col)*sizeof(double)));
     up = (double *) malloc((size_t)((col)*sizeof(double)));
     
     int *itt;
     itt = (int *) malloc((size_t)((1)*sizeof(int)));
     *itt = its1-burn1;
     
     for(i=0; i<*N; i++){
       ext_sumstat_burnin(i, its, burnin, tX, x);
       mean(itt, x, ave);
       Mean[i] = ave[0];
       stdeviation(itt, x, sd);
       SD[i] = sd[0];
       range(itt, x, low, up);
       Low[i] = low[0];
       Up[i] = up[0];
     }

     free(itt);
     free(x); free(ave); free(sd); free(low); free(up);
     
     return;
}       


/**************************************************/
// Mean, SD, low and up function together
void stats(int *n, double *x, double *ave, double *sd, double *low, double *up)
{
        int i,j,n1; 
        n1= *n;
        double tt = 0.0; 
        double vv = 0.0;
        double temp = 0.0; 

        /* total all scores */
        for (i=0; i < n1; i++) {
            tt += x[i];
            vv += x[i]*x[i];
        } 
        tt = tt/n1;
        *ave = tt;
        vv = vv/n1 - tt*tt;
        *sd = sqrt(vv);

        for (i=0; i < n1; i++) 
         for(j=i+1; j< n1; j++) {
             if(x[i]>x[j]) {
             temp=x[j];
             x[j]=x[i];
             x[i]=temp;
             }
         }
        
        *low = x[n1 * 25/1000];
        *up = x[n1 * 975/1000-1];

        return;
}

     

// Sample Variance of a variable x, with total number of observations n
void variance(int *n, double *x, double *var)
{
     int i, n1;
     n1= *n;
     
     double vv= 0;
     double *xbar;
     xbar = (double *) malloc((size_t)((n1)*sizeof(double)));

     mean(n, x, xbar);
     for(i=0; i< n1; i++) {
          vv += ((x[i]-xbar[0])*(x[i]-xbar[0]));
     }     
     vv = vv/(n1-1);
     *var = vv;

     free(xbar);
     return;
}


/**************************************************/
// Sample SD of a variable x, with total number of observations n
void stdeviation(int *n, double *x, double *sd)
{
     int i, n1;
     n1= *n;
     
     double vv= 0;
     double *xbar;
     xbar = (double *) malloc((size_t)((n1)*sizeof(double)));

     mean(n, x, xbar);
     for(i=0; i< n1; i++) {
          vv += ((x[i]-xbar[0])*(x[i]-xbar[0]));
     }     
     vv = vv/(n1-1);
     *sd = sqrt(vv);

     free(xbar);
     return;
}

/**************************************************/
// sum of a variable x, with total number of observations n
void sum( int *n, double *x, double *tot)
{
        int i,n1; 
        n1= *n;
        double tt = 0; 
        /* total all scores */
        for (i=0; i < n1; i++) {
                tt += x[i];
        } 
        *tot=tt;
        return;
} 

/**************************************************/
// sum of the integer variable x, with total number of observations n
void sumint( int *n, int *x, int *tot)
{
        int i,n1; 
        n1= *n;
        int tt = 0; 
        /* total all scores */
        for (i=0; i < n1; i++) {
                tt += x[i];
        } 
        *tot=tt;
        return;
} 

/**************************************************/
// cumulative sum of a variable x, with total number of observations n
// out shall be vector of n+1, first element is zero
void cumsum( int *n, double *x, double *out)
{
        int i,n1; 
        n1= *n;
        double tt = 0; 

        out[0] = 0.0;
        /* total all scores */
        for (i=0; i < n1; i++) {
                tt += x[i];
                out[i+1] = tt;
        } 
        return;
} 

/**************************************************/
// cumulative sum of integer variable x, with total number of observations n
// out shall be vector of n+1, first element is zero
// used in void extract_X_uneqT function

void cumsumint( int *n, int *x, int *out)
{
        int i; 
        out[0] = 0; 
        /* total all scores */
        for (i=0; i < *n; i++) {
            out[i+1] = out[i] + x[i];
        } 
        return;
} 

/*
void cumsumint( int *n, int *x, int *out)
{
        int i,n1; 
        n1= *n;
        int tt, *out1;
        out1 = (int *) malloc((size_t)((n1)*sizeof(int)));        

        out[0] = 0; 
        // total all scores 
        for (i=0; i < n1; i++) {
            tt += x[i];
            out1[i] = tt;
            out[i+1] = out1[i];
        } 
        free(out1);
        return;
} 
*/

/**************************************************/
// Mean of a variable x, with total number of observations n
void mean( int *n, double *x, double *ave)
{
        int i,n1; 
        n1= *n;
        double tt = 0; 
        /* total all scores */
        for (i=0; i < n1; i++) {
                tt += x[i];
        } 
        tt = tt/n1;
        *ave=tt;
        return;
} 

/**************************************************/
// Median of a variable x, with total number of observations n
void median(int *n, double *x, double *med)
{
     double temp;
     int i,j, n1;
     n1 = *n;
     for(i=0;i< *n;i++)
     for(j=i+1;j< *n;j++)
     {
      if(x[i]>x[j])
      {
      temp=x[j];
      x[j]=x[i];
      x[i]=temp;
      }
     }
     if(n1%2==0)
      *med = (x[n1/2]+x[(n1/2)-1])/2;
     else
      *med = x[n1/2];
      return;
}


/**************************************************/
// 2.5% and 97.5% of a variable x, with total number of observations n
void range(int *n, double *x, double *low, double *up)
{
     double temp; // intp1, *intpart;
     int i,j, n1;
     n1 = *n;
//     intpart = (double *) malloc((size_t)((1)*sizeof(double)));
          
     for(i=0;i< *n;i++)
     for(j=i+1;j< *n;j++)
     {
      if(x[i]>x[j])
      {
      temp=x[j];
      x[j]=x[i];
      x[i]=temp;
      }
     }
//     intp1 = 0.0;
//     modf((n1 * 25/1000), intpart);
//     intp1 = *intpart;             
//     *low = x[intp1];
     *low = x[n1 * 25/1000];
//     intp1 = 0.0;
//     modf((n1 * 975/1000-1), intpart);
//     intp1 = *intpart;             
//     *up = x[intp1];
     *up = x[n1 * 975/1000-1];
     
//     free(intpart);
      return;
}


/*************************************************/
// maximum of a variable x, with total number of observations n
void maximum(int *n, double *x, double *maxi)
{
     double temp;
     int i,j, n1;
     n1 = *n;
     for(i=0;i< *n;i++){
     for(j=i+1;j< *n;j++)
     {
      if(x[i]>x[j])
      {
      temp=x[j];
      x[j]=x[i];
      x[i]=temp;
      }
     }
     }
     *maxi = x[n1-1];
      return;
}

/*************************************************/
// minimum of a variable x, with total number of observations n
void minimum(int *n, double *x, double *mini)
{
     double temp;
     int i,j;
//   n1 = *n;
     for(i=0;i< *n;i++){
     for(j=i+1;j< *n;j++)
     {
      if(x[i]>x[j])
      {
      temp=x[j];
      x[j]=x[i];
      x[i]=temp;
      }
     }
     }
     
     *mini = x[0];
      return;
}


/*************************************************/
// Absolute value
double fabs(double x)
{
       double ans = 0.0;       
       if(x > 0){
          ans = x;
       }
       else if (x == 0){
          ans = x;
       }
       else {
          ans = - x;
       }
       return ans;
}          



// Absolute value of two double variables
void absol(double *a, double *b, double *out)
{
  double aa, bb, tmp, cc;
  aa = *a;
  bb = *b;
  cc = aa-bb;
  if(cc > 0){tmp = cc;}
  else{tmp = bb-aa;}
  tmp = aa - bb;
  if(tmp < 0) { cc = bb - aa; } 
  else { cc = tmp; }
  *out = cc;
  return;
}


// factorial function
//double factorial (int num)
//{
// int result=1;
// for (int i=1; i<=num; ++i)
//    result=result*=i;
// return result;
//}

int factorial ( int x ) {
    int ret = 1;
    for ( int i = 2; i <= x; i++ ) {
        ret = ret * i;
    }
    return ret;
}


// Gamma function: Spouge's approximation

double sp_gamma(double z)
{
  const int a = 12;
  static double c_space[12];
  static double *c = NULL;
  int k;
  double accm;
 
  if ( c == NULL ) {
    double k1_factrl = 1.0; // (k - 1)!*(-1)^k with 0!==1
    c = c_space;
    c[0] = sqrt(2.0*PI);
    for(k=1; k < a; k++) {
      c[k] = exp(a-k) * pow(a-k, k-0.5) / k1_factrl;
	  k1_factrl *= -k;
    }
  }
  accm = c[0];
  for(k=1; k < a; k++) {
    accm += c[k] / ( z + k );
  }
  accm *= exp(-(z+a)) * pow(z+a, z+0.5); // Gamma(z+1) 
  return accm/z;
}
 



/***********************************************/
// Convert (lat, long) to GCD (miles)
void GeoDist_miles(int *n, double *Lat, double *Long, double *out)
{
  int nn, i, j1, j; 
  double tmp, tmp1, LatI, LatJ, LongI, LongJ, u;
  
  nn = *n; 
  tmp = 0.0; tmp1 = 0.0;
  LatJ = 0.0;  LatJ = 0.0;
  LongI = 0.0;  LongJ = 0.0;
  u = 0.0;

  for(i=0;i<(nn-1);i++) 
  { 
     out[i*nn+i] = 0.0; 
     j1=i+1;

     for(j=j1;j<nn;j++)
     { 
        LatI = Lat[i];  LatJ = Lat[j];
        LongI = Long[i];  LongJ = Long[j];
        absol(&LongI, &LongJ, &tmp1);
        u =  sin(LatI*PI/180)*sin(LatJ*PI/180) + 
        cos(LatI*PI/180)*cos(LatJ*PI/180)*cos(tmp1*PI/180);
        if(u > 1.0){
        tmp = 0.0;
        }
        else if(u < -1.0){
        tmp = 0.0;
        }
        else{     
        tmp = acos(u)*180/PI*69.17334;
        }
     // we put constant 111.32387 to obtain distance in km
     // we put constant 69.17334 to obtain distance in miles
        out[i*nn+j] = tmp;
        out[j*nn+i] = out[i*nn+j]; 
     }
  }
  return;
}

// Convert (lat, long) to GCD (km)
void GeoDist_km(int *n, double *Lat, double *Long, double *out)
{
  int nn, i, j1, j; 
  double tmp, tmp1, LatI, LatJ, LongI, LongJ, u;
  
  nn = *n; 
  tmp = 0.0; tmp1 = 0.0;
  LatJ = 0.0;  LatJ = 0.0;
  LongI = 0.0;  LongJ = 0.0;
  u = 0.0;
  
  for(i=0;i<(nn-1);i++) 
  { 
     out[i*nn+i] = 0.0; 
     j1=i+1;

     for(j=j1;j<nn;j++)
     { 
        LatI = Lat[i];  LatJ = Lat[j];
        LongI = Long[i];  LongJ = Long[j];
        absol(&LongI, &LongJ, &tmp1);
        u =  sin(LatI*PI/180)*sin(LatJ*PI/180) + 
        cos(LatI*PI/180)*cos(LatJ*PI/180)*cos(tmp1*PI/180);
        if(u > 1.0){
        tmp = 0.0;
        }
        else if(u < -1.0){
        tmp = 0.0;
        }
        else{     
        tmp = acos(u)*180/PI*111.32387;
        }
     // we put constant 111.32387 to obtain distance in km
     // we put constant 69.17334 to obtain distance in miles
        out[i*nn+j] = tmp;
        out[j*nn+i] = out[i*nn+j]; 
     }
  }
  return;
}

// OR
/************************************************************/
double geodeticdistance(double *point1, double *point2)
{
  double R = 6371.0, p1rad[2], p2rad[2], d, u=0.0;
  int i;
  
  for (i=0; i<2; i++) {
    p1rad[i] = point1[i] * M_PI/180.0;
    p2rad[i] = point2[i] * M_PI/180.0;
    u += sqr(point1[i] - point2[i]); 
  }

  if (u >  0.00000001) {  
  d = sin(p1rad[1])*sin(p2rad[1])+cos(p1rad[1])*cos(p2rad[1])*
    cos(fabs(p1rad[0]-p2rad[0]));
  d = acos(d);
 /*  printf("acosd=%f distance=%f \n", d, R*d); */
  } else d =0.0; 
  return (R*d);
}

double sqr(double x)
{
  return x * x;
}
/*
// Print a matrix of X (m x n)
void matprint2(double *x, int m, int n)
{
     int i,j;
     for ( i=0; i<m; i++) { 
       printf("%d  ", i);
          for (j=0; j<n; j++)
           printf("  %4.2f ",x[i*n+j]);
           printf("\n");
      }
}
*/


// Create an Identity Matrix
void IdentityM(int *n, double *I)
{
     int i, j, n1;
     n1 =*n;
     
     for(i=0; i<n1; i++){
         for(j=0; j< n1; j++){
         if(j==i){
           I[j+n1*i] = 1.0;
         }
         else{
           I[j+n1*i] = 0.0;
         }
         }
     }
     
     return;
}            

// Create an diagonal Matrix with values
void IdentityMX(int *n, double *x, double *I)
{
     int i, j, n1;
     n1 =*n;
     
     for(i=0; i<n1; i++){
         for(j=0; j< n1; j++){
         if(j==i){
           I[j+n1*i] = x[i];
         }
         else{
           I[j+n1*i] = 0.0;
         }
         }
     }
     return;
}            
         
                       

// some of the following code are taken from GDLM software of Dou, Le, Zedak
// Calculate the matrix addition of X (m x n) and Y (m x n)
// LPACK, BLAS version is modified by Rohan Shah
void MAdd(double *x, int *xrow, int *xcol, double *y, double *out)
{
//#ifdef USE_MKL
//	vdAdd((*xrow) * (*xcol), x, y, out);
//#else
 int m, n, i, j;
 double tmp;

 m = *xrow;
 n = *xcol;

 for(i = 0; i < m; i++) {
   for(j =0; j < n; j++) {
    tmp = x[i*n+j] + y[i*n+j];
    out[i*n+j] = tmp;
   }
 }
 return;
// #endif
}

/*
void MAdd2(double *x, int *xrow, int *xcol, double *y, double *out)
{
#ifdef USE_MKL
	vdAdd((*xrow) * (*xcol), x, y, out);
#else
 int m, n, i, j;
 double tmp;

 m = *xrow;
 n = *xcol;

 for(i = 0; i < m; i++) {
   for(j =0; j < n; j++) {
    tmp = x[i*n+j] + y[i*n+j];
    out[i*n+j] = tmp;
   }
 }
 return;
#endif
}

*/

// Calculate the matrix subtraction of X (m x n) and Y (m x n)
// LPACK, BLAS version is modified by Rohan Shah
void MSub(double *x, int *xrow, int *xcol, double *y, double *out)
{
#ifdef USE_MKL
	vdSub((*xrow) * (*xcol), x, y, out);
#else      
 int nn, pp, i, j;
 double tmp;

 nn = *xrow;
 pp = *xcol;

 for(i = 0; i < nn; i++) {
   for(j =0; j < pp; j++) {
    tmp = x[i*pp+j] - y[i*pp+j];
    out[i*pp+j] = tmp;
   }
 }  
 return;
 #endif
}

/*
// Calculate the matrix product of X (q x p) and Y (p x r), 
void MProd(double *y, int *nycol, int *nyrow, double *x, 
     int *nxrow, double *out)
{ 
      
     double tmp;
     int i, j, k, r, p, q;

  r=*nycol;    
  p=*nyrow;       // nyrow = nxcol
  q=*nxrow;

  for(i = 0; i < r; i ++){ 
    for(j = 0; j < q; j ++){ 
         tmp = 0;
         for(k = 0; k < p; k++) { tmp += y[i*p+k]*x[k*q+j]; }
         out[i*q+j] = tmp; 
    }
  }
  return;
}
*/

// Calculate the matrix product of X (q x p) and Y (p x r), 
// LPACK, BLAS version is modified by Rohan Shah
void MProd(double *y, int *nycol, int *nyrow, double *x, 
     int *nxrow, double *out)
{ 
#ifdef USE_MKL
	double one = 1, zero = 0;
	dgemm("N", "N", nxrow, nycol, nyrow, &one, x, nxrow, y, nyrow, &zero, out, nxrow);
#else
     double *tmp;
     tmp = (double *) malloc((size_t)((1)*sizeof(double)));
     int i, j, k, r, p, q;

  r=*nycol;    
  p=*nyrow;       // nyrow = nxcol
  q=*nxrow;

  for(i = 0; i < r; i ++){ 
    for(j = 0; j < q; j ++){ 
         tmp[0] = 0;
         for(k = 0; k < p; k++) { tmp[0] += y[i*p+k]*x[k*q+j]; }
         out[i*q+j] = tmp[0]; 
    }
  }
  free(tmp);
  return;
#endif
}

// Transpose of a matrix
// X is nrow x ncol
void MTranspose(double *x, int *ncol, int *nrow, double *tx)
{
   
  int i , j, m, n;
  
  m = *ncol; 
  n = *nrow;

  for(i = 0; i < n; i++) { 
    for(j = 0; j < m; j++){ tx[i*m+j] = x[j*n+i]; }
  }
  return;
}

// Calculate the quadratic matrix product: xax'
void QuadMat(double *a, int *n, double *x, int *p, double *out)
{
     
  double *tmpmat, *tx;
  int nn, pp;

  nn=*n;
  pp=*p;

  tmpmat=(double *) malloc((size_t)((pp*nn)*sizeof(double)));
  tx=(double *) malloc((size_t)((nn*pp)*sizeof(double)));

  MTranspose(x, n, p, tx);
  MProd(tx, p, n, a, n, tmpmat);
  MProd(tmpmat, p, n, x, p, out);

  free(tmpmat); free(tx);
  return;
}


// Calculate the quadratic matrix product: x'ax
void QuadMat2(double *a, int *n, double *x, int *p, double *out)
{
  double *tmpmat, *tx;
  int nn, pp;

  nn=*n;
  pp=*p;

  tmpmat=(double *) malloc((size_t)((pp*nn)*sizeof(double)));
  tx=(double *) malloc((size_t)((nn*pp)*sizeof(double)));

  MProd(x, p, n, a, n, tmpmat);
  MTranspose(x, p, n, tx);
  MProd(tmpmat, p, n, tx, p, out);

  free(tmpmat); free(tx);
  return;
}
/*
// Another way to calculate x'Ax
double quad(double *x, double *A, int *p)
{
 double *z, u;
 z = (double *) malloc(p *  sizeof(double));
 mat_vec(A, x, z, p, p);
 u =  cross_vec(x, z, p);
 free (z); 
 return u;
}
 */
// Calculate x'ay  matrix
void xTay(double *x, double *A, double *y, int *p, double *xTAy)
{
      
 double u=0.0;
 int i, j, p1;

 p1 = *p;
 for (i=0; i< p1; i++) {
   for (j=0; j< p1; j++) {
     u += y[i] * A[i*p1+j] * x[j];
   }
 }
 *xTAy = u;
 return;
 }

double xTay2(double *x, double *A, double *y, int p)
{
 double u=0.0;
 int i, j;

 for (i=0; i<p; i++) {
   for (j=0; j<p; j++) {
     u += y[i] * A[i*p+j] * x[j];
   }
 }
 return u;
 }

/*
double xTay2(double x[], double A[], double y[], int p)
{
 double u=0.0;
 int i, j;

 for (i=0; i<p; i++) {
   for (j=0; j<p; j++) {
     u += y[i] * A[i*p+j] * x[j];
   }
 }
 return u;
 }
 */


/***************************************************************************/
// Determinant of a Square Matrix
void determinant_val(double mat[], int *order, double *out)
{
     double det, p;
     p = *order;
     det = determinant(mat, p);
     *out = det;
     free(mat); 
     return;
}

// Determinant
double determinant(double mat[], int order)
{    
    if(order==1)
                return mat[0];
    int i;
    double d = 0;
    for(i=0; i< order; i++)
             d += (double)pow(-1,i)*mat[i]*
             determinant(submatrix(mat, 0, i,order -1), order-1 );
    return d;
}

// Submatrix
double* submatrix(double mat[], int m, int n, int order)
{
     double* sub = (double*)malloc(order*order*sizeof(double));
     int i,j;
     for(i = 0, j=0; i<(order+1)*(order+1) && j<order*order; i++)
     {
             if(i>=(order+1)*m && i< (order+1)*(m+1))
                               i = (order+1)*(m+1);
             if(i>=n)
                     if((i-n)%(order+1)==0)        
                            continue;
             *(sub + j++) = mat[i];
     }
//    free(sub);     
    return sub;     
}     

// Calculate the Inverse of a Square Matrix
// This need the functions "sq_rt", "tinv_mat", "trans_pose" and "mat_mul"
// Adjoint Matrix [Adj(X)] is the Transpose of the Co-factor matrix
// Inverse Matrix is = Adj(X)/det(X) 

void MInv(double *S, double *inv, int *pp, double *det)
{
  double *t1; double *t; double u;
  double p;
  p = *pp;
  t1 =(double*)malloc(sizeof(double)*p*p);
  t = (double*)malloc(sizeof(double)*p*p);

  sq_rt( S, t , p , &u);
  *det = u;
  tinv_mat( t , t1 , p); 
  trans_pose(t1, t, p , p); 
  mat_mul(t, t1, inv , p+0.0 , p+0.0 , p+0.0); 
  *det = u;

free(t1);
free(t);
}

// This is used in "invMatrix"
void sq_rt(double *a, double *t, int p, double *det)
{
 int i , j , k ;
 double sum, sum1, u, v;
  for ( i = 0 ; i < p ; i ++)
      {
       for ( j = 0 ; j <p ; j++)
       t[i*p +j] = 0.0 ;  
      } 
  t[0] = sqrt(a[0]); 
 if ( p > 1 ) {
  
  v = exp(-25.0);
  t[p] = a[p] / t[0];
  u =  a[p+1]  - sqr(t[p] );
  if  ( u<v ) {
//         printf("Exited from sqrt: Matrix not pd\n");
         Rprintf("C - Error1: Exited from sqrt: Matrix is not positive definite \n");
         //;
         //exit(19);
   } else  t[p+1] = sqrt(u);
 
  for (i=2 ; i<p ; i++)
    {  
      t[i*p] =  a[i*p] /  t[0] ;
      for (j=1; j < p - 1 ; j++)
          {
            if ( j < i ) {
            sum1 = 0.0 ;
            for (k=0; k < j ; k++)
            sum1 += t[i*p + k] *  t[j*p + k]; 
            t[i*p+j] = (a[i*p + j]  - sum1) / t[j*p + j];
            }
          }
    
      sum = 0.0 ;
      for ( k = 0 ; k < i ; k ++)
      sum += sqr( t[i*p  + k] );  
      u = a[i*p+i] - sum;
         if ( u<v ) {
//                    printf("Exited from sqrt: Matrix not PD\n");
                    Rprintf("C - Error2: Failed to find determinant: Matrix not positive definite\n");
                    //;
                    //exit(19);
         } 
         else {
             t[i*p+i] = sqrt(u);
         }
    } 

 } // If loop for p > 1 
 *det = 1.0 ;
 for ( i = 0 ; i < p ; i++)
 *det = *det * sqr( t[i*p + i] ); 
}
// This is used in "invMatrix"    
void tinv_mat (double *t, double *c, int p)
{
int i , j , k;
double u;
   for ( i = 0; i <p ; i++)
       {
         for ( j = 0;  j < p ; j++)
         c[i * p  + j] = 0.0;
       }
 
  for ( i = 0 ; i < p ; i++)
  c[i * p  + i] = 1.0 / t[i * p  + i];

 if ( p > 1) {
  c[p] = -t[p]* c[0]/ t[p + 1];

 for ( i=2 ; i <p ; i ++)
 {
  for ( j = 0 ; j < p-1  ; j ++)
  {
     if ( j < i ) {
      u = 0.0 ;
      for ( k = j ; k < i ; k++)
      u = u + t[i * p  + k] * c[k * p  + j]; 
      c[i * p  + j] = - u / t[i * p  + i];
     }
   }
 }
} //  p  > 1 loop  
}
// This is used in "invMatrix"
void trans_pose (double *a, double *b, int m, int n)
{
 int i , j; 
 for ( i = 0 ; i < m ; i++)
 {
  for ( j = 0 ; j < n ; j ++)
    {  b[j * m  + i]  = a[i * n  + j] ; }
  }
}
// This is used in "invMatrix"
void mat_mul(double *ma1, double *ma2, double *mul, double no, double po, double qo)
{
  int i, j, k, n, p, q;
  n = no; p = po; q = qo;
   for( i=0; i < n; i++)
   {
      for( j=0; j < q; j++)
      {  *(mul + i * q + j) = 0.0;   }
   }
   for( i=0; i < n; i++)
   {
       for( j=0; j < q; j++)
       {
           for( k=0; k < p; k++)
             *(mul + i * q + j) = *(mul + i * q + j) + 
                   *(ma1 + i * p + k) * *(ma2 + k * q + j);
       }
   }
}

// Calculate the matrix product of X (q x p) and Y (p x r), 
void MProd2(double *ma1, double *ma2, double *mul, double *no,
     double *po, double *qo)
{ 
  int i, j, k, n, p, q;
  n = *no; p = *po; q = *qo;
   for( i=0; i < n; i++)
   {
      for( j=0; j < q; j++)
      {  *(mul + i * q + j) = 0.0;   }
   }
   for( i=0; i < n; i++)
   {
       for( j=0; j < q; j++)
       {
           for( k=0; k < p; k++)
             *(mul + i * q + j) = *(mul + i * q + j) + 
                   *(ma1 + i * p + k) * *(ma2 + k * q + j);
       }
   }
  return;
}

/*
// Inverse using LAPACK

#include <cstdio>

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

void inverse(double* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
}

int main(){

    double A [2*2] = {
        1,2,
        3,4
    };

    inverse(A, 2);

    printf("%f %f\n", A[0], A[1]);
    printf("%f %f\n", A[2], A[3]);

    return 0;
}

*/


///////////////////////////////////////////////////

/*
void mvrloop(double *x, int *xrow, int *xcol, double *y, double *out)
{
 int nn, pp, i, j;
 double tmp;

 nn = *xrow;
 pp = *xcol;

 for(i = 0; i < nn; i++) {
   for(j =0; j < pp; j++) {
    tmp = x[i*pp+j] - y[i*pp+j];
    out[i*pp+j] = tmp;
   }
 }  
 return;
}
*/
/*************************************************************************/

