//** File Name 'covariance.c' **//

#include "covariance.h"
#include "common.h"
#include "mathematics.h"
#include "randgenerator.h"


// common covariance format
void covFormat(int *cov, int *n, double *phi, double *nu, double *d, 
     double *sig_eta, double *S, double *det, double *Sinv, double *Qeta)
{
// exponential covariance
   if(cov[0] == 1){
    covExpo(n, phi, d, sig_eta, S, det, Sinv, Qeta);
   }
// gaussian covariance
   if(cov[0] == 2){
//      Rprintf("   phi: %4.4f \n",phi);
    covGaus(n, phi, d, sig_eta, S, det, Sinv, Qeta);
   }
// spherical covariance
   if(cov[0] == 3){
    covSphe(n, phi, d, sig_eta, S, det, Sinv, Qeta);
   }
// matern covariance
   if(cov[0] == 4){
    covMatern(n, phi, nu, d, sig_eta, S, det, Sinv, Qeta);
   }
   return;
}

// common covariance format
void covFormat2(int *cov, int *n, double *phi, double *nu, double *d, 
     double *sig_eta, double *det, double *Qeta)
{
   int n1;
   n1 = *n; 
   double *S;
   S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));

// exponential covariance
   if(cov[0] == 1){
    covExpo(n, phi, d, sig_eta, S, det, S, Qeta);
   }
// gaussian covariance
   if(cov[0] == 2){
    covGaus(n, phi, d, sig_eta, S, det, S, Qeta);
   }
// spherical covariance
   if(cov[0] == 3){
    covSphe(n, phi, d, sig_eta, S, det, S, Qeta);
   }
// matern covariance
   if(cov[0] == 4){
    covMatern(n, phi, nu, d, sig_eta, S, det, S, Qeta);
   }
   free(S);
   return;
}
     
// Exponential covariance 
void covExpo(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta)
{
  int i, n1;
  n1 = *n; 
  double *Q, *det1; 
  Q = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
  det1 = (double *) malloc((size_t)((1)*sizeof(double)));

  for(i = 0; i < (n1*n1); i++){
      S[i] = exp(-1.0*phi[0]*d[i]);
      Q[i] = sig2eta[0] * exp(-1.0*phi[0]*d[i]);
  }
  MInv(S, Sinv, n, det);
  MInv(Q, Qeta, n, det1);

  free(Q); free(det1);
  return;
}
// Gaussian covariance 
void covGaus(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta)
{
  int i, n1;
  n1 = *n; 
  double *Q, *det1; 
  Q = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
  det1 = (double *) malloc((size_t)((1)*sizeof(double)));

  for(i = 0; i < (n1*n1); i++){
      S[i] = exp(-1.0*phi[0]*phi[0]*d[i]*d[i]);
      Q[i] = sig2eta[0] * exp(-1.0*phi[0]*phi[0]*d[i]*d[i]);
  }
  MInv(S, Sinv, n, det);
  MInv(Q, Qeta, n, det1);

  free(Q); free(det1);
  return;
}
// Spherical covariance 
void covSphe(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta)
{
  int i, n1;
  n1 = *n; 
  double *Q, *det1; 
  Q = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
  det1 = (double *) malloc((size_t)((1)*sizeof(double)));

  for(i = 0; i < (n1*n1); i++){
    if(d[i] > 0 && d[i] <= 1.0/phi[0]){
      S[i] = 1.0-1.5*phi[0]*d[i]+0.5*(phi[0]*d[i])*(phi[0]*d[i])*(phi[0]*d[i]);
      Q[i] = sig2eta[0]*(1.0-1.5*phi[0]*d[i]+0.5*(phi[0]*d[i])*(phi[0]*d[i])*(phi[0]*d[i]));
    }
    else if(d[i] >= 1.0/phi[0]){
      S[i] = 0.0;
      Q[i] = 0.0;
    }
    else{
      S[i] = 1.0;
      Q[i] = sig2eta[0];
    }        
  }
  MInv(S, Sinv, n, det);
  MInv(Q, Qeta, n, det1);

  free(Q); free(det1);
  return;
}
// Matern covariance, nu = 3/2
void covMatern(int *n, double *phi, double *nu, double *d, double *sig2eta, 
     double *S, double *det, double *Sinv, double *Qeta)
{
  int i, n1;
  n1 = *n; 
  double *Q, *det1; 
  Q = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
  det1 = (double *) malloc((size_t)((1)*sizeof(double)));

  for(i = 0; i < (n1*n1); i++){
    if(d[i] > 0 ){
      S[i] = pow(d[i]*phi[0], nu[0])/(pow(2, nu[0]-1)*gammafn(nu[0]))*bessel_k(d[i]*phi[0], nu[0], 1.0);
      Q[i] = sig2eta[0]*(pow(d[i]*phi[0], nu[0])/(pow(2, nu[0]-1)*gammafn(nu[0]))*bessel_k(d[i]*phi[0], nu[0], 1.0));
    }  
    else{
      S[i] = 1.0;
      Q[i] = sig2eta[0];
    }        
  }
  MInv(S, Sinv, n, det);
  MInv(Q, Qeta, n, det1);

  free(Q); free(det1);
  return;
}




///////////////////////////////////////////////////////////////////////////////
// covariance format
void covF(int *cov, int *n1, int *n2, double *phi, double *nu, double *d, double *S)
{
// n1 = number of rows and n2 = numner of cols
// d = sould be n1 x n2     
// exponential covariance
   if(cov[0] == 1){
    covExp(n1, n2, phi, d, S);
   }
// gaussian covariance
   if(cov[0] == 2){
    covGau(n1, n2, phi, d, S);
   }
// spherical covariance
   if(cov[0] == 3){
    covSph(n1, n2, phi, d, S);
   }
// matern covariance
   if(cov[0] == 4){
    covMat(n1, n2, phi, nu, d, S);
   }
   return;    
}
void covExp(int *n1, int *n2, double *phi, double *d, double *S)
{
  int i,n11,n22;
  n11=*n1;
  n22=*n2;
//  int n1 = (int) round(sizeof(&d)/sizeof(double)+0.1);
//  printf("pointer: %i\n", n1);
//  printf("size: %d\n", sizeof(*d)/sizeof(double));
  for(i = 0; i < (n11*n22); i++){
      S[i] = exp(-1.0*phi[0]*d[i]);
  }
  return;
}
void covGau(int *n1, int *n2, double *phi, double *d, double *S)
{
  int i,n11,n22;
  n11=*n1;
  n22=*n2;
  for(i = 0; i < (n11*n22); i++){
      S[i] = exp(-1.0*phi[0]*phi[0]*d[i]*d[i]);
  }
  return;
}
void covSph(int *n1, int *n2, double *phi, double *d, double *S)
{
  int i,n11,n22;
  n11=*n1;
  n22=*n2;
  for(i = 0; i < (n11*n22); i++){
    if(d[i] > 0 && d[i] <= 1.0/phi[0]){
      S[i] = 1.0-1.5*phi[0]*d[i]+0.5*(phi[0]*d[i])*(phi[0]*d[i])*(phi[0]*d[i]);
    }
    else if(d[i] >= 1.0/phi[0]){
      S[i] = 0.0;
    }
    else{
      S[i] = 1.0;
    }        
  }
  return;
}
void covMat(int *n1, int *n2, double *phi, double *nu, double *d, double *S)
{
  int i,n11,n22;
  n11=*n1;
  n22=*n2;
  for(i = 0; i < (n11*n22); i++){
    if(d[i]*phi[0] > 0.0)
      S[i] = pow(d[i]*phi[0], nu[0])/(pow(2, nu[0]-1)*gammafn(nu[0]))*bessel_k(d[i]*phi[0], nu[0], 1.0);
    else
      S[i] = 1.0;
  }
}



//////////////////////////////////////////////////////////////////////


/*

// Matern covariance, nu = 3/2, for MH phi
void covMatern322(int *n, double *phi, double *d, double *det, double *Sinv)
{
  int i, n1;
  n1 = *n; 
  double *S; 
  S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));

  for(i = 0; i < (n1*n1); i++){
    if(d[i] > 0 ){
      S[i] = (1.0+phi[0]*d[i])*exp(-1.0*phi[0]*d[i]);
    }  
    else{
      S[i] = 0.0;
    }   
  }
  MInv(S, Sinv, n, det);

  free(S);
  return;
}


// Exponential covariance for MH phi 
void covExpo2(int *n, double *phi, double *d, double *det, double *Sinv)
{
  int i, n1;
  n1 = *n; 
  double *S; 
  S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));

  for(i = 0; i < (n1*n1); i++){
      S[i] = exp(-1.0*phi[0]*d[i]);
  }
  MInv(S, Sinv, n, det);

  free(S);
  return;
}

// Gaussian covariance for MH phi
void covGaus2(int *n, double *phi, double *d, double *det, double *Sinv)
{

  int i, n1;
  n1 = *n; 
  double *S;
  S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));

  for(i = 0; i < (n1*n1); i++){
      S[i] = exp(-1.0*phi[0]*phi[0]*d[i]*d[i]);
  }
  MInv(S, Sinv, n, det);

  free(S);
  return;
}

// Spherical covariance for MH phi
void covSphe2(int *n, double *phi, double *d, double *det, double *Sinv)
{
  int i, n1;
  n1 = *n; 
  double *S; 
  S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));

  for(i = 0; i < (n1*n1); i++){
    if(d[i] > 0 && d[i] <= 1.0/phi[0]){
      S[i] = 1.0-1.5*phi[0]*d[i]+0.5*(phi[0]*d[i])*(phi[0]*d[i])*(phi[0]*d[i]);
    }
    else if(d[i] >= 1.0/phi[0]){
      S[i] = 0.0;
    }
    else{
      S[i] = 1.0;
    }        
  }
  MInv(S, Sinv, n, det);

  free(S);
  return;
}

*/


/*
void testing(int int_array, double double_array)
{
    const char * format = "The %s array has %d bytes and %d elements. Size of %s is %d \n";
    printf (format, "int",
            sizeof (int_array), sizeof (int_array) / sizeof (int), "int", sizeof (int));
    printf (format, "double",
            sizeof (double_array), sizeof (double_array) / sizeof (double), "double", sizeof (double));
    // This method works even when the size of the array is not
    //   specified in the initializer. 
//    printf (format, "char",
//            sizeof (char_array), sizeof (char_array) / sizeof (char));
    return 0;
}

double len(double *in)
{
     int count=0;
     int x=0;
     while (in[count]!='\0')
     {
     count++;
     }
     printf("pointer: %d\n", count);

for(x=0;x<=count-1;x++)
{printf("\n%d\n",in[x]);}

printf("\n---------------------------------------------\n");
printf("\nnumber of elements in array ....> %d\n",x);
printf("\n-----------------------.............---------\n");

     return x;
}

*/


/*

// Covariance functions
// Exponential covariance 
void covExpo1(int *n, double *phi, double *d, double *S)
{
  int i, n1;
  n1 =*n;
//  int n1 = (int) round(sizeof(&d)/sizeof(double)+0.1);
//  printf("pointer: %i\n", n1);
//  printf("size: %d\n", sizeof(*d)/sizeof(double));
  for(i = 0; i < n1; i++){
      S[i] = exp(-1.0*phi[0]*d[i]);
  }
  return;
}

// Gaussian covariance 
void covGaus1(int *n, double *phi, double *d, double *S)
{
  int i, n1;
  n1 = *n; 
  for(i = 0; i < (n1*n1); i++){
      S[i] = exp(-1.0*phi[0]*phi[0]*d[i]*d[i]);
  }
  return;
}

// Spherical covariance 
void covSphe1(int *n, double *phi, double *d, double *S)
{
  int i, n1;
  n1 = *n; 
  for(i = 0; i < (n1*n1); i++){
    if(d[i] > 0 && d[i] <= 1.0/phi[0]){
      S[i] = 1.0-1.5*phi[0]*d[i]+0.5*(phi[0]*d[i])*(phi[0]*d[i])*(phi[0]*d[i]);
    }
    else if(d[i] >= 1.0/phi[0]){
      S[i] = 0.0;
    }
    else{
      S[i] = 1.0;
    }        
  }
  return;
}


// Matern covariance, nu = 3/2
void covMatern321(int *n, double *phi, double *d, double *S)
{
  int i, n1;
  n1 = *n; 
  for(i = 0; i < (n1*n1); i++){
    if(d[i] > 0 ){
      S[i] = (1.0+phi[0]*d[i])*exp(-1.0*phi[0]*d[i]);
    }  
    else{
      S[i] = 0.0;
    }           
  }
  return;
}

// Matern covariance
void covMatern1(int *n, double *phi, double *d, double *S)
{
  int i, n1;
  n1 = *n; 
  double tmp[1];
  double u =0.0;
  double v =0.0;
  for(i = 0; i < (n1*n1); i++){
     if (d[i]>0.00001) { 
       u = 2.0 * d[i] * phi[0];
       v =  u * bessk1(u); 
       tmp[0] = v;
     }
     else {
       v = 1.0;
       tmp[0] = v;
     }   
     S[i] = tmp[0];
  }
  return;
}




// exponential
// The "Q_eta" function, where Q_eta = Inverse of Sigma_eta
void Q_eta(double *phi, double *sig_eta, int *n, double *d, double *Qeta)
{
     double *S, *det;
     int j, n1, n2;
     n1 = *n;
     n2 = n1*n1;

     S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     det = (double *) malloc((size_t)((n2)*sizeof(double)));
     for(j=0; j< n2; j++) {
           S[j] = sig_eta[0] * exp(-(d[j]*phi[0]));
     } 
     MInv(S, Qeta, n, det);

     free(S); free(det);
     return;
}

// exponential
// The "Sinv" fnction
void S_inv(double *phi, double *d, int *n, int *constant, double *Sinv)
{
     double *S, *det, col;
     int j, n1, nn;
     n1 = *n;
     nn = n1*n1;
     col = *constant;

     S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     det = (double *) malloc((size_t)((col*col)*sizeof(double)));
     for(j=0; j< nn; j++) {
           S[j] = exp(-(d[j]*phi[0]));
     } 
     MInv(S, Sinv, n, det);

     free(S); free(det);
     return;
}

// exponential
// The "Sinv" fnction with determinant output
void S_inv_det(double *phi, double *d, int *n, int *constant, double *det, 
     double *Sinv)
{
     double *S, col;
     int j, n1, nn;
     n1 = *n;
     nn = n1*n1;
     col = *constant;

     S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     for(j=0; j< nn; j++) {
           S[j] = exp(-(d[j]*phi[0]));
     } 
     MInv(S, Sinv, n, det);

     free(S);
     return;
}

*/

////////////////////////// THE END ////////////////////////////////


/*

// common covariance format
void covFormat2(int *cov, int *n, double *phi, double *d, double *det, 
     double *Sinv)
{
     // exponential covariance
   if(cov[0] == 1){
    covExpo2(n, phi, d, det, Sinv);
   }
// gaussian covariance
   if(cov[0] == 2){
    covGaus2(n, phi, d, det, Sinv);
   }
// spherical covariance
   if(cov[0] == 3){
    covSphe2(n, phi, d, det, Sinv);
   }
// matern covariance, nu = 3/2
   if(cov[0] == 4){
    covMatern322(n, phi, d, det, Sinv);
   }
// matern covariance
   if(cov[0] == 5){
    covMatern2(n, phi, d, det, Sinv);
   }
   
   return;
}


// covariance format
void covFormat1(int *cov, int *n, double *phi, double *d, double *S)
{
// exponential covariance
   if(cov[0] == 1){
    covExpo1(n, phi, d, S);
   }
// gaussian covariance
   if(cov[0] == 2){
    covGaus1(n, phi, d, S);
   }
// spherical covariance
   if(cov[0] == 3){
    covSphe1(n, phi, d, S);
   }
// matern covariance, nu = 3/2
   if(cov[0] == 4){
    covMatern321(n, phi, d, S);
   }
// matern covariance
   if(cov[0] == 5){
    covMatern1(n, phi, d, S);
   }
   return;    
}


*/

/*
// Matern covariance
void covMatern_(int *n, double *phi, double *d, double *sig2eta, double *S, 
     double *det, double *Sinv, double *Qeta)
{
  int i, n1;
  n1 = *n; 
  double *Q, *det1, tmp[1]; 
  Q = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
  det1 = (double *) malloc((size_t)((1)*sizeof(double)));

  double u =0.0;
  double v =0.0;
  for(i = 0; i < (n1*n1); i++){
     if (d[i]>0.00001) { 
       u = 2.0 * d[i] * phi[0];
       v =  u * bessk1(u); 
       tmp[0] = v;
//       Rprintf("u=%f v=%f bess=%f\n",u, v, v/u); 
//     if (v>1.0) 
//      Rprintf("V exceeded 1! u=%f v=%f bess=%f\n",u, v, v/u); 
     }
     else {
       v = 1.0;
       tmp[0] = v;
     }   
     S[i] = tmp[0];
     Q[i] = sig2eta[0]*tmp[0];
  }
  MInv(S, Sinv, n, det);
  MInv(Q, Qeta, n, det1);

  free(Q); free(det1);
  return;
}


// Matern covariance for MH phi
void covMatern2(int *n, double *phi, double *d, double *det, double *Sinv)
{
  int i, n1;
  n1 = *n; 
  double *S, tmp[1]; 
  S = (double *) malloc((size_t)((n1*n1)*sizeof(double)));

  double u =0.0;
  double v =0.0;
  for(i = 0; i < (n1*n1); i++){
     if (d[i]>0.00001) { 
       u = 2.0 * d[i] * phi[0];
       v =  u * bessk1(u); 
       tmp[0] = v;
     }
     else {
       v = 1.0;
       tmp[0] = v;
     }   
     S[i] = tmp[0];
  }
  MInv(S, Sinv, n, det);

  free(S);
  return;
}

/////

// Bessle function for Matern covariance
double bessk1(double x)
{
	double bessi1(double x);
	double y,ans=0.0;
   
        if (x>0.000001) { 
	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
        }
	return ans;
}

double bessi1 (double x) { 
  float ax, ans; 
  double y; 

  if ((ax=fabs(x)) < 3.75) { 
    y = x/3.75; 
    y *=y; 
    ans  = ax * (0.5 + y * (0.87890594 + y* (0.51498869 + y*(0.15084934
							     + y * (0.2658733e-1 + y*(0.301532e-2+y*(0.32411e-3))))))); 
  } else { 
    y = 3.75/ax; 
    ans = 0.2282967e-1 + y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans = 0.39894228 + y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *=(exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}  

 
*/
