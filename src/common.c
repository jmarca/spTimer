//** File Name 'common.c' **//

#include "common.h"
#include "math.h"
#include "mathematics.h"
#include "randgenerator.h"
//#include "Print.h"



// ratio function for the phi sample
void ratio_fnc(double *ratio, int *constant, double *U)
{
     if( *ratio < 1.0){
         ratio[0] = ratio[0];
     }
     else {
          ratio[0] = 1.0;
     }     
     runif_val(constant, constant, U);     

     return;
}

// sort data to get 4th highest maximum
void sort_4th(double *sample, int *n, int *r, int *T, double *an4th)
{
	int i, l, t, n1, r1, T1;
	n1 =*n; r1 =*r; T1 =*T;

//	double tmp[T1];
	double *tmp;
	tmp = (double *) malloc((size_t)((T1)*sizeof(double)));
	
	for(i=0; i<n1; i++){
	    for(l=0; l<r1; l++){
	        for(t=0; t<T1; t++){
                tmp[t] = sample[t+l*T1+i*T1*r1];
		}
		qsort(tmp, T1, sizeof(double), sort_fnc);
		for(t=0; t<T1; t++){
		tmp[t] = tmp[t];
		}
		an4th[l+i*r1] = tmp[T1-4];
		//for(t=0; t<T1; t++){
		//an4th[t+l*T1+i*T1*r1] = tmp[t];
		//}
	    }
	}
	free(tmp);
	return;
}


int sort_fnc(const void *x, const void *y)
{
//	if(*x < *y) return -1;
//	if(*x == *y) return -0;
//	if(*x > *y) return 1;

	return (*(double*)x - *(double*)y); 
// This x and y checks for each pair of elements considering the elements as x and y. If x is found greater than y, then x goes before y otherwise x goes after y. If x = y then it remains on the same position. 
}


// to print in the R interface for GP models
void GP_para_printR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *beta) 
{
    int j, k;
    double phi1, sig2e1, sig2eta1, ii;
    phi1 = *phi;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;

    double num =  (iteration/report); 
    int intpart = (int)num;
  
    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f\n", 
      phi1, sig2e1, sig2eta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}    


// to print in the R interface for GP models
void GP_para_printRnu (int i, int iteration, int report, int p, double accept, 
     double *phi, double *nu, double *sig2e, double *sig2eta, double *beta) 
{
    int j, k;
    double phi1, nu1, sig2e1, sig2eta1, ii;
    phi1 = *phi;
    nu1 =*nu;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, nu: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f\n", 
      phi1, nu1, sig2e1, sig2eta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}


// to print in the R interface for GP models, for spatial beta
void GPsp_para_printR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *sig2beta, double *beta) 
{
    int j, k;
    double phi1, sig2e1, sig2eta1, sig2beta1, ii;
    phi1 = *phi;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;
    sig2beta1 = *sig2beta;    

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f, sig2beta: %4.4f\n", 
      phi1, sig2e1, sig2eta1, sig2beta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
      Rprintf(" ## Model used spatially varying parameters \n");
      Rprintf(" ## Spatially varying beta parameters are omitted in the display ");
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}    



// to print in the R interface for GP models, for spatial beta
void GPsp_para_printRnu (int i, int iteration, int report, int p, double accept, 
     double *phi, double *nu, double *sig2e, double *sig2eta, double *sig2beta, 
     double *beta) 
{
    int j, k;
    double phi1, nu1, sig2e1, sig2eta1, sig2beta1, ii;
    phi1 = *phi;
    nu1 =*nu;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;
    sig2beta1 = *sig2beta;

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, nu: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f, sig2beta: %4.4f\n", 
      phi1, nu1, sig2e1, sig2eta1, sig2beta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
      Rprintf(" ## Model used spatially varying parameters \n");
      Rprintf(" ## Spatially varying beta parameters are omitted in the display ");
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}


// to print in the R interface for GP models, for temporal beta
void GPtp_para_printR (int i, int iteration, int report, int p, int u, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *sig2delta, double *sig20,
     double *rho, double *beta) 
{
    int j, k;
    double phi1, sig2e1, sig2eta1, sig2delta1, sig201, ii;
    phi1 = *phi;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;
    sig2delta1 = *sig2delta;    
    sig201 = *sig20;    
    
    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f,\n                sig2delta: %4.4f, sig2op: %4.4f,\n", 
      phi1, sig2e1, sig2eta1, sig2delta1, sig201);
      for(k=0; k<u; k++){
      Rprintf("   rho[%d]: %4.4f", k+1, rho[k]);
      }
      Rprintf("\n");
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
      Rprintf(" ## Model used temporally varying dynamic parameters \n");
      Rprintf(" ## Dynamic beta parameters are omitted in the display ");
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}    



// to print in the R interface for GP models, for temporal beta
void GPtp_para_printRnu (int i, int iteration, int report, int p, int u, double accept, 
     double *phi, double *nu, double *sig2e, double *sig2eta, double *sig2delta, 
     double *sig20, double *rho, double *beta) 
{
    int j, k;
    double phi1, nu1, sig2e1, sig2eta1, sig2delta1, sig201, ii;
    phi1 = *phi;
    nu1 =*nu;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;
    sig2delta1 = *sig2delta;    
    sig201 = *sig20;    

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, nu: %4.4f, sig2eps: %4.4f,\n   sig2eta: %4.4f, sig2delta: %4.4f, sig2op: %4.4f,\n", 
      phi1, nu1, sig2e1, sig2eta1, sig2delta1, sig201);
      for(k=0; k<u; k++){
      Rprintf("   rho[%d]: %4.4f", k+1, rho[k]);
      }
      Rprintf("\n");
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
      Rprintf(" ## Model used temporally varying dynamic parameters \n");
      Rprintf(" ## Dynamic beta parameters are omitted in the display ");
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}




// to print in the R interface for GP models, for temporal beta
void GPsptp_para_printR (int i, int iteration, int report, int p, int u, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *sig2beta, double *sig2delta, 
     double *sig20, double *rho, double *beta) 
{
    int j, k;
    double phi1, sig2e1, sig2eta1, sig2beta1, sig2delta1, sig201, ii;
    phi1 = *phi;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;
    sig2beta1 = *sig2beta;
    sig2delta1 = *sig2delta;    
    sig201 = *sig20;    
    
    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f, sig2beta: %4.4f,\n                sig2delta: %4.4f, sig2op: %4.4f,\n", 
      phi1, sig2e1, sig2eta1, sig2beta1, sig2delta1, sig201);
      for(k=0; k<u; k++){
      Rprintf("   rho[%d]: %4.4f", k+1, rho[k]);
      }
      Rprintf("\n");
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
      Rprintf(" ## Model used spatially and temporally varying dynamic parameters \n");
      Rprintf(" ## Spatial and dynamic beta parameters are omitted in the display ");
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}    



// to print in the R interface for GP models, for temporal beta
void GPsptp_para_printRnu (int i, int iteration, int report, int p, int u, double accept, 
     double *phi, double *nu, double *sig2e, double *sig2eta, double *sig2beta, 
     double *sig2delta, double *sig20, double *rho, double *beta) 
{
    int j, k;
    double phi1, nu1, sig2e1, sig2eta1, sig2beta1, sig2delta1, sig201, ii;
    phi1 = *phi;
    nu1 =*nu;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;
    sig2beta1 = *sig2beta;
    sig2delta1 = *sig2delta;    
    sig201 = *sig20;    

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, nu: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f,\n   sig2beta: %4.4f, sig2delta: %4.4f, sig2op: %4.4f,\n", 
      phi1, nu1, sig2e1, sig2eta1, sig2beta1, sig2delta1, sig201);
      for(k=0; k<u; k++){
      Rprintf("   rho[%d]: %4.4f", k+1, rho[k]);
      }
      Rprintf("\n");
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
      Rprintf(" ## Model used spatially and temporally varying dynamic parameters \n");
      Rprintf(" ## Spatial and dynamic beta parameters are omitted in the display ");
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}



// to print in the R interface for GPP models, for spatial beta
void GPPsp_para_printRnu (int i, int iteration, int report, int p, double accept, 
     double *phi, double *nu, double *rho, double *sig2e, double *sig2eta, 
     double *sig2beta, double *beta) 
{
    int j, k;
    double phi1, nu1, rho1, sig2e1, sig2eta1, sig2beta1, ii;
    phi1 = *phi;
    nu1 =*nu;
    rho1 =*rho;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;
    sig2beta1 = *sig2beta;

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, nu: %4.4f, rho: %4.4f\n",  phi1, nu1, rho1);
      Rprintf("   sig2eps: %4.4f, sig2eta: %4.4f, sig2beta: %4.4f\n", sig2e1, sig2eta1, sig2beta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
      Rprintf(" ## Model used spatially varying parameters \n");
      Rprintf(" ## Spatially varying beta parameters are omitted in the display ");
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}
void GPPsp_para_printR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *rho, double *sig2e, double *sig2eta, 
     double *sig2beta, double *beta) 
{
    int j, k;
    double phi1, rho1, sig2e1, sig2eta1, sig2beta1, ii;
    phi1 = *phi;
    rho1 =*rho;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;
    sig2beta1 = *sig2beta;

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
  	  ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, rho: %4.4f\n",  phi1, rho1);
      Rprintf("   sig2eps: %4.4f, sig2eta: %4.4f, sig2beta: %4.4f\n", sig2e1, sig2eta1, sig2beta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
      Rprintf(" ## Model used spatially varying parameters \n");
      Rprintf(" ## Spatially varying beta parameters are omitted in the display ");
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}

// to print in the R interface for AR and GPP
void para_printR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *rho, double *sig2e, double *sig2eta, double *beta) 
{
    int j, k;
    double phi1, rho1, sig2e1, sig2eta1, ii;
    phi1 = *phi;
    rho1 = *rho;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, rho: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f\n", 
      phi1, rho1, sig2e1, sig2eta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}    

// to print in the R interface for AR and GPP
void para_printRnu (int i, int iteration, int report, int p, double accept, 
     double *phi, double *nu, double *rho, double *sig2e, double *sig2eta, double *beta) 
{
    int j, k;
    double phi1, nu1, rho1, sig2e1, sig2eta1, ii;
    phi1 = *phi;
    nu1=*nu;
    rho1 = *rho;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, nu: %4.4f, rho: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f\n", 
      phi1, nu1, rho1, sig2e1, sig2eta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}    

void para_print_spTR (int i, int iteration, int report, int p, double accept, 
     double *phi, double *sig2e, double *sig2eta, double *beta) 
{
    int j, k;
    double phi1, sig2e1, sig2eta1, ii;
    phi1 = *phi;
    sig2e1 = *sig2e;
    sig2eta1 = *sig2eta;

    double num =  (iteration/report); 
    int intpart = (int)num;

    for(j=0; j<report; j++){
    if(i==(intpart*(j+1)-1)){
      ii = (double) i;
      Rprintf("---------------------------------------------------------------\n");
      Rprintf(" Sampled: %i of %i, %3.2f%%.\n Batch Acceptance Rate (phi): %3.2f%%\n", 
      i+1, iteration, 100.0*(i+1)/iteration, 100.0*(accept/ii));
      Rprintf(" Checking Parameters: \n");
      Rprintf("   phi: %4.4f, sig2eps: %4.4f, sig2eta: %4.4f\n", 
      phi1, sig2e1, sig2eta1);
      for(k=0; k<p; k++){
      Rprintf("   beta[%d]: %4.4f", k+1, beta[k]);
      }
      Rprintf("\n---------------------------------------------------------------\n");
    }
    }
    return;
}    


// to print in the R interface
void printR (int i, int iteration) 
{
    int j;
    double intp1, *intpart;
    intpart = (double *) malloc((size_t)((1)*sizeof(double)));   
    for(j=0; j<10; j++){
    modf((iteration/10), intpart);
    intp1 = *intpart;
    if(i==(intp1*(j+1)-1)){
      Rprintf("-------------------------------------------------\n");
      Rprintf("  Sampled: %i of %i, %3.2f%%\n", i+1, iteration, 100.0*(i+1)/iteration);
      Rprintf("-------------------------------------------------\n");
    }
    }
    free(intpart);
    return;
}    


// Used in summary statistics
void ext_sumstat(int i, int *its, double *x, double *alt)
{
     int j, its1;
     its1 = *its;
     for(j=0; j<its1; j++){
        alt[j] = x[j+i*its1];
     }
     return;
}         

// Used in summary statistics
void ext_sumstat_burnin(int i, int *its, int *burnin, double *x, double *alt)
{
     int j, its1, burn1;
     its1 = *its;
     burn1 = *burnin;
     for(j=0; j<(its1-burn1); j++){
        alt[j] = x[j+burn1+i*its1];
     }
     return;
}         


// Used in predictioin_xb_ar
void extn_12(int j, int *n, double *S_12, double *S_12c)
{
     int i, row;
     row=*n;
     for(i=0; i < row; i++) {
        S_12c[i]= S_12[i+j*row];
        }
     return;
}      

// betatp => T x p ??
// Xtp => N x p
// XB => N x 1
void comb_XB_tp(int *n, int *r, int *T, int *p, double *Xtp, double *betatp, 
     int *constant, double *XB)
{
     int t, l, n1, r1, T1, p1;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     p1 =*p;

     double *X1, *beta, *XB1;
     X1 = (double *) malloc((size_t)((n1*p1)*sizeof(double)));
     beta = (double *) malloc((size_t)((p1)*sizeof(double)));
     XB1 = (double *) malloc((size_t)((n1)*sizeof(double)));     
      
     for(l=0; l<r1; l++){
         for(t=0; t<T1; t++){     
             extract_X(t, l, n, r, T, p, Xtp, X1); // nrT x p into n x p
             extract_beta_t(t, T, p, betatp, beta); // p x T into p x 1
             MProd(beta, constant, p, X1, n, XB1); // n x 1
             put_together1(l, t, n, r, T, XB, XB1);
//  int i; 
//  for(i=0; i< n1*p1; i++){
//     Rprintf("  X1: %4.4f, \n", X1[i]);      
//  }
//  for(i=0; i< p1; i++){
//     Rprintf("  beta: %4.4f, \n", beta[i]);      
//  }
//  for(i=0; i< n1; i++){
//     Rprintf("  XB: %4.4f, \n", XB1[i]);      
//  }
         }
     }

     free(X1); free(beta); free(XB1);
     return;
}

// betasp => n x q
// Xsp => N x q
// XB => N x 1
void comb_XB_sp(int *n, int *r, int *T, int *q, double *Xsp, double *betasp, 
     int *constant, double *XB)
{
     int i, j, l, t, n1, r1, T1, q1;
     n1=*n;
     r1=*r;
     T1=*T;
     q1=*q;
     
     double *XB1, *dump, *I;
     XB1 = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     dump = (double *) malloc((size_t)((n1)*sizeof(double)));
     I = (double *) malloc((size_t)((n1)*sizeof(double)));
               
     for(l=0; l<r1; l++){
        for(t=0; t<T1; t++){
           for(i=0; i<n1; i++){
             I[i] = 0.0;
           }
           for(j=0; j<q1; j++){      
             extract_X_sp2(t, l, j, n, r, T, Xsp, XB1); // n x n
             for(i=0; i<n1; i++){
               dump[i] = betasp[i+j*n1]; // n x 1
             }
             MProd(dump, constant, n, XB1, n, dump); // n x 1 
             MAdd(I, n, constant, dump, I); // n x 1
           }
           put_together1(l, t, n, r, T, XB, I);
        }
     }
     free(XB1); free(dump); free(I);
     return;
}


// betasp => m x q
// Xsp => N x q
// XB => N x 1
void comb_XB_sp_gpp(int *n, int *m, int *r, int *T, int *q, double *Xsp, 
     double *betasp, double *A, int *constant, double *XB)
{
     int i, j, l, t, n1, m1, r1, T1, q1;
     n1=*n;
     m1=*m;
     r1=*r;
     T1=*T;
     q1=*q;
     
     double *XB1, *XA, *dump, *dump1, *I;
     XB1 = (double *) malloc((size_t)((n1*n1)*sizeof(double)));
     XA = (double *) malloc((size_t)((n1*m1)*sizeof(double)));
     dump = (double *) malloc((size_t)((m1)*sizeof(double)));
     dump1 = (double *) malloc((size_t)((n1)*sizeof(double)));
     I = (double *) malloc((size_t)((n1)*sizeof(double)));
               
     for(l=0; l<r1; l++){
        for(t=0; t<T1; t++){
           for(i=0; i<n1; i++){
             I[i] = 0.0;
           }
           for(j=0; j<q1; j++){      
             extract_X_sp2(t, l, j, n, r, T, Xsp, XB1); // n x n diagonal
             MProd(A, m, n, XB1, n, XA);  // n x m  
             for(i=0; i<m1; i++){
               dump[i] = betasp[i+j*m1]; // m x 1
             }
             MProd(dump, constant, m, XA, n, dump1); // n x 1 
             MAdd(I, n, constant, dump1, I); // n x 1
           }
           put_together1(l, t, n, r, T, XB, I);
        }
     }
     free(XB1); free(dump); free(dump1); free(I);
     return;
}



// Used in dlm and spT
// extract nrT x p matrix into n x p matrix
void extract_X(int t, int l, int *n, int *r, int *T, int *p, 
     double *x, double *alt)
{
     int i, j, k, p1, n1, r1, T1;
     p1 =*p;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     
     for(i=0; i<p1; i++){
     for(j=0; j<n1; j++){         
          k = j*r1*T1 + i*n1*r1*T1 + t + l*T1;      // Not able to use in R
          alt[j+i*n1] = x[k];
     }
     }
     return;
}

// extract nrT x p matrix into n x p matrix for unequal T
void extract_X_uneqT(int t, int l, int *n, int *r, int *T, int *rT, int *p, 
     double *x, double *alt)
{
     int i, j, k, p1, n1, r1, rT1;
     p1 =*p;
     n1 =*n;
     r1 =*r;
     rT1 =*rT;

     int *T1; 
     T1 = (int *) malloc((size_t)((r1+1)*sizeof(int)));

     cumsumint(r, T, T1);
      
     for(i=0; i<p1; i++){
     for(j=0; j<n1; j++){         
          k = j*rT1 + i*n1*rT1 + t + T1[l];      // Not able to use in R
          alt[j+i*n1] = x[k];
     }
     }
     
     free(T1);
     return;
}
          
// extract nrT x p matrix into n x n diagonal matrix
void extract_X_sp2(int t, int l, int j, int *n, int *r, int *T, 
     double *x, double *alt)
{
     int i, k, n1, r1, T1;
     n1 =*n;
     r1 =*r;
     T1 =*T;
     double *I; 
     I = (double *) malloc((size_t)((n1)*sizeof(double)));   
     for(i=0; i<n1; i++){
          k = i*r1*T1 + j*n1*r1*T1 + t + l*T1;      // Not able to use in R
          I[i] = x[k];
     }
     for(i=0; i<n1; i++){
       for(k=0; k< n1; k++){
         if(k==i){
           alt[k+n1*i] = I[i];
         }
         else{
           alt[k+n1*i] = 0.0;
         }
       }
     }
     free(I);
     return;
}            


// extract nrT x p matrix into n x 1 matrix
void extract_X_sp3(int t, int l, int j, int *n, int *r, int *T, 
     double *x, double *alt)
{
     int i, k, n1, r1, T1;
     n1 =*n;
     r1 =*r;
     T1 =*T;

     for(i=0; i<n1; i++){
          k = t + l*T1 + i*T1*r1 + j*T1*n1*r1;      // Not able to use in R
          alt[i] = x[k];
     }     

     return;
}            


// extract nrT x p matrix into 1 x 1 matrix
void extract_X_sp4(int t, int l, int i, int j, int *n, int *r, int *T, 
     double *x, double *alt)
{
     int o, n1, r1, T1;
     n1 =*n;
     r1 =*r;
     T1 =*T;

          o = t + l*T1 + i*T1*r1 + j*T1*n1*r1;      // Not able to use in R
          alt[0] = x[o];

     return;
}            

// extract nrT x p matrix into 1 x 1 matrix
void extract_X_sp(int t, int l, int i, int j, int *n, int *r, int *T, 
     double *x, double *alt)
{
     int k, n1, r1, T1;
     n1 =*n;
     r1 =*r;
     T1 =*T;

          k = t + l*T1 + i*T1*r1 + j*T1*n1*r1;      // Not able to use in R
          alt[0] = x[k];

     return;
}            

// extract n x q betasp into n x 1: used in prediction
void extract_beta_sp(int j, int *n, double *betasp, double *alt)
{
     int i, k, n1;
     n1 = *n;
     
     for(i=0; i<n1; i++){
          k = i + j*n1;
          alt[i] = betasp[k];
     }
     return;
}    

// extract n x q betasp into n x (q-1): used in beta estimation
void extract_beta_sp2(int j, int *n, int *q, double *betasp, double *alt)
{
     int i, k, n1, q1;
     n1 = *n;
     q1 = *q;
     
     int array[q1];
// define array     
     for(k=0; k<q1; k++){
         array[k] = k;
     }
// delete jth element from array
//     for(k=(j-1); k<(q1-1); k++){ // starting from j = 1
     for(k=(j); k<(q1-1); k++){  // starting from j = 0
         array[k] = array[k+1];
     }
// main
     for(k=0; k<(q1-1); k++){
         for(i=0; i<n1; i++){
           alt[i+k*n1] = betasp[i+array[k]*n1];
         }
     }
     return;
}    

// Used in fm
// extract p x T into p x 1 matrix
void extract_beta_t(int t, int *T, int *p, double *beta, double *alt)
{
     int i, k, p1;
     p1=*p;
     
     for(i=0; i< p1; i++){
        k = i + t*p1;
        alt[i] = beta[k];
     }      
     return;
}

// extract p x r into p x 1 matrix
void extract_beta_l(int l, int *r, int *p, double *beta, double *alt)
{
     int i, k, p1;
     p1=*p;
     
     for(i=0; i< p1; i++){
        k = i + l*p1;
        alt[i] = beta[k];
     }      
     return;
}
    
// Used in forecast_xb_ar
// extract the nrT x P covariates into nr x p matrix
void extract_X5(int t, int *n, int *r, int *T, int *p, double *x, double *alt)
{
     int i, j, l, k, n1, p1, r1, T1;
     n1 = *n;
     p1= *p;
     r1 = *r;
     T1 = *T;
     
     for(j=0; j < p1; j++){
     for(i=0; i < n1; i++){
     for(l=0; l < r1; l++){         
         k = l*T1 + i*r1*T1 + j*n1*r1*T1 + t;      // Not able to use in R
         alt[l+i*r1+j*r1*n1] = x[k];                 // for reuse in C/C++   
     }
     }
     }
     return;
}

// Used in forecast_xb_ar
// extract the nr x p covariates into 1 x p
// 'i' = for site , 'n'
// 'l' = for year , 'r'
void extract_X4(int i, int l, int *n, int *r, int *p, double *x, double *alt)
{
     int j, k, n1, p1, r1;
     n1 = *n;
     p1= *p;
     r1 = *r;
     
     for (j=0; j < p1; j++) {
         k = j * n1 * r1 + l + i * r1;      // Not able to use in R
         alt[j] = x[k];                 // for reuse in C/C++   
     }

     return;
}

// extract the nrT x p covariates into p x 1 matrix
void extract_X41(int l, int t, int i, int *n, int *rT, int *T, int *p, 
     double *x, double *alt)
{
     int j, k, n1, rT1, T1, p1, nrT;
     n1 = *n;
     rT1 = *rT;
     T1 = *T;
     p1= *p;
     nrT= n1*rT1;
     
     for (j=0; j <p1; j++) {
         k = j*nrT + t + l*T1 + i*rT1;
         alt[j] = x[k];                 // for reuse in C/C++   
     }
     return;
}

// for unequal T
// extract the nrT x p covariates into p x 1 matrix
void extract_X41_uneqT(int l, int t, int i, int *n, int *rT, int *r, int *T, 
     int *p, double *x, double *alt)
{
     int j, k, n1, rT1, r1, p1, nrT;
     n1 = *n;
     rT1 = *rT;
     r1 =*r;
     p1= *p;
     nrT= n1*rT1;

     int *T1; 
     T1 = (int *) malloc((size_t)((r1+1)*sizeof(int)));

     cumsumint(r, T, T1);   

     for (j=0; j <p1; j++) {
         k = j*nrT + t + T1[l] + i*rT1;
//         k = j*nrT + t + l*T1 + i*rT1;
         alt[j] = x[k];                 // for reuse in C/C++   
     }
     return;
}


// extract the nrT x p covariates into p x n matrix
void extract_X21(int l, int t, int *n, int *rT, int *T, int *p, 
     double *x, double *alt)
{
     int i, j, k, n1, rT1, T1, p1, nrT;
     n1 = *n;
     rT1 = *rT;
     T1 = *T;
     p1= *p;
     nrT= n1*rT1;
     
     for (j=0; j <n1; j++) {
     for (i=0; i < p1; i++) {
         k = i*nrT + t + l*T1 + j*rT1;
         alt[i+j*p1] = x[k];                 // for reuse in C/C++   
     }
     }
     return;
}

// extract the nrT x p covariates into p x n matrix for unequal T
void extract_X21_uneqT(int l, int t, int *n, int *rT, int *r, int *T, int *p, 
     double *x, double *alt)
{
     int i, j, k, n1, rT1, p1, r1, nrT;
     n1 = *n;
     r1 =*r;
     rT1 = *rT;
     p1= *p;
     nrT= n1*rT1;

     int *T1; 
     T1 = (int *) malloc((size_t)((r1+1)*sizeof(int)));

     cumsumint(r, T, T1);   
     
     for (j=0; j <n1; j++) {
     for (i=0; i < p1; i++) {
         k = i*nrT + t + T1[l] + j*rT1;
         alt[i+j*p1] = x[k];                 // for reuse in C/C++   
     }
     }
     
     free(T1);
     return;
}


// extract the nrT x p covariates into n x p matrix
void extract_X2(int l, int t, int *n, int *rT, int *T, int *p, 
     double *x, double *alt)
{
     int i, j, k, n1, rT1, T1, p1, nrT;
     n1 = *n;
     rT1 = *rT;
     T1 = *T;
     p1= *p;
     nrT= n1*rT1;
     
     for (j=0; j <p1; j++) {
     for (i=0; i < n1; i++) {
         k = i * rT1 + l * T1 + t + j * nrT;      // Not able to use in R
         alt[i+j*n1] = x[k];                 // for reuse in C/C++   
     }
     }
     return;
}


// extract the nrT x p covariates into n x 1 matrix
void extract_X3(int l, int t, int k, int *n, int *rT, int *T, int *p, 
     double *x, double *alt)
{
     int i, k1, n1, rT1, T1, nrT;
     n1 = *n;
     rT1 = *rT;
     T1 = *T;
     nrT= n1*rT1;
     
     for (i=0; i < n1; i++) {
         k1 = i * rT1 + l * T1 + t + k *nrT;      // Not able to use in R
         alt[i] = x[k1];                 // for reuse in C/C++   
     }

     return;
}


// extract the nrT x p covariates into n x 1 matrix
void extract_X3_uneqT(int l, int t, int k, int *n, int *r, int *rT, int *T, 
     int *p, double *x, double *alt)
{
     int i, k1, n1, r1, rT1, nrT;
     n1 = *n;
     r1 = *r;
     rT1 = *rT;
     nrT= n1*rT1;

     int *T1; 
     T1 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     
     cumsumint(r, T, T1);   
     
     for (i=0; i < n1; i++) {
         k1 = i * rT1 + T1[l] + t + k *nrT;      // Not able to use in R
//         k1 = i * rT1 + l * T1 + t + k *nrT;      // Not able to use in R
         alt[i] = x[k1];                 // for reuse in C/C++   
     }

     return;
}


// Extracts the l'th year and t'th time vector
// r is the years T is the time, "rT = r x T"
// n is the number of sites

void extract_alt(int *l, int *t, int *n, int *rT, int *T, 
     double *x, double *alt)
{
     int i, k, l1, t1, n1, rT1, T1;
     l1 = *l;
     t1 = *t;
     n1 = *n;
     rT1 = *rT;
     T1 = *T;
     
     for (i=0; i < n1; i++) {
         k = i * rT1 + (l1-1) * T1 + (t1-1);      // To use it in R
         alt[i] = x[k];
     }
     return;
}
//
void extract_alt2(int l, int t, int *n, int *rT, int *T, double *x, double *alt)
{
     int i, k, n1, rT1, T1;
     n1 = *n;
     rT1 = *rT;
     T1 = *T;
     
     for (i=0; i < n1; i++) {
         k = i * rT1 + l * T1 + t;      // Not able to use in R
         alt[i] = x[k];                 // for reuse in C/C++   
     }
     return;
}

// /*
// for unequal T, ss = sum(T)
void extract_alt_uneqT(int l, int t, int *n, int *r, int *T, int *rT, 
     double *x, double *alt)
{
     int i, k, n1, r1;
     n1 = *n;
     r1 = *r;
     
     int *T1; 
     T1 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
     
     cumsumint(r, T, T1);   
        
     for(i=0; i < n1; i++){  
         k = i * rT[0] + T1[l] + t;      // Not able to use in R
         alt[i] = x[k];                 // for reuse in C/C++   
     }

     free(T1); 
    
     return;
}

//*/


// It puts together (rearrange) the "o"s
void put_together(int *n, int *r, int *T, double *x, double *out)
{
     int n1, r1, T1, rT;
     int i, j, k;
     n1 = *n;
     r1 = *r;
     T1 = *T;
     rT = r1 * T1;
     
     for(i=0; i < n1; i++) {
         for(j=0; j < rT; j ++) {
             k = j * n1 + i;
           out[j + i * rT] = x[k];
         }
      }         
      return;
}
// It puts together (rearrange) the "o"s
void put_together1(int l, int t, int *n, int *r, int *T, double *x, double *alt)
{ 
  int i, k, n1, r1, rT, T1; 
  n1 = *n;
  r1 = *r;
  T1 = *T;
  rT = r1 * T1;
  for (i=0; i<n1; i++) { 
    k = i* rT + l * T1 + t;
    x[k] = alt[i];
  }
  return;
}

// It puts together (rearrange) the "o"s for unequal T
void put_together1_uneqT(int l, int t, int *n, int *r, int *T, int *rT, 
     double *x, double *alt)
{ 
  int i, k, n1, r1, rT1; 
  n1 = *n;
  r1 = *r;
  rT1 = *rT;

  int *T1; 
  T1 = (int *) malloc((size_t)((r1+1)*sizeof(int)));
//  for(i=0; i<r1; i++){
//          T1[i] = T[i];
//  }

  cumsumint(r, T, T1);
  
  for (i=0; i<n1; i++) { 
    k = i* rT1 + T1[l] + t;
    x[k] = alt[i];
  }
  
  free(T1);
  return;
}


     
// Some other programmes used in GIBBS
void ext(double *IN, double *OUT)
{
     OUT[0] = IN[0];
     return;
}


void ext_wlt(int *m, int *r, int *T, double *wp, double *w)
{
     int j, m1, r1, T1;
     m1 = *m;
     r1 = *r;
     T1 = *T;
     for(j=0; j < m1*r1*T1; j++) {
         w[j] = wp[j];
     }
     return;
}                    


void ext_sige(double *sig_ep, double *sig_e)
{
     sig_e[0] = sig_ep[0];
     return;
}

void ext_sigeta(double *sig_etap, double *sig_eta)
{
     sig_eta[0] = sig_etap[0];
     return;
}

void ext_rho(double *rhop, double *rho)
{
     rho[0] = rhop[0];
     return;
}

void ext_beta(int *p, double *betap, double *beta)
{
     int j, p1;
     p1 = *p;
     for(j=0; j < p1; j++) {
         beta[j] = betap[j];
     }
     return;
}                    

void ext_xil(int *r, double *xi_lp, double *xi_l)
{
     int j, r1;
     r1 = *r;
     for(j=0; j < r1; j++) {
         xi_l[j] = xi_lp[j];
     }
     return;
}                    

void ext_gaml(int *n, int *r, double *gamma_lp, double *gamma_l)
{
     int j, n1, r1, nr;
     n1 = *n;
     r1 = *r;
     nr = n1 * r1;
     for(j=0; j < nr; j++) {
         gamma_l[j] = gamma_lp[j];
     }
     return;
}                    

void ext_mul(int *r, double *mu_lp, double *mu_l)
{
     int j, r1;
     r1 = *r;
     for(j=0; j < r1; j++) {
         mu_l[j] = mu_lp[j];
     }
     return;
}  

void ext_sigl(int *r, double *sig_lp, double *sig_l)
{
     int j, r1;
     r1 = *r;
     for(j=0; j < r1; j++) {
         sig_l[j] = sig_lp[j];
     }
     return;
}  

void ext_o(int *N, double *op, double *o)
{
     int j, N1;
     N1 = *N;
     for(j=0; j < N1; j++) {
         o[j] = op[j];
     }
     return;
}                    

void ext_v(int *M, double *vp, double *v)
{
     int j, M1;
     M1 = *M;
     for(j=0; j < M1; j++) {
         v[j] = vp[j];
     }
     return;
}



///////////////////////////////////////////////////////////////////////////////

////////////////

// used in dlm with only t !!
// With Day observations only
void put(int t, int *n, int *T, double *x, double *alt)
{ 
  int i, k, n1, T1; 
  n1 = *n;
  T1 = *T;
  for (i=0; i<n1; i++) { 
    k = i* T1 + t;
    x[k] = alt[i];
  }
}

// used in dlm
// Extracts the t'th time vector
// T is the time, n is the number of sites
void extract(int t, int *n, int *T, double *x, double *alt)
{
     int i, n1, T1;
     n1 = *n;
     T1 = *T;
     for (i=0; i < n1; i++) {
//         k = i * T1 + t;      // Not able to use in R
         alt[i] = x[i*T1+t];       // for reuse in C/C++   
     }
}



/*
// test reading a C function file and execute it
// not working
void test_C_read(double *input, int *T, double *out)
{
   FILE *fopen(), *fp;
   int c_fnc, i, T1;
   T1 = *T;
   fp = fopen("prog.txt","r");
   c_fnc = getc(fp);
   while (c_fnc != EOF)
   {
   		putchar(c_fnc);
		c_fnc = getc(fp);
   Rprintf("%c", c_fnc);
   }
   qsort(input, T1, sizeof(double), sort_fnc);
   for(i=0; i<T1; i++){
   out[i] = input[i];
   }         
   fclose(fp);
   return;
}     
     



// test writing in txt format
void test_RW(int *its, int *constant)
{
   int i, its1;
   its1 = *its;
   double *mu, *s2, *out, *out2;
   mu = (double *) malloc((size_t)((1)*sizeof(double)));
   s2 = (double *) malloc((size_t)((1)*sizeof(double)));
   out = (double *) malloc((size_t)((1)*sizeof(double)));
   out2 = (double *) malloc((size_t)((its1)*sizeof(double)));

   mu[0] = 0.0;
   s2[0] = 1.0;
   FILE *parafile;
   parafile = fopen("write1.txt","w");
   for(i=0; i<*its; i++){
     mvrnormal(constant, mu, s2, constant, out);
     out2[i] = out[0];
     fprintf(parafile, "%f ", out[0]);
     fprintf(parafile, "\n");
   }
   fclose(parafile);

   parafile = fopen("write2.txt","w");
   for(i=0; i<*its; i++){
     fprintf(parafile, "%f ", out2[i]);
   }
   fclose(parafile);
   
   free(mu); free(s2); free(out); free(out2);

   return;    
}


// test in R
void test(int *i, int *its, int *burnin, double *x, double *alt)
{
     int j, its1, burn1, i1;
     i1 =*i;
     its1 = *its;
     burn1 = *burnin;
     for(j=0; j<(its1-burn1); j++){
        alt[j] = x[j+burn1+i1*its1];
     }
     return;
}         

*/

