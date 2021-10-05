

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

#define ZERO 1e-10
#define TOL_INVERSE ZERO


void R_frb( double *xx, double *y, double *w, int *n, int *p, double *beta_m,
		double *scale, double *chi_res_s, double *bbetas, int *nboot,
		double *xx3, double *v2, int *bind)
{

void sampler_i(int n, int *x);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
int inverse(double **,double **, int);
void matias_vec_vec(double **, double *, double *, int);
void scalar_mat(double **, double, double **, int, int);
void scalar_vec(double *, double, double *, int);
void sum_mat(double **,double **, double **, int, int);
void sum_vec(double *, double *, double *, int);
void dif_mat(double **, double **, double **, int , int );
void dif_vec(double *, double *, double *, int);
void mat_vec(double **, double *, double *, int, int);
void mat_mat(double **, double **, double **, int, int, int);
void disp_vec(double *, int);
void disp_vec_int(int *, int);
void disp_mat(double **, int, int);

register int i,j;
double  **x, **x3,  **x2, **x4, *v_aux, *v, s=0;
int *indices;



x = (double **) malloc( sizeof(double *) * (*n) );
x2 = (double **) malloc( sizeof(double *) * (*p) );
x3 = (double **) malloc( sizeof(double *) * (*p) );
x4 = (double **) malloc( sizeof(double *) * (*p) );
indices = (int *) malloc( sizeof(int) * (*n) );
v_aux = (double *) malloc( sizeof(double) * (*p));
v = (double *) malloc( sizeof(double) * (*p));
for(i=0;i<(*n);i++)  
	x[i] = (double *) malloc( sizeof(double) * (*p) );

for(i=0;i<(*p);i++) {
	x2[i] = (double *) malloc( sizeof(double) * (*p) );
	x3[i] = (double *) malloc( sizeof(double) * (*p) );
	x4[i] = (double *) malloc( sizeof(double) * (*p) );
};

GetRNGstate(); /* get seed from R */

/* copy design matrix into x - just laziness */
for(i=0;i<(*n);i++) 
	for(j=0;j<(*p);j++)
		x[i][j]=xx[j*(*n)+i];

/* copy correction matrix into x3 - just lazy again */
for(i=0;i<(*p);i++) 
	for(j=0;j<(*p);j++)
		x3[i][j]=xx3[j*(*p)+i];

/* start bootstrap iterations */
for(i=0; i < (*nboot); i++) {
	/* get indices for bootstrap sample */
	sampler_i(*n, indices);
	R_CheckUserInterrupt();	
	reset_vec(v, (*p));
	reset_mat(x2, (*p), (*p));
	s = 0.0;
	for(j=0;j<(*n);j++) {
		scalar_vec(x[indices[j]],y[indices[j]]*w[indices[j]],v_aux,(*p));
		sum_vec(v,v_aux,v,(*p));
		matias_vec_vec(x4,x[indices[j]],x[indices[j]],(*p));
		scalar_mat(x4,w[indices[j]],x4,(*p),(*p));
		sum_mat(x2,x4,x2,(*p),(*p));
		s += chi_res_s[indices[j]]; /* Chi(res_s[indices[j]] / scale , c); */
	};
	s = s * (*scale) / .5 / (double) ((*n) - (*p))  ; 
	inverse(x2,x4,(*p));				/* x4 <- x2^-1 */
	mat_vec(x4,v,v_aux,(*p),(*p));		/* v_aux <- x4 * v */ 
	dif_vec(v_aux,beta_m,v_aux,(*p)); 	/* v_aux <- v_aux - beta_m */
	/* v has the robust bootstrapped vector, correct it */ 
	mat_vec(x3,v_aux,v,(*p),(*p));		/* v <- x3 * v_aux */ 
	scalar_vec(v2,s-(*scale),v_aux,(*p));
	sum_vec(v_aux,v,v,(*p));
	/* store the betas (R-wise!) */
	for(j=0;j<(*p);j++) {
		bbetas[j*(*nboot)+i]=v[j];
	};
	for(j=0; j< (*n); j++) {
	  bind[j*(*nboot)+i]=indices[j] + 1;
	}
};

free(indices);
for(i=0; i < (*p); i++) {
	free(x2[i]); free(x3[i]);
	free(x4[i]);
};
for(i=0; i < (*n); i++) 
	free(x[i]); 

free(x2); free(x3); free(x); 
free(x4); free(v_aux); 
free(v);

return;

}


int lu(double **a,int *P, double *x)
{
int *pp,p;
register int i,j,k;
double *kk,s;
p = *P;
if ((pp = (int *) malloc(p*sizeof(int)))==NULL)
	{ Rprintf("\nNot enough memory in LU\n");
	  Rf_error("\nNot enough memory in LU\n"); }
/* pp vector storing the permutations */
for(j=0;j<p;j++)   /* cols */
{ pp[j]=j;
  for(i=j;i<p;i++)   /* filas */
	if ( fabs( a[i][j] ) > fabs( a[pp[j]][j] ) )
		pp[j]=i;
  if ( pp[j] != j )       /* permuto las filas cambiando los punt */
	{ kk=a[j];
	  a[j]=a[pp[j]];
	  a[pp[j]]=kk;
	};
  /* salida si el sistema resulta singular (det=0)
   * se detecta si el pivote (j,j) es cero  */
/*  if ( a[j][j] == 0 ) {   free(pp);
				return(1);
				}; */
    if ( fabs(a[j][j]) < TOL_INVERSE ) {   free(pp);
				return(1);
				};
  for(k=(j+1);k<p;k++)
	a[k][j] = a[k][j] / a[j][j];
  for(k=(j+1);k<p;k++)
	for(i=(j+1);i<p;i++)
		a[k][i] = a[k][i] - a[k][j] * a[j][i];

};    /* cierra el for de j */
for(i=0;i<p;i++)
	{ s=0.0;
	  for(j=0;j<i;j++) s += a[i][j] * x[j];
	  x[i] = a[i][p] - s;          /* y[i]=a[i][p] */
	};
for(i=(p-1);i>=0;i--)
	{ s=0;
	  for(j=(i+1);j<p;j++) s += a[i][j] * x[j];
	  x[i] = (x[i] - s) / a[i][i];
	  };
free(pp);
return(0);
}


void sum_mat(double **a, double **b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++) 
		c[i][j] = a[i][j] + b[i][j];
}

void matias_vec_vec(double **a, double *v1, double *v2, int n)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		a[i][j] = v1[i] * v2[j];/* could take advantage of symmetry */
}

void scalar_mat(double **a, double b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
        for(j=0;j<m;j++)
		c[i][j]  = b * a[i][j];
}

void scalar_vec(double *a, double b, double *c, int n)
{
register int i;
for(i=0;i<n;i++)
	c[i]  = b * a[i];
}

double vecprime_vec(double *a, double *b, int n)
{
register int i;
double s = 0.0;
for(i=0;i<n;i++) s += a[i] * b[i];
return(s);
}

void sum_vec(double *a, double *b, double *c, int n)
{
register int i;
for(i=0;i<n;i++) c[i] = a[i] + b[i];
}

void dif_vec(double *a, double *b, double *c, int n)
{
register int i;
for(i=0;i<n;i++) c[i] = a[i] - b[i];
}

void dif_mat(double **a, double **b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++) 
	for(j=0;j<m;j++) c[i][j] = a[i][j] - b[i][j];
}

void mat_vec(double **a, double *b, double *c, int n, int m)
{
register int i,j; 
for(i=0;i<n;i++) 
	for(c[i]=0,j=0;j<m;j++) c[i] += a[i][j] * b[j];
}

void mat_mat(double **a, double **b, double **c, int n, 
		int m, int l)
{
register int i,j,k; 
for(i=0;i<n;i++) 
	for(j=0;j<l;j++) {
	c[i][j] = 0; 
	for(k=0;k<m;k++) c[i][j] += a[i][k] * b[k][j];
	};
}

void disp_vec(double *a, int n)
{
register int i;
Rprintf("\n");
for(i=0;i<n; i++) Rprintf("%lf ",a[i]);
Rprintf("\n");
}

void disp_vec_int(int *a, int n)
{
  register int i;
  Rprintf("\n");
  for(i=0;i<n; i++) Rprintf("%ld ",a[i]);
  Rprintf("\n");
}


int inverse(double **a, double **b, int n)
{
int lu(double **, int *, double *);
void mat_vec(double **, double *, double *, int, int);
void disp_vec(double *, int);
register int i,j,k;
double **c, *e;
c = (double **) malloc( n * sizeof(double *));
e = (double *) malloc( n * sizeof(double));
for(i=0;i<n;i++) c[i] = (double *) malloc ( (n+1) * sizeof(double) );
for(i=0;i<n;i++) {   /* i-th column */

for(j=0;j<n;j++)
	for(k=0;k<n;k++) c[j][k] = a[j][k];
for(j=0;j<i;j++) c[j][n] = 0.0;
c[i][n] = 1.0;
for(j=i+1;j<n;j++) c[j][n] = 0.0;
if( lu(c,&n,e) == 1) {
	for(i=0;i<n;i++) free(c[i]);
	free(c);free(e);
	return(1);
	};	
for(j=0;j<n;j++) b[j][i] = e[j] ;
};
for(i=0;i<n;i++) free(c[i]);
free(c);free(e);
return(0);
}

void reset_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		a[i][j] = 0.0;
}

void reset_vec(double *a, int n)
{
register int i;
for(i=0;i<n;i++) a[i] = 0.0;
}

void disp_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++) {
Rprintf("\n");
for(j=0;j<m;j++) Rprintf("%10.8f ",a[i][j]);
};
Rprintf("\n");
}

void sampler_i(int n, int *x)
{
/* function to get a random sample of
 * indices (0 to n-1)
 * *x receives the output
 * rand() returns an integer between 0 and RAND_MAX
 */
int i;
for(i=0;i<n;i++) 
	x[i] = (int) ( n * unif_rand() );
}


