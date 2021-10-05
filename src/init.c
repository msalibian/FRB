#include <R_ext/Rdynload.h>
#include "FRB.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

void R_frb( double *xx, double *y, double *w, int *n, int *p, double *beta_m,
            double *scale, double *chi_res_s, double *bbetas, int *nboot,
            double *xx3, double *v2, int *bind);
  

static R_NativePrimitiveArgType R_frb_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, 
  REALSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, INTSXP
};

static const R_CMethodDef CEntries[]  = {
  CDEF(R_frb),
  {NULL, NULL, 0}
};

void R_init_RBF(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE); 
};
