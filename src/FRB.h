/* Copied from robustbase.h */
  
#include <R.h>
#include <Rinternals.h>
#include <complex.h>
  
/* For internationalized messages */
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP : String)
#endif
  
void R_frb(double *xx, double *y, double *w, int *n, int *p, double *beta_m,
            double *scale, double *chi_res_s, double *bbetas, int *nboot,
            double *xx3, double *v2, int *bind);
