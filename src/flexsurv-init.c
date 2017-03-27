#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP flexsurv_basis_matrix(SEXP, SEXP);
extern SEXP flexsurv_basis_vector(SEXP, SEXP);
extern SEXP flexsurv_check_genf(SEXP, SEXP, SEXP, SEXP);
extern SEXP flexsurv_check_gengamma(SEXP, SEXP, SEXP);
extern SEXP flexsurv_check_gompertz(SEXP, SEXP);
extern SEXP flexsurv_check_llogis(SEXP, SEXP);
extern SEXP flexsurv_dbasis_matrix(SEXP, SEXP);
extern SEXP flexsurv_dbasis_vector(SEXP, SEXP);
extern SEXP flexsurv_dexph(SEXP);
extern SEXP flexsurv_dgenf_work(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP flexsurv_dgengamma_work(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP flexsurv_dgompertz_work(SEXP, SEXP, SEXP, SEXP);
extern SEXP flexsurv_dllogis_work(SEXP, SEXP, SEXP, SEXP);
extern SEXP flexsurv_exph(SEXP);
extern SEXP flexsurv_pgenf_work(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP flexsurv_pgengamma_work(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP flexsurv_pgompertz_work(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP flexsurv_pllogis_work(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"flexsurv_basis_matrix",   (DL_FUNC) &flexsurv_basis_matrix,   2},
    {"flexsurv_basis_vector",   (DL_FUNC) &flexsurv_basis_vector,   2},
    {"flexsurv_check_genf",     (DL_FUNC) &flexsurv_check_genf,     4},
    {"flexsurv_check_gengamma", (DL_FUNC) &flexsurv_check_gengamma, 3},
    {"flexsurv_check_gompertz", (DL_FUNC) &flexsurv_check_gompertz, 2},
    {"flexsurv_check_llogis",   (DL_FUNC) &flexsurv_check_llogis,   2},
    {"flexsurv_dbasis_matrix",  (DL_FUNC) &flexsurv_dbasis_matrix,  2},
    {"flexsurv_dbasis_vector",  (DL_FUNC) &flexsurv_dbasis_vector,  2},
    {"flexsurv_dexph",          (DL_FUNC) &flexsurv_dexph,          1},
    {"flexsurv_dgenf_work",     (DL_FUNC) &flexsurv_dgenf_work,     6},
    {"flexsurv_dgengamma_work", (DL_FUNC) &flexsurv_dgengamma_work, 5},
    {"flexsurv_dgompertz_work", (DL_FUNC) &flexsurv_dgompertz_work, 4},
    {"flexsurv_dllogis_work",   (DL_FUNC) &flexsurv_dllogis_work,   4},
    {"flexsurv_exph",           (DL_FUNC) &flexsurv_exph,           1},
    {"flexsurv_pgenf_work",     (DL_FUNC) &flexsurv_pgenf_work,     7},
    {"flexsurv_pgengamma_work", (DL_FUNC) &flexsurv_pgengamma_work, 6},
    {"flexsurv_pgompertz_work", (DL_FUNC) &flexsurv_pgompertz_work, 5},
    {"flexsurv_pllogis_work",   (DL_FUNC) &flexsurv_pllogis_work,   5},
    {NULL, NULL, 0}
};

void R_init_flexsurv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
