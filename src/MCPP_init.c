#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP bignum2func_cauchy(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bignum2func_norm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bignum2func_unif(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gethsamp(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"bignum2func_cauchy", (DL_FUNC) &bignum2func_cauchy, 5},
    {"bignum2func_norm",   (DL_FUNC) &bignum2func_norm,   5},
    {"bignum2func_unif",   (DL_FUNC) &bignum2func_unif,   5},
    {"gethsamp",           (DL_FUNC) &gethsamp,           4},
    {NULL, NULL, 0}
};

void R_init_MCPP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
