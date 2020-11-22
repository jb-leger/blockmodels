#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP dispatcher(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kmeans(SEXP, SEXP);
extern SEXP split(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"dispatcher", (DL_FUNC) &dispatcher, 5},
    {"kmeans",     (DL_FUNC) &kmeans,     2},
    {"split",      (DL_FUNC) &split,      1},
    {NULL, NULL, 0}
};

void R_init_blockmodels(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
