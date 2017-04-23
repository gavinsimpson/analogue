#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void WATpred(void *, void *, void *, void *, void *, void *, void *);
extern void WTOL(void *, void *, void *, void *, void *, void *);
extern void xx_chisq_dist(void *, void *, void *, void *, void *, void *);
extern void xx_distance(void *, void *, void *, void *, void *, void *);
extern void xx_kendall(void *, void *, void *, void *, void *, void *);
extern void xx_mixed(void *, void *, void *, void *, void *, void *, void *, void *);
extern void xy_chisq_dist(void *, void *, void *, void *, void *, void *, void *);
extern void xy_distance(void *, void *, void *, void *, void *, void *, void *);
extern void xy_kendall(void *, void *, void *, void *, void *, void *, void *);
extern void xy_mixed(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"WATpred",       (DL_FUNC) &WATpred,       7},
    {"WTOL",          (DL_FUNC) &WTOL,          6},
    {"xx_chisq_dist", (DL_FUNC) &xx_chisq_dist, 6},
    {"xx_distance",   (DL_FUNC) &xx_distance,   6},
    {"xx_kendall",    (DL_FUNC) &xx_kendall,    6},
    {"xx_mixed",      (DL_FUNC) &xx_mixed,      8},
    {"xy_chisq_dist", (DL_FUNC) &xy_chisq_dist, 7},
    {"xy_distance",   (DL_FUNC) &xy_distance,   7},
    {"xy_kendall",    (DL_FUNC) &xy_kendall,    7},
    {"xy_mixed",      (DL_FUNC) &xy_mixed,      9},
    {NULL, NULL, 0}
};

void R_init_analogue(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
