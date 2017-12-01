
#include "header.h"

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void GeoDist_km(void *, void *, void *, void *);
extern void GeoDist_miles(void *, void *, void *, void *);
extern void GIBBS_ar(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GIBBS_gp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GIBBS_sumpred_gpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GIBBS_sumpred_txt_ar(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GIBBS_sumpred_txt_gp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GIBBS_zfitsum_onephi_gpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void z_pr_its_ar(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void z_pr_its_gp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void z_pr_its_gpp1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zlt_fore_ar_its_anysite(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zlt_fore_gp_its(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zlt_fore_gpp_its(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"GeoDist_km",               (DL_FUNC) &GeoDist_km,                4},
    {"GeoDist_miles",            (DL_FUNC) &GeoDist_miles,             4},
    {"GIBBS_ar",                 (DL_FUNC) &GIBBS_ar,                 50},
    {"GIBBS_gp",                 (DL_FUNC) &GIBBS_gp,                 45},
    {"GIBBS_sumpred_gpp",        (DL_FUNC) &GIBBS_sumpred_gpp,        52},
    {"GIBBS_sumpred_txt_ar",     (DL_FUNC) &GIBBS_sumpred_txt_ar,     44},
    {"GIBBS_sumpred_txt_gp",     (DL_FUNC) &GIBBS_sumpred_txt_gp,     43},
    {"GIBBS_zfitsum_onephi_gpp", (DL_FUNC) &GIBBS_zfitsum_onephi_gpp, 59},
    {"z_pr_its_ar",              (DL_FUNC) &z_pr_its_ar,              24},
    {"z_pr_its_gp",              (DL_FUNC) &z_pr_its_gp,              22},
    {"z_pr_its_gpp1",            (DL_FUNC) &z_pr_its_gpp1,            21},
    {"zlt_fore_ar_its_anysite",  (DL_FUNC) &zlt_fore_ar_its_anysite,  24},
    {"zlt_fore_gp_its",          (DL_FUNC) &zlt_fore_gp_its,          22},
    {"zlt_fore_gpp_its",         (DL_FUNC) &zlt_fore_gpp_its,         23},
    {NULL, NULL, 0}
};

void attribute_visible R_init_spTimer(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

//
