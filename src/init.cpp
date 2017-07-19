#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spBayes.h"

static const R_CallMethodDef CallEntries[] = {
    {"adaptMetropGibbs",     (DL_FUNC) &adaptMetropGibbs,     10},
    {"idist",                (DL_FUNC) &idist,                 6},
    {"mkSpCov",              (DL_FUNC) &mkSpCov,               7},
    {"nonSpGLM_AMCMC",       (DL_FUNC) &nonSpGLM_AMCMC,       16},
    {"ptsInPoly",            (DL_FUNC) &ptsInPoly,             6},
    {"spDynLM",              (DL_FUNC) &spDynLM,              26},
    {"spGLM",                (DL_FUNC) &spGLM,                26},
    {"spGLM_AMCMC",          (DL_FUNC) &spGLM_AMCMC,          28},
    {"spGLMMisalign_AMCMC",  (DL_FUNC) &spGLMMisalign_AMCMC,  29},
    {"spLM",                 (DL_FUNC) &spLM,                 27},
    {"spLMPredict",          (DL_FUNC) &spLMPredict,          21},
    {"spLMRecover",          (DL_FUNC) &spLMRecover,          19},
    {"spMisalign",           (DL_FUNC) &spMisalign,           29},
    {"spMisalignGLMPredict", (DL_FUNC) &spMisalignGLMPredict, 17},
    {"spMisalignPredict",    (DL_FUNC) &spMisalignPredict,    19},
    {"spMisalignRecover",    (DL_FUNC) &spMisalignRecover,    16},
    {"spMPPMvDIC",           (DL_FUNC) &spMPPMvDIC,           12},
    {"spMvGLM",              (DL_FUNC) &spMvGLM,              27},
    {"spMvGLM_AMCMC",        (DL_FUNC) &spMvGLM_AMCMC,        29},
    {"spMvGLMPredict",       (DL_FUNC) &spMvGLMPredict,       16},
    {"spMvLM",               (DL_FUNC) &spMvLM,               31},
    {"spMvLMPredict",        (DL_FUNC) &spMvLMPredict,        19},
    {"spMvLMRecover",        (DL_FUNC) &spMvLMRecover,        17},
    {"spPPDynLM",            (DL_FUNC) &spPPDynLM,            28},
    {"spPPGLM",              (DL_FUNC) &spPPGLM,              28},
    {"spPPGLM_AMCMC",        (DL_FUNC) &spPPGLM_AMCMC,        30},
    {"spPPLM",               (DL_FUNC) &spPPLM,               31},
    {"spPPLMPredict",        (DL_FUNC) &spPPLMPredict,        22},
    {"spPPLMRecover",        (DL_FUNC) &spPPLMRecover,        19},
    {"spPPMvGLM",            (DL_FUNC) &spPPMvGLM,            29},
    {"spPPMvGLM_AMCMC",      (DL_FUNC) &spPPMvGLM_AMCMC,      31},
    {"spPPMvGLMPredict",     (DL_FUNC) &spPPMvGLMPredict,     17},
    {"spPPMvLM",             (DL_FUNC) &spPPMvLM,             35},
    {"spPPMvLMPredict",      (DL_FUNC) &spPPMvLMPredict,      20},
    {"spPPMvLMRecover",      (DL_FUNC) &spPPMvLMRecover,      17},
    {NULL, NULL, 0}
};

void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_sp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
