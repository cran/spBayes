#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "spBayes.h"

static R_CallMethodDef CallEntries[] = {
  {"adaptMetropGibbs", (DL_FUNC) &adaptMetropGibbs, 10},
  {"idist", (DL_FUNC) &idist, 6},
  {"mkSpCov", (DL_FUNC) &mkSpCov, 7},
  {"ptsInPoly", (DL_FUNC) &ptsInPoly, 6},
  {"spMPPMvDIC", (DL_FUNC) &spMPPMvDIC, 12},
  {"nonSpGLM_AMCMC", (DL_FUNC) &nonSpGLM_AMCMC, 16},
  {"spPPGLM", (DL_FUNC) &spPPGLM, 29},
  {"spPPGLM_AMCMC", (DL_FUNC) &spPPGLM_AMCMC, 31},
  {"spGLM", (DL_FUNC) &spGLM, 26},
  {"spGLM_AMCMC", (DL_FUNC) &spGLM_AMCMC, 28},
  {"spPPLM", (DL_FUNC) &spPPLM, 32},
  {"spLM", (DL_FUNC) &spLM, 27},
  {"spGLMMisalign_AMCMC", (DL_FUNC) &spGLMMisalign_AMCMC, 29},
  {"spMisalign", (DL_FUNC) &spMisalign, 29},
  {"spPPMvGLM", (DL_FUNC) &spPPMvGLM, 30},
  {"spPPMvGLM_AMCMC", (DL_FUNC) &spPPMvGLM_AMCMC, 32},
  {"spMvGLM", (DL_FUNC) &spMvGLM, 27},
  {"spMvGLM_AMCMC", (DL_FUNC) &spMvGLM_AMCMC, 29},
  {"spPPMvLM", (DL_FUNC) &spPPMvLM, 36},
  {"spMvLM", (DL_FUNC) &spMvLM, 31},
  {"spPPMvGLMPredict", (DL_FUNC) &spPPMvGLMPredict, 17},
  {"spMvGLMPredict", (DL_FUNC) &spMvGLMPredict, 16},
  {"spPPLMPredict", (DL_FUNC) &spPPLMPredict, 22},
  {"spLMPredict", (DL_FUNC) &spLMPredict, 21},
  {"spPPMvLMPredict", (DL_FUNC) &spPPMvLMPredict, 20},
  {"spMvLMPredict", (DL_FUNC) &spMvLMPredict, 19},
  {"spMisalignPredict", (DL_FUNC) &spMisalignPredict, 19},
  {"spMisalignGLMPredict", (DL_FUNC) &spMisalignGLMPredict, 17},
  {"spPPLMRecover", (DL_FUNC) &spPPLMRecover, 19},
  {"spLMRecover", (DL_FUNC) &spLMRecover, 19},
  {"spPPMvLMRecover", (DL_FUNC) &spPPMvLMRecover, 17},
  {"spMvLMRecover", (DL_FUNC) &spMvLMRecover, 17},
  {"spMisalignRecover", (DL_FUNC) &spMisalignRecover, 16},
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
