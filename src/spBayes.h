#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#include <R.h>
#include <Rinternals.h>

extern "C" {
  
  SEXP adaptMetropGibbs(SEXP ltd_r, SEXP starting_r, SEXP tuning_r,
			SEXP acceptRate_r,
			SEXP nBatch_r, SEXP batchLength_r,
			SEXP verbose_r, SEXP nTheta_r, 
			SEXP reportBatch, SEXP rho_r);
  
  SEXP idist(SEXP coords1_r, SEXP n1_r, SEXP coords2_r, SEXP n2_r, SEXP p_r, SEXP D_r);
  
  SEXP mkSpCov(SEXP coords_r, SEXP n_r, SEXP m_r, SEXP Psi_r, SEXP V_r, SEXP theta_r, SEXP covModel_r);
  
  SEXP ptsInPoly(SEXP verts_r, SEXP nVerts_r, SEXP pts_r, SEXP nPts_r, SEXP inPtIndx_r, SEXP nInPts_r);
  
  SEXP spMPPMvDIC(SEXP Q_r, SEXP knotsD_r, SEXP coordsKnotsD_r, SEXP n_r, SEXP m_r, SEXP q_r, 
		  SEXP Psi_r, SEXP V_r, SEXP phi_r, SEXP nu_r, SEXP covModel_r, SEXP CEps_r);
  
  SEXP nonSpGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r, SEXP weights_r,
		      SEXP betaPrior_r, SEXP betaNorm_r, SEXP betaStarting_r, SEXP betaTuning_r, 
		      SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r, SEXP amcmc_r);
  
  SEXP spPPGLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r,  SEXP family_r, SEXP weights_r,
	       SEXP m_r, SEXP knotsD_r, SEXP knotsCoordsD_r, 
	       SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	       SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP w_strStarting_r,
	       SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP nuTuning_r, SEXP betaTuning_r, SEXP w_strTuning_r,
	       SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r, SEXP weights_r,
		     SEXP m_r, SEXP knotsD_r, SEXP knotsCoordsD_r, 
		     SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
		     SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP w_strStarting_r,
		     SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP nuTuning_r, SEXP betaTuning_r, SEXP w_strTuning_r,
		     SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spGLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r, SEXP family_r, SEXP weights_r,
	     SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	     SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP wStarting_r,
	     SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP nuTuning_r, SEXP betaTuning_r, SEXP wTuning_r,
	     SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r, SEXP family_r, SEXP weights_r,
		   SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
		   SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP wStarting_r,
		   SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP nuTuning_r, SEXP betaTuning_r, SEXP wTuning_r,
		   SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP knotsD_r, SEXP knotsCoordsD_r, 
	      SEXP modPP_r, SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	      SEXP betaStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r,
	      SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP nuTuning_r, 
	      SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  
  SEXP spLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
	    SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	    SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r,
	    SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP nuTuning_r, 
	    SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spGLMMisalign_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r, SEXP family_r, SEXP weights_r, 
			   SEXP betaPrior_r, SEXP betaNorm_r, 
			   SEXP KIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
			   SEXP phiStarting_r, SEXP AStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP wStarting_r,
			   SEXP phiTuning_r, SEXP ATuning_r, SEXP nuTuning_r , SEXP betaTuning_r, SEXP wTuning_r,
			   SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMisalign(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,
		  SEXP betaPrior_r, SEXP betaNorm_r, 
		  SEXP KPrior_r, SEXP KPriorName_r, 
		  SEXP PsiPrior_r, 
		  SEXP nuUnif_r, SEXP phiUnif_r,
		  SEXP phiStarting_r, SEXP AStarting_r, SEXP PsiStarting_r, SEXP nuStarting_r, 
		  SEXP phiTuning_r, SEXP ATuning_r, SEXP PsiTuning_r, SEXP nuTuning_r, 
		  SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPMvGLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP family_r, SEXP weights_r,
		 SEXP q_r, SEXP knotsD_r, SEXP knotsCoordsD_r,
		 SEXP betaPrior_r, SEXP betaNorm_r, SEXP KIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
		 SEXP phiStarting_r, SEXP AStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP w_strStarting_r,
		 SEXP phiTuning_r, SEXP ATuning_r, SEXP nuTuning_r , SEXP betaTuning_r, SEXP w_strTuning_r,
		 SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPMvGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP family_r,  SEXP weights_r,
		       SEXP q_r, SEXP knotsD_r, SEXP knotsCoordsD_r,
		       SEXP betaPrior_r, SEXP betaNorm_r, SEXP KIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
		       SEXP phiStarting_r, SEXP AStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP w_strStarting_r,
		       SEXP phiTuning_r, SEXP ATuning_r, SEXP nuTuning_r , SEXP betaTuning_r, SEXP w_strTuning_r,
		       SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);

  SEXP spMvGLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,SEXP family_r, SEXP weights_r,
	       SEXP betaPrior_r, SEXP betaNorm_r, SEXP KIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
	       SEXP phiStarting_r, SEXP AStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP wStarting_r,
	       SEXP phiTuning_r, SEXP ATuning_r, SEXP nuTuning_r , SEXP betaTuning_r, SEXP wTuning_r,
	       SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMvGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,SEXP family_r, SEXP weights_r,
		     SEXP betaPrior_r, SEXP betaNorm_r, SEXP KIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
		     SEXP phiStarting_r, SEXP AStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP wStarting_r,
		     SEXP phiTuning_r, SEXP ATuning_r, SEXP nuTuning_r , SEXP betaTuning_r, SEXP wTuning_r,
		     SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPMvLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP g_r, SEXP knotsD_r, SEXP knotsCoordsD_r, SEXP modPP_r, 
		SEXP betaPrior_r, SEXP betaNorm_r, 
		SEXP KPrior_r, SEXP KPriorName_r, 
		SEXP PsiPrior_r, SEXP PsiPriorName_r, SEXP PsiDiag_r, 
		SEXP nuUnif_r, SEXP phiUnif_r,
		SEXP betaStarting_r, SEXP phiStarting_r, SEXP AStarting_r, SEXP LStarting_r, SEXP nuStarting_r, 
		SEXP phiTuning_r, SEXP ATuning_r, SEXP LTuning_r, SEXP nuTuning_r, 
		SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMvLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,
	      SEXP betaPrior_r, SEXP betaNorm_r, 
	      SEXP KPrior_r, SEXP KPriorName_r, 
	      SEXP PsiPrior_r, SEXP PsiPriorName_r, SEXP PsiDiag_r, 
	      SEXP nuUnif_r, SEXP phiUnif_r,
	      SEXP phiStarting_r, SEXP AStarting_r, SEXP LStarting_r, SEXP nuStarting_r, 
	      SEXP phiTuning_r, SEXP ATuning_r, SEXP LTuning_r, SEXP nuTuning_r, 
	      SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPMvGLMPredict(SEXP family_r, SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP g_r, SEXP p_r, SEXP Z_r, SEXP q_r, SEXP knotsD_r, SEXP predKnotsD_r, 
			SEXP samples_r, SEXP wSamples_r, SEXP nSamples_r, SEXP covModel_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMvGLMPredict(SEXP family_r, SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP p_r, SEXP Z_r, SEXP q_r, SEXP obsD_r, SEXP obsPredD_r,
		      SEXP samples_r, SEXP wSamples_r, SEXP nSamples_r, SEXP covModel_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPLMPredict(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP p_r, SEXP m_r, SEXP Z_r, SEXP q_r,
		     SEXP samples_r, SEXP nSamples_r, 
		     SEXP beta_r, SEXP sigmaSqIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		     SEXP nugget_r, SEXP knotsD_r, SEXP knotsObsD_r, SEXP knotsPredD_r, SEXP covModel_r, SEXP modPP_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spLMPredict(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP p_r, SEXP Z_r, SEXP q_r,
		   SEXP samples_r, SEXP nSamples_r, 
		   SEXP betaPrior_r, SEXP betaNorm_r, SEXP beta_r, SEXP sigmaSqIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		   SEXP obsD_r, SEXP obsPredD_r, SEXP covModel_r, SEXP nugget_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPMvLMPredict(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP g_r, SEXP p_r, SEXP Z_r, SEXP q_r, 
		       SEXP knotsD_r, SEXP knotsObsD_r, SEXP knotsPredD_r, 
		       SEXP samples_r, SEXP beta_r,  SEXP nSamples_r, 
		       SEXP nugget_r, SEXP PsiDiag_r, SEXP covModel_r, 
		       SEXP modPP_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMvLMPredict(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP p_r, SEXP Z_r, SEXP q_r, SEXP obsD_r, SEXP obsPredD_r, 
		     SEXP samples_r, SEXP beta_r,  SEXP nSamples_r, 
		     SEXP betaPrior_r, SEXP betaNorm_r, 
		     SEXP nugget_r, SEXP PsiDiag_r, SEXP covModel_r, 
		     SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMisalignPredict(SEXP Y_r, SEXP X_r, SEXP m_r, SEXP n_r, SEXP p_r, 
			 SEXP Z_r, SEXP nPred_r, SEXP pPred_r, 
			 SEXP obsD_r, SEXP predObsD_r,
			 SEXP samples_r, SEXP beta_r, SEXP nSamples_r, 
			 SEXP betaPrior_r, SEXP betaNorm_r, 
			 SEXP nugget_r, SEXP covModel_r, 
			 SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMisalignGLMPredict(SEXP family_r, SEXP Y_r, SEXP X_r, SEXP m_r, SEXP n_r, SEXP p_r, 
			    SEXP Z_r, SEXP nPred_r, SEXP pPred_r, 
			    SEXP obsD_r, SEXP predObsD_r,
			    SEXP samples_r, SEXP w_r, SEXP nSamples_r, 
			    SEXP covModel_r, 
			    SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPLMRecover(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP p_r, SEXP m_r,
		     SEXP samples_r, SEXP nSamples_r, 
		     SEXP beta_r, SEXP sigmaSqIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		     SEXP nugget_r, SEXP knotsD_r, SEXP knotsObsD_r, SEXP covModel_r, SEXP modPP_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spLMRecover(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
		   SEXP samples_r, SEXP nSamples_r, 
		   SEXP sigmaSqIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		   SEXP betaPrior_r, SEXP betaNorm_r, 	   
		   SEXP nugget_r, SEXP covModel_r, 
		   SEXP beta_r, SEXP w_r,
		   SEXP verbose_r, SEXP nReport_r);
  
  SEXP spPPMvLMRecover(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP g_r, SEXP p_r, 
		       SEXP knotsD_r, SEXP knotsObsD_r, 
		       SEXP samples_r, SEXP beta_r,  SEXP nSamples_r, 
		       SEXP nugget_r, SEXP PsiDiag_r, SEXP covModel_r, 
		       SEXP modPP_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMvLMRecover(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,
		     SEXP samples_r, SEXP nSamples_r, 
		     SEXP betaPrior_r, SEXP betaNorm_r, 	   
		     SEXP nugget_r, SEXP PsiDiag_r, SEXP covModel_r, 
		     SEXP beta_r, SEXP w_r,
		     SEXP verbose_r, SEXP nReport_r);
  
  SEXP spMisalignRecover(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,
			 SEXP samples_r, SEXP nSamples_r, 
			 SEXP betaPrior_r, SEXP betaNorm_r, 	   
			 SEXP nugget_r, SEXP covModel_r, 
			 SEXP beta_r, SEXP w_r,
			 SEXP verbose_r, SEXP nReport_r);

    SEXP spDynLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP Nt_r, SEXP coordsD_r,
	       SEXP beta0Norm_r, SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r, SEXP sigmaEtaIW_r, 
	       SEXP betaStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r, SEXP sigmaEtaStarting_r, 
		 SEXP phiTuning_r, SEXP nuTuning_r, SEXP covModel_r, SEXP nSamples_r, SEXP missing_r, SEXP getFitted_r, SEXP verbose_r, SEXP nReport_r);

    SEXP spPPDynLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP Nt_r, SEXP knotsD_r, SEXP coordsKnotsD_r,
		 SEXP beta0Norm_r, SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r, SEXP sigmaEtaIW_r, 
		 SEXP betaStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r, SEXP sigmaEtaStarting_r, 
		   SEXP phiTuning_r, SEXP nuTuning_r, SEXP covModel_r, SEXP nSamples_r, SEXP missing_r, SEXP getFitted_r, SEXP verbose_r, SEXP nReport_r);

}
