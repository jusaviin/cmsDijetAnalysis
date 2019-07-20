#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false false
#./processHistogramsInStepsPp.sh "data/" "data/" false false false
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false false
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false false
./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb_pfCsJets_highLeadingJetCut_xjBins_wtaAxis_noUncorr_improvisedMixing_preprocessed_2019-07-05.root" "data/dijetPbPb_pfCsJets_highLeadingJetCut_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root" true true false
./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb_pfCsJets_lowLeadingJetCut_xjBins_wtaAxis_noUncorr_improvisedMixing_preprocessed_2019-07-05.root" "data/dijetPbPb_pfCsJets_lowLeadingJetCut_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root" true true false
