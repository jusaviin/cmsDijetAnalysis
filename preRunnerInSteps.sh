#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1
./processHistogramsInSteps.sh "data/PbPbMC_GenGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_sube0_onlyGluonJets_matchedJets_JECv4_preprocessed_2019-08-09.root" "data/PbPbMC_GenGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_sube0_onlyGluonJets_matchedJets_noCorrections_JECv4_processed_2019-08-09.root" false false -1
./processHistogramsInSteps.sh "data/PbPbMC_GenGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_sube0_onlyQuarkJets_matchedJets_JECv4_preprocessed_2019-08-09.root" "data/PbPbMC_GenGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_sube0_onlyQuarkJets_matchedJets_noCorrections_JECv4_processed_2019-08-09.root" false false -1
./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_sube0_onlyGluonJets_matchedJets_JECv4_preprocessed_2019-08-09.root" "data/PbPbMC_RecoGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_sube0_onlyGluonJets_matchedJets_noCorrections_JECv4_processed_2019-08-09.root" false false -1
./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_sube0_onlyQuarkJets_matchedJets_JECv4_preprocessed_2019-08-09.root" "data/PbPbMC_RecoGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_sube0_onlyQuarkJets_matchedJets_noCorrections_JECv4_processed_2019-08-09.root" false false -1
