#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1
./processHistogramsInStepsAsymmetry.sh "data/PbPbMC_RecoGen_akFlowPuCs4PfJets_noUncorr_improvisedMixing_xjBins_jet100Trigger_subeNon0_JECv5b_preprocessed_2019-09-04.root" "data/PbPbMC_RecoGen_akFlowPuCs4PfJets_noUncorr_improvisedMixing_noCorrections_xjBins_jet100Trigger_subeNon0_JECv5b_processed_2019-09-04.root" false false -1
./processHistogramsInStepsAsymmetry.sh "data/PbPbMC_RecoGen_akFlowPuCs4PfJets_noUncorr_improvisedMixing_xjBins_sube0_JECv5b_preprocessed_2019-09-04.root" "data/PbPbMC_RecoGen_akFlowPuCs4PfJets_noUncorr_improvisedMixing_noCorrections_xjBins_sube0_JECv5b_processed_2019-09-04.root" false false -1
./processHistogramsInStepsAsymmetry.sh "data/PbPbMC_RecoGen_akFlowPuCs4PfJets_noUncorr_improvisedMixing_xjBins_subeNon0_JECv5b_preprocessed_2019-09-04.root" "data/PbPbMC_RecoGen_akFlowPuCs4PfJets_noUncorr_improvisedMixing_noCorrections_xjBins_subeNon0_JECv5b_processed_2019-09-04.root" false false -1
