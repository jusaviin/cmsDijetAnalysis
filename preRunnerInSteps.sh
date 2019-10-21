#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_JECv6_preprocessed_2019-09-26.root" "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0AtLowPt_wtaAxis_JECv6_sidebandSpilloverLooseCutUntil8_processed_2019-09-26.root" true true -1 "data/"
#./processHistogramsInStepsAtHighPt.sh "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_preprocessed_2019-09-26.root" "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0AtLowPt_wtaAxis_JECv6_sidebandSpilloverLooseCutUntil8_processed_2019-09-26.root" true true -1
./processHistogramsInSteps.sh "data/PbPbMC2018_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_centShift5_preprocessed_2019-10-12.root" "data/PbPbMC2018_GenReco_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_trackDeltaRcorrectedNoSymmetry_centShift5_processed_2019-10-12.root" false false -1 "data/"
