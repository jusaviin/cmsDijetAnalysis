#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1
./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_preprocessed_2019-08-13_fiveJobsMissing.root" "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_noStatErrorFromCorrections_lowPtResidualTrack_processed_2019-10-07_fiveJobsMissing.root" true true -1
