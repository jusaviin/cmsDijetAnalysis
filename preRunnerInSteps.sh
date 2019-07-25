#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false false
#./processHistogramsInStepsPp.sh "data/" "data/" false false false
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false false
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false false
./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_highForest_pfCsJets_noUncorr_wtaAxis_improvisedMixing_preprocessed_2019-07-22.root" "data/dijetPbPb2018_highForest_pfCsJets_noUncorr_wtaAxis_improvisedMixing_noCorrections_processed_2019-07-22.root" false false false
