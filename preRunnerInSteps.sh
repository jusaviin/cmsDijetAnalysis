#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false false
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false false
./processHistogramsInSteps.sh "data/dijetPbPb_pfCsJets_xj_noUncorr_improvisedMixing_preprocessed_2019-07-05.root" "data/dijetPbPb_pfCsJets_xj_noUncorr_improvisedMixing_allCorrections_spilloverAlsoForSubleading_processed_2019-07-05.root" true true false
