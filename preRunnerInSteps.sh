#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false false
#./processHistogramsInStepsPp.sh "data/" "data/" false false false
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false false
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false false
./processHistogramsInStepsAsymmetryPp.sh "data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_negativeVz_wtaAxis_preprocessed_2019-08-13.root" "data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_negativeVz_wtaAxis_allCorrections_processed_2019-08-13.root" true false false
./processHistogramsInStepsAsymmetryPp.sh "data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_positiveVz_wtaAxis_preprocessed_2019-08-13.root" "data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_positiveVz_wtaAxis_allCorrections_processed_2019-08-13.root" true false false
