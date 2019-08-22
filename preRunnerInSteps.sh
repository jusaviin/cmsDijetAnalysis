#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1
./processHistogramsInStepsAsymmetryPp.sh "data/ppData2017_highForest_pfJets_improvisedMixing_JECv2_wtaAxisPpreprocessed_2019-08-13.root" "data/ppData2017_highForest_pfJets_improvisedMixing_JECv2_wtaAxis_allCorrections_processed_2019-08-13.root" true false -1
