#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1

./processHistogramsInStepsPp.sh "data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_xjBins_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-10-17.root" "data/ppData2017_JECv4_jetSpectra_processed_2019-10-17.root" false false -1
