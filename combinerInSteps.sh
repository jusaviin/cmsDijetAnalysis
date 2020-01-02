#!/bin/bash

#./combineHistogramsInSteps.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsPp.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsAsymmetry.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsAsymmetryPp.sh "data/" "data/" "data/" 1 1

./combineHistogramsInSteps.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_jet100trig_xjBins_wtaAxis_improvisedMixing_highJetCut_preprocessed_2019-10-13.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_jet100trig_xjBins_wtaAxis_improvisedMixing_lowJetCut_preprocessed_2019-10-13.root" "rofl.root" 1 1
