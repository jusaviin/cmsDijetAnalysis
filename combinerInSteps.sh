#!/bin/bash

#./combineHistogramsInSteps.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsPp.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsAsymmetry.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsAsymmetryPp.sh "data/" "data/" "data/" 1 1

./combineHistogramsInStepsAsymmetry.sh "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentSh_5eveMix_onlyQuark_jet100Tri_preprocessed_2019-12-20.root" "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentSh_5eveMix_xjBins_onlyGluon_jet100Tri_preprocessed_2019-12-20.root" "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentSh_5eveMix_xjBins_quarkGluonCombine_25pMoreQuark_jet100Tri_preprocessed_2019-12-20.root" 1.25 0.75
