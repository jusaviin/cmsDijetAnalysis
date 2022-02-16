#!/bin/bash

#./combineHistogramsInSteps.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsPp.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsAsymmetry.sh "data/" "data/" "data/" 1 1
#./combineHistogramsInStepsAsymmetryPp.sh "data/" "data/" "data/" 1 1

#./combineHistogramsInStepsAsymmetry.sh "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentSh_5eveMix_onlyQuark_jet100Tri_preprocessed_2019-12-20.root" "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentSh_5eveMix_xjBins_onlyGluon_jet100Tri_preprocessed_2019-12-20.root" "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentSh_5eveMix_xjBins_quarkGluonCombine_25pMoreQuark_jet100Tri_preprocessed_2019-12-20.root" 1.25 0.75

./combineHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_multWeight_subeNon0_jetEta1v3_onlyQuarkJets_improvisedMixing_preprocessed_2022-01-25_part0.root" "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_multWeight_subeNon0_jetEta1v3_onlyGluonJets_improvisedMixing_preprocessed_2022-01-25_part0.root" "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_multWeight_subeNon0_jetEta1v3_25pMoreQuark_improvisedMixing_preprocessed_2022-01-25_part0.root" 1.25 0.75
./combineHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_multWeight_jetEta1v3_onlyQuarkJets_improvisedMixing_preprocessed_2022-01-25_part0.root" "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_multWeight_jetEta1v3_onlyGluonJets_improvisedMixing_preprocessed_2022-01-25_part0.root" "data/PbPbMC2018_RecoGen_akCaloJet_dihadron_multWeight_jetEta1v3_25pMoreQuark_improvisedMixing_preprocessed_2022-01-25_part0.root" 1.25 0.75
