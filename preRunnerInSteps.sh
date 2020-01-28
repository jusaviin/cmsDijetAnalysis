#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsManual.sh "data/" "data/" false false -1 "data/"

#./processHistogramsInStepsAsymmetryPp.sh "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20eveMix_quarkGluonCombine_25pMoreQuark_dijetWeight_JECv4_preprocessed_2020-01-10.root" "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20eveMix_quarkGluonCombine_25pMoreQuark_dijetWeight_JECv4_onlySeagull_processed_2020-01-10.root" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_genJetTrackCorrection_wtaAxis_processed_2020-01-24.root" true true -1 "data/"
./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6_allCorrectionsExectTrackingDeltaR_wtaAxis_processed_2020-01-24.root" true true -1 "data/"
