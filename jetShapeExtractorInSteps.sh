#!/bin/bash

#./extractJetShapeInSteps.sh "data/" "data/"
#./extractJetShapeInStepsPp.sh "data/" "data/"
#./extractJetShapeInStepsAsymmetry.sh "data/" "data/"
#./extractJetShapeInStepsAsymmetryPp.sh "data/" "data/"

./extractJetShapeInStepsAsymmetry.sh "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_allCorrectionsForClosure_processed_2020-09-15.root" "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_allCorrectionsForClosure_onlyFinalResults_processed_2020-09-15.root"
./extractJetShapeInStepsAsymmetry.sh "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_processed_2019-10-12.root" "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_onlyFinalResults_processed_2019-10-12.root"

#./extractJetShapeInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_fixSeagull_allCorrections_processed_2020-02-17.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_finalTuning_onlyJetShapes_processed_2020-02-17.root"

#./extractJetShapeInSteps.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_allPtTrackDeltaR_processed_2019-10-16.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_allPtTrackDeltaR_onlyJetShape_processed_2019-10-16.root"

#./extractJetShapeInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_tuningForSeagull_wtaAxis_processed_2020-02-04.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_onlyJetShapa_manualTuning_wtaAxis_processed_2020-02-04.root"
