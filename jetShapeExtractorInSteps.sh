#!/bin/bash

#./extractJetShapeInSteps.sh "data/" "data/"
#./extractJetShapeInStepsPp.sh "data/" "data/"
#./extractJetShapeInStepsAsymmetry.sh "data/" "data/"
#./extractJetShapeInStepsAsymmetryPp.sh "data/" "data/"

#./extractJetShapeInSteps.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_allPtTrackDeltaR_processed_2019-10-16.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_allPtTrackDeltaR_onlyJetShape_processed_2019-10-16.root"

#./extractJetShapeInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_tuningForSeagull_wtaAxis_processed_2020-02-04.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_onlyJetShapa_manualTuning_wtaAxis_processed_2020-02-04.root"
./extractJetShapeInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_tuningForSeagull_wtaAxis_processed_2020-02-04.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_onlyJetShapa_manualTuning_wtaAxis_processed_2020-02-04.root"