#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsManual.sh "data/" "data/" false false -1 "data/"

./processHistogramsInSteps.sh "rofl.root" "rofl_preprocessed.root" false false 2 "data/"

#./processHistogramsInSteps.sh "data/dihadronPbPb2018_leadingTrigger3-4_subleadingTrigger4-8_preprocessed_2020-03-16.root" "data/dihadronPbPb2018_leadingTrigger3-4_subleadingTrigger4-8_noCorrections_processed_2020-03-16.root" false false -1 "data/"

#./processHistogramsInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_fixSeagull_allCorrections_processed_2020-02-17.root" true true -1 "data/"

#./processHistogramsInStepsManual.sh "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_onlySeagull_processed_2019-10-10.root" "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_tuning_processed_2019-10-21.root" true true -1 "data/"

#./processHistogramsInStepsManual.sh "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_preprocessed_2019-09-28.root" "data/ppMC2017_RecoReco_Pythia8_pfJets_xjBins_wtaAxis_noUncorr_20EventsMixed_JECv4_tuning_processed_2019-12-04.root" false false -1 "data/"

#./processHistogramsInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "oneBinWithJff.root" true true -1 "data/"

#./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6minus_wtaAxis_preprocessed_2019-12-16.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6minus_wtaAxis_allCorrections_processed_2020-02-14.root" true true -1 "data/"

#./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6plus_wtaAxis_preprocessed_2019-12-29.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6plus_wtaAxis_allCorrections_processed_2020-02-14.root" true true -1 "data/"
