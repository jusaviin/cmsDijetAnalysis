#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsManual.sh "data/" "data/" false false -1 "data/"

#./processHistogramsInStepsAsymmetryPp.sh "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20eveMix_quarkGluonCombine_25pMoreQuark_dijetWeight_JECv4_preprocessed_2020-01-10.root" "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20eveMix_quarkGluonCombine_25pMoreQuark_dijetWeight_JECv4_onlySeagull_processed_2020-01-10.root" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrectionsExceptSpillover_wtaAxis_processed_2020-02-06.root" true false -1 "data/"
./processHistogramsInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_tuningForSeagull_wtaAxis_processed_2020-02-04.root" true true -1 "data/"

#./processHistogramsInStepsManual.sh "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_xjBins_quarkGluonCombined_wta_subeNon0_centShift5_preprocessed_2019-12-13.root" "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5eveMix_quarkGluonCombined_wta_subeNon0_centShift5_seagullTuning_processed_2019-12-13.root" false false -1 "data/PbPbMC2018_RecoGen_akFlowJet_noUncOrInc_5pCentSh_5eveMix_quarkGluon_jet100Tri_preprocessed_2019-12-20.root"


#./processHistogramsInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "oneBinToRuleThemAll.root" true true -1 "data/"

#./processHistogramsInStepsAsymmetryPp.sh "data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_preprocessed_2020-01-31.root" "data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-04-02.root" true false -1
