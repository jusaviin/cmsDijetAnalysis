#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1 "data/"
./processHistogramsInStepsPp.sh "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noCorrelations_dijetWeight_JECv4_2019-10-31.root" "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noCorrelations_dijetWeight_JECv4_processed_2019-10-31.root" false false -1

#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_preprocessed_2019-10-07.root" "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlySeagullAndTracking_fineTunedRange_mixingScale16_processed_2019-10-07.root" false false -1 "data/"

#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_preprocessed_2019-10-07.root" "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_onlySeagullUptoHigh_mixingScale16_processed_2019-10-07.root" false false -1 "data/"

#./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_preprocessed_2019-09-26.root" "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSuboNon0_wtaAxis_onlySpilloverAlsoTunedSubleading_JECv6_processed_2019-09-26.root" false true -1 "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_subeNon0_wtaAxis_JECv6_preprocessed_2019-09-26.root"

#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5pShiftedCent_5eveMix_jet100Trigger_preprocessed_2019-10-10.root" "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_5pShiftedCent_5eveMix_jet100Trigger_noCorrections_processed_2019-10-24.root" false false -1 "data/"

#./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_preprocessed_2019-09-26.root" "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_mixingFromSubeNon0AtLowPt_newTry_wtaAxis_JECv6_noCorrections_processed_2019-09-26.root" false false -1 "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_subeNon0_wtaAxis_JECv6_preprocessed_2019-09-26.root"
#./processHistogramsInStepsAtHighPt.sh "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_preprocessed_2019-09-26.root" "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_mixingFromSubeNon0AtLowPt_newTry_wtaAxis_JECv6_noCorrections_processed_2019-09-26.root" false false -1
