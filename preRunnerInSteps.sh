#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsAsymmetry.sh "data/" "data/" false false -1 "data/"
#./processHistogramsInStepsAsymmetryPp.sh "data/" "data/" false false -1
#./processHistogramsInStepsManual.sh "data/" "data/" false false -1 "data/"

./processHistogramsInStepsManual.sh "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_onlyJets_jetWeight_JECv6_2021-05-11.root" "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_onlyJets_jetWeight_JECv6_processed_2021-05-11.root" false false -1 "data/"
./processHistogramsInStepsManual.sh "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_triggered_jetWeight_JECv6_2021-05-11.root" "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_triggered_jetWeight_JECv6_processed_2021-05-11.root" false false -1 "data/"
./processHistogramsInStepsManual.sh "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_triggered_noJetWeight_JECv6_2021-05-11.root" "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_triggered_noJetWeight_JECv6_processed_2021-05-11.root" false false -1 "data/"

#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_noQcut_doubleDijet_preprocessed_2021-04-28.root" "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_noQcut_doubleDijet_noCorrections_processed_2021-04-28.root" false false -1 "data/"
#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove1p8_doubleDijet_preprocessed_2021-04-28.root" "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove1p8_doubleDijet_noCorrections_processed_2021-04-28.root" false false -1 "data/"
#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove2_doubleDijet_preprocessed_2021-04-28.root" "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove2_doubleDijet_noCorrections_processed_2021-04-28.root" false false -1 "data/"
#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove2p2_doubleDijet_preprocessed_2021-04-28.root" "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove2p2_doubleDijet_noCorrections_processed_2021-04-28.root" false false -1 "data/"
#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove2p5_doubleDijet_preprocessed_2021-04-28.root" "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove2p5_doubleDijet_noCorrections_processed_2021-04-28.root" false false -1 "data/"
#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove2p8_doubleDijet_preprocessed_2021-04-28.root" "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove2p8_doubleDijet_noCorrections_processed_2021-04-28.root" false false -1 "data/"
#./processHistogramsInSteps.sh "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove3p3_doubleDijet_preprocessed_2021-04-28.root" "data/PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_improvisedMixing_qVectorAbove3p3_doubleDijet_noCorrections_processed_2021-04-28.root" false false -1 "data/"


#./processHistogramsInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_10eveMix_20pSmear_wtaAxis_preprocessed_2020-04-04.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_10eveMix_20pSmear_wtaAxis_allCorrectionsUnsmeared_processed_2020-04-04.root" true true -1 "data/"

#./processHistogramsInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_fixSeagull_allCorrections_processed_2020-02-17.root" true true -1 "data/"

#./processHistogramsInStepsManual.sh "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_4pShiftedCent_5eveMix_jet100Trigger_onlySeagull_processed_2019-10-10.root" "data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_4pShiftedCent_5eveMix_jet100Trigger_allCorrections_tuning_processed_2019-10-21.root" true true -1 "data/"

#./processHistogramsInStepsManual.sh "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_preprocessed_2019-09-28.root" "data/ppMC2017_RecoReco_Pythia8_pfJets_xjBins_wtaAxis_noUncorr_20EventsMixed_JECv4_tuning_processed_2019-12-04.root" false false -1 "data/"

#./processHistogramsInStepsManual.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_preprocessed_2019-10-06.root" "oneBinWithJff.root" true true -1 "data/"

#./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6minus_wtaAxis_preprocessed_2019-12-16.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6minus_wtaAxis_allCorrections_processed_2020-02-14.root" true true -1 "data/"

#./processHistogramsInStepsAsymmetry.sh "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6plus_wtaAxis_preprocessed_2019-12-29.root" "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6plus_wtaAxis_allCorrections_processed_2020-02-14.root" true true -1 "data/"
