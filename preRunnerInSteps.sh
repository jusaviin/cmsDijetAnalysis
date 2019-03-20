#!/bin/bash

#./processHistogramsInSteps.sh "data/" "data/" false false
#./processHistogramsInSteps.sh "data/PbPbMC_GenGen_skims_pfJets_noUncorr_noMixing_jetMatching_noJetLimit_2019-01-31.root" "data/PbPbMC_GenGen_skims_pfJets_noUncorr_noMixing_jetMatching_noJetLimit_noCorrelations_processed_2019-01-31.root" false false
#./processHistogramsInSteps.sh "data/PbPbMC_GenReco_skims_pfJets_noUncorr_noMixing_noJetLimit_2019-01-29.root" "data/PbPbMC_GenReco_skims_pfJets_noUncorr_noMixing_noJetLimit_noCorrelations_processed_2019-01-29.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_GenGen_skims_pfJets_sube0_noUncorr_matchedJets_noMixing_2019-03-18.root" "data/PbPbMC_GenGen_skims_pfJets_sube0_noUncorr_matchedJets_improvisedMixing_processed_2019-03-18.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_skims_pfJets_sube0_noUncorr_matchedJets_noMixing_2019-03-18.root" "data/PbPbMC_RecoGen_skims_pfJets_sube0_noUncorr_matchedJets_improvisedMixing_processed_2019-03-18.root" false false
