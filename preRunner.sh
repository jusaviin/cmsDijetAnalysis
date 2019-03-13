#!/bin/bash

#./processHistograms.sh "data/" "data/" false false
#./processHistograms.sh "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_noJetLimit_2019-01-15.root" "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_noJetLimit_smoothedMixing_adjustedBackground_processed_2019-01-15.root" false false
#./processHistograms.sh "data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_noJetLimit_2019-01-15.root" "data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_noJetLimit_smoothedMixing_adjustedBackground_processed_`2019-01-15.root" false false
./processHistograms.sh "data/dijet_ppMC_GenGen_Pythia6_pfJets_noUncorr_matchedJets_WTAaxis_2019-03-08.root" "data/dijet_ppMC_GenGen_Pythia6_pfJets_noUncorr_matchedJets_WTAaxis_processed_2019-03-08.root" false false
./processHistograms.sh "data/dijet_ppMC_RecoGen_Pythia6_pfJets_noUncorr_matchedJets_WTAaxis_2019-03-08.root" "data/dijet_ppMC_RecoGen_Pythia6_pfJets_noUncorr_matchedJets_WTAaxis_processed_2019-03-08.root" false false
