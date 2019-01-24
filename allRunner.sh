#!/bin/bash

#./processHistograms.sh "data/dijet_pp_highForest_pfJets_2018-09-14.root" "data/dijet_pp_highForest_pfJets_processed_2018-09-14.root" true true
#./processHistograms.sh "data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_2018-09-15.root" "data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root" true true
#./processHistograms.sh "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_2018-09-15.root" "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_newProcessing_processed_2018-09-15.root" false false
#./processHistograms.sh "data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_2019-01-14.root" "data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_noJffCorrection_processed_2019-01-14.root" false false
./processHistograms.sh "data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_2019-01-14.root" "data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_noCorrections_adjustedBackground_processed_2019-01-14.root" false false
