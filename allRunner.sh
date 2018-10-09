#!/bin/bash

./processHistograms.sh "data/dijet_pp_highForest_pfJets_2018-09-14.root" "data/dijet_pp_highForest_pfJets_processed_2018-09-14.root" true
./processHistograms.sh "data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_2018-09-15.root" "data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root" true
#./processHistograms.sh "data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_2018-09-15.root" "data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root" false
#./processHistograms.sh "data/dijet_ppMC_GenReco_mergedSkims_Pythia6_2018-07-27.root" "data/dijet_ppMC_GenReco_mergedSkims_Pythia6_processed_2018-08-16.root" false
#./processHistograms.sh "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_2018-09-15.root" "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root" false
