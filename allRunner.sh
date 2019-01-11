#!/bin/bash

#./processHistograms.sh "data/dijet_pp_highForest_pfJets_2018-09-14.root" "data/dijet_pp_highForest_pfJets_processed_2018-09-14.root" true true
#./processHistograms.sh "data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_2018-09-15.root" "data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root" true true
./processHistograms.sh "data/dijetPbPb_skims_pfJets_noUncorrected_10mixedEvents_2018-11-11.root" "data/dijetPbPb_skims_pfJets_noUncorrected_10mixedEvents_smoothedMixing_noCorrections_2019-01-07.root" false false
./processHistograms.sh "data/dijet_pp_highForest_pfJets_2018-12-03.root" "data/dijet_pp_highForest_pfJets_smoothedMixing_noCorrections_2019-01-07.root" false false
#./processHistograms.sh "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_2018-09-15.root" "data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root" false false
