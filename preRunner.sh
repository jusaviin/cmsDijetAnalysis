#!/bin/bash

./processHistograms.sh "data/dijet_pp_highForest_pfJets_pfCandAxis_2018-12-07.root" "data/dijet_pp_highForest_pfJets_pfCandAxis_smoothedMixing_noCorrections_processed_2018-12-07.root" false false
./processHistograms.sh "data/dijet_pp_highForest_pfJets_pfCandAxis_2018-12-07.root" "data/dijet_pp_highForest_pfJets_pfCandAxis_smoothedMixing_processed_2018-12-07.root" true false
