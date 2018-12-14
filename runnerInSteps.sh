#!/bin/bash

./processHistogramsInSteps.sh "data/PbPbMC_RecoReco_skims_pfJets_pfCandAxis_noMixing_2018-12-06.root" "data/PbPbMC_RecoReco_skims_pfJets_pfCandAxis_noMixing_processed_2018-12-06.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_skims_pfJets_pfCandAxis_noMixing_2018-12-06.root" "data/PbPbMC_RecoGen_skims_pfJets_pfCandAxis_noMixing_processed_2018-12-06.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_GenReco_skims_pfJets_pfCandAxis_noMixing_2018-12-06.root" "data/PbPbMC_GenReco_skims_pfJets_pfCandAxis_noMixing_processed_2018-12-06.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_GenGen_skims_pfJets_pfCandAxis_noMixing_2018-12-06.root" "data/PbPbMC_GenGen_skims_pfJets_pfCandAxis_noMixing_processed_2018-12-06.root" false false
