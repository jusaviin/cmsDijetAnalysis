#!/bin/bash

./processHistogramsInSteps.sh "data/PbPbMC_GenGen_skims_pfJets_noMixing_noJetLimit_2019-01-10.root" "data/PbPbMC_GenGen_skims_pfJets_noMixing_noJetLimit_noCorrelations_processed_2019-01-10.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_GenReco_skims_pfJets_noMixing_noJetLimit_2019-01-10.root" "data/PbPbMC_GenReco_skims_pfJets_noMixing_noJetLimit_noCorrelations_processed_2019-01-10.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_skims_pfJets_noMixing_noJetLimit_2019-01-10.root" "data/PbPbMC_RecoGen_skims_pfJets_noMixing_noJetLimit_noCorrelations_processed_2019-01-10.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_RecoReco_skims_pfJets_noMixing_noJetLimit_2019-01-10.root" "data/PbPbMC_RecoReco_skims_pfJets_noMixing_noJetLimit_noCorrelations_processed_2019-01-10.root" false false

