#!/bin/bash

#./processHistogramsInSteps.sh "data/PbPbMC_GenGen_skims_pfJets_noUncorrected_3eventsMixed_sube0_2018-10-09.root" "data/PbPbMC_GenGen_skims_pfJets_noInclusiveOrUncorrected_3eventsMixed_sube0_smoothedMixing_processed_2018-10-30.root" false false
#./processHistogramsInSteps.sh "data/PbPbMC_GenGen_skims_pfJets_noInclusiveOrUncorrected_3eventsMixed_subeNon0_2018-10-01.root" "data/PbPbMC_GenGen_skims_pfJets_noInclusiveOrUncorrected_3eventsMixed_subeNon0_processed_2018-10-01.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_skims_pfJets_noUncorrected_3eventsMixed_sube0_2018-10-09.root" "data/PbPbMC_RecoGen_skims_pfJets_noInclusiveOrUncorrected_3eventsMixed_sube0_smoothedMixing_processed_2018-10-30.root" false false
./processHistogramsInSteps.sh "data/PbPbMC_RecoGen_skims_pfJets_noUncorrected_3eventsMixed_subeNon0_2018-10-09.root" "data/PbPbMC_RecoGen_skims_pfJets_noInclusiveOrUncorrected_3eventsMixed_subeNon0_smoothedMixing_processed_2018-10-30.root" false false