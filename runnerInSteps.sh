#!/bin/bash

#./processHistogramsInSteps.sh "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_2018-10-19.root" "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_smoothedMixing_processed_2018-10-19.root" true true
#./processHistogramsInSteps.sh "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_2018-10-19.root" "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_noCorrections_smoothedMixing_processed_2018-10-19.root" false false
./processHistogramsInSteps.sh "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_2018-10-19.root" "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_noCorrections_smoothedMixing_processed_2018-11-19.root" false false
./processHistogramsInSteps.sh "data/dijetPbPb_skims_pfJets_noUncorrected_10mixedEvents_2018-11-11.root" "data/dijetPbPb_pfJets_skims_noUncorrected_10mixedEvents_noCorrections_smoothedMixing_processed_2018-11-19.root" false false
