#!/bin/bash

#./processHistogramsInSteps.sh "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_2018-10-19.root" "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_smoothedMixing_processed_2018-10-19.root" true true
#./processHistogramsInSteps.sh "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_2018-10-19.root" "data/dijetPbPb_pfJets_noInclusiveOrUncorrected_noCorrections_smoothedMixing_processed_2018-10-19.root" false false
./processHistogramsInSteps.sh "data/dijetPbPb_pfJets_3eventsMixed_noUncorrected_2018-10-02.root" "data/dijetPbPb_pfJets_3eventsMixed_noUncorrected_noCorrections_smoothedMixing_processed_2018-10-02.root" false false
