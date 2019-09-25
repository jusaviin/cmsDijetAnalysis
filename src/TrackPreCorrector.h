#ifndef TRACKPRECORRECTOR_H
#define TRACKPRECORRECTOR_H


#include "TH2D.h"
#include "TFile.h"
#include <iostream>

class TrackPreCorrector{
  
public:
  TrackPreCorrector(TString filePath);
  virtual ~TrackPreCorrector();
  double GetTrackWeight(Double_t pt, Double_t eta, Double_t phi, Double_t centrality);
private:
  Int_t BinarySearch(Double_t key, Double_t *array, Int_t iMax, Int_t iMin );
  TH2D* fMinimumBiasRatios[4][7];
  Int_t fNTrackPtBins = 7;
  Int_t fNCentralityBins = 4;
  Double_t fTrackPtBins[8] = {0.7,1,2,3,4,8,12,999};
  Double_t fCentralityBins[5] = {0,10,30,50,100};
};

#endif
