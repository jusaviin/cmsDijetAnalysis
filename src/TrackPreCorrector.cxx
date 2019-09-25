#include "TrackPreCorrector.h"

/*
 * Constructor for the precorrector.
 *
 *  TString fileName = File from which the the precorrector is initialized
 */
TrackPreCorrector::TrackPreCorrector(TString fileName){
  TFile *file = TFile::Open(fileName);
  for(Int_t iCentrality = 0; iCentrality < fNCentralityBins; iCentrality++){
    for(Int_t iTrackPt = 0; iTrackPt < fNTrackPtBins; iTrackPt++){
      fMinimumBiasRatios[iCentrality][iTrackPt] = (TH2D*) file->Get(Form("minimumBiasTrackRatio_C%dT%d", iCentrality, iTrackPt));
    }
  }
}

/*
 * Destructor
 */
TrackPreCorrector::~TrackPreCorrector()
{
  // Destructor
}

/*
 * Find the correct bin from the given values. Unnecessarily complicated function from Xiao, but seems to work.
 *
 *  Double_t key = Value for which the correct bin is needed
 *  Double_t *array = Array holding the bin boundaries for the searhed value
 *  Int_t iMax = Maximum index used for the searching
 *  Int_t iMin = Minimum index used for the searching
 *
 *  return: Index in the array correcponding to key
 */
Int_t TrackPreCorrector::BinarySearch(Double_t key, Double_t *array, Int_t iMax, Int_t iMin ){
  if(key > array[iMax] ) return -1;
  if(key < array[iMin] ) return -1;
  Int_t mid = floor(Double_t(iMax + iMin)/2);
  //	cout<<mid<<endl;
  if(mid == iMin ) return mid;
  if(array[mid] > key) return BinarySearch(key, array, mid, iMin);
  else if( array[mid] < key) return BinarySearch(key, array, iMax, mid);
  else return mid;
}

/*
 * Get the weight for the track to match 2018 to 2015. Used in temporary correction
 *
 *  Double_t pt = Track pT
 *  Double_t eta = Track eta
 *  Double_t phi = Track phi
 *  Int_t centrality = Event centrality
 *
 *  return: Track weigth matching the 2018 distribution to 2015 distribution
 */
Double_t TrackPreCorrector::GetTrackWeight(Double_t pt, Double_t eta, Double_t phi, Double_t centrality){
  int iTrackPt = BinarySearch(pt, fTrackPtBins, fNTrackPtBins, 0);
  int iCentrality = BinarySearch(centrality, fCentralityBins, fNCentralityBins, 0);
  if(iTrackPt < 0 || iCentrality < 0) {
    std::cout << "iTrackPt = " << iTrackPt << ", iCentrality = " << iCentrality << std::endl;
    std::cout << "Track pT = " << pt << ", centrality = " << centrality << std::endl;
    return 0;
  }
  return fMinimumBiasRatios[iCentrality][iTrackPt]->GetBinContent(fMinimumBiasRatios[iCentrality][iTrackPt]->FindBin(phi,eta));
}

