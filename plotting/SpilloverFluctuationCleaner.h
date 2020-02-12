#ifndef SPILLOVERFLUCTUATIONCLEANER_H
#define SPILLOVERFLUCTUATIONCLEANER_H

// Own includes
#include "DijetHistogramManager.h"

/*
 * SpilloverFluctuationCleaner class
 *
 * Class whose purpose is to have te configuration to clean fluctuating bins caused by the spillover correction
 */
class SpilloverFluctuationCleaner {
  
public:
 
  SpilloverFluctuationCleaner();   // Contructor
  ~SpilloverFluctuationCleaner();  // Destructor
  
  void CleanSpilloverFluctuationDeltaR(TH1* jetShapeHistogram, const int iJetTrack, const int iAsymmetry, const int iCentrality, const int iTrackPt) const;

};

#endif
