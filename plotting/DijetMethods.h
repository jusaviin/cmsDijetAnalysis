#ifndef DIJETMETHODS_H
#define DIJETMETHODS_H

/*
 * This class is a collection of methods that are used to process the
 * results produced by the dijet analysis
 */

// C++ includes
#include <iostream>

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>

class DijetMethods{
  
public:
  
  DijetMethods();   // Constructor
  DijetMethods(const DijetMethods& in); // Copy constructor
  ~DijetMethods();  // Destructor
  DijetMethods& operator=(const DijetMethods& in); // Equal sign operator
  
  TH2D* MixedEventCorrect(TH2D *sameEventHistogram, TH2D *leadingMixedEventHistogram, TH2D *subleadingMixedEventHistogram); // Mixed event correction for a two-dimensional histogram
  TH2D* SubtractBackground(TH2D *leadingHistogramWithBackground, TH2D *subleadingHistogramWithBackground); // Subtract background from a two-dimensional leading histogram
  
  // Getter for the most recent background distribution
  TH2D* GetBackground() const; // Getter for the most recent background distribution used to subtract the background
  
  // Setters for mixed event configuration
  void SetMixedEventFitRegion(const double etaRange);  // Setter for deltaEta range used for normalizing the mixed event
  
  // Setters for background subtraction configuration
  void SetBackgroundDeltaEtaRegion(const double minDeltaEta, const double maxDeltaEta); // Setter for background deltaEta region
  
private:
  
  // =============================================
  // =========== Mixed event correction ==========
  // =============================================
  
  double fMixedEventFitRegion;  // Region in deltaEta in which a constant fit is done to normalize mixed event distributions
  
  // =============================================
  // =========== Background subtraction ==========
  // =============================================
  
  TH2D *fBackgroundDistribution;  // Remember the background distribution from the most recent background subtraction
  double fMinBackgroundDeltaEta;  // Minimum deltaEta for background subtraction region
  double fMaxBackgroundDeltaEta;  // Maximum deltaEta for background subtraction region
  
  // Private methods
  double GetMixedEventScale(TH2D* mixedEventHistogram); // Find the normalization scale for the mixed event histogram
  TH1D* ProjectBackgroundDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram); // Project deltaPhi distribution out of a two-dimensional deltaPhi-deltaEta distribution
  
};

#endif
