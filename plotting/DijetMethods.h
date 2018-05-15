#ifndef DIJETMETHODS_H
#define DIJETMETHODS_H

/*
 * This class is a collection of methods that are used to process the
 * results produced by the dijet analysis
 */

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>

class DijetMethods{
  
public:
  
  DijetMethods();   // Constructor
  ~DijetMethods();  // Destructor
  
  TH2D* MixedEventCorrect(TH2D *sameEventHistogram, TH2D *mixedEventHistogram); // Mixed event correction for a two dimensional histogram
  
  // Setters for mixed event configuration
  void SetMixedEventFitRegion(const double etaRange);  // Setter for deltaEta range used for normalizing the mixed event
  
private:
  
  // =============================================
  // =========== Mixed event correction ==========
  // =============================================
  
  double fMixedEventFitRegion;  // Region in deltaEta in which a constant fit is done to normalize mixed event distributions
  
};

#endif
