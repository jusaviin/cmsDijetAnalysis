/*
 * Implementation of the methods
 */

// Own includes
#include "DijetMethods.h"

/*
 *  Constructor
 */
DijetMethods::DijetMethods() :
  fMixedEventFitRegion(0.2)
{
  // Constructor
}

/*
 *  Destructor
 */
DijetMethods::~DijetMethods()
{
  // Destructor
}

/*
 * Do the mixed event correction and return corrected TH2D.
 * The idea here is, that we divide the same event histogram with a normalized mixed event histogram.
 * The normalization is obtained by first projecting the deltaEta distribution out of the
 * two-dimensional histogram and then fitting a constant to the central region in deltaEta.
 * This gives the value of the highest bins in the two dimensional histogram while suppressing
 * single bin fluctuations.
 *
 *  TODO: In Hallie's code the scaling factor is the average of the leading and subleading factors.
 *        Decide if this is necessary or if we can do the correction separately.
 *        Also there is some smoothening of the mixed event distribution by taking average of
 *        different sides of deltaEta to suppress fluctuations on edges. See if something like
 *        this needs to be implemented here.
 *
 * Arguments:
 *  TH2D* sameEventHistogram = Histogram with correlation from the same event
 *  TH2D* mixedEventHistogram = Histogram with correlation from different events
 *
 *  return: Corrected same event histogram
 */
TH2D* DijetMethods::MixedEventCorrect(TH2D *sameEventHistogram, TH2D *mixedEventHistogram){
  
  // Clone the same event histogram for correction
  char newName[100];
  sprintf(newName,"%sCorrected",sameEventHistogram->GetName());
  TH2D* correctedHistogram = (TH2D*) sameEventHistogram->Clone(newName);
  
  // In the 2D histograms deltaPhi is x-axis and deltaEta y-axis. We need deltaEta for the correction
  TH1D *hDeltaEtaMixed = mixedEventHistogram->ProjectionY("MixedDeltaEtaProjection",1,mixedEventHistogram->GetNbinsX());
  
  // Use a constant fit function to fit the projected histogram
  hDeltaEtaMixed->Fit("pol0","0Q","",-fMixedEventFitRegion,fMixedEventFitRegion);
  
  // The normalization scale is the fit result divided by the number of deltaPhi bins integrated for one deltaEta bin
  double scale = hDeltaEtaMixed->GetFunction("pol0")->GetParameter(0) / mixedEventHistogram->GetNbinsX();
  
  // Normalize the mixed event histogram and do the correction
  mixedEventHistogram->Scale(1.0/scale);
  correctedHistogram->Divide(mixedEventHistogram);
  
  // Return the corrected histogram
  return correctedHistogram;
}

// Setter for deltaEta range used for normalizing the mixed event
void DijetMethods::SetMixedEventFitRegion(const double etaRange){
  fMixedEventFitRegion = etaRange;
}
