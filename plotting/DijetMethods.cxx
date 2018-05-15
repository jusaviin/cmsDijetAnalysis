/*
 * Implementation of the methods
 */

// Own includes
#include "DijetMethods.h"

/*
 *  Constructor
 */
DijetMethods::DijetMethods() :
  fMixedEventFitRegion(0.2),
  fMinBackgroundDeltaEta(1.5),
  fMaxBackgroundDeltaEta(2.5)
{
  // Constructor
  fBackgroundDistribution = NULL;
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

/*
 * Subtract the background from the given TH2D distribution.
 * It is assumed that the large deltaEta region defined by variables fMinBackgroundDeltaEta
 * and fMaxBackgroundDeltaEta is not correlated with the jet. Then the deltaPhi distribution
 * is projected out from this deltaEta region and a new two-dimensional distribution is
 * generated where every deltaEta value in a deltaPhi strip has the same value, given by the
 * distribution projected in the bakcground region. This background is then subtracted from
 * the distribution to get the signal distribution.
 *
 * Arguments:
 *  TH2D *histogramWithBackground = Two-dimensional histogram with deltaPhi as x-axis and deltaEta as y-axis
 *
 *  return: Background subtracted histogram
 */
TH2D* DijetMethods::SubtractBackground(TH2D *histogramWithBackground){
  
  // TODO: Away side should not be used for background subtraction because of eta swing
  
  // Start by finding the bin indices for the defined deltaEta region
  // Apply a little offset for the defined eta region borders to avoid bin border effects
  int lowNegativeDeltaEtaBin  = histogramWithBackground->GetYaxis()->FindBin(-fMaxBackgroundDeltaEta+0.0001);
  int highNegativeDeltaEtaBin = histogramWithBackground->GetYaxis()->FindBin(-fMinBackgroundDeltaEta-0.0001);
  int lowPositiveDeltaEtaBin  = histogramWithBackground->GetYaxis()->FindBin(fMinBackgroundDeltaEta+0.0001);
  int highPositiveDeltaEtaBin = histogramWithBackground->GetYaxis()->FindBin(fMaxBackgroundDeltaEta-0.0001);
 
  // Calculate the number of deltaEta bins in the background region
  int nBinsBackgroundRegion = highNegativeDeltaEtaBin-lowNegativeDeltaEtaBin+highPositiveDeltaEtaBin-lowPositiveDeltaEtaBin+2;
  
  // Project out the deltaPhi distribution in the background region in near side
  char histogramName[200];
  sprintf(histogramName,"%sBackgroundDeltaPhi",histogramWithBackground->GetName());
  TH1D *backgroundDeltaPhi = histogramWithBackground->ProjectionX(histogramName,lowNegativeDeltaEtaBin,highNegativeDeltaEtaBin);
  backgroundDeltaPhi->Add(histogramWithBackground->ProjectionX("dummyName",lowPositiveDeltaEtaBin,highPositiveDeltaEtaBin));
  
  // Scale the projected deltaPhi distribution with the number of deltaEta bins projected over to retain normalization
  backgroundDeltaPhi->Scale(1.0/nBinsBackgroundRegion);
  
  // Construct the two-dimensional background histogram by populating the whole deltaEta region from background deltaPhi histogram
  double deltaPhiValue;   // Value in the deltaPhi histogram bin
  double deltaPhiError;   // Error of the deltaPhi histogram bin
  sprintf(histogramName,"%sBackground",histogramWithBackground->GetName());
  fBackgroundDistribution = (TH2D*)histogramWithBackground->Clone(histogramName);
  
  // Loop over deltaPhi bins
  for(int iDeltaPhi = 1; iDeltaPhi <= fBackgroundDistribution->GetNbinsX(); iDeltaPhi++){
    
    // Read the values from the deltaPhi background histogram
    deltaPhiValue = backgroundDeltaPhi->GetBinContent(iDeltaPhi);
    deltaPhiError = backgroundDeltaPhi->GetBinError(iDeltaPhi);
    
    // Loop over deltaEta bins
    for(int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
      
      // Insert the values to the two-dimensional background histogram
      fBackgroundDistribution->SetBinContent(iDeltaPhi,iDeltaEta,deltaPhiValue);
      fBackgroundDistribution->SetBinError(iDeltaPhi,iDeltaEta,deltaPhiError);
      
      /*
       * Note: In Hallie's code there is a multiplication by sqrt(bins) for the bin error here
       *       I do not think it is needed, since the error should shrink as we calculate average
       *       and then assign this average to each bin. The sqrt(bins) gives the average error in
       *       each bin you are likely to see, but I do not think this is what we want here. If
       *       we just fill the histogram with the average, we should use the error of the average
       *       which comes without sqrt(bins).
       */
      
    } // DeltaEta loop
  } // DeltaPhi loop
  
  // Subtract the generated background from the histogram including the background
  sprintf(histogramName,"%sBackgroundSubtracted",histogramWithBackground->GetName());
  TH2D *backgroundSubtractedHistogram = (TH2D*)histogramWithBackground->Clone();
  backgroundSubtractedHistogram->Add(fBackgroundDistribution,-1);
  
  // Set all the negative bins to zero
  for(int iDeltaPhi = 1; iDeltaPhi <= backgroundSubtractedHistogram->GetNbinsX(); iDeltaPhi++){
    for(int iDeltaEta = 1; iDeltaEta <= backgroundSubtractedHistogram->GetNbinsY(); iDeltaEta++){
      if(backgroundSubtractedHistogram->GetBinContent(iDeltaPhi,iDeltaEta) < 0){
        backgroundSubtractedHistogram->SetBinContent(iDeltaPhi,iDeltaEta,0);
      }
    } // DeltaEta loop
  } // DeltaPhi loop
  
  // Return the background subtracted histogram with negative values removed
  return backgroundSubtractedHistogram;
}

// Getter for the most recent background distribution used to subtract the background
TH2D* DijetMethods::GetBackground() const{
  return fBackgroundDistribution;
}

// Setter for deltaEta range used for normalizing the mixed event
void DijetMethods::SetMixedEventFitRegion(const double etaRange){
  fMixedEventFitRegion = etaRange;
}

// Setter for background deltaEta region
void DijetMethods::SetBackgroundDeltaEtaRegion(const double minDeltaEta, const double maxDeltaEta){
  fMinBackgroundDeltaEta = minDeltaEta;
  fMaxBackgroundDeltaEta = maxDeltaEta;
  
  // Check that minimum value is smaller tham maximum value
  if(fMinBackgroundDeltaEta > fMaxBackgroundDeltaEta){
    fMinBackgroundDeltaEta = maxDeltaEta;
    fMaxBackgroundDeltaEta = minDeltaEta;
  }
}
