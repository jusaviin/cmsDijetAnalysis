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
 * This same is done for leading and subleading jet-track mixed event histograms and an average
 * is taken from these. This is done to ensure that the normalization is the same, because both
 * leading and subleading distributions are needed for the background subrtaction.
 *
 *  TODO: There is some smoothening of the mixed event distribution by taking average of
 *        different sides of deltaEta to suppress fluctuations on edges in Hallie's code.
 *        See if something like this needs to be implemented here.
 *
 * Arguments:
 *  TH2D* sameEventHistogram = Histogram with correlation from the same event
 *  TH2D* leadingMixedEventHistogram = Leading jet-mixed event track histogram
 *  TH2D* subleadingMixedEventHistogram = Subleading jet-mixed event track histogram
 *
 *  return: Corrected same event histogram
 */
TH2D* DijetMethods::MixedEventCorrect(TH2D *sameEventHistogram, TH2D *leadingMixedEventHistogram, TH2D *subleadingMixedEventHistogram){
  
  // Clone the same event histogram for correction
  char newName[100];
  sprintf(newName,"%sCorrected",sameEventHistogram->GetName());
  TH2D* correctedHistogram = (TH2D*) sameEventHistogram->Clone(newName);
  
  // Calculate the average scale from leading and subleading event scales
  double leadingScale = GetMixedEventScale(leadingMixedEventHistogram);
  double subleadingScale = GetMixedEventScale(subleadingMixedEventHistogram);
  double averageScale = (leadingScale+subleadingScale)/2.0;
  
  // Normalize the mixed event histogram and do the correction
  leadingMixedEventHistogram->Scale(1.0/averageScale);
  correctedHistogram->Divide(leadingMixedEventHistogram);
  
  // We need to rescale the mixed event histogram in the end, because it might later be used for another correction
  leadingMixedEventHistogram->Scale(averageScale);
  
  // Return the corrected histogram
  return correctedHistogram;
}

/*
 * Find the scale to normalize a mixed event histogram
 *
 *  Arguments:
 *   TH2D* mixedEventHistogram = Mixed event histogram from which the scale is seeked
 *
 *   return: Scale to be used for normalization
 */
double DijetMethods::GetMixedEventScale(TH2D* mixedEventHistogram){
  
  // In the 2D histograms deltaPhi is x-axis and deltaEta y-axis. We need deltaEta for the correction
  TH1D *hDeltaEtaMixed = mixedEventHistogram->ProjectionY("MixedDeltaEtaProjection",1,mixedEventHistogram->GetNbinsX());
  
  // Use a constant fit function to fit the projected histogram
  hDeltaEtaMixed->Fit("pol0","0Q","",-fMixedEventFitRegion,fMixedEventFitRegion);
  
  // The normalization scale is the fit result divided by the number of deltaPhi bins integrated for one deltaEta bin
  return (hDeltaEtaMixed->GetFunction("pol0")->GetParameter(0) / mixedEventHistogram->GetNbinsX());
}

/*
 * Subtract the background from the given leading TH2D distribution.
 * It is assumed that the large deltaEta region defined by variables fMinBackgroundDeltaEta
 * and fMaxBackgroundDeltaEta is not correlated with the jet on the near side. The deltaPhi
 * distribution is projected out from this deltaEta region from two histograms, one for the
 * leading jet-track correlations and one for the subleading-jet track correlations. A new
 * distribution is generated where every deltaEta value in a deltaPhi strip has the same value,
 * given by the leading distribution projected in the background region for the near side and
 * and the near side of the subleading distribution for the away side. This way the wide away
 * side jet peak is avoided and we can generate the background for the whole deltaPhi region.
 * This background is then subtracted from the distribution to get the signal distribution.
 *
 * Arguments:
 *  TH2D *leadingHistogramWithBackground = Leading-jet track correlation histogram with deltaPhi as x-axis and deltaEta as y-axis
 *  TH2D *subleadingHistogramWithBackground = Subleading-jet track correlation histogram with deltaPhi as x-axis and deltaEta as y-axis
 *
 *  return: Background subtracted leading jet-track correlation histogram
 */
TH2D* DijetMethods::SubtractBackground(TH2D *leadingHistogramWithBackground, TH2D *subleadingHistogramWithBackground){
  
  // Start by projecting the deltaPhi distribution from the leading and subleading histograms in the background region
  TH1D *backgroundDeltaPhiLeading = ProjectBackgroundDeltaPhi(leadingHistogramWithBackground);
  TH1D *backgroundDeltaPhiSubleading = ProjectBackgroundDeltaPhi(subleadingHistogramWithBackground);
  
  // Construct the two-dimensional background histogram by populating the whole deltaEta region from background deltaPhi histogram
  char histogramName[200]; // Helper variable for histogram naming
  double deltaPhiValueNear;    // Value in the leading deltaPhi histogram bin
  double deltaPhiErrorNear;    // Error of the leading deltaPhi histogram bin
  double deltaPhiValueAway;    // Value in the subleading deltaPhi histogram bin
  double deltaPhiErrorAway;    // Error of the subleading deltaPhi histogram bin
  
  sprintf(histogramName,"%sBackground",leadingHistogramWithBackground->GetName());
  fBackgroundDistribution = (TH2D*)leadingHistogramWithBackground->Clone(histogramName);
  int nDeltaPhiBins = fBackgroundDistribution->GetNbinsX(); // Number of deltaPhi bins

  // Loop over deltaPhi bins and fill the leading jet-track correlation result in the near side
  // and the subleading set-track correlation result in the away side
  for(int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins/2; iDeltaPhi++){
    
    // Read the values from the deltaPhi background histogram
    deltaPhiValueNear = backgroundDeltaPhiLeading->GetBinContent(iDeltaPhi);
    deltaPhiErrorNear = backgroundDeltaPhiLeading->GetBinError(iDeltaPhi);
    deltaPhiValueAway = backgroundDeltaPhiSubleading->GetBinContent(iDeltaPhi);
    deltaPhiErrorAway = backgroundDeltaPhiSubleading->GetBinError(iDeltaPhi);
    
    // Loop over deltaEta bins
    for(int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
      
      // Insert the values to the two-dimensional background histogram
      fBackgroundDistribution->SetBinContent(iDeltaPhi,iDeltaEta,deltaPhiValueNear);
      fBackgroundDistribution->SetBinError(iDeltaPhi,iDeltaEta,deltaPhiErrorNear);
      fBackgroundDistribution->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,deltaPhiValueAway);
      fBackgroundDistribution->SetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,deltaPhiErrorAway);
      
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
  sprintf(histogramName,"%sBackgroundSubtracted",leadingHistogramWithBackground->GetName());
  TH2D *backgroundSubtractedHistogram = (TH2D*)leadingHistogramWithBackground->Clone();
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

/*
 * Project the deltaPhi distribution out of a two-dimensional deltaPhi-deltaEta distribution
 *
 *  Arguments:
 *   TH2D* deltaPhiDeltaEtaHistogram = two-dimensional deltaPhi-deltaEta distribution
 *
 *   return: DeltaPhi distribution projected from the background deltaEta region
 */
TH1D* DijetMethods::ProjectBackgroundDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram){
  
  // Start by finding the bin indices for the defined deltaEta region
  // Apply a little offset for the defined eta region borders to avoid bin border effects
  int lowNegativeDeltaEtaBin  = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(-fMaxBackgroundDeltaEta+0.0001);
  int highNegativeDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(-fMinBackgroundDeltaEta-0.0001);
  int lowPositiveDeltaEtaBin  = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(fMinBackgroundDeltaEta+0.0001);
  int highPositiveDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(fMaxBackgroundDeltaEta-0.0001);
  
  // Calculate the number of deltaEta bins in the background region
  int nBinsBackgroundRegion = highNegativeDeltaEtaBin-lowNegativeDeltaEtaBin+highPositiveDeltaEtaBin-lowPositiveDeltaEtaBin+2;
  
  // Project out the deltaPhi distribution in the background
  char histogramName[200];
  sprintf(histogramName,"%sBackgroundDeltaPhi",deltaPhiDeltaEtaHistogram->GetName());
  TH1D *backgroundDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName,lowNegativeDeltaEtaBin,highNegativeDeltaEtaBin);
  backgroundDeltaPhi->Add(deltaPhiDeltaEtaHistogram->ProjectionX("dummyName",lowPositiveDeltaEtaBin,highPositiveDeltaEtaBin));
  
  // Scale the projected deltaPhi distribution with the number of deltaEta bins projected over to retain normalization
  backgroundDeltaPhi->Scale(1.0/nBinsBackgroundRegion);
  
  // Return the scaled projection
  return backgroundDeltaPhi;
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
