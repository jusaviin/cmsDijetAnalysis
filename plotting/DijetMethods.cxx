/*
 * Implementation of the methods
 */

// Own includes
#include "DijetMethods.h"

/*
 * Combination of zeroth and first order polynomial for seagull correction
 * Functional form: f(x) = c, if |x| < 0.5  <=> f(x) = a*x^2+b*x+c if |x| > 0.5
 */
double seagullPoly2(double *x, double *par){
  if(x[0] < 0.5 && x[0] > -0.5) return par[0];
  return par[0]+par[1]*x[0]*x[0]+par[2]*x[0];
}

/*
 * Combination of zeroth and first order polynomial for seagull correction
 * Functional form: f(x) = a + b*e^(c*x)
 */
double seagullExp(double *x, double *par){
  return par[0]+par[1]*TMath::Exp(x[0]*par[2]);
  //return par[0]+x[0]*par[1]+x[0]*x[0]*par[2]+x[0]*x[0]*x[3]*par[3];
}

/*
 * Combination of two Gauss functions of different widths without overlap
 */
double doubleGaussNoOverlap(double *x, double *par){
  if(x[0] < 0.4 && x[0] > -0.4) return par[0]/(TMath::Sqrt(2*TMath::Pi())*par[1])*TMath::Exp(-0.5*TMath::Power(x[0]/par[1],2));
  return par[2]/(TMath::Sqrt(2*TMath::Pi())*par[3])*TMath::Exp(-0.5*TMath::Power(x[0]/par[3],2));
}

/*
 * Combination of two Gauss functions of different widths without overlap
 */
double doubleGaussNoOverlap2D(double *x, double *par){
  if(TMath::Sqrt(x[0]*x[0]+x[1]*x[1]) < 0.4) return par[0]/(2*TMath::Pi()*par[1]*par[2])*TMath::Exp(-0.5*TMath::Power(x[0]/par[1],2))*TMath::Exp(-0.5*TMath::Power(x[1]/par[2],2));
  return par[3]/(2*TMath::Pi()*par[4]*par[5])*TMath::Exp(-0.5*TMath::Power(x[0]/par[4],2))*TMath::Exp(-0.5*TMath::Power(x[1]/par[5],2));
}

/*
 *  Constructor
 */
DijetMethods::DijetMethods() :
  fMixedEventFitRegionLow(-0.2),
  fMixedEventFitRegionHigh(0.2),
  fMixedEventNormalizationMethod(kSingle),
  fSmoothMixing(false),
  fMaximumDeltaEta(4),
  fMinBackgroundDeltaEta(1.5),
  fMaxBackgroundDeltaEta(2.5),
  fAdjustBackground(false),
  fnOverlapBins(3),
  fBackgroundErrorScalingFactor(1),
  fMinBackgroundDeltaPhi(1.4),
  fMaxBackgroundDeltaPhi(1.7),
  fSeagullRebin(4),
  fSeagullChi2Limit(1.15),
  fnBinsProjectedOver(0),
  fMaxSignalDeltaEta(1.0),
  fSpilloverYieldError(0),
  fJetShapeNormalizationMethod(kBinWidth),
  fnRebinDeltaEta(0),
  fnRebinDeltaPhi(0),
  fPairAcceptancePositiveLevel(0),
  fPairAcceptanceNegativeLevel(0),
  fBackgroundSubtractionInnerMean(0),
  fBackgroundSubtractionOuterMean(0)
{
  // Constructor
  fNormalizedMixedEventHistogram = NULL;
  fBackgroundDistribution = NULL;
  fBackgroundOverlap = NULL;
  fSpilloverDeltaEta = NULL;
  fSpilloverDeltaPhi = NULL;
  fSpilloverFitDeltaEta = NULL;
  fSpilloverFitDeltaPhi = NULL;
  fhJetShapeCounts = NULL;
  fhJetShapeBinMap = NULL;
  fBackgroundEtaProjection = NULL;
  fSeagullFit = NULL;
  fRebinDeltaEta = nullptr;
  fRBins = nullptr;
  fRebinDeltaPhi = nullptr;
  
  // Set default RBins for jet shape
  const int nRBins = 16;  // Number of R-bins for jet shape histograms
  double rBins[nRBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0,1.25,1.5}; // R-bin boundaries for jet shape histogram
  SetJetShapeBinEdges(nRBins,rBins);
}

/*
 *  Copy constructor
 */
DijetMethods::DijetMethods(const DijetMethods& in) :
  fNormalizedMixedEventHistogram(in.fNormalizedMixedEventHistogram),
  fMixedEventFitRegionLow(in.fMixedEventFitRegionLow),
  fMixedEventFitRegionHigh(in.fMixedEventFitRegionHigh),
  fMixedEventNormalizationMethod(in.fMixedEventNormalizationMethod),
  fSmoothMixing(in.fSmoothMixing),
  fMaximumDeltaEta(in.fMaximumDeltaEta),
  fBackgroundDistribution(in.fBackgroundDistribution),
  fBackgroundOverlap(in.fBackgroundOverlap),
  fMinBackgroundDeltaEta(in.fMinBackgroundDeltaEta),
  fMaxBackgroundDeltaEta(in.fMaxBackgroundDeltaEta),
  fAdjustBackground(in.fAdjustBackground),
  fnOverlapBins(in.fnOverlapBins),
  fBackgroundErrorScalingFactor(in.fBackgroundErrorScalingFactor),
  fBackgroundEtaProjection(in.fBackgroundEtaProjection),
  fSeagullFit(in.fSeagullFit),
  fMinBackgroundDeltaPhi(in.fMinBackgroundDeltaPhi),
  fMaxBackgroundDeltaPhi(in.fMaxBackgroundDeltaPhi),
  fSeagullRebin(in.fSeagullRebin),
  fSeagullChi2Limit(in.fSeagullChi2Limit),
  fnBinsProjectedOver(in.fnBinsProjectedOver),
  fMaxSignalDeltaEta(in.fMaxSignalDeltaEta),
  fSpilloverDeltaEta(in.fSpilloverDeltaEta),
  fSpilloverDeltaPhi(in.fSpilloverDeltaPhi),
  fSpilloverFitDeltaEta(in.fSpilloverFitDeltaEta),
  fSpilloverFitDeltaPhi(in.fSpilloverFitDeltaPhi),
  fSpilloverYieldError(in.fSpilloverYieldError),
  fJetShapeNormalizationMethod(in.fJetShapeNormalizationMethod),
  fhJetShapeCounts(in.fhJetShapeCounts),
  fhJetShapeBinMap(in.fhJetShapeBinMap),
  fPairAcceptancePositiveLevel(in.fPairAcceptancePositiveLevel),
  fPairAcceptanceNegativeLevel(in.fPairAcceptanceNegativeLevel),
  fBackgroundSubtractionInnerMean(in.fBackgroundSubtractionInnerMean),
  fBackgroundSubtractionOuterMean(in.fBackgroundSubtractionOuterMean)
{
  // Copy constructor
  SetBinBoundaries(in.fnRBins,in.fRBins,fnRBins,&fRBins);
  SetBinBoundaries(in.fnRebinDeltaPhi,in.fRebinDeltaPhi,fnRebinDeltaPhi,&fRebinDeltaPhi);
  SetBinBoundaries(in.fnRebinDeltaEta,in.fRebinDeltaEta,fnRebinDeltaEta,&fRebinDeltaEta);
}

/*
 * Assingment operator
 */
DijetMethods& DijetMethods::operator=(const DijetMethods& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fNormalizedMixedEventHistogram = in.fNormalizedMixedEventHistogram;
  fMixedEventFitRegionLow = in.fMixedEventFitRegionLow;
  fMixedEventFitRegionHigh = in.fMixedEventFitRegionHigh;
  fMixedEventNormalizationMethod = in.fMixedEventNormalizationMethod;
  fSmoothMixing = in.fSmoothMixing;
  fMaximumDeltaEta = in.fMaximumDeltaEta;
  fBackgroundDistribution = in.fBackgroundDistribution;
  fBackgroundOverlap = in.fBackgroundOverlap;
  fMinBackgroundDeltaEta = in.fMinBackgroundDeltaEta;
  fMaxBackgroundDeltaEta = in.fMaxBackgroundDeltaEta;
  fAdjustBackground = in.fAdjustBackground;
  fnOverlapBins = in.fnOverlapBins;
  fBackgroundErrorScalingFactor = in.fBackgroundErrorScalingFactor;
  fBackgroundEtaProjection = in.fBackgroundEtaProjection;
  fSeagullFit = in.fSeagullFit;
  fMinBackgroundDeltaPhi = in.fMinBackgroundDeltaPhi;
  fMaxBackgroundDeltaPhi = in.fMaxBackgroundDeltaPhi;
  fSeagullRebin = in.fSeagullRebin;
  fSeagullChi2Limit = in.fSeagullChi2Limit;
  fnBinsProjectedOver = in.fnBinsProjectedOver;
  fMaxSignalDeltaEta = in.fMaxSignalDeltaEta;
  fJetShapeNormalizationMethod = in.fJetShapeNormalizationMethod;
  fSpilloverDeltaEta = in.fSpilloverDeltaEta;
  fSpilloverDeltaPhi = in.fSpilloverDeltaPhi;
  fSpilloverFitDeltaEta = in.fSpilloverFitDeltaEta;
  fSpilloverFitDeltaPhi = in.fSpilloverFitDeltaPhi;
  fSpilloverYieldError = in.fSpilloverYieldError;
  fhJetShapeCounts = in.fhJetShapeCounts;
  fhJetShapeBinMap = in.fhJetShapeBinMap;
  fPairAcceptancePositiveLevel = in.fPairAcceptancePositiveLevel;
  fPairAcceptanceNegativeLevel = in.fPairAcceptanceNegativeLevel;
  fBackgroundSubtractionInnerMean = in.fBackgroundSubtractionInnerMean;
  fBackgroundSubtractionOuterMean = in.fBackgroundSubtractionOuterMean;
  
  SetBinBoundaries(in.fnRBins,in.fRBins,fnRBins,&fRBins);
  SetBinBoundaries(in.fnRebinDeltaPhi,in.fRebinDeltaPhi,fnRebinDeltaPhi,&fRebinDeltaPhi);
  SetBinBoundaries(in.fnRebinDeltaEta,in.fRebinDeltaEta,fnRebinDeltaEta,&fRebinDeltaEta);
  
  return *this;
}

/*
 *  Destructor
 */
DijetMethods::~DijetMethods()
{
  // Destructor
  if(fRBins) delete [] fRBins;
  if(fRebinDeltaEta) delete [] fRebinDeltaEta;
  if(fRebinDeltaPhi) delete [] fRebinDeltaPhi;
}

/*
 * Improvise a mixed event distribution by looking at the same event distribution between the deltaPhi peaks and
 * expanding the distribution in this region over the whole deltaPhi space.
 *
 * Arguments:
 *  const TH2D *sameEventHistogram = Histogram with correlation from the same event
 *
 *  return: Improvised mixed event distribution
 */
TH2D* DijetMethods::ImproviseMixedEvent(const TH2D *sameEventHistogram){
  
  // First, project the deltaEta from the background deltaPhi region
  fBackgroundEtaProjection = ProjectRegionDeltaEta(sameEventHistogram,fMinBackgroundDeltaPhi,fMaxBackgroundDeltaPhi,"improvisedMixingDeltaEta");
  
  // After the projection is done, propagate the eta values back to the length of whole two-dimensional distribution
  char histogramName[200];     // Helper variable for histogram naming
  sprintf(histogramName,"%sImprovisedBackground",sameEventHistogram->GetName());
  TH2D *improvisedMixedEventHistogram = (TH2D*)sameEventHistogram->Clone(histogramName);
  int nDeltaPhiBins = sameEventHistogram->GetNbinsX(); // Number of deltaPhi bins
  int nDeltaEtaBins = sameEventHistogram->GetNbinsY(); // Number of deltaEta bins
  double binWidthDeltaPhi = sameEventHistogram->GetXaxis()->GetBinWidth(1); // Bin width of a deltaPhi bin
  
  double deltaEtaValue, deltaEtaError;
  
  for(int iDeltaEta = 1; iDeltaEta <= nDeltaEtaBins; iDeltaEta++){
    
    // Read the values from the deltaEta background histogram
    deltaEtaValue = fBackgroundEtaProjection->GetBinContent(iDeltaEta);
    deltaEtaError = fBackgroundEtaProjection->GetBinError(iDeltaEta);
    
    // Repopulate the deltaPhi axis
    for(int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins; iDeltaPhi++){
      
      // Insert the values to the two-dimensional improvised mixed event histogram
      improvisedMixedEventHistogram->SetBinContent(iDeltaPhi,iDeltaEta,deltaEtaValue/binWidthDeltaPhi);
      improvisedMixedEventHistogram->SetBinError(iDeltaPhi,iDeltaEta,deltaEtaError*TMath::Sqrt(fnBinsProjectedOver)/binWidthDeltaPhi);
      
    } // deltaPhi loop
  } // deltaEta loop
  
  return improvisedMixedEventHistogram;
  
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
 * Arguments:
 *  const TH2D* sameEventHistogram = Histogram with correlation from the same event
 *  const TH2D* leadingMixedEventHistogram = Leading jet-mixed event track histogram
 *  const TH2D* subleadingMixedEventHistogram = Subleading jet-mixed event track histogram
 *
 *  return: Corrected same event histogram
 */
TH2D* DijetMethods::MixedEventCorrect(const TH2D *sameEventHistogram, const TH2D *leadingMixedEventHistogram, const TH2D *subleadingMixedEventHistogram, const bool avoidPeaks){

  // Clone the same event histogram for correction
  char newName[100];
  sprintf(newName,"%sCorrected",sameEventHistogram->GetName());
  TH2D* correctedHistogram = (TH2D*) sameEventHistogram->Clone(newName);
  
  // Set bins above maximum deltaEta to zero. This helps to suppress fluctuations at high deltaEta
  int minimumDeltaEtaBin = correctedHistogram->GetYaxis()->FindBin(-fMaximumDeltaEta-0.001);
  int maximumDeltaEtaBin = correctedHistogram->GetYaxis()->FindBin(fMaximumDeltaEta+0.001);
  for(int iDeltaPhi = 1; iDeltaPhi <= correctedHistogram->GetNbinsX(); iDeltaPhi++){
    for(int iDeltaEta = 1; iDeltaEta <= minimumDeltaEtaBin; iDeltaEta++){
      correctedHistogram->SetBinContent(iDeltaPhi,iDeltaEta,0);
      correctedHistogram->SetBinError(iDeltaPhi,iDeltaEta,0);
    }
    for(int iDeltaEta = maximumDeltaEtaBin; iDeltaEta <= correctedHistogram->GetNbinsY(); iDeltaEta++){
      correctedHistogram->SetBinContent(iDeltaPhi,iDeltaEta,0);
      correctedHistogram->SetBinError(iDeltaPhi,iDeltaEta,0);
    }
  }
  
  // Calculate the average scale from leading and subleading event scales
  double leadingScale = GetMixedEventScale(leadingMixedEventHistogram,avoidPeaks);
  double subleadingScale = 0;
  if(fMixedEventNormalizationMethod == kAverage) subleadingScale = GetMixedEventScale(subleadingMixedEventHistogram,avoidPeaks);
  double averageScale = (leadingScale+subleadingScale)/2.0;
  double normalizationScale = leadingScale;
  if(fMixedEventNormalizationMethod == kAverage) normalizationScale = averageScale;

  // Normalize the mixed event histogram and do the correction
  sprintf(newName,"%sNormalized",leadingMixedEventHistogram->GetName());
  fNormalizedMixedEventHistogram = (TH2D*)leadingMixedEventHistogram->Clone(newName);
  fNormalizedMixedEventHistogram->Scale(1.0/normalizationScale);
  
  // If configured to smoothen the mixing distribution, take an average in each phi slice
  // Do not do this if there is chance to expand peak to flat region
  // TODO: For 2018, we might want to smoothen high deltaEta and leave low deltaEta be.
  // Will have to see how the mixed events look like with new jet corrections
  if(fSmoothMixing && !avoidPeaks){
    
    // First project out phi and normalize the projection to one in the central region
    TH1D *etaProjection = leadingMixedEventHistogram->ProjectionY("_eta");
    int lowFitBin = etaProjection->FindBin(fMixedEventFitRegionLow);
    int highFitBin = etaProjection->FindBin(fMixedEventFitRegionHigh);
    double mean = 0;
    for(int i = lowFitBin; i <= highFitBin; i++) mean += etaProjection->GetBinContent(i);
    mean /= ((highFitBin-lowFitBin+1)*leadingMixedEventHistogram->GetNbinsX());
    etaProjection->Scale(1.0/mean);
    
    // Then fill to the mixing histogram to each phi bin the average value over all phi
    for(int iEta = 1; iEta <= leadingMixedEventHistogram->GetNbinsY(); iEta++){
      for(int iPhi = 1; iPhi <= leadingMixedEventHistogram->GetNbinsX(); iPhi++){
        if(iEta<=highFitBin && iEta>=lowFitBin){
          fNormalizedMixedEventHistogram->SetBinContent(iPhi,iEta,1);
          fNormalizedMixedEventHistogram->SetBinError(iPhi,iEta,0);
        } else {
          fNormalizedMixedEventHistogram->SetBinContent(iPhi, iEta, etaProjection->GetBinContent(iEta)/leadingMixedEventHistogram->GetNbinsX());
          fNormalizedMixedEventHistogram->SetBinError(iPhi, iEta, etaProjection->GetBinError(iEta)/sqrt(leadingMixedEventHistogram->GetNbinsX()));
        }
      }
    }
  }
  
  correctedHistogram->Divide(fNormalizedMixedEventHistogram);
  
  // Return the corrected histogram
  return correctedHistogram;
}

/*
 * Find the scale to normalize a mixed event histogram
 *
 *  Arguments:
 *   const TH2D* mixedEventHistogram = Mixed event histogram from which the scale is seeked
 *   const bool findPeak = Find the highest 2x2 peak on the histogram and use that as normalization scale
 *
 *   return: Scale to be used for normalization
 */
double DijetMethods::GetMixedEventScale(const TH2D* mixedEventHistogram, const bool findPeak){
  
  // In the 2D histograms deltaPhi is x-axis and deltaEta y-axis. We need deltaEta for the correction
  TH1D *hDeltaEtaMixed;
  int nBinsProjectedOver;
  double binSum;
  double maxSum = 0;
  
  if(findPeak){
    
    // Find the highest 2x2 bin peak value from the histogram and use that as the normalization scale
    for(int iDeltaPhi = 1; iDeltaPhi < mixedEventHistogram->GetNbinsX(); iDeltaPhi++){
      for(int iDeltaEta = 1; iDeltaEta < mixedEventHistogram->GetNbinsY(); iDeltaEta++){
        binSum = mixedEventHistogram->GetBinContent(iDeltaPhi,iDeltaEta);
        binSum += mixedEventHistogram->GetBinContent(iDeltaPhi+1,iDeltaEta);
        binSum += mixedEventHistogram->GetBinContent(iDeltaPhi,iDeltaEta+1);
        binSum += mixedEventHistogram->GetBinContent(iDeltaPhi+1,iDeltaEta+1);
        if(binSum > maxSum) maxSum = binSum;
      }
    }
    binSum = maxSum;
    nBinsProjectedOver = 4;
    
  } else {
    
    // Use the whole DeltaPhi and central region in DeltaEta to get the normalization for mixed event
    hDeltaEtaMixed = mixedEventHistogram->ProjectionY("MixedDeltaEtaProjection",1,mixedEventHistogram->GetNbinsX());
    nBinsProjectedOver = mixedEventHistogram->GetNbinsX();
    
    // Use a constant fit function to fit the projected histogram
    hDeltaEtaMixed->Fit("pol0","0Q","",fMixedEventFitRegionLow,fMixedEventFitRegionHigh);
    binSum = hDeltaEtaMixed->GetFunction("pol0")->GetParameter(0);
  }
  
  // The normalization scale is the fit result divided by the number of deltaPhi bins integrated for one deltaEta bin
  return (binSum / nBinsProjectedOver);
}

/*
 * Apply a seagull correction to the histogram
 *
 * After an ideal mixed event correction, the corrected histogram will have a flat deltaEta distribution at large deltaPhi.
 * This is not the case in real life, but the so called eta wings emerge, meaning that there is more yield at large deltaEta.
 * The purpose of the seagull correction is to correct for this fact and make the background eta distribution flat.
 * The correction assumes that the deltaEta shape does not depend on deltaPhi.
 *
 * Algorithm used to make the correction:
 *  1) Project out the deltaEta distribution at background deltaPhi region
 *  2) Fit a combination of a zeroth and second order polynomial to the distribution
 *      -> f(x) = c, if |x| < 0.5  <=> f(x) = a*x^2+b*x+c if |x| > 0.5
 *  3) Interpret c in f(x) as the true level of deltaEta (the level at zero)
 *  4) Use the ratio c/f(deltaEta) as a correction factor for each point in the histogram
 *
 * Arguments:
 *  const TH2D *mixedEventCorrectedHistogram = Two-dimensional deltaEta-deltaPhi distribution, that is corrected by the mixed event
 *  const int normalizationMethod = 0: Assume no dip in the middle. 1: Assume dip in the middle and symmetrize background deltaEta.
 *  const int vetoFlag = 0: Regular correction. 1: Skip constant check. 2: Skip correction
 *
 *  return: The two dimensional distribution corrected for the seagull effect
 */
TH2D* DijetMethods::DoSeagullCorrection(const TH2D *mixedEventCorrectedHistogram, const int normalizationMethod, const int vetoFlag){
  
  // Project out the deltaEta distribution at background deltaPhi region
  int minDeltaPhiBin = mixedEventCorrectedHistogram->GetXaxis()->FindBin(fMinBackgroundDeltaPhi+0.0001);
  int maxDeltaPhiBin = mixedEventCorrectedHistogram->GetXaxis()->FindBin(fMaxBackgroundDeltaPhi+0.0001);
  fBackgroundEtaProjection = mixedEventCorrectedHistogram->ProjectionY("SeagullEtaProjection",minDeltaPhiBin,maxDeltaPhiBin);
  
  // Scale the projected deltaEta distribution with the number of deltaPhi bins projected over and deltaPhi bin width to retain normalization
  fBackgroundEtaProjection->Scale(1.0/(maxDeltaPhiBin-minDeltaPhiBin+1));
  fBackgroundEtaProjection->Scale(mixedEventCorrectedHistogram->GetXaxis()->GetBinWidth(1));
  fBackgroundEtaProjection->Rebin(fSeagullRebin);
  fBackgroundEtaProjection->Scale(1.0/fSeagullRebin);
  
  // Clone the original histogram
  TH2D *seagullCorrectedHistogram = (TH2D*) mixedEventCorrectedHistogram->Clone();
  
  // Check if a constant fit gives a good description of the background deltaEta distribution
  // In this case, do not apply the correction
  if(vetoFlag != 1){
    fBackgroundEtaProjection->Fit("pol0","","",-3,3);
    fSeagullFit = fBackgroundEtaProjection->GetFunction("pol0");
    double chi2PerNdf = fSeagullFit->GetChisquare() / fSeagullFit->GetNDF();
    if(chi2PerNdf < fSeagullChi2Limit || vetoFlag == 2) return seagullCorrectedHistogram;
    
    // Remove the constant fit from the list of functions if it is deemed bad
    fBackgroundEtaProjection->RecursiveRemove(fSeagullFit);
  }
  
  // Symmetrize the positive and negative sides of the distribution
  if(normalizationMethod == 1){
    double negativeContent, negativeError;
    double positiveContent, positiveError;
    double averageContent, averageError;
    int nDeltaEtaBins = fBackgroundEtaProjection->GetNbinsX();
    for(int iBin = 1; iBin <= nDeltaEtaBins/2; iBin++){
      
      // Average the contents from the positive and negative sides of the histogram
      negativeContent = fBackgroundEtaProjection->GetBinContent(iBin);
      negativeError   = fBackgroundEtaProjection->GetBinError(iBin);
      positiveContent = fBackgroundEtaProjection->GetBinContent(nDeltaEtaBins-iBin+1);
      positiveError   = fBackgroundEtaProjection->GetBinError(nDeltaEtaBins-iBin+1);
      averageContent  = (positiveContent+negativeContent)/2.0;
      averageError    = (positiveError+negativeError)/2.0;
      
      // Set the average as the new content for both sides
      fBackgroundEtaProjection->SetBinContent(iBin,averageContent);
      fBackgroundEtaProjection->SetBinError(iBin,averageError);
      fBackgroundEtaProjection->SetBinContent(nDeltaEtaBins-iBin+1,averageContent);
      fBackgroundEtaProjection->SetBinError(nDeltaEtaBins-iBin+1,averageError);
    }
  }
  
  // Prepare the fit function
  fSeagullFit = new TF1("seagullFit",seagullPoly2,-3,3,3);
  double initialLevel = fBackgroundEtaProjection->GetBinContent(fBackgroundEtaProjection->FindBin(0));
  fSeagullFit->SetParameters(initialLevel,0,0);
  fBackgroundEtaProjection->Fit(fSeagullFit,"","",-3,3);
  
  // Fit the projected distribution
  double backgroundLevel = fSeagullFit->GetParameter(0);
  
  // If we want to use exponential in some bins, redifine the seagull fit
  // Note that we want to have the same background level estimation in both cases
  if(normalizationMethod == 1){
    fBackgroundEtaProjection->RecursiveRemove(fSeagullFit);
    fSeagullFit = new TF1("seagullFitExp",seagullExp,-3,3,3);
    initialLevel = fBackgroundEtaProjection->GetBinContent(fBackgroundEtaProjection->FindBin(0));
    fSeagullFit->SetParameters(initialLevel,-1,-1);
    fBackgroundEtaProjection->Fit(fSeagullFit,"","",0,3);
    backgroundLevel = fSeagullFit->GetParameter(0); // TODO: Check if this or constant fit is beta
  }
  
  // Apply the correction to the input 2D histogram
  double binEta;
  double seagullCorrection;
  double binContent;
  double binError;
  
  for(int iEta = 1; iEta <= seagullCorrectedHistogram->GetNbinsY(); iEta++){
    
    // Calculate the correction for this deltaEta bin
    binEta = seagullCorrectedHistogram->GetYaxis()->GetBinCenter(iEta);
    if(normalizationMethod == 1) binEta = TMath::Abs(binEta);
    seagullCorrection = backgroundLevel/fSeagullFit->Eval(binEta);
    
    // Apply the correction for all deltaPhi bins in the deltaEta strip
    for(int iPhi = 1; iPhi <= seagullCorrectedHistogram->GetNbinsX(); iPhi++){
      binContent = seagullCorrectedHistogram->GetBinContent(iPhi,iEta);
      binError = seagullCorrectedHistogram->GetBinError(iPhi,iEta);
      seagullCorrectedHistogram->SetBinContent(iPhi,iEta,binContent*seagullCorrection);
      seagullCorrectedHistogram->SetBinError(iPhi,iEta,binError*seagullCorrection);
    }
  }
    
  // Return the corrected histogram
  return seagullCorrectedHistogram;
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
 *  bool isInclusive = Flag for inclusive correlations. True = inclusive. False = dijet.
 *
 *  return: Background subtracted leading jet-track correlation histogram
 */
TH2D* DijetMethods::SubtractBackground(TH2D *leadingHistogramWithBackground, TH2D *subleadingHistogramWithBackground, double maxDeltaEta, bool isInclusive){
  
  // Start by projecting the deltaPhi distribution from the leading and subleading histograms in the background region
  // For inclusive distribution, same inclusive distribution should be given as both input histograms
  TH1D *backgroundDeltaPhiLeading = ProjectBackgroundDeltaPhi(leadingHistogramWithBackground);
  TH1D *backgroundDeltaPhiSubleading = ProjectBackgroundDeltaPhi(subleadingHistogramWithBackground);
  
  // Construct the two-dimensional background histogram by populating the whole deltaEta region from background deltaPhi histogram
  char histogramName[200];     // Helper variable for histogram naming
  double deltaPhiValueNear;    // Value in the leading deltaPhi histogram bin
  double deltaPhiErrorNear;    // Error of the leading deltaPhi histogram bin
  double deltaPhiValueAway;    // Value in the subleading deltaPhi histogram bin
  double deltaPhiErrorAway;    // Error of the subleading deltaPhi histogram bin
  
  // Variables needed in the background level adjustment
  double leadingOverlap = 0;
  double leadingOverlapError = 0;
  double subleadingOverlap = 0;
  double subleadingOverlapError = 0;
  
  sprintf(histogramName,"%sBackground",leadingHistogramWithBackground->GetName());
  fBackgroundDistribution = (TH2D*)leadingHistogramWithBackground->Clone(histogramName);
  int nDeltaPhiBins = fBackgroundDistribution->GetNbinsX(); // Number of deltaPhi bins
  
  // Initialize also a distribution for background overlap for normalization level check
  sprintf(histogramName,"%sBackgroundOverlap",leadingHistogramWithBackground->GetName());
  fBackgroundOverlap = (TH2D*)leadingHistogramWithBackground->Clone(histogramName);
  
  // We need the deltaEta bin width for normalization purposes
  double binWidthDeltaEta = leadingHistogramWithBackground->GetYaxis()->GetBinWidth(1);
  int offset = isInclusive ? nDeltaPhiBins/2 : 0;  // Apply offset for inclusive histograms to scan over whole deltaPhi space

  // Do not apply background to region where there is no content
  int minFilledDeltaEtaBin = leadingHistogramWithBackground->GetYaxis()->FindBin(-maxDeltaEta+0.001);
  int maxFilledDeltaEtaBin = leadingHistogramWithBackground->GetYaxis()->FindBin(maxDeltaEta-0.001);
  
  // Loop over deltaPhi bins and fill the leading jet-track correlation result in the near side
  // and the subleading set-track correlation result in the away side
  for(int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins/2; iDeltaPhi++){
    
    // Read the values from the deltaPhi background histogram
    deltaPhiValueNear = backgroundDeltaPhiLeading->GetBinContent(iDeltaPhi);
    deltaPhiErrorNear = backgroundDeltaPhiLeading->GetBinError(iDeltaPhi);
    deltaPhiValueAway = backgroundDeltaPhiSubleading->GetBinContent(iDeltaPhi+offset);
    deltaPhiErrorAway = backgroundDeltaPhiSubleading->GetBinError(iDeltaPhi+offset);
    
    // Calculate the sum of leading and subleading backgrounds in the overlap region
    if(iDeltaPhi <= fnOverlapBins){
      subleadingOverlap += deltaPhiValueAway;
      subleadingOverlapError += deltaPhiErrorAway;
    }
    if(iDeltaPhi > (nDeltaPhiBins/2 - fnOverlapBins)){
      leadingOverlap += deltaPhiValueNear;
      leadingOverlapError += deltaPhiErrorNear;
    }
    
    // Repopulate the deltaEta axis
    for(int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
      
      // If there is no content in the histogram, do not subtract background from it
      if(iDeltaEta >= minFilledDeltaEtaBin && iDeltaEta <= maxFilledDeltaEtaBin){

        // Insert the values to the two-dimensional background histogram
        fBackgroundDistribution->SetBinContent(iDeltaPhi,iDeltaEta,deltaPhiValueNear/binWidthDeltaEta);
        fBackgroundDistribution->SetBinError(iDeltaPhi,iDeltaEta,deltaPhiErrorNear*TMath::Sqrt(fnBinsProjectedOver)/binWidthDeltaEta);
        fBackgroundDistribution->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,deltaPhiValueAway/binWidthDeltaEta);
        fBackgroundDistribution->SetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,deltaPhiErrorAway*TMath::Sqrt(fnBinsProjectedOver)/binWidthDeltaEta);
        
      }
      
      // Set the contants of the background overlap to zero
      fBackgroundOverlap->SetBinContent(iDeltaPhi,iDeltaEta,0);
      fBackgroundOverlap->SetBinError(iDeltaPhi,iDeltaEta,0);
      fBackgroundOverlap->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,0);
      fBackgroundOverlap->SetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,0);
      
      /*
       * Note: There is a multiplication by sqrt(bins) for the bin error to get the average error
       *       in each bin. We want the average error rather than error of the average here, because
       *       otherwise the errors become too small. You can see this the most clearly by plotting
       *       the distribution and see that adjacent bins fluctuate much more than would be expected
       *       from the error bars. Thus the sqrt(bins) to properly treat the errors.
       */
      
    } // DeltaEta loop
  } // DeltaPhi loop
  
  // For the background overlap distribution, fill a few bins on each side of the gluing point
  for(int iDeltaPhi = nDeltaPhiBins/2 + 1; iDeltaPhi <= nDeltaPhiBins/2 + fnOverlapBins; iDeltaPhi++){
    
    // Read the values from the deltaPhi background histogram
    deltaPhiValueNear = backgroundDeltaPhiLeading->GetBinContent(iDeltaPhi);
    deltaPhiErrorNear = backgroundDeltaPhiLeading->GetBinError(iDeltaPhi);
    
    // For the away side we need to continue to the left, which rotates to the few last bins in deltaPhi
    deltaPhiValueAway = backgroundDeltaPhiSubleading->GetBinContent(iDeltaPhi + nDeltaPhiBins/2 - fnOverlapBins);
    deltaPhiErrorAway = backgroundDeltaPhiSubleading->GetBinError(iDeltaPhi + nDeltaPhiBins/2 - fnOverlapBins);
    
    // Calculate the sum of leading and subleading backgrounds in the overlap region
    leadingOverlap += deltaPhiValueNear;
    leadingOverlapError += deltaPhiErrorNear;
    subleadingOverlap += deltaPhiValueAway;
    subleadingOverlapError += deltaPhiErrorAway;
    
    // Repopulate the deltaEta axis
    for(int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
      
      // Do not set anything for the overlap if there is no content in the original histogram
      if(iDeltaEta >= minFilledDeltaEtaBin && iDeltaEta <= maxFilledDeltaEtaBin){
        
        // Set the contants of the background overlap
        fBackgroundOverlap->SetBinContent(iDeltaPhi,iDeltaEta,deltaPhiValueNear/binWidthDeltaEta);
        fBackgroundOverlap->SetBinError(iDeltaPhi,iDeltaEta,deltaPhiErrorNear*TMath::Sqrt(fnBinsProjectedOver)/binWidthDeltaEta);
        fBackgroundOverlap->SetBinContent(iDeltaPhi-fnOverlapBins,iDeltaEta,deltaPhiValueAway/binWidthDeltaEta);
        fBackgroundOverlap->SetBinError(iDeltaPhi-fnOverlapBins,iDeltaEta,deltaPhiErrorAway*TMath::Sqrt(fnBinsProjectedOver)/binWidthDeltaEta);
        
      }
    } // DeltaEta loop
  } // DeltaPhi loop
  
  // In the optimal case the ratio of overlapping bins from leading and subleading backgrounds is 1.
  // We can adjust the subleading background level based on the sum of background level in the overlap
  // region to make this happen. Nothing to adjust for inclusive jet-track correlation
  double leadingScalingFactor = 1;
  double subleadingScalingFactor = 1;
  if(fAdjustBackground && !isInclusive){
    double leadingToSubleadingRatio = leadingOverlap/subleadingOverlap;
    double subleadingToLeadingRatio = subleadingOverlap/leadingOverlap;
    double leadingToSubleadingRatioError = TMath::Sqrt(TMath::Power(leadingOverlapError/subleadingOverlap,2)+TMath::Power(leadingOverlap*subleadingOverlapError/TMath::Power(subleadingOverlap,2),2));
    double scaledLeadingContent, scaledLeadingError;
    double scaledSubleadingContent, scaledSubleadingError;
    
    // Only do the adjustment if the error of the ratio is somewhat reasonable
    if(leadingToSubleadingRatioError < 1){
      
      // Scale both sides to match in the middle
      if(leadingToSubleadingRatio > 1){
        leadingScalingFactor = 1 + ((leadingToSubleadingRatio-1)/2);
        subleadingScalingFactor = 1 - ((1-subleadingToLeadingRatio)/2);
      } else {
        leadingScalingFactor = 1 - ((1-leadingToSubleadingRatio)/2);
        subleadingScalingFactor = 1 + ((subleadingToLeadingRatio-1)/2);
      }
      
      for(int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins/2; iDeltaPhi++){
        for(int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
          
          // If there is no content in the histogram, do not do adjustments
          if(iDeltaEta >= minFilledDeltaEtaBin && iDeltaEta <= maxFilledDeltaEtaBin){
            
            // Scale the actual distribution
            scaledLeadingContent = fBackgroundDistribution->GetBinContent(iDeltaPhi,iDeltaEta)/leadingScalingFactor;
            scaledLeadingError = fBackgroundDistribution->GetBinError(iDeltaPhi,iDeltaEta)/leadingScalingFactor;
            scaledSubleadingContent = fBackgroundDistribution->GetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta)/subleadingScalingFactor;
            scaledSubleadingError = fBackgroundDistribution->GetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta)/subleadingScalingFactor;
            fBackgroundDistribution->SetBinContent(iDeltaPhi,iDeltaEta,scaledLeadingContent);
            fBackgroundDistribution->SetBinError(iDeltaPhi,iDeltaEta,scaledLeadingError);
            fBackgroundDistribution->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,scaledSubleadingContent);
            fBackgroundDistribution->SetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,scaledSubleadingError);
            
            // Scale also the overlap histogram for debugging purposes
            // Note that for overlap leading and subleading sides change with respect to pi/2
            scaledLeadingContent = fBackgroundOverlap->GetBinContent(iDeltaPhi,iDeltaEta)/subleadingScalingFactor;
            scaledLeadingError = fBackgroundOverlap->GetBinError(iDeltaPhi,iDeltaEta)/subleadingScalingFactor;
            scaledSubleadingContent = fBackgroundOverlap->GetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta)/leadingScalingFactor;
            scaledSubleadingError = fBackgroundOverlap->GetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta)/leadingScalingFactor;
            fBackgroundOverlap->SetBinContent(iDeltaPhi,iDeltaEta,scaledLeadingContent);
            fBackgroundOverlap->SetBinError(iDeltaPhi,iDeltaEta,scaledLeadingError);
            fBackgroundOverlap->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,scaledSubleadingContent);
            fBackgroundOverlap->SetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,scaledSubleadingError);
            
          } // Bins that are actually filled
          
        } // DeltaEta loop
      } // DeltaPhi loop
    } // Reasonable error if
  } // Adjusting the background if
  
  // In the optimal case the difference in yields of overlapping bins from leading and subleading backgrounds is 0.
  // We can adjust the background levels to make this happen. Nothing to adjust for inclusive jet-track correlation
  /*if(fAdjustBackground && !isInclusive){
    double leadingToSubleadingRatio = leadingOverlap/subleadingOverlap;
    double subleadingToLeadingRatio = subleadingOverlap/leadingOverlap;
    double leadingToSubleadingRatioError = TMath::Sqrt(TMath::Power(leadingOverlapError/subleadingOverlap,2)+TMath::Power(leadingOverlap*subleadingOverlapError/TMath::Power(subleadingOverlap,2),2));
    leadingAdjustment = (leadingOverlap-subleadingOverlap)/2;
    subleadingAdjustment = (subleadingOverlap-leadingOverlap)/2;
    double errorAverage = (leadingOverlapError-subleadingOverlapError)/2;
    double scaledLeadingContent, scaledSubleadingContent;
    
    // Only do the adjustment if the error of the ratio is somewhat reasonable
    if(leadingToSubleadingRatioError < 1){
      
      // Scale both sides to match in the middle
      if(leadingToSubleadingRatio > 1){
        leadingScalingFactor = 1 + ((leadingToSubleadingRatio-1)/2);
        subleadingScalingFactor = 1 - ((1-subleadingToLeadingRatio)/2);
      } else {
        leadingScalingFactor = 1 - ((1-leadingToSubleadingRatio)/2);
        subleadingScalingFactor = 1 + ((subleadingToLeadingRatio-1)/2);
      }
      
      leadingAdjustment = (leadingOverlap/(fnOverlapBins*2))*leadingScalingFactor;
      subleadingAdjustment = (subleadingOverlap/(fnOverlapBins*2))*subleadingScalingFactor;
      
      for(int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins/2; iDeltaPhi++){
        for(int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
          
          // If there is no content in the histogram, do not do adjustments
          if(iDeltaEta >= minFilledDeltaEtaBin && iDeltaEta <= maxFilledDeltaEtaBin){
            
            // Scale the actual distribution
            scaledLeadingContent = fBackgroundDistribution->GetBinContent(iDeltaPhi,iDeltaEta)-leadingAdjustment;
            scaledSubleadingContent = fBackgroundDistribution->GetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta)-subleadingAdjustment;
            fBackgroundDistribution->SetBinContent(iDeltaPhi,iDeltaEta,scaledLeadingContent);
            fBackgroundDistribution->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,scaledSubleadingContent);
            
            // Scale also the overlap histogram for debugging purposes
            // Note that for overlap leading and subleading sides change with respect to pi/2
            scaledLeadingContent = fBackgroundOverlap->GetBinContent(iDeltaPhi,iDeltaEta)-subleadingAdjustment;
            scaledSubleadingContent = fBackgroundOverlap->GetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta)-leadingAdjustment;
            fBackgroundOverlap->SetBinContent(iDeltaPhi,iDeltaEta,scaledLeadingContent);
            fBackgroundOverlap->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,scaledSubleadingContent);
            
          } // Bins that are actually filled
          
        } // DeltaEta loop
      } // DeltaPhi loop
    } // Reasonable error if
    
  } // Adjusting the background if*/
  
  // Subtract the generated background from the histogram including the background
  sprintf(histogramName,"%sBackgroundSubtracted",leadingHistogramWithBackground->GetName());
  TH2D *backgroundSubtractedHistogram = (TH2D*)leadingHistogramWithBackground->Clone(histogramName);
  backgroundSubtractedHistogram->Add(fBackgroundDistribution,-1);
  
  // Define a scaling factor for the errors of the background distribution, if it is projected for the full deltaEta.
  // The scaling factor for the errors is needed, because when taking a projection, root will treat the error as
  // error = sum(Errors) / sqrt(nBins). This is the correct treatment for the error if the errors on the bins do not
  // depend on the errors of other bins. In this function the background is projected from the region
  // fMinBackgroundDeltaEta < |deltaEta| < fMaxBackgroundDeltaEta. Thus the errors should be scaled with the square
  // root of the number of bins in this area, not in the whole distribution. Thus without scaling the errors will
  // be underestimated. The scaling factor is the square root of the times there are more bins in the distribution
  // than what is used for the background estimation.
  fBackgroundErrorScalingFactor = TMath::Sqrt(maxDeltaEta/(fMaxBackgroundDeltaEta-fMinBackgroundDeltaEta));
  
  // Set all the negative bins to zero.
  // Note: This seems to be cause bias due to fluctuations arising from small bin size. Thus the part is commented out.
  /*for(int iDeltaPhi = 1; iDeltaPhi <= backgroundSubtractedHistogram->GetNbinsX(); iDeltaPhi++){
    for(int iDeltaEta = 1; iDeltaEta <= backgroundSubtractedHistogram->GetNbinsY(); iDeltaEta++){
      if(backgroundSubtractedHistogram->GetBinContent(iDeltaPhi,iDeltaEta) < 0){
        backgroundSubtractedHistogram->SetBinContent(iDeltaPhi,iDeltaEta,0);
      }
    } // DeltaEta loop
  } // DeltaPhi loop*/
  
  // Return the background subtracted histogram with negative values removed
  return backgroundSubtractedHistogram;
}

/*
 * Combine deltaPhi distribution by taking near side from the two input histograms
 *
 *  Arguments:
 *   const TH2D *leadingHistogramWithBackground = Histogram from which the near side of deltaPhi is taken
 *   const TH2D *subleadingHistogramWithBackground = Histogram from which the away side of the deltaPhi is taken
 *   const double minDeltaEta = Minimum deltaEta from wihch deltaPhi is projected
 *   const double maxDeltaEta = Maximum deltaEta from which deltaPhi is projected
 *   const char* newName = Name to be given to the new histogram
 *   const bool oneSide = Do not do symmetric deltaEta projection
 */
TH1D* DijetMethods::CombineDeltaPhi(const TH2D *leadingHistogramWithBackground, const TH2D *subleadingHistogramWithBackground, const double minDeltaEta, const double maxDeltaEta, const char* newName, const bool oneSide){
  
  // Project the given region from the two-dimensional histograms
  TH1D *nearSideProjection = ProjectRegionDeltaPhi(leadingHistogramWithBackground, minDeltaEta, maxDeltaEta, newName, oneSide);
  TH1D *awaySideProjection = ProjectRegionDeltaPhi(subleadingHistogramWithBackground, minDeltaEta, maxDeltaEta, Form("%s2",newName), oneSide);
  
  // Set the near side of the away side projection as the away side of the near side projection
  int nDeltaPhiBins = nearSideProjection->GetNbinsX();
  double deltaPhiValueAway, deltaPhiErrorAway;
  for(int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins/2; iDeltaPhi++){
    
    deltaPhiValueAway = awaySideProjection->GetBinContent(iDeltaPhi);
    deltaPhiErrorAway = awaySideProjection->GetBinError(iDeltaPhi);
    
    nearSideProjection->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,deltaPhiValueAway);
    nearSideProjection->SetBinError(iDeltaPhi+nDeltaPhiBins/2,deltaPhiErrorAway);
    
  } // Loop over half of the deltaPhi bins
  
  // Return the combined deltaPhi distribution
  return nearSideProjection;
  
}

/*
 * Get the adjustment factors such that if leading and subleading backgrounds are multiplied by them, they match nicely in the middle
 *
 *  Arguments:
 *   TH1 *leadingBackground = Background histogram for the leading side
 *   TH1 *leadingBackgroundOverlap = Overlapping bins from leading and subleading side
 *
 *  return: std::tuple<double,double> Where the first number is factor for leading side and the second for subleading side
 */
std::tuple<double,double> DijetMethods::GetBackgroundAdjustmentFactors(TH1 *leadingBackground, TH1 *leadingBackgroundOverlap) const{
  
  // Set the yield in overlapping bins to zero
  double leadingOverlap = 0;
  double leadingOverlapError = 0;
  double subleadingOverlap = 0;
  double subleadingOverlapError = 0;
  
  // Get the number of deltaPhi bins in the background histogram
  int nDeltaPhiBins = leadingBackground->GetNbinsX();
  
  // Calculate the yield in the overlapping bins
  for(int iDeltaPhi = nDeltaPhiBins/2 - fnOverlapBins + 1; iDeltaPhi <= nDeltaPhiBins/2; iDeltaPhi++){
    leadingOverlap += leadingBackground->GetBinContent(iDeltaPhi);
    leadingOverlapError += leadingBackground->GetBinError(iDeltaPhi);
    subleadingOverlap += leadingBackgroundOverlap->GetBinContent(iDeltaPhi);
    subleadingOverlapError += leadingBackgroundOverlap->GetBinError(iDeltaPhi);
  }
  
  for(int iDeltaPhi = nDeltaPhiBins/2 + 1; iDeltaPhi <= nDeltaPhiBins/2 + fnOverlapBins; iDeltaPhi++){
    subleadingOverlap += leadingBackground->GetBinContent(iDeltaPhi);
    subleadingOverlapError += leadingBackground->GetBinError(iDeltaPhi);
    leadingOverlap += leadingBackgroundOverlap->GetBinContent(iDeltaPhi);
    leadingOverlapError += leadingBackgroundOverlap->GetBinError(iDeltaPhi);
  }
  
  // Calculate the scaling factor from the yields
  double leadingScalingFactor = 1;
  double subleadingScalingFactor = 1;
  
  double leadingToSubleadingRatio = leadingOverlap/subleadingOverlap;
  double subleadingToLeadingRatio = subleadingOverlap/leadingOverlap;
  double leadingToSubleadingRatioError = TMath::Sqrt(TMath::Power(leadingOverlapError/subleadingOverlap,2)+TMath::Power(leadingOverlap*subleadingOverlapError/TMath::Power(subleadingOverlap,2),2));
  
  if(leadingToSubleadingRatioError < 1){
    
    // Scale both sides to match in the middle
    if(leadingToSubleadingRatio > 1){
      leadingScalingFactor = 1 + ((leadingToSubleadingRatio-1)/2);
      subleadingScalingFactor = 1 - ((1-subleadingToLeadingRatio)/2);
    } else {
      leadingScalingFactor = 1 - ((1-leadingToSubleadingRatio)/2);
      subleadingScalingFactor = 1 + ((subleadingToLeadingRatio-1)/2);
    }
    
  }
  
  return std::make_tuple(leadingScalingFactor,subleadingScalingFactor);
}

/*
 * Symmetrize histogram up to the given limits. Set everything to zero above that.
 *
 * Arguments:
 *  const TH2D *histogramToBeSymmetrized = Histogram that is symmetrized
 *  const double maxR = Maximum radius until which the histogram is symmetrized around zero
 *
 *  return: Symmetrized two-dimensional histogram
 */
TH2D* DijetMethods::SymmetrizeHistogram(const TH2D *histogramToBeSymmetrized, const double maxR){

  // Clone the given histogram to get the original binning
  TH2D *symmetrizedHistogram = (TH2D*) histogramToBeSymmetrized->Clone(Form("%sSymmetrized",histogramToBeSymmetrized->GetName()));
  
  // Find the important bin indices to properly symmetrize the distribution
  int deltaPhiZeroPlus = symmetrizedHistogram->GetXaxis()->FindBin(0.001);
  int deltaPhiZeroMinus = symmetrizedHistogram->GetXaxis()->FindBin(-0.001);
  int deltaEtaZeroPlus = symmetrizedHistogram->GetYaxis()->FindBin(0.001);
  int deltaEtaZeroMinus = symmetrizedHistogram->GetYaxis()->FindBin(-0.001);
  
  // Symmetrize the distribution with respect to zero
  double averageContent, averageError;
  int phiBinNegative, phiBinPositive, etaBinNegative, etaBinPositive;
  double binCenterEta, binCenterPhi;
  
  // Loop over deltaPhi bins
  for(int iPhi = 1; iPhi <= symmetrizedHistogram->GetNbinsX(); iPhi++){
    
    // Find symmetric indices around zero for deltaPhi
    if(iPhi < deltaPhiZeroPlus){
      phiBinNegative = iPhi;
      phiBinPositive = deltaPhiZeroPlus + (deltaPhiZeroMinus - iPhi);
    } else {
      phiBinPositive = iPhi;
      phiBinNegative = deltaPhiZeroMinus - (iPhi - deltaPhiZeroPlus);
    }
    binCenterPhi = symmetrizedHistogram->GetXaxis()->GetBinCenter(iPhi);
    
    for(int iEta = 1; iEta <= symmetrizedHistogram->GetNbinsY(); iEta++){
      
      binCenterEta = symmetrizedHistogram->GetYaxis()->GetBinCenter(iEta);
      
      // If the bins are out of specified range, set the content to zero
      if(TMath::Sqrt(binCenterEta*binCenterEta+binCenterPhi*binCenterPhi) > maxR){
        symmetrizedHistogram->SetBinContent(iPhi,iEta,0);
        symmetrizedHistogram->SetBinError(iPhi,iEta,0);
        continue;
      }
      
      // Find symmetric indices around zero for deltaEta
      if(iEta < deltaEtaZeroPlus){
        etaBinNegative = iEta;
        etaBinPositive = deltaEtaZeroPlus + (deltaEtaZeroMinus - iEta);
      } else {
        etaBinPositive = iEta;
        etaBinNegative = deltaEtaZeroMinus - (iEta - deltaEtaZeroPlus);
      }
      
      // Now that we have the symmetric indices around zero, we can calculate the new content for symmetrized histogram
      averageContent = histogramToBeSymmetrized->GetBinContent(phiBinNegative,etaBinNegative);
      averageContent += histogramToBeSymmetrized->GetBinContent(phiBinNegative,etaBinPositive);
      averageContent += histogramToBeSymmetrized->GetBinContent(phiBinPositive,etaBinPositive);
      averageContent += histogramToBeSymmetrized->GetBinContent(phiBinPositive,etaBinNegative);
      averageContent /= 4.0;
      
      averageError = TMath::Power(histogramToBeSymmetrized->GetBinError(phiBinNegative,etaBinNegative),2);
      averageError += TMath::Power(histogramToBeSymmetrized->GetBinError(phiBinNegative,etaBinPositive),2);
      averageError += TMath::Power(histogramToBeSymmetrized->GetBinError(phiBinPositive,etaBinPositive),2);
      averageError += TMath::Power(histogramToBeSymmetrized->GetBinError(phiBinPositive,etaBinNegative),2);
      averageError = TMath::Sqrt(averageError)/4.0;
      
      // After we have calculated the value and error, we can set these values to the symmetrized histogram
      symmetrizedHistogram->SetBinContent(iPhi,iEta,averageContent);
      symmetrizedHistogram->SetBinError(iPhi,iEta,averageError);
      
    } // deltaEta loop
  } // deltaPhi loop
  
  // Return the symmetrized histogram
  return symmetrizedHistogram;
  
}

/*
 * Method for finding the spillover yield from the mixed event corrected distribution. This method fits a constant line
 * to the tails of the DeltaEta projection of the distribution distribution and uses bin counting to get the histogram yield
 * on top of the fitted line. The idea is that in later stage the yield of the fitted gaussians can be fixed to this number.
 * If maxEtaNormalizationRange is set to a smaller number than minEtaNormalizationRange, no line is fitted and the yield is
 * obtained directly from bin counting.
 *
 * Arguments:
 *  TH2D *onlyHydjetHistogram = Mixed event corrected histogram from the HYDJET simulation
 *  double minEtaNormalizationRange = Minimum eta range to which the fitted constant is normalized
 *  double maxEtaNormalizationRange = Maximum eta range to which the fitted constant is normalized
 *
 *  return: Integral of the histogram over the fitted line
 */
double DijetMethods::GetSpilloverYield(TH2D *onlyHydjetHistogram, double minEtaNormalizationRange, double maxEtaNormalizationRange){
  
  // Define the deltaPhi range over which the deltaEta histogram is projected
  double spilloverPhiRange = 1.5;
  
  // Project the deltaEta histogram from the two-dimensional histogram
  fSpilloverDeltaEta = ProjectRegionDeltaEta(onlyHydjetHistogram,-spilloverPhiRange,spilloverPhiRange,"spilloverNormalizationDeltaEta");
  fSpilloverDeltaEta->Scale(fnBinsProjectedOver);
  
  // If maximum fit range is larger than minimum fit range, fit a constant to the given range and determine the background level
  // Otherwise the background level is set to zero.
  double negativeLevel = 0;
  double negativeLevelError = 0;
  double positiveLevel = 0;
  double positiveLevelError = 0;
  if(maxEtaNormalizationRange > minEtaNormalizationRange){
    fSpilloverDeltaEta->Fit("pol0","0","",-maxEtaNormalizationRange,-minEtaNormalizationRange);
    negativeLevel = fSpilloverDeltaEta->GetFunction("pol0")->GetParameter(0);
    negativeLevelError = fSpilloverDeltaEta->GetFunction("pol0")->GetParError(0);
    fSpilloverDeltaEta->Fit("pol0","0","",minEtaNormalizationRange,maxEtaNormalizationRange);
    positiveLevel = fSpilloverDeltaEta->GetFunction("pol0")->GetParameter(0);
    positiveLevelError = fSpilloverDeltaEta->GetFunction("pol0")->GetParError(0);
  }
  double averageLevel = (negativeLevel+positiveLevel)/2;
  double averageLevelError = TMath::Sqrt(negativeLevelError*negativeLevelError+positiveLevelError*positiveLevelError)/2;
  
  // Calculate the bin count from the histogram over the background level
  int minimumYieldBin = fSpilloverDeltaEta->FindBin(-spilloverPhiRange+0.001);
  int maximumYieldBin = fSpilloverDeltaEta->FindBin(spilloverPhiRange-0.001);
  double yieldError = 0;
  double yield = fSpilloverDeltaEta->IntegralAndError(minimumYieldBin,maximumYieldBin,yieldError,"width") - 2*spilloverPhiRange*averageLevel;
  fSpilloverYieldError = TMath::Sqrt(yieldError*yieldError+4*spilloverPhiRange*spilloverPhiRange*averageLevelError*averageLevelError);
  
  // Return the obtained yield
  return yield;
}

/*
 * Get the spillover correction from only HYDJET histogram
 *
 * In heavy ion collisions fluctuations of the background can sometimes be erronously recinstructed as jets.
 * To suppress these contributions, a spillover correction is applied. This effect is estimated by reconstructing
 * jets from only the soft background generetad by the HYDJET. To suppress fluctuations in HYDJET, the deltaEta
 * and deltaPhi projections are fitted with Gaussian functions and a two dimensional Gaussian function is
 * constructed from these. Then this two dimensional function is transformed into a histogram which can be
 * subtracted from the actual distribution to correct for the spillover effect.
 *
 * Arguments:
 *  TH2D *onlyHydjetHistogram = Two-dimensional deltaEta-deltaPhi histogram from soft HYDJET event
 *  int fitMethod = Method used for fitting deltaEta and deltaPhi distributions. 0 = Gauss, 1 = Gauss+constant, 2 = Gauss+Gauss
 *  double spilloverEtaFitRange = Range in deltaEta that is used to fit the spillover deltaEta distribution
 *  double spilloverPhiFitRange = Range in deltaPhi that is used to fit the spillover deltaPhi distribution
 *  double lowConstantRange = Low range to which the constant is fitted in constant + Gauss fit
 *  double highConstantRange = High range to which the constant is fitted in constant + Gauss fit
 *  double fixedYield = If number > 0, fix the yield of the Gaussian fit to this number
 *  double fixedEtaWidth = If number > 0, fix the width of the Gaussian in deltaEta fit to this number
 *  double fixedPhiWidth = If number > 0, fix the width of the Gaussian in deltaPhi fit to this number
 *
 *  return: Two-dimensional histogram with the estimation of the spillover effect
 */
TH2D* DijetMethods::GetSpilloverCorrection(TH2D *onlyHydjetHistogram, int fitMethod, double spilloverEtaFitRange, double spilloverPhiFitRange, double lowConstantRange, double highConstantRange, double fixedYield, double fixedEtaWidth, double fixedPhiWidth){
  
  // Define the range over which the spillover correction is calculated
  double spilloverEtaRange = 1.5;
  double spilloverPhiRange = 1.5;
  
  // First get the projections for the deltaEta and deltaPhi distributions
  fSpilloverDeltaEta = ProjectRegionDeltaEta(onlyHydjetHistogram,-spilloverPhiRange,spilloverPhiRange,"spilloverDeltaEta");
  fSpilloverDeltaEta->Scale(fnBinsProjectedOver);
  fSpilloverDeltaPhi = ProjectRegionDeltaPhi(onlyHydjetHistogram,0,spilloverEtaRange,"spilloverDeltaPhi");
  fSpilloverDeltaPhi->Scale(fnBinsProjectedOver);
  
  // Do some rebin
  int etaRebin = 1; // Old: 4
  int phiRebin = 1; // Old: 2
  if(etaRebin > 1){
    fSpilloverDeltaEta->Rebin(etaRebin);
    fSpilloverDeltaEta->Scale(1.0/etaRebin);
  }
  if(phiRebin > 1){
    fSpilloverDeltaPhi->Rebin(phiRebin);
    fSpilloverDeltaPhi->Scale(1.0/phiRebin);
  }
  
  // Fit one dimensional Gaussians to the projected deltaEta and deltaPhi distributions
  if(fitMethod == 0){
    fSpilloverFitDeltaEta = FitGauss(fSpilloverDeltaEta, spilloverEtaFitRange, spilloverPhiRange, fixedYield, fixedEtaWidth);
    fSpilloverFitDeltaPhi = FitGauss(fSpilloverDeltaPhi, spilloverPhiFitRange, spilloverEtaRange, fixedYield, fixedPhiWidth);
  } else if(fitMethod == 1){
    fSpilloverFitDeltaEta = FitGaussAndConstant(fSpilloverDeltaEta, spilloverEtaFitRange, spilloverPhiRange, lowConstantRange, highConstantRange, fixedYield, fixedEtaWidth);
    fSpilloverFitDeltaPhi = FitGaussAndConstant(fSpilloverDeltaPhi, spilloverPhiFitRange, spilloverEtaRange, 1, 2, fixedYield, fixedPhiWidth);
  }else {
    fSpilloverFitDeltaEta = FitDoubleGauss(fSpilloverDeltaEta, spilloverEtaFitRange, spilloverPhiRange);
    fSpilloverFitDeltaPhi = FitDoubleGauss(fSpilloverDeltaPhi, spilloverPhiFitRange, spilloverEtaRange);
  }
  
  // Combine the one-dimensional fits to a two-dimensional Gaussian distribution
  TF2 *gauss2D;
  
  if(fitMethod < 2){
    gauss2D = new TF2("gauss2D", "[0]/(2*TMath::Pi()*[1]*[2])*TMath::Exp(-0.5*TMath::Power(x/[1],2))*TMath::Exp(-0.5*TMath::Power(y/[2],2))",-TMath::Pi()/2,3*TMath::Pi()/2,-5,5);
    gauss2D->SetParameter(0,(fSpilloverFitDeltaEta->GetParameter(0)+fSpilloverFitDeltaEta->GetParameter(0))/2.0);
    gauss2D->SetParameter(1,fSpilloverFitDeltaPhi->GetParameter(1));
    gauss2D->SetParameter(2,fSpilloverFitDeltaEta->GetParameter(1));
  } else {
    gauss2D = new TF2("gauss2D",doubleGaussNoOverlap2D,-TMath::Pi()/2,3*TMath::Pi()/2,-5,5,6);
    gauss2D->SetParameter(0,(fSpilloverFitDeltaEta->GetParameter(0)+fSpilloverFitDeltaEta->GetParameter(0))/2.0);
    gauss2D->SetParameter(1,fSpilloverFitDeltaPhi->GetParameter(1));
    gauss2D->SetParameter(2,fSpilloverFitDeltaEta->GetParameter(1));
    gauss2D->SetParameter(3,(fSpilloverFitDeltaEta->GetParameter(2)+fSpilloverFitDeltaEta->GetParameter(2))/2.0);
    gauss2D->SetParameter(4,fSpilloverFitDeltaPhi->GetParameter(3));
    gauss2D->SetParameter(5,fSpilloverFitDeltaEta->GetParameter(3));
  }
  
  // Create a new two-dimensional histogram from the two-dimensional Gaussian function
  char histogramName[150];
  sprintf(histogramName,"%sSpillover",onlyHydjetHistogram->GetName());
  TH2D *spilloverCorrection = (TH2D*) onlyHydjetHistogram->Clone(histogramName);
  
  // Set all the bins to zero before filling the histogram with function contents
  for(int iPhiBin = 0; iPhiBin < spilloverCorrection->GetNbinsX(); iPhiBin++){
    for(int iEtaBin = 0; iEtaBin < spilloverCorrection->GetNbinsY(); iEtaBin++){
      spilloverCorrection->SetBinContent(iPhiBin,iEtaBin,0);
      spilloverCorrection->SetBinError(iPhiBin,iEtaBin,0);
    }
  }
  
  spilloverCorrection->Eval(gauss2D,"A");
  
  // Set the errors to zero after filling the histogram from the two-dimensional function
  // The systematic error estimation should be done separately, errors would be overestimated otherwise
  // By default, root just assigns the square root of the bin content as an error for each bin
  for(int iPhiBin = 1; iPhiBin <= spilloverCorrection->GetNbinsX(); iPhiBin++){
    for(int iEtaBin = 1; iEtaBin <= spilloverCorrection->GetNbinsY(); iEtaBin++){
      spilloverCorrection->SetBinError(iPhiBin,iEtaBin,0);
    }
  }
    
  // Return the spillover correction
  return spilloverCorrection;
}

/*
 * Fit a Gaussian function to a histogram and return the fit function
 *
 * Arguments:
 *  TH1D *fittedHistogram = Histogram to be fitted with the Gaussian
 *  double fitRange = Range for the fit
 *  double normalizationRange = Range for yield normalization
 *  double fixedYield = If > 0, fix the Gaussian fit yield to this number
 *  double fixedWidth = If > 0, fix the width of the Gaussian fit to this number
 *
 *  return: Gaussian function fitted to the histogram
 */
TF1* DijetMethods::FitGauss(TH1D *fittedHistogram, double fitRange, double normalizationRange, double fixedYield, double fixedWidth){
  
  // Create a Gaussian function for the fit
  TF1 *gaussFit = new TF1("gaussFit","[0]/(TMath::Sqrt(2*TMath::Pi())*[1])*TMath::Exp(-0.5*TMath::Power(x/[1],2))+[2]",-fitRange,fitRange);
  
  // Calculate the integral of the fitted histogram over the normalization range
  double fitYield = fittedHistogram->Integral(fittedHistogram->FindBin(-normalizationRange+0.001),fittedHistogram->FindBin(normalizationRange-0.001),"width");
  
  // Fix the yield of the Gaussian function to the integral and give some rought estimate for the width
  gaussFit->SetParameter(0,fitYield);
  if(fixedYield > 0) gaussFit->FixParameter(0,fixedYield);
  gaussFit->SetParameter(1,0.3);
  if(fixedWidth > 0) gaussFit->FixParameter(1,fixedWidth);
  gaussFit->FixParameter(2,0); // Fix the constant to zero. Add zero here to make this compatible with Gauss+constant fit elwhere in the code.
  
  // Fit the histogram over the defined range
  fittedHistogram->Fit("gaussFit","","",-fitRange,fitRange);
  
  // Return the fit function
  return gaussFit;
  
}

/*
 * Fit a Gaussian function together with a constant to a histogram and return the fit function
 *
 * Arguments:
 *  TH1D *fittedHistogram = Histogram to be fitted with the Gaussian
 *  double fitRange = Range for the fit
 *  double normalizationRange = Range for yield normalization
 *  double lowConstantRange = Low range to which the constant is fitted
 *  double highConstantRange = High range to which the constant is fitted
 *  double fixedYield = If > 0, fix the Gaussian fit yield to this number
 *  double fixedWidth = If > 0, fix the width of the Gaussian fit to this number
 *
 *  return: Gaussian function fitted to the histogram
 */
TF1* DijetMethods::FitGaussAndConstant(TH1D *fittedHistogram, double fitRange, double normalizationRange, double lowConstantRange, double highConstantRange, double fixedYield, double fixedWidth){
  
  // Create a Gaussian function for the fit
  TF1 *gaussFit = new TF1("gaussFit","[0]/(TMath::Sqrt(2*TMath::Pi())*[1])*TMath::Exp(-0.5*TMath::Power(x/[1],2))+[2]",-fitRange,fitRange);
  
  // Calculate the integral of the fitted histogram over the normalization range
  double fitYield = fittedHistogram->Integral(fittedHistogram->FindBin(-normalizationRange+0.001),fittedHistogram->FindBin(normalizationRange-0.001),"width");
  
  // First fit the constant to the specified range
  fittedHistogram->Fit("pol0","","",-highConstantRange,-lowConstantRange);
  double negativeLevel = fittedHistogram->GetFunction("pol0")->GetParameter(0);
  fittedHistogram->Fit("pol0","","",lowConstantRange,highConstantRange);
  double positiveLevel = fittedHistogram->GetFunction("pol0")->GetParameter(0);
  double averageLevel = (negativeLevel+positiveLevel)/2;
  
  // Give initial estimates for the parameters
  gaussFit->SetParameter(0,fitYield-(averageLevel*3));
  if(fixedYield > 0) gaussFit->FixParameter(0,fixedYield);
  gaussFit->SetParameter(1,0.3);
  if(fixedWidth > 0) gaussFit->FixParameter(1,fixedWidth);
  gaussFit->FixParameter(2,averageLevel);
  
  // Fit the histogram over the defined range
  fittedHistogram->Fit("gaussFit","","",-fitRange,fitRange);
  
  // Return the fit function
  return gaussFit;
  
}

/*
 * Fit narrow and wide Gaussians function to a histogram and return the fit function
 *
 * Arguments:
 *  TH1D *fittedHistogram = Histogram to be fitted with two Gaussians
 *  double fitRange = Range for the fit
 *  double normalizationRange = Range for yield normalization
 *
 *  return: Double Gaussian function fitted to the histogram
 */
TF1* DijetMethods::FitDoubleGauss(TH1D *fittedHistogram, double fitRange, double normalizationRange){
  
  // Create a Gaussian function for the fit
  //TF1 *gaussFit = new TF1("gaussFit","[0]/(TMath::Sqrt(2*TMath::Pi())*[1])*TMath::Exp(-0.5*TMath::Power(x/[1],2))+[2]/(TMath::Sqrt(2*TMath::Pi())*[3])*TMath::Exp(-0.5*TMath::Power(x/[3],2))",-fitRange,fitRange);
  TF1 *gaussFit = new TF1("gaussFit",doubleGaussNoOverlap,-fitRange,fitRange,4);
  
  // Calculate the integral of the fitted histogram over the normalization range
  double fitYield = fittedHistogram->Integral(fittedHistogram->FindBin(-normalizationRange+0.001),fittedHistogram->FindBin(normalizationRange-0.001),"width");
  
  // Initialize the yield of the Gaussian function to the integral and give some rought estimate for the width
  gaussFit->SetParameter(0,fitYield*0.8);
  gaussFit->SetParameter(1,0.3);
  gaussFit->SetParameter(2,fitYield*0.2);
  gaussFit->SetParameter(3,0.8);
  gaussFit->SetParLimits(0,0,fitYield);
  gaussFit->SetParLimits(1,0.05,0.38);
  gaussFit->SetParLimits(2,0,fitYield);
  gaussFit->SetParLimits(3,0.38,10);
  
  // Fit the histogram over the defined range
  fittedHistogram->Fit("gaussFit","","",-fitRange,fitRange);
  
  // Return the fit function
  return gaussFit;
  
}

/*
 * Project the deltaEta distribution out of a two-dimensional deltaPhi-deltaEta distribution
 *
 *  Arguments:
 *   TH2D* deltaPhiDeltaEtaHistogram = two-dimensional deltaPhi-deltaEta distribution
 *   const double minDeltaPhi = Minimum deltaPhi value for the projected region
 *   const double maxDeltaPhi = Maximum deltaPhi value for the projected region
 *   const char* newName = Name that is added to the projected histogram
 *
 *   return: DeltaPhi distribution projected from the input deltaEta region
 */
TH1D* DijetMethods::ProjectRegionDeltaEta(const TH2D* deltaPhiDeltaEtaHistogram, const double minDeltaPhi, const double maxDeltaPhi, const char* newName){
  
  // Start by finding the bin indices for the defined deltaEta region
  // Apply a little offset for the defined eta region borders to avoid bin border effects
  int lowDeltaPhiBin  = deltaPhiDeltaEtaHistogram->GetXaxis()->FindBin(minDeltaPhi+0.0001);
  int highDeltaPhiBin = deltaPhiDeltaEtaHistogram->GetXaxis()->FindBin(maxDeltaPhi-0.0001);
  
  // Calculate the number of deltaPhi bins in the projected region
  fnBinsProjectedOver = highDeltaPhiBin-lowDeltaPhiBin+1;
  
  // Project deltaEta distribution from the defined deltaPhi region
  char histogramName[200];
  sprintf(histogramName,"%s%s",deltaPhiDeltaEtaHistogram->GetName(),newName);
  TH1D *projectedDeltaEta;
  projectedDeltaEta = deltaPhiDeltaEtaHistogram->ProjectionY(histogramName,lowDeltaPhiBin,highDeltaPhiBin);
  
  // Scale the projected deltaEta distribution with the number of deltaPhi bins projected over and deltaPhi bin width to retain normalization
  projectedDeltaEta->Scale(1.0/fnBinsProjectedOver);
  projectedDeltaEta->Scale(deltaPhiDeltaEtaHistogram->GetXaxis()->GetBinWidth(1));
  
  // Return the scaled projection
  return projectedDeltaEta;
}

/*
 * Project the deltaPhi distribution out of a two-dimensional deltaPhi-deltaEta distribution
 *
 *  Arguments:
 *   TH2D* deltaPhiDeltaEtaHistogram = two-dimensional deltaPhi-deltaEta distribution
 *   const double minDeltaEta = Minimum deltaEta value for the projected region
 *   const double maxDeltaEta = Maximum deltaEta value for the projected region
 *   const char* newName = Name that is added to the projected histogram
 *   const bool oneSide = Only look at the given side of deltaEta, do not take symmetric region from the opposite side
 *
 *   return: DeltaPhi distribution projected from the input deltaEta region
 */
TH1D* DijetMethods::ProjectRegionDeltaPhi(const TH2D* deltaPhiDeltaEtaHistogram, const double minDeltaEta, const double maxDeltaEta, const char* newName, const bool oneSide){
  
  // If the minimum is zero, do the projection from -maxDeltaEta to maxDeltaEta
  bool oneRegion = (TMath::Abs(minDeltaEta) < 0.01);
  
  // Start by finding the bin indices for the defined deltaEta region
  // Apply a little offset for the defined eta region borders to avoid bin border effects
  int lowNegativeDeltaEtaBin  = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(-maxDeltaEta+0.0001);
  int highNegativeDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(-minDeltaEta-0.0001);
  int lowPositiveDeltaEtaBin  = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(minDeltaEta+0.0001);
  int highPositiveDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(maxDeltaEta-0.0001);
  
  // Calculate the number of deltaEta bins in the projected region
  fnBinsProjectedOver = highNegativeDeltaEtaBin-lowNegativeDeltaEtaBin+highPositiveDeltaEtaBin-lowPositiveDeltaEtaBin+2;
  if(oneRegion) fnBinsProjectedOver = fnBinsProjectedOver-highNegativeDeltaEtaBin+lowPositiveDeltaEtaBin-1;
  if(oneSide) fnBinsProjectedOver = highPositiveDeltaEtaBin-lowPositiveDeltaEtaBin-1;
  
  // Project out the deltaPhi distribution in the defined region
  char histogramName[200];
  sprintf(histogramName,"%s%s",deltaPhiDeltaEtaHistogram->GetName(),newName);
  TH1D *projectedDeltaPhi;
  if(oneRegion) {
    projectedDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName,lowNegativeDeltaEtaBin,highPositiveDeltaEtaBin);
  } else if (oneSide){
    projectedDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName,lowPositiveDeltaEtaBin,highPositiveDeltaEtaBin);
  } else {
    projectedDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName,lowNegativeDeltaEtaBin,highNegativeDeltaEtaBin);
    projectedDeltaPhi->Add(deltaPhiDeltaEtaHistogram->ProjectionX("dummyName",lowPositiveDeltaEtaBin,highPositiveDeltaEtaBin));
  }
  
  // Scale the projected deltaPhi distribution with the number of deltaEta bins projected over and deltaEta bin width to retain normalization
  projectedDeltaPhi->Scale(1.0/fnBinsProjectedOver);
  projectedDeltaPhi->Scale(deltaPhiDeltaEtaHistogram->GetYaxis()->GetBinWidth(1));
  
  // Return the scaled projection
  return projectedDeltaPhi;
}

/*
 * Project the background deltaPhi distribution out of a two-dimensional deltaPhi-deltaEta distribution
 *
 *  Arguments:
 *   TH2D* deltaPhiDeltaEtaHistogram = two-dimensional deltaPhi-deltaEta distribution
 *
 *   return: DeltaPhi distribution projected from the background deltaEta region
 */
TH1D* DijetMethods::ProjectBackgroundDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram){
  return ProjectRegionDeltaPhi(deltaPhiDeltaEtaHistogram,fMinBackgroundDeltaEta,fMaxBackgroundDeltaEta,"BackgroundDeltaPhi");
}

/*
 * Project the signal deltaPhi distribution out of a two-dimensional deltaPhi-deltaEta distribution
 *
 *  Arguments:
 *   TH2D* deltaPhiDeltaEtaHistogram = two-dimensional deltaPhi-deltaEta distribution
 *
 *   return: DeltaPhi distribution projected from the signal deltaEta region
 */
TH1D* DijetMethods::ProjectSignalDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram){
  return ProjectRegionDeltaPhi(deltaPhiDeltaEtaHistogram,0,fMaxSignalDeltaEta,"SignalDeltaPhi");
}

/*
 * Extract the jet shape from two-dimensional histogram
 *
 *  Arguments:
 *   TH2D *backgroundSubtractedHistogram = Background subtracted deltaPhi-deltaEta histogram
 *
 *   return: TH1D histogram with extracted jet shape
 */
TH1D* DijetMethods::GetJetShape(TH2D *backgroundSubtractedHistogram){
  
  // Create a new TH1D for the jet shape
  char histogramName[200];
  sprintf(histogramName,"%sJetShape",backgroundSubtractedHistogram->GetName());
  TH1D *jetShapeHistogram = new TH1D(histogramName,histogramName,fnRBins,fRBins);
  
  // Create another histogram to count how many bins in two-dimensional histogram correspond to filled bins in jet shape histogram
  sprintf(histogramName,"%sJetShapeCounts",backgroundSubtractedHistogram->GetName());
  fhJetShapeCounts = (TH1D*) jetShapeHistogram->Clone(histogramName);
  
  // Create a two dimensional histogram to collect information on which Rbin each deltaEta-deltaPhi phi is assigned
  sprintf(histogramName,"%sRbinMap",backgroundSubtractedHistogram->GetName());
  fhJetShapeBinMap = (TH2D*) backgroundSubtractedHistogram->Clone(histogramName);
  
  // Set all the bins to zero for the RBin to deltaPhi-deltaEta mapping histogram
  for(int iPhiBin = 0; iPhiBin < fhJetShapeBinMap->GetNbinsX(); iPhiBin++){
    for(int iEtaBin = 0; iEtaBin < fhJetShapeBinMap->GetNbinsY(); iEtaBin++){
      fhJetShapeBinMap->SetBinContent(iPhiBin,iEtaBin,0);
      fhJetShapeBinMap->SetBinError(iPhiBin,iEtaBin,0);
    }
  }
  
  // To get the yield, bin area normalization must be multiplied out
  double binWidthDeltaPhi = backgroundSubtractedHistogram->GetXaxis()->GetBinWidth(1);
  double binWidthDeltaEta = backgroundSubtractedHistogram->GetYaxis()->GetBinWidth(1);
  double oneBinArea = binWidthDeltaPhi*binWidthDeltaEta;
  
  // The jet shape is calculated from the region close to near side peak, so we need to find bin indices for that region
  // Apply small offset to avoid problems with bin boundaries
  int minBinPhi = backgroundSubtractedHistogram->GetXaxis()->FindBin((-TMath::Pi()/2.0)+0.0001);
  int maxBinPhi = backgroundSubtractedHistogram->GetXaxis()->FindBin((TMath::Pi()/2.0)-0.0001);
  int minBinEta = backgroundSubtractedHistogram->GetYaxis()->FindBin((-TMath::Pi()/2.0)+0.0001);
  int maxBinEta = backgroundSubtractedHistogram->GetYaxis()->FindBin((TMath::Pi()/2.0)-0.0001);
  
  // Helper variables for calcalations inside the loop
  double deltaEta;   // deltaEta
  double deltaPhi;   // deltaPhi
  double R;          // sqrt(deltaEta^2+deltaPhi^2)
  int Rbin;          // Bin in jet shape histogram corresponding to calculated R
  double newContent; // Content to be added to the jet shape histogram
  double newError;   // Error of the new content to be added
  double oldContent; // Content in the jet shape histogram before the new addition
  double oldError;   // Error of the content already in the jet shape histogram
  double oldCounts;  // Number pf counts in the count histogram
  
  // Loop over the selected bins, calculate R, and fill the bin content to histogram
  for(int iPhiBin = minBinPhi; iPhiBin <= maxBinPhi; iPhiBin++){
    deltaPhi = backgroundSubtractedHistogram->GetXaxis()->GetBinCenter(iPhiBin); // For deltaPhi value, use bin center
    
    for(int iEtaBin = minBinEta; iEtaBin <= maxBinEta; iEtaBin++){
      deltaEta = backgroundSubtractedHistogram->GetYaxis()->GetBinCenter(iEtaBin); // For deltaEta value, use bin center
      
      // Calculate R = sqrt(deltaEta^2+deltaPhi^2)
      R = TMath::Sqrt(TMath::Power(deltaPhi,2)+TMath::Power(deltaEta,2));
      
      // Find the correct bin in the jet shape histogram to fill for this R
      Rbin = jetShapeHistogram->FindBin(R);
      
      // Get the content to be added from the two-dimensional histogram, multiplying out the bin area normalization
      newContent = backgroundSubtractedHistogram->GetBinContent(iPhiBin,iEtaBin)*oneBinArea;
      newError = backgroundSubtractedHistogram->GetBinError(iPhiBin,iEtaBin)*oneBinArea;
      
      // Get the content already in the jet shape histogram
      oldContent = jetShapeHistogram->GetBinContent(Rbin);
      oldError = jetShapeHistogram->GetBinError(Rbin);
      oldCounts = fhJetShapeCounts->GetBinContent(Rbin);
      
      // Add the new content on top of the old content. Add errors in quadrature
      jetShapeHistogram->SetBinContent(Rbin,newContent+oldContent);
      jetShapeHistogram->SetBinError(Rbin,TMath::Sqrt(TMath::Power(newError,2)+TMath::Power(oldError,2)));
      
      // Increment the counts by one. There is no error for counts, since we know the histogram binning exactly
      fhJetShapeCounts->SetBinContent(Rbin,oldCounts+1);
      fhJetShapeCounts->SetBinError(Rbin,0);
      
      // Set the Rbin mapping to deltaEta-deltaPhi bins
      fhJetShapeBinMap->SetBinContent(iPhiBin,iEtaBin,Rbin);
      
    } // deltaEta loop
  } // deltaPhi loop
  
  // Normalize each bin in the jet shape histogram by the number of bins in two-dimensional histogram corresponding to that bin
  if(fJetShapeNormalizationMethod == kBinWidth){
    jetShapeHistogram->Scale(1.0,"width");

    // To correct for the fact the use rectangles to integrate a circular area, we need to weight the bins on how accurately
    // the rectangular area reflects the actual circles.
    double totalBinArea;
    double ringArea;
    
    for(int iRBin = 1; iRBin <= jetShapeHistogram->GetNbinsX(); iRBin++){
      totalBinArea = fhJetShapeCounts->GetBinContent(iRBin)*oneBinArea;
      ringArea = TMath::Pi()*fRBins[iRBin]*fRBins[iRBin] - TMath::Pi()*fRBins[iRBin-1]*fRBins[iRBin-1];
      jetShapeHistogram->SetBinContent(iRBin,jetShapeHistogram->GetBinContent(iRBin)*(ringArea/totalBinArea));
      jetShapeHistogram->SetBinError(iRBin,jetShapeHistogram->GetBinError(iRBin)*(ringArea/totalBinArea));
    }
    
  } else if (fJetShapeNormalizationMethod == kBinArea){
    jetShapeHistogram->Divide(fhJetShapeCounts);
  }
  
  // Return the calculated jet shape histogram
  return jetShapeHistogram;
  
}

/*
 * Do a Fourier fit for the background deltaPhi distribution
 *
 *  TH1D* backgroundDeltaPhi = Background deltaPhi histogram
 *  const int maxVn = Largest vn included in the Fourier fit
 *
 *  return = The Fourier fit function
 */
TF1* DijetMethods::FourierFit(TH1D* backgroundDeltaPhi, const int maxVn){
  
  // Define the fourier fit. Use fit parameters up to maxVn
  char fourierFormula[300] = "[0]*(1";
  for(int vn = 1; vn <= maxVn; vn++){
    sprintf(fourierFormula,"%s+2.0*[%d]*TMath::Cos(%d.0*x)",fourierFormula,vn,vn);
  }
  sprintf(fourierFormula,"%s)",fourierFormula);
  
  TF1 *fourier = new TF1("fourier",fourierFormula,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);

  // Set names for the parameters
  fourier->SetParName(0,"BkgLevel");
  for(int vn = 1; vn <= maxVn; vn++){
    sprintf(fourierFormula,"v%d",vn);
    fourier->SetParName(vn,fourierFormula);
  }
  
  // Set initial values for the parameters
  fourier->SetParameter(0,backgroundDeltaPhi->GetBinContent(backgroundDeltaPhi->FindBin(TMath::Pi()/2)));
  fourier->SetParameter(1,0);
  fourier->SetParameter(2,0.03);
  for(int vn = 3; vn <= maxVn; vn++){
    fourier->SetParameter(vn,0);
  }
  
  // Set limits such that the parameters must remain sensible
  fourier->SetParLimits(1,-1.0,1.0);
  for(int vn = 2; vn <= maxVn; vn++){
    fourier->SetParLimits(vn,0,1.0);
  }
  
  // Do the fit!
  backgroundDeltaPhi->Fit("fourier","","",-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
  
  // Return the fitted distribution
  return fourier;
}

/*
 * Systematic uncertainty estimation for the pair acceptance correction.
 * To estimate the uncertainty, the mixed event corrected deltaEta histogram is fitted with
 * a constant to the regions -2.5 < deltaEta < -1.5 and 1.5 < deltaEta < 2.5 separately.
 * The difference between these two fits is assigned as the systematic uncertainty coming
 * from the pair acceptance correction
 */
double DijetMethods::EstimateSystematicsForPairAcceptanceCorrection(const TH1* deltaEtaHistogram){
  
  TH1D *fitHistogram = (TH1D*) deltaEtaHistogram->Clone("temporaryFitHistogram");
  double pairAcceptanceLowLimit = 1.5;
  double pairAcceptanceHighLimit = 2.5;
  
  // Fit a constant to the specified ranges
  fitHistogram->Fit("pol0","","",-pairAcceptanceHighLimit,-pairAcceptanceLowLimit);
  fPairAcceptanceNegativeLevel = 0;
  if(fitHistogram->GetFunction("pol0")) fPairAcceptanceNegativeLevel = fitHistogram->GetFunction("pol0")->GetParameter(0);
  
  fitHistogram->Fit("pol0","","",pairAcceptanceLowLimit,pairAcceptanceHighLimit);
  fPairAcceptancePositiveLevel = 0;
  if(fitHistogram->GetFunction("pol0")) fPairAcceptancePositiveLevel = fitHistogram->GetFunction("pol0")->GetParameter(0);
  
  // Return the absolute value of the diffence between the two levels.
  return TMath::Abs(fPairAcceptancePositiveLevel-fPairAcceptanceNegativeLevel);
  
}

/*
 * Systematic uncertainty estimation for the background subtraction.
 * To estimate the uncertainty, the background subtracted deltaEta histogram is fitted in
 * four different regions, -2.5 < deltaEta < -2.0, -2.0 < deltaEta < -1.5, 1.5 < deltaEta < 2.0
 * and 2.0 < deltaEta < 2.5. We then calculate the average for the inner part 1.5 < |deltaEta| < 2.0
 * and the outer part 2.0 < |deltaEta| < 2.5 and assign larger for these deviation from zero as the
 * systematic uncertainty coming from the background subtraction.
 */
double DijetMethods::EstimateSystematicsForBackgroundSubtraction(const TH1* deltaEtaHistogram){
  
  TH1D *fitHistogram = (TH1D*) deltaEtaHistogram->Clone("temporaryFitHistogram");
  double backgroundLowLimit = 1.5;
  double backgroundMidPoint = 2.0;
  double backgroundHighLimit = 2.5;
  
  // Fit a constant to the specified ranges
  fitHistogram->Fit("pol0","","",-backgroundMidPoint,-backgroundLowLimit);
  double innerNegativeLevel = 0;
  if(fitHistogram->GetFunction("pol0")) innerNegativeLevel = fitHistogram->GetFunction("pol0")->GetParameter(0);
  
  fitHistogram->Fit("pol0","","",-backgroundHighLimit,-backgroundMidPoint);
  double outerNegativeLevel = 0;
  if(fitHistogram->GetFunction("pol0")) outerNegativeLevel = fitHistogram->GetFunction("pol0")->GetParameter(0);
  
  fitHistogram->Fit("pol0","","",backgroundLowLimit,backgroundMidPoint);
  double innerPositiveLevel = 0;
  if(fitHistogram->GetFunction("pol0")) innerPositiveLevel = fitHistogram->GetFunction("pol0")->GetParameter(0);
  
  fitHistogram->Fit("pol0","","",backgroundMidPoint,backgroundHighLimit);
  double outerPositiveLevel = 0;
  if(fitHistogram->GetFunction("pol0")) outerPositiveLevel = fitHistogram->GetFunction("pol0")->GetParameter(0);
  
  // Calculate means for outer and inner levels
  fBackgroundSubtractionInnerMean = TMath::Abs((innerNegativeLevel+innerPositiveLevel)/2.0);
  fBackgroundSubtractionOuterMean = TMath::Abs((outerNegativeLevel+outerPositiveLevel)/2.0);
  
  // Return the larger deviation from zero
  if(fBackgroundSubtractionOuterMean > fBackgroundSubtractionInnerMean) return fBackgroundSubtractionOuterMean;
  return fBackgroundSubtractionInnerMean;
  
}

/*
 * Some systematic uncertainty estimations are calculated from deltaEta histograms.
 * We cannot just straighforwardly propagate those to deltaR, but we need to take
 * into account the ring structure. The ring is taken into account by multiplying by
 * the area of the current ring and then normalizing by the deltaR bin width.
 */
void DijetMethods::PropagateDeltaEtaToDeltaR(TH1* errorHistogram, const double deltaEtaError){
  
  double ring;
  for(int iDeltaR = 1; iDeltaR <= errorHistogram->GetNbinsX(); iDeltaR++){
    ring = TMath::Pi()*(pow(errorHistogram->GetBinLowEdge(iDeltaR)+errorHistogram->GetBinWidth(iDeltaR),2) - pow(errorHistogram->GetBinLowEdge(iDeltaR),2))/errorHistogram->GetBinWidth(iDeltaR);
    errorHistogram->SetBinContent(iDeltaR,ring*deltaEtaError);
  }
}

/*
 * Rebin one dimensional histogram with asymmetric bin edges
 *
 * Arguments:
 *  TH1D* histogramInNeedOfRebinning = Histogram to be rebinned
 *  const int nBins = Number of bins for the rebinned histogram
 *  const double* binEdges = Bin edges for the rebinned histogram
 */
TH1D* DijetMethods::RebinAsymmetric(TH1D* histogramInNeedOfRebinning, const int nBins, const double* binEdges){
  
  // First, check that the new bin boundaries are also bin boundaries in the original histogram
  bool binsGood = CheckBinBoundaries(nBins,binEdges,histogramInNeedOfRebinning->GetXaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because given bin borders do not match with the bin borders of the original histogram!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  // Clone the original histogram
  TH1D *clonedHistogram = (TH1D*) histogramInNeedOfRebinning->Clone("_rebinned");
  
  // Set the new binning for histogram (destroys content in each bin)
  clonedHistogram->SetBins(nBins,binEdges);
  
  // Make sure that each bin is set to zero
  for(int iBin = 1; iBin <= clonedHistogram->GetNbinsX(); iBin++){
    clonedHistogram->SetBinContent(iBin,0);
    clonedHistogram->SetBinError(iBin,0);
  }
  
  // Add the contents back to the histogram that was rebinned
  double binContent, binError, binCenter, oldContent, oldError;
  int newBin;
  double binWidth;
  for(int iBin = 1; iBin <= histogramInNeedOfRebinning->GetNbinsX(); iBin++){
    
    // Read the contents from the non-rebinned histogram
    binWidth = histogramInNeedOfRebinning->GetBinWidth(iBin);
    binContent = histogramInNeedOfRebinning->GetBinContent(iBin)*binWidth;  // Remove previous bin width normalization
    binError = histogramInNeedOfRebinning->GetBinError(iBin)*binWidth;      // Remove previous bin width normalization
    binCenter = histogramInNeedOfRebinning->GetBinCenter(iBin);
    
    // Add the contents to the rebinned histgram
    newBin = clonedHistogram->FindBin(binCenter);
    oldContent = clonedHistogram->GetBinContent(newBin);
    oldError = clonedHistogram->GetBinError(newBin);
    clonedHistogram->SetBinContent(newBin,binContent+oldContent);
    clonedHistogram->SetBinError(newBin,TMath::Sqrt(binError*binError+oldError*oldError));
  }
  
  // Normalize the bin contents to bin width
  for(int iBin = 1; iBin <= clonedHistogram->GetNbinsX(); iBin++){
    binWidth = clonedHistogram->GetBinWidth(iBin);
    binContent = clonedHistogram->GetBinContent(iBin);
    binError = clonedHistogram->GetBinError(iBin);
    clonedHistogram->SetBinContent(iBin,binContent/binWidth);
    clonedHistogram->SetBinError(iBin,binError/binWidth);
  }
  
  // Return the rebinned histogram
  return clonedHistogram;
}


/*
 * Rebin a two-dimensional deltaPhi-deltaEta histogram
 *
 *  Arguments:
 *   TH2D *histogramInNeedOfRebinning = DeltaPhi-deltaEta histogram to be rebinned
 *
 *   return: Rebinned histogram
 */
TH2D* DijetMethods::RebinHistogram(TH2D *histogramInNeedOfRebinning){
  
  // First, check that the new bin boundaries are also bin boundaries in the original histogram
  bool binsGood = CheckBinBoundaries(fnRebinDeltaPhi,fRebinDeltaPhi,histogramInNeedOfRebinning->GetXaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because of a bin edge problem in deltaPhi axis!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  binsGood = CheckBinBoundaries(fnRebinDeltaEta,fRebinDeltaEta,histogramInNeedOfRebinning->GetYaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because of a bin edge problem in deltaEta axis!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  // Root does not offer a method to directly rebin a two-dimensional histogram, so I have implemented my own
  // Helper variables for rebinning
  double currentBinContent;
  double currentBinError;
  double nonRebinnedContent;
  double nonRebinnedError;
  int newHistogramIndex;
  double deltaPhiValue;
  double deltaEtaValue;
  
  // Create the histogram with new binning
  char newName[200];
  sprintf(newName,"%sRebinned",histogramInNeedOfRebinning->GetName());
  TH2D *rebinnedHistogram = new TH2D(newName,newName,fnRebinDeltaPhi,fRebinDeltaPhi,fnRebinDeltaEta,fRebinDeltaEta);
  
  // Loop over all the bins in the old histogram and insert the content to the new histogram
  for(int iDeltaPhi = 1; iDeltaPhi <= histogramInNeedOfRebinning->GetNbinsX(); iDeltaPhi++){
    deltaPhiValue = histogramInNeedOfRebinning->GetXaxis()->GetBinCenter(iDeltaPhi);
    for(int iDeltaEta = 1; iDeltaEta <= histogramInNeedOfRebinning->GetNbinsY(); iDeltaEta++){
      deltaEtaValue = histogramInNeedOfRebinning->GetYaxis()->GetBinCenter(iDeltaEta);
      
      // Find the global bin index from the new histogram correcponding to the bin in the old histogram
      newHistogramIndex = rebinnedHistogram->FindBin(deltaPhiValue,deltaEtaValue);
      
      // Add the bin content from the old histogram to the new histogram, adding errors in quadrature
      nonRebinnedContent = histogramInNeedOfRebinning->GetBinContent(iDeltaPhi,iDeltaEta);
      nonRebinnedError = histogramInNeedOfRebinning->GetBinError(iDeltaPhi,iDeltaEta);  
      currentBinContent = rebinnedHistogram->GetBinContent(newHistogramIndex);
      currentBinError = rebinnedHistogram->GetBinError(newHistogramIndex);
      rebinnedHistogram->SetBinContent(newHistogramIndex,currentBinContent+nonRebinnedContent);
      rebinnedHistogram->SetBinError(newHistogramIndex,TMath::Sqrt(TMath::Power(currentBinError,2)+TMath::Power(nonRebinnedError,2)));
    } // DeltaEta loop
  } // DeltaPhi loop
  
  // After rebinning, we need to do bin area normalization for each bin
  double binWidthDeltaPhi;
  double binWidthDeltaEta;
  for(int iDeltaPhi = 1; iDeltaPhi <= rebinnedHistogram->GetNbinsX(); iDeltaPhi++){
    binWidthDeltaPhi = rebinnedHistogram->GetXaxis()->GetBinWidth(iDeltaPhi);
    for(int iDeltaEta = 1; iDeltaEta <= rebinnedHistogram->GetNbinsY(); iDeltaEta++){
      binWidthDeltaEta = rebinnedHistogram->GetYaxis()->GetBinWidth(iDeltaEta);
      currentBinContent = rebinnedHistogram->GetBinContent(iDeltaPhi,iDeltaEta);
      currentBinError = rebinnedHistogram->GetBinError(iDeltaPhi,iDeltaEta);
      rebinnedHistogram->SetBinContent(iDeltaPhi,iDeltaEta,currentBinContent/(binWidthDeltaEta*binWidthDeltaPhi));  // TODO: If already normalized by bin area, need to
      rebinnedHistogram->SetBinError(iDeltaPhi,iDeltaEta,currentBinError/(binWidthDeltaEta*binWidthDeltaPhi));      // undo that here before renormalizing
    }
  }
  
  return rebinnedHistogram;
}

/*
 * Checker that new bin boundaries correspond to old ones
 *
 *  Arguments:
 *   const int nCheckedBins = Number of new bins to be checked
 *   const double *checkedBins = Bin boundaries to be checked
 *   const TAxis *originalAxis = Original axis against with the new bins are checked
 *
 *   return: True, if all the new bin boundaries can be found from the original axis. False if not.
 */
bool DijetMethods::CheckBinBoundaries(const int nCheckedBins, const double *checkedBins, const TAxis *originalAxis){
  
  // Flag, if the current bin is a bin boundary in the histogram to be rebinned
  bool binOK = false;
  
  // First, check that the bin boundaries for the rebinned histogram match with the old bin boundaries
  for(int iCheckedBin = 0; iCheckedBin < nCheckedBins + 1; iCheckedBin++){
    binOK = false;
    for(int iOldBin = 1; iOldBin <= originalAxis->GetNbins()+1; iOldBin++){
      
      // We the bin edge is close enough to one original bin, accept the bin
      if(TMath::Abs(originalAxis->GetBinLowEdge(iOldBin)-checkedBins[iCheckedBin]) < 1e-4){
        binOK = true;
        break;
      }
      
    } // Loop over bins in the original axis
    if(!binOK){ // If the bin is not in original histogram, print error message and return false
      std::cout << "The bin boundary " << checkedBins[iCheckedBin] << " is not a bin boundary in the original histogram!" << std::endl;
      return false;
    }
  } // Loop over bins to be checked
  
  // If all is good, return true
  return true;
}

// Getter for most recent two-dimensional jet shape count histogram
TH2D* DijetMethods::GetJetShapeBinMap() const{
  return fhJetShapeBinMap;
}

// Getter for most recent jet shape count histogram
TH1D* DijetMethods::GetJetShapeCounts() const{
  return fhJetShapeCounts;
}

// Getter for the most recent normalized mixed event histogram
TH2D* DijetMethods::GetNormalizedMixedEvent() const{
  return fNormalizedMixedEventHistogram;
}

// Getter for the most recent background distribution used to subtract the background
TH2D* DijetMethods::GetBackground() const{
  return fBackgroundDistribution;
}

// Getter for the most recent background overlap distribution for normalization check
TH2D* DijetMethods::GetBackgroundOverlap() const{
  return fBackgroundOverlap;
}

// Getter for deltaEta distribution on background deltaPhi region used for seagull fit
TH1D* DijetMethods::GetBackgroundEta() const{
  return fBackgroundEtaProjection;
}

// Getter for the most recent seagull fit
TF1* DijetMethods::GetSeagullFit() const{
  return fSeagullFit;
}

// Getter for the most recent spillover deltaEta distribution
TH1D* DijetMethods::GetSpilloverDeltaEta() const{
  return fSpilloverDeltaEta;
}

// Getter for the most recent spillover deltaPhi distribution
TH1D* DijetMethods::GetSpilloverDeltaPhi() const{
  return fSpilloverDeltaPhi;
}

// Getter for the most recent fit to spillover deltaEta distribution
TF1* DijetMethods::GetSpilloverDeltaEtaFit() const{
  return fSpilloverFitDeltaEta;
}

// Getter for the most recent fit to spillover deltaPhi distribution
TF1* DijetMethods::GetSpilloverDeltaPhiFit() const{
  return fSpilloverFitDeltaPhi;
}

// Getter for the most recent spillover yield error
double DijetMethods::GetSpilloverYieldError() const{
  return fSpilloverYieldError;
}

// Getter for background error scaling factor
double DijetMethods::GetBackgroundErrorScalingFactor() const{
  return fBackgroundErrorScalingFactor;
}

// Get the number of bins projected over in the previously done projection
int DijetMethods::GetNBinsProjectedOver() const{
  return fnBinsProjectedOver;
}

// Setter for deltaEta range used for normalizing the mixed event
void DijetMethods::SetMixedEventFitRegion(const double etaRangeLow, const double etaRangeHigh){
  
  if(etaRangeLow > etaRangeHigh){
    fMixedEventFitRegionLow = -etaRangeLow;
    fMixedEventFitRegionHigh = etaRangeLow;
  } else {
    fMixedEventFitRegionLow = etaRangeLow;
    fMixedEventFitRegionHigh = etaRangeHigh;
  }
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

// Setter for signal deltaEta region
void DijetMethods::SetSignalDeltaEtaRegion(const double maxDeltaEta){
  fMaxSignalDeltaEta = maxDeltaEta;
}

// Setter for background adjustment
void DijetMethods::SetBackgroundAdjustment(const bool adjust, const int overlapBins){
  fAdjustBackground = adjust;
  fnOverlapBins = overlapBins;
}

// Setter for new bin borders
void DijetMethods::SetBinBoundaries(const int nBins, double *binBorders, int& copyNbins, double *copyBinBorders[]){
  
  // Copy the number of bins and bin contents
  copyNbins = nBins;
  *copyBinBorders = new double[nBins+1];
  for(int iBin = 0; iBin < nBins+1; iBin++){
    (*copyBinBorders)[iBin] = binBorders[iBin];
  }
}

// Setter for R-binning for jet shape histograms
void DijetMethods::SetJetShapeBinEdges(const int nBins, double *binBorders){
  if(fRBins) delete [] fRBins;  // Delete the memory allocation before allocating new memory
  SetBinBoundaries(nBins,binBorders,fnRBins,&fRBins);
}

// Setter for two-dimensional histogram rebinning information
void DijetMethods::SetRebinBoundaries(const int nRebinDeltaEta, double *deltaEtaBorders, const int nRebinDeltaPhi, double *deltaPhiBorders){
  if(fRebinDeltaEta) delete [] fRebinDeltaEta;  // Delete the memory allocation before allocating new memory
  SetBinBoundaries(nRebinDeltaEta,deltaEtaBorders,fnRebinDeltaEta,&fRebinDeltaEta);
  if(fRebinDeltaPhi) delete [] fRebinDeltaPhi;  // Delete the memory allocation before allocating new memory
  SetBinBoundaries(nRebinDeltaPhi,deltaPhiBorders,fnRebinDeltaPhi,&fRebinDeltaPhi);
}

/*
 * Check the sanity of a given index for normalization methods
 *
 *  Arguments:
 *   const int normalizationType = Index for which the saniry check is made
 *   maxIndex = First index out of bounds
 *
 *   return: 0 for negative input, maxIndex-1 for too large input, just the input for valid input
 */
int DijetMethods::CheckNormalizationSanity(const int normalizationType, const int maxIndex){
  if(normalizationType < 0) return 0;
  if(normalizationType >= maxIndex) return maxIndex-1;
  return normalizationType;
}

// Setter for jet shape normalization method
void DijetMethods::SetJetShapeNormalization(const int normalizationType){
  fJetShapeNormalizationMethod = CheckNormalizationSanity(normalizationType,knJetShapeNormalizations);
}

// Setter for mixed event normalization method
void DijetMethods::SetMixedEventNormalization(const int normalizationType, const bool smoothenMixing){
  fMixedEventNormalizationMethod = CheckNormalizationSanity(normalizationType,knMixedEventNormalizations);
  fSmoothMixing = smoothenMixing;
}

// Setter for deltaPhi region considered as background in seagull correction and mixed event improvising
void DijetMethods::SetBackgroundDeltaPhiRegion(const double minDeltaPhi, const double maxDeltaPhi){
  fMinBackgroundDeltaPhi = minDeltaPhi;
  fMaxBackgroundDeltaPhi = maxDeltaPhi;
}

// Setter for the amount of rebin applied to deltaEta histogram before fit in seagull correction
void DijetMethods::SetSeagullRebin(const int nRebin){
  fSeagullRebin = nRebin;
}

// Getter for positive level in pair acceptance systematic uncertainty estimation
double DijetMethods::GetPairAcceptancePositiveLevel(){
  return fPairAcceptancePositiveLevel;
}

// Getter for negative level in pair acceptance systematic uncertainty estimation
double DijetMethods::GetPairAcceptanceNegativeLevel(){
  return fPairAcceptanceNegativeLevel;
}

// Getter for inner mean in background subtraction systematic uncertainty estimation
double DijetMethods::GetBackgroundSubtractionInnerMean(){
  return fBackgroundSubtractionInnerMean;
}

// Getter for outer mean in background subtraction systematic uncertainty estimation
double DijetMethods::GetBackgroundSubtractionOuterMean(){
  return fBackgroundSubtractionOuterMean;
}
