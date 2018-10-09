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
 *  Constructor
 */
DijetMethods::DijetMethods() :
  fMixedEventFitRegion(0.2),
  fMixedEventNormalizationMethod(kSingle),
  fMinBackgroundDeltaEta(1.5),
  fMaxBackgroundDeltaEta(2.5),
  fMinBackgroundDeltaPhi(1.4),
  fMaxBackgroundDeltaPhi(1.7),
  fSeagullRebin(4),
  fMaxSignalDeltaEta(1.0),
  fJetShapeNormalizationMethod(kBinWidth),
  fnRebinDeltaEta(0),
  fnRebinDeltaPhi(0)
{
  // Constructor
  fBackgroundDistribution = NULL;
  fBackgroundOverlap = NULL;
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
  fMixedEventFitRegion(in.fMixedEventFitRegion),
  fMixedEventNormalizationMethod(in.fMixedEventNormalizationMethod),
  fBackgroundDistribution(in.fBackgroundDistribution),
  fBackgroundOverlap(in.fBackgroundOverlap),
  fMinBackgroundDeltaEta(in.fMinBackgroundDeltaEta),
  fMaxBackgroundDeltaEta(in.fMaxBackgroundDeltaEta),
  fBackgroundEtaProjection(in.fBackgroundEtaProjection),
  fSeagullFit(in.fSeagullFit),
  fMinBackgroundDeltaPhi(in.fMinBackgroundDeltaPhi),
  fMaxBackgroundDeltaPhi(in.fMaxBackgroundDeltaPhi),
  fSeagullRebin(in.fSeagullRebin),
  fMaxSignalDeltaEta(in.fMaxSignalDeltaEta),
  fJetShapeNormalizationMethod(in.fJetShapeNormalizationMethod),
  fhJetShapeCounts(in.fhJetShapeCounts),
  fhJetShapeBinMap(in.fhJetShapeBinMap)
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
  
  fMixedEventFitRegion = in.fMixedEventFitRegion;
  fMixedEventNormalizationMethod = in.fMixedEventNormalizationMethod;
  fBackgroundDistribution = in.fBackgroundDistribution;
  fBackgroundOverlap = in.fBackgroundOverlap;
  fMinBackgroundDeltaEta = in.fMinBackgroundDeltaEta;
  fMaxBackgroundDeltaEta = in.fMaxBackgroundDeltaEta;
  fBackgroundEtaProjection = in.fBackgroundEtaProjection;
  fSeagullFit = in.fSeagullFit;
  fMinBackgroundDeltaPhi = in.fMinBackgroundDeltaPhi;
  fMaxBackgroundDeltaPhi = in.fMaxBackgroundDeltaPhi;
  fSeagullRebin = in.fSeagullRebin;
  fMaxSignalDeltaEta = in.fMaxSignalDeltaEta;
  fJetShapeNormalizationMethod = in.fJetShapeNormalizationMethod;
  fhJetShapeCounts = in.fhJetShapeCounts;
  fhJetShapeBinMap = in.fhJetShapeBinMap;
  
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
 * Do the mixed event correction and return corrected TH2D.
 * The idea here is, that we divide the same event histogram with a normalized mixed event histogram.
 * The normalization is obtained by first projecting the deltaEta distribution out of the
 * two-dimensional histogram and then fitting a constant to the central region in deltaEta.
 * This same is done for leading and subleading jet-track mixed event histograms and an average
 * is taken from these. This is done to ensure that the normalization is the same, because both
 * leading and subleading distributions are needed for the background subrtaction.
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
  double subleadingScale = 0;
  if(fMixedEventNormalizationMethod == kAverage) subleadingScale = GetMixedEventScale(subleadingMixedEventHistogram);
  double averageScale = (leadingScale+subleadingScale)/2.0;
  double normalizationScale = leadingScale;
  if(fMixedEventNormalizationMethod == kAverage) normalizationScale = averageScale;

  // Normalize the mixed event histogram and do the correction
  leadingMixedEventHistogram->Scale(1.0/normalizationScale);
  correctedHistogram->Divide(leadingMixedEventHistogram);
  
  // We need to rescale the mixed event histogram in the end, because it might later be used for another correction
  leadingMixedEventHistogram->Scale(normalizationScale);
  
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
 *  TH2D *mixedEventCorrectedHistogram = Two-dimensional deltaEta-deltaPhi distribution, that is corrected by the mixed event
 *
 *  return: The two dimensional distribution corrected for the seagull effect
 */
TH2D* DijetMethods::DoSeagullCorrection(TH2D *mixedEventCorrectedHistogram){
  
  // Project out the deltaEta distribution at background deltaPhi region
  int minDeltaPhiBin = mixedEventCorrectedHistogram->GetXaxis()->FindBin(fMinBackgroundDeltaPhi+0.0001);
  int maxDeltaPhiBin = mixedEventCorrectedHistogram->GetXaxis()->FindBin(fMaxBackgroundDeltaPhi+0.0001);
  fBackgroundEtaProjection = mixedEventCorrectedHistogram->ProjectionY("SeagullEtaProjection",minDeltaPhiBin,maxDeltaPhiBin);
  
  // Scale the projected deltaEta distribution with the number of deltaPhi bins projected over and deltaPhi bin width to retain normalization
  fBackgroundEtaProjection->Scale(1.0/(maxDeltaPhiBin-minDeltaPhiBin+1));
  fBackgroundEtaProjection->Scale(mixedEventCorrectedHistogram->GetXaxis()->GetBinWidth(1));
  fBackgroundEtaProjection->Rebin(fSeagullRebin);
  fBackgroundEtaProjection->Scale(1.0/fSeagullRebin);
  
  // Prepare the fit function
  fSeagullFit = new TF1("seagullFit",seagullPoly2,-3,3,3);
  double initialLevel = fBackgroundEtaProjection->GetBinContent(fBackgroundEtaProjection->FindBin(0));
  fSeagullFit->SetParameters(initialLevel,0,0);
  
  // Fit the projected distribution
  fBackgroundEtaProjection->Fit(fSeagullFit,"","",-3,3);
  double backgroundLevel = fSeagullFit->GetParameter(0);
  
  // Apply the correction to the input 2D histogram
  TH2D *seagullCorrectedHistogram = (TH2D*) mixedEventCorrectedHistogram->Clone();
  double binEta;
  double seagullCorrection;
  double binContent;
  
  for(int iEta = 1; iEta <= seagullCorrectedHistogram->GetNbinsY(); iEta++){
    
    // Calculate the correction for this deltaEta bin
    binEta = seagullCorrectedHistogram->GetYaxis()->GetBinCenter(iEta);
    seagullCorrection = backgroundLevel/fSeagullFit->Eval(binEta);
    
    // Apply the correction for all deltaPhi bins in the deltaEta strip
    for(int iPhi = 1; iPhi <= seagullCorrectedHistogram->GetNbinsX(); iPhi++){
      binContent = seagullCorrectedHistogram->GetBinContent(iPhi,iEta);
      seagullCorrectedHistogram->SetBinContent(iPhi,iEta,binContent*seagullCorrection);
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
TH2D* DijetMethods::SubtractBackground(TH2D *leadingHistogramWithBackground, TH2D *subleadingHistogramWithBackground, bool isInclusive){
  
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
  
  sprintf(histogramName,"%sBackground",leadingHistogramWithBackground->GetName());
  fBackgroundDistribution = (TH2D*)leadingHistogramWithBackground->Clone(histogramName);
  int nDeltaPhiBins = fBackgroundDistribution->GetNbinsX(); // Number of deltaPhi bins
  
  // Initialize also a distribution for background overlap for normalization level check
  sprintf(histogramName,"%sBackgroundOverlap",leadingHistogramWithBackground->GetName());
  fBackgroundOverlap = (TH2D*)leadingHistogramWithBackground->Clone(histogramName);
  
  // We need the deltaEta bin width for normalization purposes
  double binWidthDeltaEta = leadingHistogramWithBackground->GetYaxis()->GetBinWidth(1);
  int offset = isInclusive ? nDeltaPhiBins/2 : 0;  // Apply offset for inclusive histograms to scan over whole deltaPhi space

  // Loop over deltaPhi bins and fill the leading jet-track correlation result in the near side
  // and the subleading set-track correlation result in the away side
  for(int iDeltaPhi = 1; iDeltaPhi <= nDeltaPhiBins/2; iDeltaPhi++){
    
    // Read the values from the deltaPhi background histogram
    deltaPhiValueNear = backgroundDeltaPhiLeading->GetBinContent(iDeltaPhi);
    deltaPhiErrorNear = backgroundDeltaPhiLeading->GetBinError(iDeltaPhi);
    deltaPhiValueAway = backgroundDeltaPhiSubleading->GetBinContent(iDeltaPhi+offset);
    deltaPhiErrorAway = backgroundDeltaPhiSubleading->GetBinError(iDeltaPhi+offset);
    
    // Repopulate the deltaEta axis
    for(int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
      
      // Insert the values to the two-dimensional background histogram
      fBackgroundDistribution->SetBinContent(iDeltaPhi,iDeltaEta,deltaPhiValueNear/binWidthDeltaEta);
      fBackgroundDistribution->SetBinError(iDeltaPhi,iDeltaEta,deltaPhiErrorNear/binWidthDeltaEta);
      fBackgroundDistribution->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,deltaPhiValueAway/binWidthDeltaEta);
      fBackgroundDistribution->SetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,deltaPhiErrorAway/binWidthDeltaEta);
      
      // Set the contants of the background overlap to zero
      fBackgroundOverlap->SetBinContent(iDeltaPhi,iDeltaEta,0);
      fBackgroundOverlap->SetBinError(iDeltaPhi,iDeltaEta,0);
      fBackgroundOverlap->SetBinContent(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,0);
      fBackgroundOverlap->SetBinError(iDeltaPhi+nDeltaPhiBins/2,iDeltaEta,0);
      
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
  
  // For the background overlap distribution, fill three bins on each side of the gluing point
  for(int iDeltaPhi = nDeltaPhiBins/2 + 1; iDeltaPhi <= nDeltaPhiBins/2 + 3; iDeltaPhi++){
    
    // Read the values from the deltaPhi background histogram
    deltaPhiValueNear = backgroundDeltaPhiLeading->GetBinContent(iDeltaPhi);
    deltaPhiErrorNear = backgroundDeltaPhiLeading->GetBinError(iDeltaPhi);
    
    // For the away side we need to continue to the left, which rotates to the three last bins in deltaPhi
    deltaPhiValueAway = backgroundDeltaPhiSubleading->GetBinContent(iDeltaPhi + nDeltaPhiBins/2 - 3);
    deltaPhiErrorAway = backgroundDeltaPhiSubleading->GetBinError(iDeltaPhi + nDeltaPhiBins/2 - 3);
    
    // Repopulate the deltaEta axis
    for(int iDeltaEta = 1; iDeltaEta <= fBackgroundDistribution->GetNbinsY(); iDeltaEta++){
      // Set the contants of the background overlap to zero
      fBackgroundOverlap->SetBinContent(iDeltaPhi,iDeltaEta,deltaPhiValueNear/binWidthDeltaEta);
      fBackgroundOverlap->SetBinError(iDeltaPhi,iDeltaEta,deltaPhiErrorNear/binWidthDeltaEta);
      fBackgroundOverlap->SetBinContent(iDeltaPhi-3,iDeltaEta,deltaPhiValueAway/binWidthDeltaEta);
      fBackgroundOverlap->SetBinError(iDeltaPhi-3,iDeltaEta,deltaPhiErrorAway/binWidthDeltaEta);
    } // DeltaEta loop
  } // DeltaPhi loop
  
  // Subtract the generated background from the histogram including the background
  sprintf(histogramName,"%sBackgroundSubtracted",leadingHistogramWithBackground->GetName());
  TH2D *backgroundSubtractedHistogram = (TH2D*)leadingHistogramWithBackground->Clone(histogramName);
  backgroundSubtractedHistogram->Add(fBackgroundDistribution,-1);
  
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
 *
 *  return: Two-dimensional histogram with the estimation of the spillover effect
 */
TH2D* DijetMethods::GetSpilloverCorrection(TH2D *onlyHydjetHistogram){
  
  // Define the range over which the spillover correction is calculated
  double spilloverEtaRange = 1.5;
  double spilloverPhiRange = 1.5;
  
  // First get the projections for the deltaEta and deltaPhi distributions
  TH1D *deltaEtaProjection = ProjectRegionDeltaEta(onlyHydjetHistogram,-spilloverPhiRange,spilloverPhiRange,"spillover");
  TH1D *deltaPhiProjection = ProjectRegionDeltaPhi(onlyHydjetHistogram,0,spilloverEtaRange,"spillover");
  
  // Fit one dimensional Gaussians to the projected deltaEta and deltaPhi distributions
  TF1 *gaussForEta = FitGauss(deltaEtaProjection,spilloverEtaRange);
  TF1 *gaussForPhi = FitGauss(deltaPhiProjection,spilloverPhiRange);
  
  // Combine the one-dimensional fits to a two-dimensional Gaussian distribution
  TF2 *gauss2D = new TF2("gauss2D", "[0]/(2*TMath::Pi()*[1]*[2])*TMath::Exp(-0.5*TMath::Power(x/[1],2))*TMath::Exp(-0.5*TMath::Power(y/[2],2))",-TMath::Pi()/2,3*TMath::Pi()/2,-5,5);
  gauss2D->SetParameter(0,gaussForEta->GetParameter(0));
  gauss2D->SetParameter(1,gaussForPhi->GetParameter(1));
  gauss2D->SetParameter(2,gaussForEta->GetParameter(1));
  
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
  
  // Return the spillover correction
  return spilloverCorrection;
}

/*
 * Fit a Gaussian function to a histogram and return the fit function
 *
 * Arguments:
 *  TH1D *fittedHistogram = Histogram to be fitted with the Gaussian
 *  double fitRange = Range for the fit
 *
 *  return: Gaussian function fitted to the histogram
 */
TF1* DijetMethods::FitGauss(TH1D *fittedHistogram, double fitRange){
  
  // Create a Gaussian function for the fit
  TF1 *gaussFit = new TF1("gaussFit","[0]/(TMath::Sqrt(2*TMath::Pi())*[1])*TMath::Exp(-0.5*TMath::Power(x/[1],2))",-fitRange,fitRange);
  
  // Calculate the integral of the fitted histogram over the fit range
  double fitYield = fittedHistogram->Integral(fittedHistogram->FindBin(-fitRange+0.001),fittedHistogram->FindBin(fitRange-0.001),"width");
  
  // Fix the yield of the Gaussian function to the integral and give some rought estimate for the width
  gaussFit->FixParameter(0,fitYield);
  gaussFit->SetParameter(1,0.3);
  
  // Fit the histogram over the defined range
  fittedHistogram->Fit("gaussFit","","",-fitRange,fitRange);
  
  // Return the fit function
  return gaussFit;
  
}

/*
 * Project the deltaPhi distribution out of a two-dimensional deltaPhi-deltaEta distribution
 *
 *  Arguments:
 *   TH2D* deltaPhiDeltaEtaHistogram = two-dimensional deltaPhi-deltaEta distribution
 *   const double minDeltaEta = Minimum deltaEta value for the projected region
 *   const double maxDeltaEta = Maximum deltaEta value for the projected region
 *   const char* newName = Name that is added to the projected histogram
 *
 *   return: DeltaPhi distribution projected from the input deltaEta region
 */
TH1D* DijetMethods::ProjectRegionDeltaEta(TH2D* deltaPhiDeltaEtaHistogram, const double minDeltaPhi, const double maxDeltaPhi, const char* newName){
  
  // Start by finding the bin indices for the defined deltaEta region
  // Apply a little offset for the defined eta region borders to avoid bin border effects
  int lowDeltaPhiBin  = deltaPhiDeltaEtaHistogram->GetXaxis()->FindBin(minDeltaPhi+0.0001);
  int highDeltaPhiBin = deltaPhiDeltaEtaHistogram->GetXaxis()->FindBin(maxDeltaPhi-0.0001);
  
  // Calculate the number of deltaPhi bins in the projected region
  int nBinsProjectedRegion = highDeltaPhiBin-lowDeltaPhiBin+1;
  
  // Project out the deltaPhi distribution in the background
  char histogramName[200];
  sprintf(histogramName,"%s%s",deltaPhiDeltaEtaHistogram->GetName(),newName);
  TH1D *projectedDeltaEta;
  projectedDeltaEta = deltaPhiDeltaEtaHistogram->ProjectionY(histogramName,lowDeltaPhiBin,highDeltaPhiBin);
  
  // Scale the projected deltaEta distribution with the number of deltaPhi bins projected over and deltaPhi bin width to retain normalization
  projectedDeltaEta->Scale(1.0/nBinsProjectedRegion);
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
 *
 *   return: DeltaPhi distribution projected from the input deltaEta region
 */
TH1D* DijetMethods::ProjectRegionDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram, const double minDeltaEta, const double maxDeltaEta, const char* newName){
  
  // If the minimum is zero, do the projection from -maxDeltaEta to maxDeltaEta
  bool oneRegion = (minDeltaEta < 0.01);
  
  // Start by finding the bin indices for the defined deltaEta region
  // Apply a little offset for the defined eta region borders to avoid bin border effects
  int lowNegativeDeltaEtaBin  = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(-maxDeltaEta+0.0001);
  int highNegativeDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(-minDeltaEta-0.0001);
  int lowPositiveDeltaEtaBin  = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(minDeltaEta+0.0001);
  int highPositiveDeltaEtaBin = deltaPhiDeltaEtaHistogram->GetYaxis()->FindBin(maxDeltaEta-0.0001);
  
  // Calculate the number of deltaEta bins in the projected region
  int nBinsProjectedRegion = highNegativeDeltaEtaBin-lowNegativeDeltaEtaBin+highPositiveDeltaEtaBin-lowPositiveDeltaEtaBin+2;
  if(oneRegion) nBinsProjectedRegion = nBinsProjectedRegion-highNegativeDeltaEtaBin+lowPositiveDeltaEtaBin-1;
  
  // Project out the deltaPhi distribution in the defined region
  char histogramName[200];
  sprintf(histogramName,"%s%s",deltaPhiDeltaEtaHistogram->GetName(),newName);
  TH1D *projectedDeltaPhi;
  if(oneRegion) {
    projectedDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName,lowNegativeDeltaEtaBin,highPositiveDeltaEtaBin);
  } else {
    projectedDeltaPhi = deltaPhiDeltaEtaHistogram->ProjectionX(histogramName,lowNegativeDeltaEtaBin,highNegativeDeltaEtaBin);
    projectedDeltaPhi->Add(deltaPhiDeltaEtaHistogram->ProjectionX("dummyName",lowPositiveDeltaEtaBin,highPositiveDeltaEtaBin));
  }
  
  // Scale the projected deltaPhi distribution with the number of deltaEta bins projected over and deltaEta bin width to retain normalization
  projectedDeltaPhi->Scale(1.0/nBinsProjectedRegion);
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
    }
    
  } else if (fJetShapeNormalizationMethod == kBinArea){
    jetShapeHistogram->Divide(fhJetShapeCounts);
  }
  
  // Return the calculated jet shape histogram
  return jetShapeHistogram;
  
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
 *   TAxis *originalAxis = Original axis against with the new bins are checked
 *
 *   return: True, if all the new bin boundaries can be found from the original axis. False if not.
 */
bool DijetMethods::CheckBinBoundaries(const int nCheckedBins, const double *checkedBins, TAxis *originalAxis){
  
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

// Setter for signal deltaEta region
void DijetMethods::SetSignalDeltaEtaRegion(const double maxDeltaEta){
  fMaxSignalDeltaEta = maxDeltaEta;
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
void DijetMethods::SetMixedEventNormalization(const int normalizationType){
  fMixedEventNormalizationMethod = CheckNormalizationSanity(normalizationType,knMixedEventNormalizations);
}

// Setter for deltaPhi region considered as background in seagull correction
void DijetMethods::SetBackgroundDeltaPhiRegionSeagull(const double minDeltaPhi, const double maxDeltaPhi){
  fMinBackgroundDeltaPhi = minDeltaPhi;
  fMaxBackgroundDeltaPhi = maxDeltaPhi;
}

// Setter for the amount of rebin applied to deltaEta histogram before fit in seagull correction
void DijetMethods::SetSeagullRebin(const int nRebin){
  fSeagullRebin = nRebin;
}
