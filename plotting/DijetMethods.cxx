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
  fMixedEventNormalizationMethod(kSingle),
  fMinBackgroundDeltaEta(1.5),
  fMaxBackgroundDeltaEta(2.5),
  fJetShapeNormalizationMethod(kBinWidth),
  fnRBins(0),
  fnRebinDeltaEta(0),
  fnRebinDeltaPhi(0)
{
  // Constructor
  fBackgroundDistribution = NULL;
  fhJetShapeCounts = NULL;
  fhJetShapeBinMap = NULL;
  fRBins = nullptr;
  fRebinDeltaEta = nullptr;
  fRebinDeltaPhi = nullptr;
}

/*
 *  Copy constructor
 */
DijetMethods::DijetMethods(const DijetMethods& in) :
  fMixedEventFitRegion(in.fMixedEventFitRegion),
  fMixedEventNormalizationMethod(in.fMixedEventNormalizationMethod),
  fBackgroundDistribution(in.fBackgroundDistribution),
  fMinBackgroundDeltaEta(in.fMinBackgroundDeltaEta),
  fMaxBackgroundDeltaEta(in.fMaxBackgroundDeltaEta),
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
  fMinBackgroundDeltaEta = in.fMinBackgroundDeltaEta;
  fMaxBackgroundDeltaEta = in.fMaxBackgroundDeltaEta;
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
  char histogramName[200];     // Helper variable for histogram naming
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
  TH2D *backgroundSubtractedHistogram = (TH2D*)leadingHistogramWithBackground->Clone(histogramName);
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
      
      // Get the content to be added from the two-dimensional histogram
      newContent = backgroundSubtractedHistogram->GetBinContent(iPhiBin,iEtaBin);
      newError = backgroundSubtractedHistogram->GetBinError(iPhiBin,iEtaBin);
      
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
