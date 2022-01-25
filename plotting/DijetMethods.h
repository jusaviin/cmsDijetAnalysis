#ifndef DIJETMETHODS_H
#define DIJETMETHODS_H

/*
 * This class is a collection of methods that are used to process the
 * results produced by the dijet analysis
 */

// C++ includes
#include <iostream>
#include <tuple>

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TAxis.h>

class DijetMethods{
  
public:
  
  enum enumJetShapeNormalization{kBinWidth,kBinArea,knJetShapeNormalizations};   // Enumeration for used normalization for jet shape
  enum enumMixedEventNormalization{kSingle,kAverage,knMixedEventNormalizations}; // Enumeration for used normalization for mixed events
  
  DijetMethods();   // Constructor
  DijetMethods(const DijetMethods& in); // Copy constructor
  ~DijetMethods();  // Destructor
  DijetMethods& operator=(const DijetMethods& in); // Equal sign operator
  
  TH2D* MixedEventCorrect(const TH2D *sameEventHistogram, const TH2D *leadingMixedEventHistogram, const TH2D *subleadingMixedEventHistogram, const bool avoidPeaks); // Mixed event correction for a two-dimensional histogram
  TH2D* DoSeagullCorrection(const TH2D *mixedEventCorrectedHistogram, const int normalizationMethod = 0, const int vetoFlag = 0, const double middleArea = 0.5);  // Apply a seagull correction to the histogram
  TH2D* SubtractBackground(TH2D *leadingHistogramWithBackground, TH2D *subleadingHistogramWithBackground, double maxDeltaEta, bool isInclusive = false); // Subtract background from a two-dimensional leading histogram
  TH2D* ImproviseMixedEvent(const TH2D *sameEventHistogram); // Improvise mixed event distribution from background deltaPhi region of the same event histogram
  double GetSpilloverYield(TH2D *onlyHydjetHistogram, double minEtaNormalizationRange, double maxEtaNormalizationRange); // Getter for the dpillover yield from the mixed event corrected distribution
  double GetSpilloverYieldError() const; // Getter for the most recent spillover yield error
  TH2D* GetSpilloverCorrection(TH2D *onlyHydjetHistogram, int fitMethod, double spilloverEtaFitRange = 1.5, double spilloverPhiFitRange = 1.5, double lowConstantRange = 1, double highConstantRange = 2, double fixedYield = 0, double fixedEtaWidth = 0, double fixedPhiWidth = 0);  // Get the spillover correction from only hydjet histogram
  TH1D* GetJetShape(TH2D *backgroundSubtractedHistogram); // Extract the jet shape from the two-dimensional histogram
  TH1D* ProjectAnalysisYieldDeltaEta(TH2D *backgroundSubtractedHistogram, const double lowPtEdge, const double highPtEdge, const bool symmetrize = false, const bool splitZero = false);  // Project the deltaEta yield in final analysis binning
  TH1D* RebinAsymmetric(TH1D *histogramInNeedOfRebinning, const int nBins, const double* binEdges); // Asymmetric rebinning for one-dimensional histograms
  TH2D* RebinHistogram(TH2D *histogramInNeedOfRebinning); // Rebin a two-dimensional deltaPhi-deltaEta histogram
  TH2D* RebinHistogram(TH2D *histogramInNeedOfRebinning, const int nBinsX, const double* binBordersX, const int nBinsY, const double* binBordersY, const bool undoBinArea, const bool normalizeBinArea); // Rebin a two-dimensional histogram with given bin borders
  void SymmetrizeDeltaEta(TH1D *histogramInNeedOfSymmetrization); // Symmetrize a given deltaEta histogram
  TH2D* SymmetrizeHistogram(const TH2D *histogramToBeSymmetrized, const double maxR, const double fillValue = 0, const bool onlyCut = false); // Symmetrize eta and phi in a histogram up to given maximum radius
  TH1D* ProjectSignalDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram); // Project deltaPhi distribution in the signal region in eta out of a two-dimensional deltaPhi-deltaEta distribution
  TH1D* ProjectBackgroundDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram); // Project deltaPhi distribution in the background region out of a two-dimensional deltaPhi-deltaEta distribution
  TF1* FourierFit(TH1D* backgroundDeltaPhi, const int maxVn, const bool onlyNearSideFit = false); // Do a Fourier fit for the background deltaPhi distribution
  TH1D* CombineDeltaPhi(const TH2D *leadingHistogramWithBackground, const TH2D *subleadingHistogramWithBackground, const double minDeltaEta, const double maxDeltaEta, const char* newName, const bool oneSide = false);
  void NormalizeMatrix(TH2D *histogramInNeedOfNormalization, const double value = 1, const int direction = 1);  // Normalize rows or columns of a 2D histogram to a given value
  void NormalizeColumns(TH2D *histogramInNeedOfNormalization, const double value = 1);  // Normalize all the columns of a 2-D histogram to a given value
  void NormalizeRows(TH2D *histogramInNeedOfNormalization, const double value = 1);  // Normalize all the rows of a 2-D histogram to a given value
  TH2D* RotateHistogram(TH2D *originalHistogram); // Rotate two dimensional histogram 90 degrees
  
  // Project a region in one direction from a two-dimensional histogram
  TH1D* ProjectRegionDeltaPhi(const TH2D* deltaPhiDeltaEtaHistogram, const double minDeltaEta, const double maxDeltaEta, const char* newName, const bool oneSide = false);  // Project deltaPhi distribution out of a two-dimensional deltaPhi-deltaEta distribution
  TH1D* ProjectRegionDeltaEta(const TH2D* deltaPhiDeltaEtaHistogram, const double minDeltaPhi, const double maxDeltaPhi, const char* newName);  // Project deltaEta distribution out of a two-dimensional deltaPhi-deltaEta distribution
  int GetNBinsProjectedOver() const; // Get the number of bins projected over in the previously done projection
  
  // Methods for estimating the systematic uncertainties
  double EstimateSystematicsForPairAcceptanceCorrection(const TH1* deltaEtaHistogram);
  double EstimateSystematicsForBackgroundSubtraction(const TH1* deltaEtaHistogram);
  void PropagateDeltaEtaToDeltaR(TH1* errorHistogram, const double deltaEtaError);
  
  // Getters for QA values for systematic uncertaintie
  double GetPairAcceptancePositiveLevel();
  double GetPairAcceptanceNegativeLevel();
  double GetBackgroundSubtractionInnerMean();
  double GetBackgroundSubtractionOuterMean();
  
  // Getters for produced distributions
  TH2D* GetNormalizedMixedEvent() const; // Getter for the most recent normalized mixed event histogram
  TH2D* GetBackground() const;           // Getter for the most recent background distribution used to subtract the background
  TH2D* GetBackgroundOverlap() const;    // Getter for the most recent background overlap distribution for normalization check
  TH1D* GetJetShapeCounts() const;       // Getter for the jet shape count distribution
  TH2D* GetJetShapeBinMap() const;       // Getter for the map between R bins and deltaEta-deltaPhi bins
  TH1D* GetBackgroundEta() const;        // Getter for deltaEta distribution on background deltaPhi region used for seagull fit
  TF1* GetSeagullFit() const;            // Getter for the most recent seagull fit
  TH1D* GetSpilloverDeltaEta() const;    // Getter for projected eta distribution for spillover calculation
  TH1D* GetSpilloverDeltaPhi() const;    // Getter for projected phi distribution for spillover calculation
  TF1* GetSpilloverDeltaEtaFit() const;  // Getter for fit to projected eta distribution for spillover calculation
  TF1* GetSpilloverDeltaPhiFit() const;  // Getter for fit to projected phi distribution for spillover calculation
  
  // Getter for background error scaling factor
  double GetBackgroundErrorScalingFactor() const; // Getter for background error scaling factor
  std::tuple<double,double> GetBackgroundAdjustmentFactors(TH1 *leadingBackground, TH1 *leadingBackgroundOverlap) const; // Getter for background adjustment factors
  
  // Setters for mixed event configuration
  void SetMixedEventFitRegion(const double etaRangeLow, const double etaRangeHigh = 0);  // Setter for deltaEta range used for normalizing the mixed event
  void SetMixedEventNormalization(const int normalizationType, const bool smoothenMixing); // Setter for normalization method used for mixed event distributions
  
  // Setters for background subtraction configuration
  void SetBackgroundDeltaEtaRegion(const double minDeltaEta, const double maxDeltaEta, const bool oneSide = false); // Setter for background deltaEta region
  void SetSignalDeltaEtaRegion(const double maxDeltaEta);                               // Setter for signal deltaEta region
  void SetBackgroundAdjustment(const bool adjust, const int overlapBins);               // Setter for background adjustment
  
  // Setters for jet shape calculation configuration
  void SetJetShapeBinEdges(const int nBins, double *binBorders); // Setter for R-binning for jet shape histograms
  void SetJetShapeNormalization(const int normalizationType);    // Setter for jet shape normalization method
  
  // Setter for two-dimensional histogram rebinning information
  void SetRebinBoundaries(const int nRebinDeltaEta, double *deltaEtaBorders, const int nRebinDeltaPhi, double *deltaPhiBorders); // Setter for deltaEta and deltaPhi rebin borders
  
  // Setter for parameters for seagull correction
  void SetBackgroundDeltaPhiRegion(const double minDeltaPhi, const double maxDeltaPhi);  // Setter for deltaPhi region considered as background in seagull correction and mixed event improvising
  void SetSeagullRebin(const int nRebin);  // Setter for the amount of rebin applied to deltaEta histogram before fit in seagull correction
  
private:
  
  // =============================================
  // =========== Mixed event correction ==========
  // =============================================
  
  TH2D *fNormalizedMixedEventHistogram; // Mixed event histogram after normalization
  double fMixedEventFitRegionLow;  // Low bound of the region in deltaEta in which a constant fit is done to normalize mixed event distributions
  double fMixedEventFitRegionHigh;  // High bound of the region in deltaEta in which a constant fit is done to normalize mixed event distributions
  int fMixedEventNormalizationMethod; // Normalization method used for mixed event distributions
  bool fSmoothMixing; // Smoothen the mixing distribution in phi
  double fMaximumDeltaEta;  // Maximum allowed deltaEta in the corrected distribution
  
  // =============================================
  // =========== Background subtraction ==========
  // =============================================
  
  TH2D *fBackgroundDistribution;        // Remember the background distribution from the most recent background subtraction
  TH2D *fBackgroundOverlap;             // Fill a few points over the gluing point to see how well the gluing works
  double fMinBackgroundDeltaEta;        // Minimum deltaEta for background subtraction region
  double fMaxBackgroundDeltaEta;        // Maximum deltaEta for background subtraction region
  bool fOneBackgroundRegion;            // Only use positive or negative eta for background estimation
  bool fAdjustBackground;               // Adjust the level of background based on the ratio of leading and subleading sides in overlapping bins
  int fnOverlapBins;                    // Number of overlapping bins to be used in background adjustment and overlap
  double fBackgroundErrorScalingFactor; // Scaling factor to be used for errors when projecting deltaPhi out of whole distribution
  
  // =============================================
  // ============ Seagull correction =============
  // =============================================
  
  TH1D *fBackgroundEtaProjection;  // Projected deltaEta in the background region of deltaPhi
  TF1 *fSeagullFit;                // Function used to do the seagull fit
  double fMinBackgroundDeltaPhi;   // Minimum deltaPhi for seagull fit region
  double fMaxBackgroundDeltaPhi;   // Maximum deltaPhi for seagull fit region
  int fSeagullRebin;               // Rebin applied to deltaEta histogram before fitting
  double fSeagullChi2Limit;        // Maximum chi2 allowed for constant fit such that correction is not applied
  
  // =============================================
  // ============ DeltaPhi projections ===========
  // =============================================
  int fnBinsProjectedOver;        // Number of bins projected over when projection a region of histogram
  double fMaxSignalDeltaEta;      // Maximum deltaEta value accepted for the signal region
  
  // =============================================
  // ============ Spillover correction ===========
  // =============================================
  TH1D *fSpilloverDeltaEta;       // DeltaEta projection for spillover calculation
  TH1D *fSpilloverDeltaPhi;       // DeltaPhi projection for spillover calculation
  TF1 *fSpilloverFitDeltaEta;     // Fit to deltaEta projection during spillover calculation
  TF1 *fSpilloverFitDeltaPhi;     // Fit to deltaPhi projection during spillover calculation
  double fSpilloverYieldError;    // Error in the latest spillover yield calculation
  
  // =============================================
  // =========== Jet shape calculation ===========
  // =============================================
  
  int fJetShapeNormalizationMethod; // Normalization method used for jet shape distribution, either bin area or bin width
  TH1D *fhJetShapeCounts;   // How many bins from the two-dimensional histogram correspond to one jet shape bin
  TH2D *fhJetShapeBinMap;   // Information to which jet shape bin each deltaPhi-deltaEta bin is assigned
  int fnRBins;              // Number of R-bins for jet shape histograms
  double *fRBins;           // R-bin boundaries for jet shape histograms
  
  // ==============================================
  // ==== Rebinning two-dimensional histograms ====
  // ==============================================
  int fnRebinDeltaEta;    // Number of new deltaEta bins
  double* fRebinDeltaEta; // Bin boundaries for the new deltaEta bins
  int fnRebinDeltaPhi;    // Number of new deltaPhi bins
  double* fRebinDeltaPhi; // Bin boundaries for the new deltaPhi bins
  
  // ==============================================
  // ===== Systematic uncertainty estimation ======
  // ==============================================
  double fPairAcceptancePositiveLevel;
  double fPairAcceptanceNegativeLevel;
  double fBackgroundSubtractionInnerMean;
  double fBackgroundSubtractionOuterMean;
  
  // Private methods
  double GetMixedEventScale(const TH2D* mixedEventHistogram, const bool findPeak); // Find the normalization scale for the mixed event histogram
  TF1* FitGauss(TH1D* fittedHistogram, double fitRange, double normalizationRange, double fixedYield, double fixedWidth, double fixedConstant = 0);  // Fit a Gaussian function to a histogram and return the fit function
  TF1* FitGaussAndConstant(TH1D* fittedHistogram, double fitRange, double normalizationRange, double lowConstantRange, double highConstantRange, double fixedYield, double fixedWidth);  // Fit a Gaussian function together with a constant to a histogram and return the fit function
  TF1* FitDoubleGauss(TH1D* fittedHistogram, double fitRange, double normalizationRange);  // Fit a narrow and wide Gauss functions to a histogram and return the fit function
  void SetBinBoundaries(const int nBins, double *binBorders, int& copyNbins, double *copyBinBorders[]); // Setter for bin boundaries
  bool CheckBinBoundaries(const int nCheckedBins, const double *checkedBins, const TAxis *originalAxis); // Checker that new bin boundaries correspond to old ones
  int CheckNormalizationSanity(const int normalizationType, const int maxIndex); // Sanity check for input normalizations
  
};

#endif
