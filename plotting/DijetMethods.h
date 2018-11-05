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
  
  TH2D* MixedEventCorrect(TH2D *sameEventHistogram, TH2D *leadingMixedEventHistogram, TH2D *subleadingMixedEventHistogram); // Mixed event correction for a two-dimensional histogram
  TH2D* DoSeagullCorrection(TH2D *mixedEventCorrectedHistogram);  // Apply a seagull correction to the histogram
  TH2D* SubtractBackground(TH2D *leadingHistogramWithBackground, TH2D *subleadingHistogramWithBackground, bool isInclusive = false); // Subtract background from a two-dimensional leading histogram
  TH2D* GetSpilloverCorrection(TH2D *onlyHydjetHistogram);  // Get the spillover correction from only hydjet histogram
  TH1D* GetJetShape(TH2D *backgroundSubtractedHistogram); // Extract the jet shape from the two-dimensional histogram
  TH2D* RebinHistogram(TH2D *histogramInNeedOfRebinning); // Rebin a two-dimensional deltaPhi-deltaEta histogram
  TH1D* ProjectSignalDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram); // Project deltaPhi distribution in the signal region in eta out of a two-dimensional deltaPhi-deltaEta distribution
  TH1D* ProjectBackgroundDeltaPhi(TH2D* deltaPhiDeltaEtaHistogram); // Project deltaPhi distribution in the background region out of a two-dimensional deltaPhi-deltaEta distribution
  
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
  
  // Setters for mixed event configuration
  void SetMixedEventFitRegion(const double etaRange);  // Setter for deltaEta range used for normalizing the mixed event
  void SetMixedEventNormalization(const int normalizationType, const bool smoothenMixing); // Setter for normalization method used for mixed event distributions
  
  // Setters for background subtraction configuration
  void SetBackgroundDeltaEtaRegion(const double minDeltaEta, const double maxDeltaEta); // Setter for background deltaEta region
  void SetSignalDeltaEtaRegion(const double maxDeltaEta);                               // Setter for signal deltaEta region
  
  // Setters for jet shape calculation configuration
  void SetJetShapeBinEdges(const int nBins, double *binBorders); // Setter for R-binning for jet shape histograms
  void SetJetShapeNormalization(const int normalizationType);    // Setter for jet shape normalization method
  
  // Setter for two-dimensional histogram rebinning information
  void SetRebinBoundaries(const int nRebinDeltaEta, double *deltaEtaBorders, const int nRebinDeltaPhi, double *deltaPhiBorders); // Setter for deltaEta and deltaPhi rebin borders
  
  // Setter for parameters for seagull correction
  void SetBackgroundDeltaPhiRegionSeagull(const double minDeltaPhi, const double maxDeltaPhi);  // Setter for deltaPhi region considered as background in seagull correction
  void SetSeagullRebin(const int nRebin);  // Setter for the amount of rebin applied to deltaEta histogram before fit in seagull correction
  
private:
  
  // =============================================
  // =========== Mixed event correction ==========
  // =============================================
  
  TH2D *fNormalizedMixedEventHistogram; // Mixed event histogram after normalization
  double fMixedEventFitRegion;  // Region in deltaEta in which a constant fit is done to normalize mixed event distributions
  int fMixedEventNormalizationMethod; // Normalization method used for mixed event distributions
  bool fSmoothMixing; // Smoothen the mixing distribution in phi
  
  // =============================================
  // =========== Background subtraction ==========
  // =============================================
  
  TH2D *fBackgroundDistribution;  // Remember the background distribution from the most recent background subtraction
  TH2D *fBackgroundOverlap;       // Fill a few points over the gluing point to see how well the gluing works
  double fMinBackgroundDeltaEta;  // Minimum deltaEta for background subtraction region
  double fMaxBackgroundDeltaEta;  // Maximum deltaEta for background subtraction region
  
  // =============================================
  // ============ Seagull correction =============
  // =============================================
  
  TH1D *fBackgroundEtaProjection;  // Projected deltaEta in the background region of deltaPhi
  TF1 *fSeagullFit;                // Function used to do the seagull fit
  double fMinBackgroundDeltaPhi;   // Minimum deltaPhi for seagull fit region
  double fMaxBackgroundDeltaPhi;   // Maximum deltaPhi for seagull fit region
  int fSeagullRebin;               // Rebin applied to deltaEta histogram before fitting
  
  // =============================================
  // ============ DeltaPhi projections ===========
  // =============================================
  double fMaxSignalDeltaEta;      // Maximum deltaEta value accepted for the signal region
  
  // =============================================
  // ============ Spillover correction ===========
  // =============================================
  TH1D *fSpilloverDeltaEta;       // DeltaEta projection for spillover calculation
  TH1D *fSpilloverDeltaPhi;       // DeltaPhi projection for spillover calculation
  TF1 *fSpilloverFitDeltaEta;     // Fit to deltaEta projection during spillover calculation
  TF1 *fSpilloverFitDeltaPhi;     // Fir to deltaPhi projection during spillover calculation
  
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
  
  // Private methods
  double GetMixedEventScale(TH2D* mixedEventHistogram); // Find the normalization scale for the mixed event histogram
  TF1* FitGauss(TH1D* fittedHistogram, double fitRange);  // Fit a Gaussian function to a histogram and return the fit function
  TH1D* ProjectRegionDeltaPhi(const TH2D* deltaPhiDeltaEtaHistogram, const double minDeltaEta, const double maxDeltaEta, const char* newName);  // Project deltaPhi distribution out of a two-dimensional deltaPhi-deltaEta distribution
  TH1D* ProjectRegionDeltaEta(const TH2D* deltaPhiDeltaEtaHistogram, const double minDeltaPhi, const double maxDeltaPhi, const char* newName);  // Project deltaEta distribution out of a two-dimensional deltaPhi-deltaEta distribution
  void SetBinBoundaries(const int nBins, double *binBorders, int& copyNbins, double *copyBinBorders[]); // Setter for bin boundaries
  bool CheckBinBoundaries(const int nCheckedBins, const double *checkedBins, TAxis *originalAxis); // Checker that new bin boundaries correspond to old ones
  int CheckNormalizationSanity(const int normalizationType, const int maxIndex); // Sanity check for input normalizations
  
};

#endif
