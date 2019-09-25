#ifndef DIJETCOMPARINGDRAWER_H
#define DIJETCOMPARINGDRAWER_H

// C++ includes
#include <iostream>
#include <tuple>      // For returning several arguments in a transparent manner

// Root includes
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

// Own includes
#include "DijetCard.h"
#include "JDrawer.h"
#include "DijetMethods.h"
#include "DijetHistogramManager.h"

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class DijetComparingDrawer {
  
private:
  static const int knMaxRatios = 10;  // Maximum number or ratio plots per one canvas
  int fColors[knMaxRatios] = {kRed,kBlue,kMagenta,kCyan,kGreen+4,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};
  
  void SetHistogramStyle(TH1 *histogram, double rangeX, const char* xTitle); // Style setting for histograms in big canvas
  
public:
  
  DijetComparingDrawer(DijetHistogramManager *baseHistograms);  // Constructor
  ~DijetComparingDrawer();                 // Destructor
  
  void DrawHistograms();          // Draw the histograms
  
  // Add histograms to draw together with base histograms
  void AddHistogramToDraw(DijetHistogramManager *additionalHistogram);
  void AddLegendComment(TString comment);
  
  // Setter for dijets
  void SetDrawDijetHistograms(const bool drawOrNot);  // Setter for drawing dijet histograms
  
  // Setters for single jets
  void SetDrawLeadingJetHistograms(const bool drawOrNot);    // Setter for drawing leading jet histograms
  void SetDrawSubleadingJetHistograms(const bool drawOrNot); // Setter for drawing subleading jet histograms
  void SetDrawAnyJetHistograms(const bool drawOrNot);        // Setter for drawing all jet histograms
  void SetDrawAnyLeadingJetHistograms(const bool drawOrNot); // Setter for drawing all leading jet histograms
  void SetDrawAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny, const bool drawAnyLeading);   // Setter for drawing jet histograms
  
  // Setters for tracks
  void SetDrawTracks(const bool drawOrNot);            // Setter for drawing tracks
  void SetDrawTracksUncorrected(const bool drawOrNot); // Setter for drawing uncorrected tracks
  void SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected); // Setter for drawing all track histograms
  void SetDrawInclusiveTracks(const bool drawOrNot);            // Setter for drawing tracks
  void SetDrawInclusiveTracksUncorrected(const bool drawOrNot); // Setter for drawing uncorrected tracks
  void SetDrawAllInclusiveTracks(const bool drawTracks, const bool drawUncorrected); // Setter for drawing all track histograms
  
  // Setters for leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelations(const bool drawOrNot);            // Setter for drawing leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelationsUncorrected(const bool drawOrNot); // Setter for drawing uncorrected leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelationsPtWeighted(const bool drawOrNot);  // Setter for drawing pT weighted leading jet-track correlations
  void SetDrawAllTrackLeadingJetCorrelations(const bool drawLeading, const bool drawUncorrected, const bool drawPtWeighted); // Setter for drawing all correlations related to tracks and leading jets
  
  // Setters for drawing subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelations(const bool drawOrNot);            // Setter for drawing subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelationsUncorrected(const bool drawOrNot); // Setter for drawing uncorrected subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelationsPtWeighted(const bool drawOrNot);  // Setter for drawing pT weighted subleading jet-track correlations
  void SetDrawAllTrackSubleadingJetCorrelations(const bool drawSubleading, const bool drawUncorrected, const bool drawPtWeighted); // Setter for drawing all correlations related to tracks and subleading jets
  
  // Setters for drawing inclusive jet-track correlations
  void SetDrawTrackInclusiveJetCorrelations(const bool drawOrNot);            // Setter for drawing inclusive jet-track correlations
  void SetDrawTrackInclusiveJetCorrelationsPtWeighted(const bool drawOrNot);  // Setter for drawing pT weighted inclusive jet-track correlations
  void SetDrawAllTrackInclusiveJetCorrelations(const bool drawInclusive, const bool drawPtWeighted); // Setter for drawing all correlations related to tracks and inclusive jets
  
  // Setters for drawing different histograms derived from jet-track correlations
  void SetDrawJetTrackDeltaPhi(const bool drawOrNot);                 // Setter for drawing jet-track deltaPhi correlations
  void SetDrawJetTrackDeltaEta(const bool drawOrNot);                 // Setter for drawing jet-track deltaEta correlations
  void SetDrawJetTrackDeltaEtaDeltaPhi(const bool drawOrNot);         // Setter for drawing jet-track deltaEta-deltaPhi correlations
  void SetDrawJetTrackDeltas(const bool deltaPhi, const bool deltaEta, const bool deltaEtaDeltaPhi); // Setter for drawing all the jet-track deltaEta/Phi correlations
  void SetDrawJetShape(const bool drawOrNot);                         // Setter for drawing jet shapes
  void SetDrawJetShapeCounts(const bool drawOrNot);                   // Setter for drawing jet shape counts
  void SetDrawAllJetShapes(const bool jetShape, const bool counts);   // Setter for drawing all different jet shape histograms
  void SetDrawJetShapeBinMap(const bool drawOrNot);                   // Setter for drawing bin mapping between Rbins and deltaEta-deltaPhi bins
  
  // Setters for drawing different correlation types (same event, mixed event, corrected)
  void SetDrawSameEvent(const bool drawOrNot);              // Setter for drawing same event correlation distributions
  void SetDrawMixedEvent(const bool drawOrNot);             // Setter for drawing mixed event correlation distributions
  void SetDrawNormalizedMixedEvent(const bool drawOrNot);   // Setter for drawing normalized mixed event correlation distributions
  void SetDrawCorrectedCorrelations(const bool drawOrNot);  // Setter for drawing corrected correlation distributions
  void SetDrawCorrelationTypes(const bool sameEvent, const bool mixedEvent, const bool normalizedMixedEvent, const bool corrected); // Setter for drawing different correlation types
  void SetDrawBackgroundSubtracted(const bool drawOrNot);   // Setter for drawing background subtracted jet-track correlation histograms
  void SetDrawBackground(const bool drawOrNot);             // Setter for drawing the generated background distributions
  
  // Setter for drawing the event mixing check
  void SetDrawEventMixingCheck(const bool drawOrNot, const bool zoom, const int distributionType);  // Setter for drawing the event mixing check
  
  // Setters for figure saving and logarithmic axes
  void SetSaveFigures(const bool saveOrNot, const char *format, const char *comment = "");  // Setter for saving the figures to a file
  void SetLogPt(const bool isLog);          // Setter for logarithmic pT axis
  void SetLogCorrelation(const bool isLog); // Setter for logarithmic z axis for correlation plots
  void SetLogJetShape(const bool isLog);    // Setter for logarithmic jet shape drawing
  void SetLogAxes(const bool pt, const bool correlation, const bool jetShape); // Setter for logarithmic axes
  
  // Setter for scaling and rebinning the histograms
  void SetApplyScaling(const bool applyScaling); // Set if we should scale the histograms with their integral before comparing them
  void SetJetPtRebin(const bool doRebin);  // Tell if we want to rebin the jet pT histograms
  
  // Setters for ratio plots
  void SetUseDifferenceInRatioPlot(const bool useDifference);  // Setter for plotting difference instead of ratio to lower pad
  void SetRatioZoomMin(const double minValue);  // Setter for minimum value of y-axis in ratio plots
  void SetRatioZoomMax(const double maxValue);  // Setter for maximum value of y-axis in ratio plots
  void SetRatioZoom(const double minValue, const double maxValue);  // Setter for y-axis values in ratio plots
  void SetRatioLabel(TString label);        // Setter for y-axis label in ratio plots
  
  // Setters for drawing style and colors
  void SetColorPalette(const int color);     // Setter for color palette
  void SetDrawingStyle2D(const char* style); // Setter for 2D drawing style
  void SetDrawingStyle3D(const char* style); // Setter for 3D drawing style
  void SetDrawingStyles(const int color, const char* style2D, const char* style3D); // Setter for drawing styles
  
  // Setters for drawing ranges for different bins
  void SetCentralityBinRange(const int first, const int last); // Setter for drawn centrality bins
  void SetTrackPtBinRange(const int first, const int last);    // Setter for drawn track pT bins
  void SetAsymmetryBin(const int asymmetry);                   // Setter for the selected asymmetry bin
  
private:
  
  // Data members
  JDrawer *fDrawer;                       // JDrawer for drawing the histograms
  DijetHistogramManager *fBaseHistograms; // Histograms with respect to which ratios are takes
  DijetHistogramManager *fAddedHistograms[knMaxRatios];  // Histograms drawn together with the base histogram
  int fnAddedHistograms;                  // Number of histograms added for drawing
  
  // ==============================================================
  // == Helper variables to store histograms before drawing them ==
  // ==============================================================
  TH1D *fMainHistogram;
  TH1D *fComparisonHistogram[knMaxRatios];
  TH1D *fRatioHistogram[knMaxRatios];
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
  bool fDrawDijets;                                                               // Draw the dijet histograms
  bool fDrawSingleJets[DijetHistogramManager::knSingleJetCategories];             // Draw the single jet histograms
  bool fDrawTracks[DijetHistogramManager::knTrackCategories];                     // Draw the track histograms
  bool fDrawJetTrackCorrelations[DijetHistogramManager::knJetTrackCorrelations];  // Draw the jet-track correlation histograms
  
  // ==============================================
  // ============== Drawing settings ==============
  // ==============================================
  
  // Flags for drawing individual histograms within histogram groups
  bool fDrawCorrelationType[DijetHistogramManager::knCorrelationTypes];  // Draw different event correalation types
  bool fDrawJetTrackDeltaPhi;                     // Draw the jet-track deltaPhi correlation
  bool fDrawJetTrackDeltaEta;                     // Draw the jet-track deltaEta correlation
  bool fDrawJetTrackDeltaEtaDeltaPhi;             // Draw the jet-track deltaEta-deltaPhi correlation
  bool fDrawJetShape[DijetHistogramManager::knJetShapeTypes];            // Draw the jet shape histograms
  bool fDrawEventMixingCheck;                     // Draw the event mixing check histograms
  
  // Choose if you want to write the figures to pdf file
  bool fSaveFigures;           // Flag for saving the figures to file
  const char* fFigureFormat;   // Format in which the figures are saved
  const char* fFigureComment;  // Comment added to figure name
  
  // Choose if we should scale the histograms before comparing them
  bool fApplyScaling;
  double fScalingFactors[knMaxRatios+1];
  
  // Logarithmic scales for figures
  bool fLogPt;          // pT distributions
  bool fLogCorrelation; // track-jet deltaPhi-deltaEta distributions
  bool fLogJetShape;    // Jet shape distributions
  
  // Zooming for ratio plots
  bool fUseDifferenceInsteadOfRatio;  // Instead of ratio, draw the difference of the two distributions
  double fRatioZoomMin;               // Lower y-axis boundary in ratio plots
  double fRatioZoomMax;               // Upper y-axis boundary in ratio plots
  TString fRatioLabel;                // Label given to ratio plots y-axes
  bool fEventMixingZoom;              // Zoom closer to background region for event mixing check figures
  int fEventMixingDistribution;       // Choose which type of distribution is used in the event mixing check
  
  // Plotting style for 2D and 3D plots
  int fColorPalette;      // Used color palatte for drawing
  const char* fStyle2D;   // Used option for two-dimensional drawing style
  const char* fStyle3D;   // Used option for three-dimensional drawing style
  
  // Rebinning
  bool fRebinJetPt;       // Rebin the single jet pT distributions
  
  // Drawn centrality bins
  int fFirstDrawnCentralityBin;  // First centrality bin that is drawn
  int fLastDrawnCentralityBin;   // Last centrality bin that is drawn
  int fFirstDrawnTrackPtBin;     // First track pT bin that is drawn
  int fLastDrawnTrackPtBin;      // Last track pT bin that is drawn
  int fAsymmetryBin;             // Asymmetry bin selected to be drawn
  
  // Comments given for legends
  TString fLegendComment[knMaxRatios+1];
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int first, int last);  // Sanity check for binning
  
  // Methods for drawing
  void DrawDijetHistograms();     // Draw dijet histograms
  void DrawSingleJetHistograms(); // Draw single jet histograms
  void DrawTrackHistograms();     // Draw track histograms
  void DrawJetTrackCorrelationHistograms(); // Draw jet-track correlation histograms
  void DrawJetShapeHistograms();  // Draw jet shape histograms
  void DrawJetShapeMCComparison();  // Draw jet shape Monte Carlo comparison histograms
  void DrawEventMixingCheck();  // Draw different deltaEta regions in deltaPhi histograms to the same figure to ensure event mixing gives good background
  void SetupLegend(TLegend *legend, TString centralityString = "", TString trackString = "", TString asymmetryString = ""); // Common legend style setup for figures
  void SaveFigure(TString figureName, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""); // Save the figure from current canvas to file
  void PrepareRatio(TString name, int rebin, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0, int bin6 = 0, double minRange = 1000, double maxRange = -1000); // Prepare the ratio histograms out of input histograms
  void DrawToUpperPad(const char* xTitle, const char* yTitle, bool logAxis = false); // Draw the histograms to the same figure in the upper pad of JDrawer
  void DrawToLowerPad(const char* xTitle, const char* yTitle, const double zoomMin, const double zoomMax); // Draw the ratios to the lower pad of the JDrawer
  void ZoomToRegion(const double maxZoomValue, const int nZoomBins, const double scaleFactor, const bool bothSides, const bool asymmetricZoom);  // Zoom the y-axis scale to the specified region of the distribution
  std::tuple<double,double> GetHistogramAverageAndDifferenceInRegion(TH1D *histogram, const double maxZoomValue, const int nZoomBins, const bool bothSides); // Get the average and absolute difference from a histogram in the specific area
  
  // Find the per jet scaling factor
  void FindScalingFactors(const char* histogramName, int iJetCategory, int iCentrality, int iAsymmetry);
  
  
};

#endif
