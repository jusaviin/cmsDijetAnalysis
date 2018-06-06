#ifndef DIJETCOMPARINGDRAWER_H
#define DIJETCOMPARINGDRAWER_H

// C++ includes
#include <iostream>

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
  
public:
  
  DijetComparingDrawer(DijetHistogramManager *baseHistograms);  // Constructor
  ~DijetComparingDrawer();                 // Destructor
  
  void DrawHistograms();          // Draw the histograms
  
  // Add histograms to draw together with base histograms
  void AddHistogramToDraw(DijetHistogramManager *additionalHistogram);
  
  // Setters for single jets
  void SetDrawLeadingJetHistograms(const bool drawOrNot);    // Setter for drawing leading jet histograms
  void SetDrawSubleadingJetHistograms(const bool drawOrNot); // Setter for drawing subleading jet histograms
  void SetDrawAnyJetHistograms(const bool drawOrNot);        // Setter for drawing all jet histograms
  void SetDrawAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny);   // Setter for drawing jet histograms
  
  // Setters for tracks
  void SetDrawTracks(const bool drawOrNot);            // Setter for drawing tracks
  void SetDrawTracksUncorrected(const bool drawOrNot); // Setter for drawing uncorrected tracks
  void SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected); // Setter for drawing all track histograms
  
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
  void SetDrawCorrectedCorrelations(const bool drawOrNot);  // Setter for drawing corrected correlation distributions
  void SetDrawCorrelationTypes(const bool sameEvent, const bool mixedEvent, const bool corrected); // Setter for drawing different correlation types
  void SetDrawBackgroundSubtracted(const bool drawOrNot);   // Setter for drawing background subtracted jet-track correlation histograms
  void SetDrawBackground(const bool drawOrNot);             // Setter for drawing the generated background distributions
  
  // Setters for figure saving and logarithmic axes
  void SetSaveFigures(const bool saveOrNot, const char *format);  // Setter for saving the figures to a file
  void SetLogPt(const bool isLog);          // Setter for logarithmic pT axis
  void SetLogCorrelation(const bool isLog); // Setter for logarithmic z axis for correlation plots
  void SetLogJetShape(const bool isLog);    // Setter for logarithmic jet shape drawing
  void SetLogAxes(const bool pt, const bool correlation, const bool jetShape); // Setter for logarithmic axes
  
  // Setters for drawing style and colors
  void SetColorPalette(const int color);     // Setter for color palette
  void SetDrawingStyle2D(const char* style); // Setter for 2D drawing style
  void SetDrawingStyle3D(const char* style); // Setter for 3D drawing style
  void SetDrawingStyles(const int color, const char* style2D, const char* style3D); // Setter for drawing styles
  
  // Setters for drawing ranges for different bins
  void SetCentralityBinRange(const int first, const int last); // Setter for drawn centrality bins
  void SetTrackPtBinRange(const int first, const int last);    // Setter for drawn track pT bins
  
private:
  
  // Data members
  JDrawer *fDrawer;                       // JDrawer for drawing the histograms
  DijetHistogramManager *fBaseHistograms; // Histograms with respect to which ratios are takes
  DijetHistogramManager *fAddedHistograms[knMaxRatios];  // Histograms drawn together with the base histogram
  int fnAddedHistograms;                  // Number of histograms added for drawing
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
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
  
  // Choose if you want to write the figures to pdf file
  bool fSaveFigures;           // Flag for saving the figures to file
  const char* fFigureFormat;   // Format in which the figures are saved
  
  // Logarithmic scales for figures
  bool fLogPt;          // pT distributions
  bool fLogCorrelation; // track-jet deltaPhi-deltaEta distributions
  bool fLogJetShape;    // Jet shape distributions
  
  // Plotting style for 2D and 3D plots
  int fColorPalette;      // Used color palatte for drawing
  const char* fStyle2D;   // Used option for two-dimensional drawing style
  const char* fStyle3D;   // Used option for three-dimensional drawing style
  
  // Drawn centrality bins
  int fFirstDrawnCentralityBin;  // First centrality bin that is drawn
  int fLastDrawnCentralityBin;   // Last centrality bin that is drawn
  int fFirstDrawnTrackPtBin;     // First track pT bin that is drawn
  int fLastDrawnTrackPtBin;      // Last track pT bin that is drawn
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int first, int last);  // Sanity check for binning
  
  // Methods for drawing
  void DrawSingleJetHistograms(); // Draw single jet histograms
  void DrawTrackHistograms();     // Draw track histograms
  void DrawJetTrackCorrelationHistograms(); // Draw jet-track correlation histograms
  void DrawJetShapeHistograms();  // Draw jet shape histograms
  void SetupLegend(TLegend *legend, TH1D *mainHistogram, TH1D *additionalHistogram[knMaxRatios], TString centralityString = "", TString trackString = ""); // Common legend style setup for figures
  void SaveFigure(TString figureName, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""); // Save the figure from current canvas to file
  TH1D* PrepareRatio(TH1D *mainHistogram, TH1D *additionalHistogram[knMaxRatios], TH1D *hRatio[knMaxRatios], TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0); // Prepare the ratio histograms out of input histograms
  void DrawToUpperPad(TH1D *mainHistogram, TH1D *additionalHistogram[knMaxRatios], const char* xTitle, const char* yTitle, bool logAxis); // Draw the histograms to the same figure in the upper pad of JDrawer
  
  
};

#endif
