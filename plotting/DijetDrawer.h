#ifndef DIJETDRAWER_H
#define DIJETDRAWER_H

// C++ includes
#include <bitset>

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

// Own includes
#include "JDrawer.h"
#include "DijetHistogramManager.h"
#include "stackHist.h"   // Xiao's histogram stacking class

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class DijetDrawer {
  
public:
  
  // Indices for different background drawing styles
  enum enumBackgroundStyles{kDrawOverlap,kOverlapZoom,kDrawFit,kDrawFitComposition,knBackgroundStyles};
  
  DijetDrawer(DijetHistogramManager *inputHistograms);  // Constructor
  ~DijetDrawer();                 // Destructor
  
  void DrawHistograms();          // Draw the histograms
  void DrawJetShapeStack();       // Draw stack figures combining all pT bins for jet shape histograms
  
  // Setters for event information and dijets
  void SetDrawEventInformation(const bool drawOrNot); // Setter for drawing event information
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
  
  // Setters for drawing inclusive jet-track correlation
  void SetDrawTrackInclusiveJetCorrelations(const bool drawOrNot);            // Setter for drawing inclusive jet-track correlations
  void SetDrawTrackInclusiveJetCorrelationsPtWeighted(const bool drawOrNot);  // Setter for drawing pT weighted inclusive jet-track correlations
  void SetDrawAllTrackInclusiveJetCorrelations(const bool drawInclusive, const bool drawPtWeighted); // Setter for drawing all inclusive jet-track correlations
  
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
  void SetDrawSameMixedDeltaEtaRatio(const bool drawOrNot); // Setter for drawing same and mixed event ratio for deltaEta plots in the UE region
  
  // Setters for figure saving and logarithmic axes
  void SetSaveFigures(const bool saveOrNot, const char *format, const TString suffix);  // Setter for saving the figures to a file
  void SetLogPt(const bool isLog);          // Setter for logarithmic pT axis
  void SetLogCorrelation(const bool isLog); // Setter for logarithmic z axis for correlation plots
  void SetLogJetShape(const bool isLog);    // Setter for logarithmic jet shape drawing
  void SetLogAxes(const bool pt, const bool correlation, const bool jetShape); // Setter for logarithmic axes
  void SetNormalizeJetShape(const bool normalization);  // Setter for normalization of jet shape
  
  // Setters for drawing style and colors
  void SetColorPalette(const int color);     // Setter for color palette
  void SetDrawingStyle2D(const char* style); // Setter for 2D drawing style
  void SetDrawingStyle3D(const char* style); // Setter for 3D drawing style
  void SetDrawingStyles(const int color, const char* style2D, const char* style3D); // Setter for drawing styles
  void SetBackgroundDrawStyle(const int style); // Setter for background deltaPhi draw styles
  
  
private:
  
  // Data members
  DijetHistogramManager *fHistograms; // Manager for all the drawn histograms
  TString fSystemAndEnergy;           // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy;    // Same a before but without white spaces and dots
  TString fFigureSaveNameAppend;      // Text that can be appended to standard figure naming scheme
  JDrawer *fDrawer;                   // JDrawer for drawing the histograms
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
  bool fDrawEventInformation;                                                     // Draw the event information histograms
  bool fDrawDijetHistograms;                                                      // Draw the dijet histograms
  bool fDrawSingleJets[DijetHistogramManager::knSingleJetCategories];             // Draw the single jet histograms
  bool fDrawTracks[DijetHistogramManager::knTrackCategories];                     // Draw the track histograms
  bool fDrawJetTrackCorrelations[DijetHistogramManager::knJetTrackCorrelations];  // Draw the jet-track correlation histograms
  
  // ==============================================
  // ============== Drawing settings ==============
  // ==============================================
  
  // Flags for drawing individual histograms within histogram groups
  bool fDrawCorrelationType[DijetHistogramManager::knCorrelationTypes];  // Draw different event correalation types
  bool fDrawSameMixedDeltaEtaRatio;                                      // Draw the ratio of same event between peaks and mixed event
  bool fDrawJetTrackDeltaPhi;                                            // Draw the jet-track deltaPhi correlation
  bool fDrawJetTrackDeltaEta;                                            // Draw the jet-track deltaEta correlation
  bool fDrawJetTrackDeltaEtaDeltaPhi;                                    // Draw the jet-track deltaEta-deltaPhi correlation
  bool fDrawJetShape[DijetHistogramManager::knJetShapeTypes];            // Draw the jet shape histograms
  bool fBackgroundDrawStyle[knBackgroundStyles];                         // Bacroung drawing style settings
  
  // Choose if you want to write the figures to pdf file
  bool fSaveFigures;           // Flag for saving the figures to file
  const char* fFigureFormat;   // Format in which the figures are saved
  
  // Normalization settings
  bool fNormalizeJetShape;     // Normalize jet shape to one, meaning draw rho instead of P
  
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
  
  // Methods for drawing
  void DrawEventInformation();    // Draw the event information histograms
  void DrawDijetHistograms();     // Draw the dijet histograms
  void DrawSingleJetHistograms(); // Draw single jet histograms
  void DrawTrackHistograms();     // Draw track histograms
  void DrawJetTrackCorrelationHistograms(); // Draw jet-track correlation histograms
  void DrawJetShapeHistograms();  // Draw jet shape histograms
  void SetupLegend(TLegend *legend, TString centralityString = "", TString trackString = ""); // Common legend style setup for figures
  void SaveFigure(TString figureName, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""); // Save the figure from current canvas to file
  
};

#endif
