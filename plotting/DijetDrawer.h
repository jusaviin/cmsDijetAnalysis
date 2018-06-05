#ifndef DIJETDRAWER_H
#define DIJETDRAWER_H

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

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class DijetDrawer {
 
private:
  // Dimensions for histogram arrays
  static const int knCentralityBins = 4;   // Number of centrality bins
  static const int knTrackPtBins = 6;      // Number of track pT bins
  
  // Indices for different correlation types
  enum enumCorrelationTypes{kSameEvent,kMixedEvent,kCorrected,kBackgroundSubtracted,kBackground,kBackgroundOverlap,kJetShapeBinMap,knCorrelationTypes};
  
  // Indices for different jet-track correlation categories
  enum enumJetTrackCorrelation {kTrackLeadingJet, kUncorrectedTrackLeadingJet, kPtWeightedTrackLeadingJet, kTrackSubleadingJet, kUncorrectedTrackSubleadingJet, kPtWeightedTrackSubleadingJet, knJetTrackCorrelations};
  
  // Indices for different track histogram categories
  enum enumTrackHistograms{kTrack, kUncorrectedTrack, knTrackCategories};
  
  // Indices for different single jet histogram categories
  enum enumSingleJet{kLeadingJet, kSubleadingJet, kAnyJet, knSingleJetCategories};
  
  // Indices for different deltaPhi bins
  enum enumDeltaPhiBins{kWholePhi, kNearSide, kAwaySide, kBetweenPeaks, knDeltaPhiBins};
  
  // Indices for different jet shape histograms
  enum enumJetShape{kJetShape, kJetShapeBinCount, knJetShapeTypes};
  
  // Naming for different correlation types
  TString fCorrelationTypeString[knCorrelationTypes] = {"Same Event","Mixed Event","Corrected","Background subtracted","Background","Background overlap","Jet Shape Bin Map"};
  TString fCompactCorrelationTypeString[knCorrelationTypes] = {"_SameEvent","_MixedEvent","_Corrected","_BackgroundSubtracted","_Background","_BackgroundOverlap","_JetShapeBinMap"};
  
  // Naming for jet-track correlation histograms
  const char* fJetTrackHistogramNames[knJetTrackCorrelations] = {"trackLeadingJet","trackLeadingJetUncorrected","trackLeadingJetPtWeighted","trackSubleadingJet","trackSubleadingJetUncorrected","trackSubleadingJetPtWeighted"}; // Names that different histograms have in the input file
  const char* fJetTrackAxisNames[knJetTrackCorrelations] = {"Track-LJet","UC Track-LJet","p_{T}w Track-LJet","Track-SJet","UC Track-SJet","p_{T}w Track-SJet"}; // Names attached to the figure axes
  
  // Naming for track histograms
  const char* fTrackHistogramNames[knTrackCategories] = {"track","trackUncorrected"}; // Names that different track histograms have in the input file
  const char* fTrackAxisNames[knTrackCategories] = {"Track","Uncorrected track"}; // Names attached to the figure axes
  
  // Naming for single jet histograms
  const char* fSingleJetHistogramName[knSingleJetCategories] = {"leadingJet","subleadingJet","anyJet"}; // Names that different single jet histograms have in the input file
  const char* fSingleJetAxisNames[knSingleJetCategories] = {"Leading jet","Subleading jet","Any jet"}; // Names attached to the figure axes
  
  // Naming for jet shape histograms
  const char* fJetShapeHistogramName[knJetShapeTypes] = {"JetShape","JetShapeCounts"};
  const char* fJetShapeYAxisNames[knJetShapeTypes] = {"#rho(#DeltaR)","counts"};
  
public:
  
  DijetDrawer(TFile *inputFile);  // Constructor
  ~DijetDrawer();                 // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void ProcessHistograms();       // Do the mixed event correction, subtract the background and calculate jet shape
  void DrawHistograms();          // Draw the histograms
  
  // Setters for binning information
  void SetCentralityBins(double *binBorders); // Set up centrality bin indices according to provided bin borders
  void SetTrackPtBins(double *binBorders);    // Set up track pT bin indices according to provided bin borders
  void SetDeltaPhiBins(double *lowBinBorders, double *highBinBorders, TString deltaPhiStrings[knDeltaPhiBins], TString compactDeltaPhiStrings[knDeltaPhiBins]); //  Set up deltaPhi bin indices according to provided bin borders and bin names
  
  // Setters for event information and dijets
  void SetDrawEventInformation(const bool drawOrNot); // Setter for drawing event information
  void SetDrawDijetHistograms(const bool drawOrNot);  // Setter for drawing dijet histograms
  
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
  void SetDrawSameMixedDeltaEtaRatio(const bool drawOrNot); // Setter for drawing same and mixed event ratio for deltaEta plots in the UE region
  
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
  
  // Setters for mixed event configuration
  void SetMixedEventFitRegion(const double etaRange);  // Setter for deltaEta range used for normalizing the mixed event
  
  // Setter for used DijetMethods
  void SetDijetMethods(DijetMethods* newMethods); // Setter for used DijetMethods
  
private:
  
  // Data members
  TFile *fInputFile;               // File from which the histograms are read
  DijetCard *fCard;                // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;        // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy; // Same a before but without white spaces and dots
  JDrawer *fDrawer;                // JDrawer for drawing the histograms
  DijetMethods *fMethods;          // DijetMethods for processing the loaded histograms
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
  bool fDrawEventInformation;                              // Draw the event information histograms
  bool fDrawDijetHistograms;                               // Draw the dijet histograms
  bool fDrawSingleJets[knSingleJetCategories];             // Draw the single jet histograms
  bool fDrawTracks[knTrackCategories];                     // Draw the track histograms
  bool fDrawJetTrackCorrelations[knJetTrackCorrelations];  // Draw the jet-track correlation histograms
  
  // ==============================================
  // ======== Flags for histograms to load ========
  // ==============================================
  
  bool fLoadJetTrackCorrelations[knJetTrackCorrelations];  // Load the jet-track correlation histograms
  
  // ==============================================
  // ============== Drawing settings ==============
  // ==============================================
  
  // Flags for drawing individual histograms within histogram groups
  bool fDrawCorrelationType[knCorrelationTypes];  // Draw different event correalation types
  bool fDrawSameMixedDeltaEtaRatio;               // Draw the ratio of same event between peaks and mixed event
  bool fDrawJetTrackDeltaPhi;                     // Draw the jet-track deltaPhi correlation
  bool fDrawJetTrackDeltaEta;                     // Draw the jet-track deltaEta correlation
  bool fDrawJetTrackDeltaEtaDeltaPhi;             // Draw the jet-track deltaEta-deltaPhi correlation
  bool fDrawJetShape[knJetShapeTypes];            // Draw the jet shape histograms
  
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
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[knCentralityBins+1];     // Indices for centrality bins in centrality binned histograms
  double fCentralityBinBorders[knCentralityBins+1];  // Centrality bin borders, from which bin indices are obtained
  int fTrackPtBinIndices[knTrackPtBins+1];           // Indices for track pT bins in track pT binned histograms
  int fFineTrackPtBinIndices[knTrackPtBins+1];       // Indices for track pT bins in fine track pT binned histograms
  double fTrackPtBinBorders[knTrackPtBins+1];        // Track pT bin borders, from which bin indices are obtained
  int fLowDeltaPhiBinIndices[knDeltaPhiBins];        // Indices for low bin borders in deltaPhi binned histograms
  int fHighDeltaPhiBinIndices[knDeltaPhiBins];       // Indices for high bin borders in deltaPhi binned histograms
  TString fDeltaPhiString[knDeltaPhiBins];           // Names for different deltaPhi bins
  TString fCompactDeltaPhiString[knDeltaPhiBins];    // Names added to figure names for deltaPhi bins
  
  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================
  
  // Event information histograms
  TH1D *fhVertexZ;         // Vertex z position
  TH1D *fhEvents;          // Number of events surviving different event cuts
  TH1D *fhTrackCuts;       // Number of tracks surviving different track cuts
  TH1D *fhCentrality;      // Centrality of all events
  TH1D *fhCentralityDijet; // Centrality of dijet events
  
  // Histograms for single jets
  TH1D *fhJetPt[knSingleJetCategories][knCentralityBins];      // Jet pT histograms
  TH1D *fhJetPhi[knSingleJetCategories][knCentralityBins];     // Jet phi histograms
  TH1D *fhJetEta[knSingleJetCategories][knCentralityBins];     // Jet eta histograms
  TH2D *fhJetEtaPhi[knSingleJetCategories][knCentralityBins];  // 2D eta-phi histogram for jets
  
  // Histograms for dijets
  TH1D *fhDijetDphi[knCentralityBins];                  // Dijet deltaPhi histograms
  TH1D *fhDijetAsymmetry[knCentralityBins];             // Dijet asymmetry histograms
  TH2D *fhDijetLeadingVsSubleadingPt[knCentralityBins]; // Leading versus subleading jet pT 2D histograms
  
  // Histograms for tracks in dijet events
  TH1D *fhTrackPt[knTrackCategories][knCorrelationTypes][knCentralityBins];                      // Track pT histograms
  TH1D *fhTrackPhi[knTrackCategories][knCorrelationTypes][knCentralityBins][knTrackPtBins+1];    // Track phi histograms
  TH1D *fhTrackEta[knTrackCategories][knCorrelationTypes][knCentralityBins][knTrackPtBins+1];    // Track eta histograms
  TH2D *fhTrackEtaPhi[knTrackCategories][knCorrelationTypes][knCentralityBins][knTrackPtBins+1]; // 2D eta-phi histogram for track
  
  // Histograms for track-leading jet correlations
  TH1D *fhJetTrackDeltaPhi[knJetTrackCorrelations][knCorrelationTypes][knCentralityBins][knTrackPtBins];                 // DeltaPhi between jet and track
  TH1D *fhJetTrackDeltaEta[knJetTrackCorrelations][knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // DeltaEta between jet and track
  TH2D *fhJetTrackDeltaEtaDeltaPhi[knJetTrackCorrelations][knCorrelationTypes][knCentralityBins][knTrackPtBins];         // DeltaEta and deltaPhi between jet and track
  
  // Jet shape histograms
  TH1D *fhJetShape[knJetShapeTypes][knJetTrackCorrelations][knCentralityBins][knTrackPtBins];  // Jet shape histograms
  
  // Private methods
  void SetBinIndices(const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinIndices(const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetDrawJetTrackCorrelations(const bool drawOrNot, const int primaryIndex, const int connectedIndex); // Setter for drawing and loading the jet-track correlation histograms
  int GetConnectedIndex(const int jetTrackIndex) const;
  
  // Finders for histograms with different amount of restrictions
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex); // Extract a histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0); // Extract a histogram using given axis restrictions from THnSparseD
  
  // Loaders for different groups of histograms
  void LoadSingleJetHistograms(); // Loader for single jet histograms
  void LoadDijetHistograms(); // Loader for dijet histograms
  void LoadTrackHistograms(); // Loader for track histograms
  void LoadJetTrackCorrelationHistograms(); // Loader for jet-track correlation histograms
  
  // Methods for drawing
  void DrawEventInformation();    // Draw the event information histograms
  void DrawDijetHistograms();     // Draw the dijet histograms
  void DrawSingleJetHistograms(); // Draw single jet histograms
  void DrawTrackHistograms();     // Draw track histograms
  void DrawJetTrackCorrelationHistograms(); // Draw jet-track correlation histograms
  void DrawJetShapeHistograms();  // Draw jet shape histograms
  void SetupLegend(TLegend *legend, TString centralityString = "", TString trackString = ""); // Common legend style setup for figures
  void SaveFigure(TString figureName, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""); // Save the figure from current canvas to file
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int first, int last); // Sanity check for given binning
  
  // Mixed event correction, background subtraction and jet shape calculation
  void DoMixedEventCorrection();  // Apply mixed event correction for jet-track correlation histograms
  void SubtractBackgroundAndCalculateJetShape();  // Subtract the background from the distributions and use these histograms to calculate jet shape
  
};

#endif
