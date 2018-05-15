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
  static const int knDeltaPhiBins = 4;     // Number of delta phi slices (whole phi/near side/away side/between peaks)
  
  // Indices for different correlation types
  enum enumCorrelationTypes{kSameEvent,kMixedEvent,kCorrected,kBackgroundSubtracted,knCorrelationTypes};
  
  // Indices for different jet-track correlation categories
  enum enumJetTrackCorrelation {kTrackLeadingJet, kUncorrectedTrackLeadingJet, kPtWeightedTrackLeadingJet, kTrackSubleadingJet, kUncorrectedTrackSubleadingJet, kPtWeightedTrackSubleadingJet, knJetTrackCorrelations};
  
  // Indices for different track histogram categories
  enum enumTrackHistograms{kTrack, kUncorrectedTrack, knTrackCategories};
  
  // Indices for different single jet histogram categories
  enum enumSingleJet{kLeadingJet, kSubleadingJet, kAnyJet, knSingleJetCategories};
  
  // Naming for different correlation types
  TString fCorrelationTypeString[knCorrelationTypes] = {"Same Event","Mixed Event","Corrected","Background subtracted"};
  TString fCompactCorrelationTypeString[knCorrelationTypes] = {"_SameEvent","_MixedEvent","_Corrected","_Background subtracted"};
  
  // Naming for jet-track correlation histograms
  const char* fJetTrackHistogramNames[knJetTrackCorrelations] = {"trackLeadingJet","trackLeadingJetUncorrected","trackLeadingJetPtWeighted","trackSubleadingJet","trackSubleadingJetUncorrected","trackSubleadingJetPtWeighted"}; // Names that different histograms have in the input file
  const char* fJetTrackAxisNames[knJetTrackCorrelations] = {"Track-LJet","UC Track-LJet","p_{T}w Track-LJet","Track-SJet","UC Track-SJet","p_{T}w Track-SJet"}; // Names attached to the figure axes
  
  // Naming for track histograms
  const char* fTrackHistogramNames[knTrackCategories] = {"track","trackUncorrected"}; // Names that different track histograms have in the input file
  const char* fTrackAxisNames[knTrackCategories] = {"Track","Uncorrected track"}; // Names attached to the figure axes
  
  // Naming for single jet histograms
  const char* fSingleJetHistogramName[knSingleJetCategories] = {"leadingJet","subleadingJet","anyJet"}; // Names that different single jet histograms have in the input file
  const char* fSingleJetAxisNames[knSingleJetCategories] = {"Leading jet","Subeading jet","Any jet"}; // Names attached to the figure axes
  
public:
  
  DijetDrawer(TFile *inputFile);  // Constructor
  ~DijetDrawer();                 // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void DoMixedEventCorrection();  // Apply mixed event correction for jet-track correlation histograms
  void DrawHistograms();          // Draw the histograms
  
  // Setters for binning information
  void SetCentralityBins(double *binBorders); // Set up centrality bin indices according to provided bin borders
  void SetTrackPtBins(double *binBorders);    // Set up track pT bin indices according to provided bin borders
  void SetDeltaPhiBins(double *lowBinBorders, double *highBinBorders, TString deltaPhiStrings[knDeltaPhiBins], TString compactDeltaPhiStrings[knDeltaPhiBins]); //  Set up deltaPhi bin indices according to provided bin borders and bin names
  
  // Setters for event information and dijets
  void SetDrawEventInformation(bool drawOrNot); // Setter for drawing event information
  void SetDrawDijetHistograms(bool drawOrNot);  // Setter for drawing dijet histograms
  
  // Setters for single jets
  void SetDrawLeadingJetHistograms(bool drawOrNot);    // Setter for drawing leading jet histograms
  void SetDrawSubleadingJetHistograms(bool drawOrNot); // Setter for drawing subleading jet histograms
  void SetDrawAnyJetHistograms(bool drawOrNot);        // Setter for drawing all jet histograms
  void SetDrawAllJets(bool drawLeading, bool drawSubleading, bool drawAny);   // Setter for drawing jet histograms
  
  // Setters for tracks
  void SetDrawTracks(bool drawOrNot);            // Setter for drawing tracks
  void SetDrawTracksUncorrected(bool drawOrNot); // Setter for drawing uncorrected tracks
  void SetDrawAllTracks(bool drawTracks, bool drawUncorrected); // Setter for drawing all track histograms
  
  // Setters for leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelations(bool drawOrNot);            // Setter for drawing leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelationsUncorrected(bool drawOrNot); // Setter for drawing uncorrected leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelationsPtWeighted(bool drawOrNot);  // Setter for drawing pT weighted leading jet-track correlations
  void SetDrawAllTrackLeadingJetCorrelations(bool drawLeading, bool drawUncorrected, bool drawPtWeighted); // Setter for drawing all correlations related to tracks and leading jets
  
  // Setters for drawing subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelations(bool drawOrNot);            // Setter for drawing subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelationsUncorrected(bool drawOrNot); // Setter for drawing uncorrected subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelationsPtWeighted(bool drawOrNot);  // Setter for drawing pT weighted subleading jet-track correlations
  void SetDrawAllTrackSubleadingJetCorrelations(bool drawSubleading, bool drawUncorrected, bool drawPtWeighted); // Setter for drawing all correlations related to tracks and subleading jets
  
  // Setters for drawing different correlation types (same event, mixed event, corrected)
  void SetDrawSameEvent(bool drawOrNot);             // Setter for drawing same event correlation distributions
  void SetDrawMixedEvent(bool drawOrNot);            // Setter for drawing mixed event correlation distributions
  void SetDrawCorrectedCorrelations(bool drawOrNot); // Setter for drawing corrected correlation distributions
  void SetDrawCorrelationTypes(bool sameEvent, bool mixedEvent, bool corrected); // Setter for drawing different correlation types
  void SetDrawBackgroundSubtracted(bool drawOrNot);  // Setter for drawing background subtracted jet-track correlation histograms
  void SetDrawSameMixedDeltaEtaRatio(bool drawOrNot); // Setter for drawing same and mixed event ratio for deltaEta plots in the UE region
  
  // Setters for figure saving and logarithmic axes
  void SetSaveFigures(bool saveOrNot, const char *format);  // Setter for saving the figures to a file
  void SetLogPt(bool isLog); // Setter for logarithmic pT axis
  void SetLogCorrelation(bool isLog); // Setter for logarithmic z axis for correlation plots
  void SetLogAxes(bool pt, bool correlation); // Setter for logarithmic axes
  
  // Setters for drawing style and colors
  void SetColorPalette(int color); // Setter for color palette
  void SetDrawingStyle2D(const char* style); // Setter for 2D drawing style
  void SetDrawingStyle3D(const char* style); // Setter for 3D drawing style
  void SetDrawingStyles(int color, const char* style2D, const char* style3D); // Setter for drawing styles
  
  // Setters for drawing ranges for different bins
  void SetCentralityBinRange(const int first, const int last); // Setter for drawn centrality bins
  void SetTrackPtBinRange(const int first, const int last);    // Setter for drawn track pT bins
  
  // Setters for mixed event configuration
  void SetMixedEventFitRegion(const double etaRange);  // Setter for deltaEta range used for normalizing the mixed event
  
private:
  
  // Data members
  TFile *fInputFile;  // File from which the histograms are read
  DijetCard *fCard;   // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;  // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy; // Same a before but without white spaces and dots
  JDrawer *fDrawer;   // JDrawer for drawing the histograms
  DijetMethods *fMethods;  // DijetMethods for processing the loaded histograms
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
  bool fDrawEventInformation;
  bool fDrawDijetHistograms;
  bool fDrawSingleJets[knSingleJetCategories];
  bool fDrawTracks[knTrackCategories];
  bool fDrawJetTrackCorrelations[knJetTrackCorrelations];
  
  // ==============================================
  // ============== Drawing settings ==============
  // ==============================================
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool fDrawCorrelationType[knCorrelationTypes];
  bool fDrawSameMixedDeltaEtaRatio;
  
  // Choose if you want to write the figures to pdf file
  bool fSaveFigures;
  const char* fFigureFormat;
  
  // Logarithmic scales for figures for pT distributions
  bool fLogPt;          // pT distributions
  bool fLogCorrelation; // track-jet deltaPhi-deltaEta distributions
  
  // Plotting style for 2D and 3D plots
  int fColorPalette;
  const char* fStyle2D;
  const char* fStyle3D;
  
  // Drawn centrality bins
  int fFirstDrawnCentralityBin;
  int fLastDrawnCentralityBin;
  int fFirstDrawnTrackPtBin;
  int fLastDrawnTrackPtBin;
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[knCentralityBins+1];
  double fCentralityBinBorders[knCentralityBins+1];
  int fTrackPtBinIndices[knTrackPtBins+1];
  int fFineTrackPtBinIndices[knTrackPtBins+1];
  double fTrackPtBinBorders[knTrackPtBins+1];
  int fLowDeltaPhiBinIndices[knDeltaPhiBins];
  int fHighDeltaPhiBinIndices[knDeltaPhiBins];
  TString fDeltaPhiString[knDeltaPhiBins];
  TString fCompactDeltaPhiString[knDeltaPhiBins];
  
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
  
  // Private methods
  void SetBinIndices(const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinIndices(const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis); // Read the bin indices for given bin borders
  
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
  void SetupLegend(TLegend *legend, TString centralityString = "", TString trackString = ""); // Common legend style setup for figures
  void SaveFigure(TString figureName, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""); // Save the figure from current canvas to file
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int first, int last); // Sanity check for given binning
  
};

#endif
