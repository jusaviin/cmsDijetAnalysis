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

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class DijetDrawer {
 
private:
  // Dimensions for histogram arrays
  static const int knCentralityBins = 4;   // Number of centrality bins
  static const int knTrackPtBins = 6;      // Number of track pT bins
  static const int knCorrelationTypes = 3; // Number of correlation type bins (same event/mixed event/corrected same event)
  static const int knDeltaPhiBins = 4;     // Number of delta phi slices (whole phi/near side/away side/between peaks)
  
public:
  
  DijetDrawer(TFile *inputFile);  // Constructor
  ~DijetDrawer();                 // Destructor
  
  void LoadHistograms();         // Load the histograms from the inputfile
  void DoMixedEventCorrection(); // Apply mixed event correction for jet-track correlation histograms
  void DrawHistograms();         // Draw the histograms
  
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
  
private:
  
  // Data members
  TFile *fInputFile;  // File from which the histograms are read
  DijetCard *fCard;   // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;  // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy; // Same a before but without white spaces and dots
  JDrawer *fDrawer;   // JDrawer for drawing the histograms
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
  bool fDrawEventInformation;
  bool fDrawDijetHistograms;
  bool fDrawLeadingJetHistograms;
  bool fDrawSubleadingJetHistograms;
  bool fDrawAnyJetHistograms;
  bool fDrawTracks;
  bool fDrawUncorrectedTracks;
  bool fDrawTrackLeadingJetCorrelations;
  bool fDrawUncorrectedTrackLeadingJetCorrelations;
  bool fDrawPtWeightedTrackLeadingJetCorrelations;
  bool fDrawTrackSubleadingJetCorrelations;
  bool fDrawUncorrectedTrackSubleadingJetCorrelations;
  bool fDrawPtWeightedTrackSubleadingJetCorrelations;
  
  // ==============================================
  // ============== Drawing settings ==============
  // ==============================================
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool fDrawSameEvent;
  bool fDrawMixedEvent;
  bool fDrawCorrected;
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
  double fTrackPtBinBorders[knTrackPtBins+1];
  int fLowDeltaPhiBinIndices[knDeltaPhiBins];
  int fHighDeltaPhiBinIndices[knDeltaPhiBins];
  TString fDeltaPhiString[knDeltaPhiBins];
  TString fCompactDeltaPhiString[knDeltaPhiBins];
  TString fCorrelationTypeString[knCorrelationTypes];
  TString fCompactCorrelationTypeString[knCorrelationTypes];
  
  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================
  
  // Vertex z position
  TH1D *fhVertexZ;
  
  // Number of events surviving different event cuts
  TH1D *fhEvents;
  
  // Number of tracks surviving different track cuts
  TH1D *fhTrackCuts;
  
  // Centrality of all and dijet events
  TH1D *fhCentrality;
  TH1D *fhCentralityDijet;
  
  // Histograms for leading jets
  TH1D *fhLeadingJetPt[knCentralityBins];               // Leading jet pT histograms
  TH1D *fhLeadingJetPhi[knCentralityBins];              // Leading jet phi histograms
  TH1D *fhLeadingJetEta[knCentralityBins];              // Leading jet eta histograms
  TH2D *fhLeadingJetEtaPhi[knCentralityBins];           // 2D eta-phi histogram for leading jet
  
  // Histograms for subleading jets
  TH1D *fhSubleadingJetPt[knCentralityBins];            // Subleading jet pT histograms
  TH1D *fhSubleadingJetPhi[knCentralityBins];           // Subleading jet phi histograms
  TH1D *fhSubleadingJetEta[knCentralityBins];           // Subleading jet eta histograms
  TH2D *fhSubleadingJetEtaPhi[knCentralityBins];        // 2D eta-phi histogram for subleading jet
  
  // Histograms for all jets
  TH1D *fhAnyJetPt[knCentralityBins] ;                  // Any jet pT histograms
  TH1D *fhAnyJetPhi[knCentralityBins];                  // Any jet phi histograms
  TH1D *fhAnyJetEta[knCentralityBins];                  // Any jet eta histograms
  TH2D *fhAnyJetEtaPhi[knCentralityBins];               // 2D eta-phi histogram for all jets
  
  // Histograms for dijets
  TH1D *fhDijetDphi[knCentralityBins];                  // Dijet deltaPhi histograms
  TH1D *fhDijetAsymmetry[knCentralityBins];             // Dijet asymmetry histograms
  TH2D *fhDijetLeadingVsSubleadingPt[knCentralityBins]; // Leading versus subleading jet pT 2D histograms
  
  // Histograms for tracks in dijet events
  TH1D *fhTrackPt[knCorrelationTypes][knCentralityBins];                    // Track pT histograms
  TH1D *fhTrackPhi[knCorrelationTypes][knCentralityBins];                   // Track phi histograms
  TH1D *fhTrackEta[knCorrelationTypes][knCentralityBins];                   // Track eta histograms
  TH2D *fhTrackEtaPhi[knCorrelationTypes][knCentralityBins];                // 2D eta-phi histogram for track
  
  // Histograms for uncorrected tracks in dijet events
  TH1D *fhTrackPtUncorrected[knCorrelationTypes][knCentralityBins];         // Uncorrected track pT histograms
  TH1D *fhTrackPhiUncorrected[knCorrelationTypes][knCentralityBins];        // Uncorrected track phi histograms
  TH1D *fhTrackEtaUncorrected[knCorrelationTypes][knCentralityBins];        // Uncorrected track eta histograms
  TH2D *fhTrackEtaPhiUncorrected[knCorrelationTypes][knCentralityBins];     // 2D eta-phi histogram for uncorrected tracks
  
  // Histograms for track-leading jet correlations
  TH1D *fhTrackLeadingJetDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // DeltaPhi between track and leading jet
  TH1D *fhTrackLeadingJetDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // DeltaEta between track and leading jet
  TH2D *fhTrackLeadingJetDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for uncorrected track-leading jet correlations
  TH1D *fhTrackLeadingJetDeltaPhiUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // Uncorrected deltaPhi between track and leading jet
  TH1D *fhTrackLeadingJetDeltaEtaUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // Uncorrected deltaEta between track and leading jet
  TH2D *fhTrackLeadingJetDeltaEtaDeltaPhiUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // Uncorrected deltaEta and deltaPhi between track and leading jet
  
  // Histograms for pT weighted track-leading jet correlations
  TH1D *fhTrackLeadingJetDeltaPhiPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // pT weighted deltaPhi between track and leading jet
  TH1D *fhTrackLeadingJetDeltaEtaPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // pT weighted deltaEta between track and leading jet
  TH2D *fhTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // pT weighted deltaEta and deltaPhi between track and leading jet
  
  // Histograms for track-subleading jet correlations
  TH1D *fhTrackSubleadingJetDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // DeltaPhi between track and subleading jet
  TH1D *fhTrackSubleadingJetDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // DeltaEta between track and subleading jet
  TH2D *fhTrackSubleadingJetDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // DeltaEta and deltaPhi between track and subleading jet
  
  // Histograms for uncorrected track-subleading jet correlations
  TH1D *fhTrackSubleadingJetDeltaPhiUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // Uncorrected deltaPhi between track and subleading jet
  TH1D *fhTrackSubleadingJetDeltaEtaUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // Uncorrected deltaEta between track and subleading jet
  TH2D *fhTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // Uncorrected deltaEta and deltaPhi between track and subleading jet
  
  // Histograms for pT weighted track-subleading jet correlations
  TH1D *fhTrackSubleadingJetDeltaPhiPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // pT weighted deltaPhi between track and subleading jet
  TH1D *fhTrackSubleadingJetDeltaEtaPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // pT weighted deltaEta between track and subleading jet
  TH2D *fhTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // pT weighted deltaEta and deltaPhi between track and subleading jet
  
  // Private methods
  void SetBinIndices(const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinIndices(const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis); // Read the bin indices for given bin borders
  
  // Finders for histograms with different amount of restrictions
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex); // Extract a histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0); // Extract a histogram using given axis restrictions from THnSparseD
  
  // Loaders for different groups of histograms
  void LoadSingleJetHistograms(TH1D *hJetPt[knCentralityBins], TH1D* hJetPhi[knCentralityBins], TH1D* hJetEta[knCentralityBins], TH2D* hJetEtaPhi[knCentralityBins], const char* name, const int iCentralityAxis); // Loader for single jet histograms
  void LoadDijetHistograms(TH1D *hDeltaPhi[knCentralityBins], TH1D* hAsymmetry[knCentralityBins], TH2D* hLeadingSubleadingPt[knCentralityBins], const char* name); // Loader for dijet histograms
  void LoadTrackHistograms(TH1D *hTrackPt[knCorrelationTypes][knCentralityBins], TH1D *hTrackPhi[knCorrelationTypes][knCentralityBins], TH1D *hTrackEta[knCorrelationTypes][knCentralityBins], TH2D *hTrackEtaPhi[knCorrelationTypes][knCentralityBins], const char* name); // Loader for track histograms
  void LoadJetTrackCorrelationHistograms(TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins], TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], const char* name); // Loader for jet-track correlation histograms
  
  // Methods for drawing
  void DrawEventInformation(); // Draw the event information histograms
  void DrawDijetHistograms();  // Draw the dijet histograms
  void DrawSingleJetHistograms(TH1D *hJetPt[knCentralityBins], TH1D* hJetPhi[knCentralityBins], TH1D* hJetEta[knCentralityBins], TH2D* hJetEtaPhi[knCentralityBins], const char* nameForAxis, const char* nameForSave); // Draw single jet histograms
  void DrawTrackHistograms(TH1D *hTrackPt[knCorrelationTypes][knCentralityBins], TH1D *hTrackPhi[knCorrelationTypes][knCentralityBins], TH1D *hTrackEta[knCorrelationTypes][knCentralityBins], TH2D *hTrackEtaPhi[knCorrelationTypes][knCentralityBins], const char* nameForAxis, const char* nameForSave); // Draw track histograms
  void DrawJetTrackCorrelationHistograms(TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins], TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], const char* nameForAxis, const char* nameForSave); // Draw jet-track correlation histograms
  void SetupLegend(TLegend *legend, TString centralityString = "", TString trackString = ""); // Common legend style setup for figures
  void SaveFigure(TString figureName, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""); // Save the figure from current canvas to file
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int first, int last); // Sanity check for given binning
  
  // Methods for mixed event correction
  void ApplyMixedEventCorrection(TH1D *hDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins], TH1D *hDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins], TH2D *hDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins]); // Mixed event correction for jet-track correlation histograms
  TH2D* MixedEventCorrect(TH2D *sameEventHistogram, TH2D *mixedEventHistogram); // Mixed event correction for one two dimensional histogram

  
};

#endif
