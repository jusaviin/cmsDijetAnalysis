#ifndef DIJETHISTOGRAMMANAGER_H
#define DIJETHISTOGRAMMANAGER_H

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
class DijetHistogramManager {
 
public:
  
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
  
private:
  // Dimensions for histogram arrays
  static const int knCentralityBins = 4;   // Number of centrality bins
  static const int knTrackPtBins = 6;      // Number of track pT bins
  
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
  
  DijetHistogramManager(TFile *inputFile);                // Constructor
  DijetHistogramManager(const DijetHistogramManager& in); // Copy constructor
  ~DijetHistogramManager();                               // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void ProcessHistograms();       // Do the mixed event correction, subtract the background and calculate jet shape
  
  // Setters for binning information
  void SetCentralityBins(double *binBorders); // Set up centrality bin indices according to provided bin borders
  void SetTrackPtBins(double *binBorders);    // Set up track pT bin indices according to provided bin borders
  void SetDeltaPhiBins(double *lowBinBorders, double *highBinBorders, TString deltaPhiStrings[knDeltaPhiBins], TString compactDeltaPhiStrings[knDeltaPhiBins]); //  Set up deltaPhi bin indices according to provided bin borders and bin names
  
  // Setters for event information and dijets
  void SetLoadEventInformation(const bool loadOrNot); // Setter for drawing event information
  void SetLoadDijetHistograms(const bool loadOrNot);  // Setter for drawing dijet histograms
  
  // Setters for single jets
  void SetLoadLeadingJetHistograms(const bool loadOrNot);    // Setter for drawing leading jet histograms
  void SetLoadSubleadingJetHistograms(const bool loadOrNot); // Setter for drawing subleading jet histograms
  void SetLoadAnyJetHistograms(const bool loadOrNot);        // Setter for drawing all jet histograms
  void SetLoadAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny);   // Setter for drawing jet histograms
  
  // Setters for tracks
  void SetLoadTracks(const bool loadOrNot);            // Setter for drawing tracks
  void SetLoadTracksUncorrected(const bool loadOrNot); // Setter for drawing uncorrected tracks
  void SetLoadAllTracks(const bool drawTracks, const bool drawUncorrected); // Setter for drawing all track histograms
  
  // Setters for leading jet-track correlations
  void SetLoadTrackLeadingJetCorrelations(const bool loadOrNot);            // Setter for drawing leading jet-track correlations
  void SetLoadTrackLeadingJetCorrelationsUncorrected(const bool loadOrNot); // Setter for drawing uncorrected leading jet-track correlations
  void SetLoadTrackLeadingJetCorrelationsPtWeighted(const bool loadOrNot);  // Setter for drawing pT weighted leading jet-track correlations
  void SetLoadAllTrackLeadingJetCorrelations(const bool drawLeading, const bool drawUncorrected, const bool drawPtWeighted); // Setter for drawing all correlations related to tracks and leading jets
  
  // Setters for drawing subleading jet-track correlations
  void SetLoadTrackSubleadingJetCorrelations(const bool loadOrNot);            // Setter for drawing subleading jet-track correlations
  void SetLoadTrackSubleadingJetCorrelationsUncorrected(const bool loadOrNot); // Setter for drawing uncorrected subleading jet-track correlations
  void SetLoadTrackSubleadingJetCorrelationsPtWeighted(const bool loadOrNot);  // Setter for drawing pT weighted subleading jet-track correlations
  void SetLoadAllTrackSubleadingJetCorrelations(const bool drawSubleading, const bool drawUncorrected, const bool drawPtWeighted); // Setter for drawing all correlations related to tracks and subleading jets
  
  // Setters for drawing ranges for different bins
  void SetCentralityBinRange(const int first, const int last); // Setter for drawn centrality bins
  void SetTrackPtBinRange(const int first, const int last);    // Setter for drawn track pT bins
  
  // Setter for used DijetMethods
  void SetDijetMethods(DijetMethods* newMethods); // Setter for used DijetMethods
  
  // Getters for number of bins in histograms
  int GetNCentralityBins() const; // Getter for the number of centrality bins
  int GetNTrackPtBins() const; // Getter for the number of track pT bins
  double GetCentralityBinBorder(const int iCentrality) const;  // Getter for i:th centrality bin border
  double GetTrackPtBinBorder(const int iTrackPt) const;        // Getter for i:th track pT bin border
  
  // Getter for histogram and axis naming
  TString GetCorrelationTypeString(int iCorrelationType) const;        // Getter for correlation type string
  TString GetCompactCorrelationTypeString(int iCorrelationType) const; // Getter for compact correlation type string
  
  TString GetDeltaPhiString(int iDeltaPhiRegion) const;        // Getter for deltaPhi string
  TString GetCompactDeltaPhiString(int iDeltaPhiRegion) const; // Getter for compact deltaPhi string
  
  const char* GetJetTrackHistogramName(int iJetTrackCorrelation) const; // Getter for jet-track correlation histogram name
  const char* GetJetTrackAxisName(int iJetTrackCorrelation) const;      // Getter for name suitable for x-axis in a given jet-track correlation histogram
  
  const char* GetTrackHistogramName(int iTrackType) const; // Getter for track histogram name
  const char* GetTrackAxisName(int iTrackType) const;      // Getter for name suitable for x-axis in a given track histogram
  
  const char* GetSingleJetHistogramName(int iJetType) const; // Getter for single jet histogram name
  const char* GetSingleJetAxisName(int iJetType) const;      // Getter for name suitable for x-axis in a given single jet histogram
  
  const char* GetJetShapeHistogramName(int iJetShapeType) const; // Getter for jet shape histogram name
  const char* GetJetShapeAxisName(int iJetShapeType) const;      // Getter for name suitable for x-axis in a given jet shape histogram
  
  TString GetSystem() const;  // Getter for collision system
  
  // Getters for event information histograms
  TH1D* GetHistogramVertexZ() const;         // Getter for z-vertex histogram
  TH1D* GetHistogramEvents() const;          // Getter for histogram for number of events surviving different event cuts
  TH1D* GetHistogramTrackCuts() const;       // Getter for histogram for number of tracks surviving different track cuts
  TH1D* GetHistogramCentrality() const;      // Getter for centrality histogram in all events
  TH1D* GetHistogramCentralityDijet() const; // Getter for centrality histogram in dijet events
  
  // Getters for single jet histograms
  TH1D* GetHistogramJetPt(const int iJetType, const int iCentrality) const;      // Jet pT histograms
  TH1D* GetHistogramJetPhi(const int iJetType, const int iCentrality) const;     // Jet phi histograms
  TH1D* GetHistogramJetEta(const int iJetType, const int iCentrality) const;     // Jet eta histograms
  TH2D* GetHistogramJetEtaPhi(const int iJetType, const int iCentrality) const;  // 2D eta-phi histogram for jets
  
  // Getters for dijet histograms
  TH1D* GetHistogramDijetDeltaPhi(const int iCentrality) const;              // Dijet deltaPhi histograms
  TH1D* GetHistogramDijetAsymmetry(const int iCentrality) const;             // Dijet asymmetry histograms
  TH2D* GetHistogramDijetLeadingVsSubleadingPt(const int iCentrality) const; // Leading versus subleading jet pT 2D histograms
  
  // Getters for histograms for tracks in dijet events
  TH1D* GetHistogramTrackPt(const int iTrackType, const int iCorrelationType, const int iCentrality) const;                      // Track pT histograms
  TH1D* GetHistogramTrackPhi(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const;    // Track phi histograms
  TH1D* GetHistogramTrackEta(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const;    // Track eta histograms
  TH2D* GetHistogramTrackEtaPhi(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const; // 2D eta-phi histogram for track
  
  // Getters for track-leading jet correlation histograms
  TH1D* GetHistogramJetTrackDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, const int iCentrality, const int iTrackPt) const;                 // DeltaPhi between jet and track
  TH1D* GetHistogramJetTrackDeltaEta(const int iJetTrackCorrelation, const int iCorrelationType, const int iCentrality, const int iTrackPt, const int iDeltaPhiRegion) const; // DeltaEta between jet and track
  TH2D* GetHistogramJetTrackDeltaEtaDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, const int iCentrality, const int iTrackPt) const;         // DeltaEta and deltaPhi between jet and track
  
  // Getters for jet shape histograms
  TH1D* GetHistogramJetShape(const int iJetShapeType, const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt) const;  // Jet shape histograms
  
  TH1D* GetOneDimensionalHistogram(TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0) const; // Getter for any one-dimensional histogram based on input string
  TH2D* GetTwoDimensionalHistogram(TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0) const; // Getter for any two-dimensional histogram based on input string
  
private:
  
  // Data members
  TFile *fInputFile;               // File from which the histograms are read
  DijetCard *fCard;                // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;        // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy; // Same a before but without white spaces and dots
  DijetMethods *fMethods;          // DijetMethods for processing the loaded histograms
  
  // ==============================================
  // ======== Flags for histograms to load ========
  // ==============================================
  
  bool fLoadEventInformation;                              // Draw the event information histograms
  bool fLoadDijetHistograms;                               // Draw the dijet histograms
  bool fLoadSingleJets[knSingleJetCategories];             // Draw the single jet histograms
  bool fLoadTracks[knTrackCategories];                     // Draw the track histograms
  bool fLoadJetTrackCorrelations[knJetTrackCorrelations];  // Draw the jet-track correlation histograms
  bool fLoad2DHistograms;                                  // Load also two-dimensional (eta,phi) and (deltaEta,deltaPhi) histograms
  
  // ==============================================
  // ======== Ranges of histograms to load ========
  // ==============================================
  
  // Drawn centrality bins
  int fFirstLoadedCentralityBin;  // First centrality bin that is drawn
  int fLastLoadedCentralityBin;   // Last centrality bin that is drawn
  int fFirstLoadedTrackPtBin;     // First track pT bin that is drawn
  int fLastLoadedTrackPtBin;      // Last track pT bin that is drawn
  
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
  void SetLoadJetTrackCorrelations(const bool loadOrNot, const int primaryIndex, const int connectedIndex); // Setter for drawing and loading the jet-track correlation histograms
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
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int first, int last) const; // Sanity check for given binning
  int BinIndexCheck(const int nBins, const int binIndex) const; // Check that given index is in defined range
  
  // Mixed event correction, background subtraction and jet shape calculation
  void DoMixedEventCorrection();  // Apply mixed event correction for jet-track correlation histograms
  void SubtractBackgroundAndCalculateJetShape();  // Subtract the background from the distributions and use these histograms to calculate jet shape
  
};

#endif
