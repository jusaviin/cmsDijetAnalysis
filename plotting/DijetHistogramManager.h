#ifndef DIJETHISTOGRAMMANAGER_H
#define DIJETHISTOGRAMMANAGER_H

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

// Own includes
#include "DijetCard.h"
#include "JDrawer.h"
#include "DijetMethods.h"
#include "../src/DijetHistograms.h"

class JffCorrector;  // Need forward declaration of the JffCorrector as this class depends on DijetHistogramManager
class SeagullConfiguration; // Also SeagullConfiguration depends on DijetHistogramManager

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class DijetHistogramManager {
 
public:
  
  // Indices for different correlation types
  enum enumCorrelationTypes{kSameEvent, kMixedEvent, kMixedEventNormalized, kCorrected, kBackgroundSubtracted, kBackground, kBackgroundOverlap, kJetShapeBinMap, knCorrelationTypes};
  
  // Indices for different jet-track correlation categories
  enum enumJetTrackCorrelation {kTrackLeadingJet, kUncorrectedTrackLeadingJet, kPtWeightedTrackLeadingJet, kTrackSubleadingJet, kUncorrectedTrackSubleadingJet, kPtWeightedTrackSubleadingJet, kTrackInclusiveJet, kPtWeightedTrackInclusiveJet, knJetTrackCorrelations};
  
  // Indices for different track histogram categories
  enum enumTrackHistograms{kTrack, kUncorrectedTrack, kInclusiveTrack, kUncorrectedInclusiveTrack, knTrackCategories};
  
  // Indices for different single jet histogram categories
  enum enumSingleJet{kLeadingJet, kSubleadingJet, kAnyJet, kAnyLeadingJet, knSingleJetCategories};
  
  // Indices for different deltaPhi bins
  enum enumDeltaPhiBins{kWholePhi, kNearSide, kAwaySide, kBetweenPeaks, knDeltaPhiBins};
  
  // Indices for different deltaEta bins
  enum enumDeltaEtaBins{kWholeEta, kSignalEtaRegion, kBackgroundEtaRegion, knDeltaEtaBins};
  
  // Indices for different jet shape histograms
  enum enumJetShape{kJetShape, kJetShapeBinCount, knJetShapeTypes};
  
  // Indices for processing steps for jet-track correlation histograms
  enum enumProcessingStep{kMixedEventCorrection, kSeagullCorrection, kTrackDeltaRCorrection, kSpilloverCorrection, kBackgroundSubtraction, kJffCorrection, knProcessingSteps};
  
  // Dimensions for histogram arrays
  static const int kMaxCentralityBins = 5;       // Maximum allowed number of centrality bins
  static const int kMaxTrackPtBins = 10;         // Maximum allowed number of track pT bins
  static const int knFittedFlowComponents = 4;   // Number of fitted flow components
  static const int kMaxAsymmetryBins = 5;        // Maximum allowed number of dijet asymmetry bins
  static const int knGenJetPtBins = 45;          // Number of generator level jet pT bins for jet pT closures
  static const int knJetPtBins = 3;              // Number of leading jet pT bins for asymmetry histograms
  static const int knJetEtaBins = 50;            // Number of jet eta bins for jet pT closures
  
private:
  
  // Naming for different correlation types
  TString fCorrelationTypeString[knCorrelationTypes] = {"Same Event","Mixed Event","Normalized ME","Corrected","Background subtracted","Background","Background overlap","Jet Shape Bin Map"};
  TString fCompactCorrelationTypeString[knCorrelationTypes] = {"_SameEvent","_MixedEvent","_NormalizedME","_Corrected","_BackgroundSubtracted","_Background","_BackgroundOverlap","_JetShapeBinMap"};
  
  // Naming for jet-track correlation histograms
  const char* fJetTrackHistogramNames[knJetTrackCorrelations] = {"trackLeadingJet","trackLeadingJetUncorrected","trackLeadingJetPtWeighted","trackSubleadingJet","trackSubleadingJetUncorrected","trackSubleadingJetPtWeighted","trackJetInclusive","trackJetInclusivePtWeighted"}; // Names that different histograms have in the input file
  const char* fJetTrackAxisNames[knJetTrackCorrelations] = {"Track-LJet","UC Track-LJet","p_{T}w Track-LJet","Track-SJet","UC Track-SJet","p_{T}w Track-SJet","Trk-IncJet","p_{T}w Trk-IncJet"}; // Names attached to the figure axes
  
  // Naming for track histograms
  const char* fTrackHistogramNames[knTrackCategories] = {"track","trackUncorrected","trackInclusive","trackInclusiveUncorrected"}; // Names that different track histograms have in the input file
  const char* fTrackAxisNames[knTrackCategories] = {"Track","Uncorrected track", "Inclusive track", "UC Inclusive Track"}; // Names attached to the figure axes
  
  // Naming for single jet histograms
  const char* fSingleJetHistogramName[knSingleJetCategories] = {"leadingJet","subleadingJet","anyJet","anyLeadingJet"}; // Names that different single jet histograms have in the input file
  const char* fSingleJetAxisNames[knSingleJetCategories] = {"Leading jet","Subleading jet","Inclusive jet","Leading jet"}; // Names attached to the figure axes
  
  // Naming for jet shape histograms
  const char* fJetShapeHistogramName[knJetShapeTypes] = {"JetShape","JetShapeCounts"};
  const char* fJetShapeYAxisNames[knJetShapeTypes] = {"P(#Deltar)","counts"};
  
  // Naming for deltaEta bins
  const char* fDeltaEtaString[knDeltaEtaBins] = {"","Signal #Delta#eta","Background #Delta#eta"};
  const char* fCompactDeltaEtaString[knDeltaEtaBins] = {"","_SignalDeltaEta","_BackgroundDeltaEta"};
  
  // Naming for closure particle
  const char* fClosureParticleName[DijetHistograms::knClosureParticleTypes+1] = {"_quark","_gluon",""};
  
public:
  
  DijetHistogramManager();                                    // Default constructor
  DijetHistogramManager(TFile *inputFile);                    // Constructor
  DijetHistogramManager(TFile *inputFile, TFile* mixingFile); // Constructor with mixing file
  DijetHistogramManager(TFile *inputFile, DijetCard *card);   // Constructor with card
  DijetHistogramManager(TFile *inputFile, TFile* mixingFile, DijetCard *card); // Constructor with mixing file and a card
  DijetHistogramManager(const DijetHistogramManager& in);     // Copy constructor
  ~DijetHistogramManager();                                   // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void ProcessHistograms();       // Do the mixed event correction, subtract the background and calculate jet shape
  void Write(const char* fileName, const char* fileOption);          // Write all the loaded histograms into a file
  void WriteJetShape(const char* fileName, const char* fileOption);  // Write only the jet shape histograms into a file
  void WriteSkim(const char* fileName, const char* fileOption);      // Write jet shapes and final deltaEta histograms into a file
  void LoadProcessedHistograms(); // Load processed histograms from the inputfile
  void ApplyJffCorrection(JffCorrector *jffCorrectionFinder);  // Apply the JFF correction to relevant histograms
  void CalculateJetShape();       // Calculate the jet shape from the background subtracted histograms
  void NormalizeJetShape();       // Normalize the jet shape histograms
  void ProjectFinalDeltaEta(const bool splitZero = false);    // Project the final deltaEta yield results from two dimensional distributions
  
  // Setters for binning information
  void SetCentralityBins(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true); // Set up centrality bin indices according to provided bin borders
  void SetTrackPtBins(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true);    // Set up track pT bin indices according to provided bin borders
  void SetDeltaPhiBins(const bool readBinsFromFile, const double *lowBinBorders, const double *highBinBorders, TString deltaPhiStrings[knDeltaPhiBins], TString compactDeltaPhiStrings[knDeltaPhiBins], bool setIndices = true); //  Set up deltaPhi bin indices according to provided bin borders and bin names
  
  // Setters for event information and dijets
  void SetLoadEventInformation(const bool loadOrNot); // Setter for loading event information
  void SetLoadDijetHistograms(const bool loadOrNot);  // Setter for loading dijet histograms
  
  // Setters for single jets
  void SetLoadLeadingJetHistograms(const bool loadOrNot);    // Setter for loading leading jet histograms
  void SetLoadSubleadingJetHistograms(const bool loadOrNot); // Setter for loading subleading jet histograms
  void SetLoadAnyJetHistograms(const bool loadOrNot);        // Setter for loading all jet histograms
  void SetLoadAnyLeadingJetHistograms(const bool loadOrNot); // Setter for loading all leading jet histograms
  void SetLoadAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny, const bool drawAnyLeading);   // Setter for loading jet histograms
  
  // Setters for tracks
  void SetLoadTracks(const bool loadOrNot);            // Setter for loading tracks
  void SetLoadTracksUncorrected(const bool loadOrNot); // Setter for loading uncorrected tracks
  void SetLoadAllTracks(const bool drawTracks, const bool drawUncorrected); // Setter for loading all track histograms
  void SetLoadInclusiveTracks(const bool loadOrNot);            // Setter for loading tracks
  void SetLoadInclusiveTracksUncorrected(const bool loadOrNot); // Setter for loading uncorrected tracks
  void SetLoadAllInclusiveTracks(const bool drawTracks, const bool drawUncorrected); // Setter for loading all track histograms
  
  // Setters for loading leading jet-track correlations
  void SetLoadTrackLeadingJetCorrelations(const bool loadOrNot);            // Setter for loading leading jet-track correlations
  void SetLoadTrackLeadingJetCorrelationsUncorrected(const bool loadOrNot); // Setter for loading uncorrected leading jet-track correlations
  void SetLoadTrackLeadingJetCorrelationsPtWeighted(const bool loadOrNot);  // Setter for loading pT weighted leading jet-track correlations
  void SetLoadAllTrackLeadingJetCorrelations(const bool drawLeading, const bool drawUncorrected, const bool drawPtWeighted); // Setter for loading all correlations related to tracks and leading jets
  
  // Setters for loading subleading jet-track correlations
  void SetLoadTrackSubleadingJetCorrelations(const bool loadOrNot);            // Setter for loading subleading jet-track correlations
  void SetLoadTrackSubleadingJetCorrelationsUncorrected(const bool loadOrNot); // Setter for loading uncorrected subleading jet-track correlations
  void SetLoadTrackSubleadingJetCorrelationsPtWeighted(const bool loadOrNot);  // Setter for loading pT weighted subleading jet-track correlations
  void SetLoadAllTrackSubleadingJetCorrelations(const bool drawSubleading, const bool drawUncorrected, const bool drawPtWeighted); // Setter for loading all correlations related to tracks and subleading jets
  
  // Setters for loading inclusive jet-track correlations
  void SetLoadTrackInclusiveJetCorrelations(const bool loadOrNot);           // Setter for loading inclusive jet-track correlations
  void SetLoadTrackInclusiveJetCorrelationsPtWeighted(const bool loadOrNot); // Setter for leading pT weighted inclusive jet-track correlations
  void SetLoadAllTrackInclusiveJetCorrelations(const bool loadInclusive, const bool loadPtWeighted);        // Setter for loading all correlations related to tracks and inclusive jets
  void SetCorrelationJetFlavor(const int iFlavor);  // For Monte Carlo, can select if we are looking for quark or gluon initiated jets
  
  // Setter for loading additional histograms
  void SetLoad2DHistograms(const bool loadOrNot);           // Setter for loading two-dimensional histograms
  void SetLoadJetPtClosureHistograms(const bool loadOrNot); // Setter for loading jet pT closure histograms
  
  // Setters for ranges for different bins
  void SetCentralityBinRange(const int first, const int last); // Setter for centrality bin range
  void SetTrackPtBinRange(const int first, const int last);    // Setter for track pT bin range
  void SetAsymmetryBinRange(const int first, const int last);  // Setter for processing asymmetry bins
  void SetPreprocess(const int preprocess);                    // Setter for preprocessing (only load and write same and/or mixed event)
  
  // Setter for used DijetMethods
  void SetDijetMethods(DijetMethods* newMethods); // Setter for used DijetMethods
  
  // Setters for corrections
  void SetJffCorrection(TFile *jffFile, const bool applyCorrection);              // Setter for JFF corrector and flag
  void SetSpilloverCorrection(TFile* spilloverFile, const bool applyCorrection);  // Setter for spillover corrector and flag
  void SetTrackDeltaRCorrection(TFile* trackingFile, const bool applyCorrection, const bool applyResidualScale); // Setter for residual tracking corrector and flag
  void SetSeagullCorrection(const bool applyCorrection);                          // Setter for seagull correction flag
  void SetManualSpilloverCleaning(const bool applyCleaning);                      // Setter for manually clearing fluctuations from spillover correction
  
  // Getters for number of bins in histograms
  int GetNCentralityBins() const; // Getter for the number of centrality bins
  int GetNTrackPtBins() const;    // Getter for the number of track pT bins
  int GetNJetPtBins() const;      // Getter for the number of jet pT bins
  int GetNAsymmetryBins() const;  // Getter for the number of dijet asymmetry AJ bins
  double GetCentralityBinBorder(const int iCentrality) const;  // Getter for i:th centrality bin border
  double GetTrackPtBinBorder(const int iTrackPt) const;        // Getter for i:th track pT bin border
  double GetJetPtBinBorder(const int iJetPt) const;            // Getter for i:th jet pT bin border
  double GetDeltaPhiBorderLow(const int iDeltaPhi) const;      // Getter for i:th low deltaPhi border
  double GetDeltaPhiBorderHigh(const int iDeltaPhi) const;     // Getter for i:th high deltaPhi border
  
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
  TH1D* GetHistogramVertexZ() const;            // Getter for z-vertex histogram
  TH1D* GetHistogramVertexZWeighted() const;    // Getter for weighted z-vertex histogram
  TH1D* GetHistogramVertexZDijet() const;       // Getter for z-vertex histogram in dijet events
  TH1D* GetHistogramEvents() const;             // Getter for histogram for number of events surviving different event cuts
  TH1D* GetHistogramTrackCuts() const;          // Getter for histogram for number of tracks surviving different track cuts
  TH1D* GetHistogramCentrality() const;         // Getter for centrality histogram in all events
  TH1D* GetHistogramCentralityWeighted() const; // Getter for weighted centrality histogram in all events
  TH1D* GetHistogramCentralityDijet() const;    // Getter for centrality histogram in dijet events
  
  TH1D* GetHistogramMultiplicity(int iCentrality) const;               // Getter for multiplicity from all events
  TH1D* GetHistogramMultiplicityWeighted(int iCentrality) const;       // Getter for efficiency weighted multiplicity from all events
  TH1D* GetHistogramMultiplicityDijet(int iCentrality) const;          // Getter for multiplicity from dijet events
  TH1D* GetHistogramMultiplicityDijetWeighted(int iCentrality) const;  // Getter for efficiency weighted multiplicity from dijet events
  TH2D* GetHistogramMultiplicityMap() const;                           // Getter for multiplicity vs. centrality map
  TH2D* GetHistogramMultiplicityMapDijet() const;                      // Getter for multiplicity vs. centrality map in dijet events
  TH2D* GetHistogramWeightedMultiplicityMap() const;                   // Getter for efficiency weighted multiplicity vs. centrality map
  TH2D* GetHistogramWeightedMultiplicityMapDijet() const;              // Getter for efficiency weighted multiplicity vs. centrality map in dijet events
  
  // Getters for single jet histograms
  TH1D* GetHistogramJetPt(const int iJetType, int iCentrality, int iAsymmetry = kMaxAsymmetryBins) const;     // Jet pT histograms
  TH1D* GetHistogramJetPhi(const int iJetType, int iCentrality, int iAsymmetry = kMaxAsymmetryBins) const;    // Jet phi histograms
  TH1D* GetHistogramJetEta(const int iJetType, int iCentrality, int iAsymmetry = kMaxAsymmetryBins) const;    // Jet eta histograms
  TH2D* GetHistogramJetEtaPhi(const int iJetType, int iCentrality, int iAsymmetry = kMaxAsymmetryBins) const; // 2D eta-phi histogram for jets
  
  // Getters for dijet histograms
  TH1D* GetHistogramDijetDeltaPhi(const int iCentrality) const;                                                  // Dijet deltaPhi histograms
  TH1D* GetHistogramDijetAsymmetry(const int iCentrality, const int iJetPt = knJetPtBins) const;                 // Dijet asymmetry AJ histograms
  TH1D* GetHistogramDijetXj(const int iCentrality, const int iJetPt = knJetPtBins) const;                        // Dijet asymmetry xJ histograms
  TH2D* GetHistogramDijetLeadingVsSubleadingPt(const int iCentrality, int iAsymmetry = kMaxAsymmetryBins) const; // Leading versus subleading jet pT 2D histograms
  TH2D* GetHistogramDijetXjMatrix(const int iCentrality) const;  // Matrix between reconstructed and generated dijet xj values
  
  // Getters for histograms for tracks in dijet events
  TH1D* GetHistogramTrackPt(const int iTrackType, const int iCorrelationType, const int iCentrality) const;                      // Track pT histograms
  TH1D* GetHistogramTrackPhi(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const;    // Track phi histograms
  TH1D* GetHistogramTrackEta(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const;    // Track eta histograms
  TH2D* GetHistogramTrackEtaPhi(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const; // 2D eta-phi histogram for track
  
  // Getters for track-leading jet correlation histograms
  TH1D* GetHistogramJetTrackDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, int iAsymmetry, const int iCentrality, const int iTrackPt, const int iDeltaEta) const;  // DeltaPhi between jet and track
  TH1D* GetHistogramJetTrackDeltaEta(const int iJetTrackCorrelation, const int iCorrelationType, int iAsymmetry, const int iCentrality, const int iTrackPt, const int iDeltaPhiRegion) const; // DeltaEta between jet and track
  TH2D* GetHistogramJetTrackDeltaEtaDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, int iAsymmetry, const int iCentrality, const int iTrackPt) const;         // DeltaEta and deltaPhi between jet and track
  TH1D* GetHistogramJetTrackDeltaEtaFinal(const int iJetTrackCorrelation, int iAsymmetry, const int iCentrality, const int iTrackPt) const; // DeltaEta between jet and track using final binning
  
  // Getters for jet shape histograms
  TH1D* GetHistogramJetShape(const int iJetShapeType, const int iJetTrackCorrelation, int iAsymmetry, int iCentrality, const int iTrackPt) const;  // Jet shape histograms
  
  // Getter for jet pT closure histograms
  TH1D* GetHistogramJetPtClosure(const int iClosureType, const int iGenPtBin, const int iEtaBin, const int iCentrality, int iAsymmetry, const int iClosureParticle) const; // Jet pT closure
  TH2D* GetHistogramJetPtClosureReactionPlane(const int iClosureType, const int iGenPtBin, const int iEtaBin, const int iCentrality, int iAsymmetry, const int iClosureParticle) const; // Jet pT closure with respect to the reaction plane
  
  TH1D* GetOneDimensionalHistogram(TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0, int bin6 = 0) const; // Getter for any one-dimensional histogram based on input string
  TH2D* GetTwoDimensionalHistogram(TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0) const; // Getter for any two-dimensional histogram based on input string
  
  // Check if any of the loaded histograms is NULL
  bool CheckNull(const bool preprocess) const;   // Check if any of the loaded histograms is NULL
  
  // Getters for the loaded centrality and track pT bins
  int GetFirstCentralityBin() const;  // Get the first loaded centrality bin
  int GetLastCentralityBin() const;   // Get the last loaded centrality bin
  int GetFirstTrackPtBin() const;     // Get the first loaded track pT bin
  int GetLastTrackPtBin() const;      // Get the last loaded track pT bin
  int GetFirstAsymmetryBin() const;   // Get the first loaded asymmetry bin
  int GetLastAsymmetryBin() const;    // Get the last loaded asymmetry bin
  
  // Getters for normalization information
  int GetNEvents() const;                      // Getter for the number of events passing the cuts
  int GetNDijets() const;                      // Getter for the number of dijets
  double GetPtIntegral(const int iCentrality, int iAsymmetry = kMaxAsymmetryBins) const; // Getter for integral over leading jet pT in a given centrality and dijet asymmetry bin
  double GetAnyLeadingJetPtIntegral(int iCentrality) const; // Getter for integral over all leading jets with pT > 120 GeV in a given centrality bin
  double GetInclusiveJetPtIntegral(int iCentrality, const double minPt = 120) const; // Getter for integral over inclusive jet pT above minPt in a given centrality bin
  double GetJetShapeNormalizationFactor(const int iJetTrack, const int iCentrality, int iAsymmetry = kMaxAsymmetryBins) const;       // Factor for normalizing capital rho to lower case rho
  
  // Setter for avoiding possible peaks in mixed event distribution
  void SetAvoidMixingPeak(const bool avoid);     // Avoid peak region in mixed event distribution
  void SetImproviseMixing(const bool improvise); // Create mixed event distributions from deltaPhi sideband region
  void SetDefaultMixingDeltaEtaFitRange(const double fitRange);  // Default fit range used to normalize the mixed event
  
  // Setter for histogram processing options
  void SetProcessingStartLevel(const int processingLevel);  // Select from which level we start processing the histograms
  
  // Getter for the card
  DijetCard* GetCard() const;  // Getter for the JCard
  
  // Histogram combining
  void CombineHistograms(DijetHistogramManager *combinedHistograms, double factorThis = 1, double factorCombine = 1);   // Combine histograms between this histogram manager and provided histogram manager
  
private:
  
  // Data members
  TFile *fInputFile;                  // File from which the histograms are read
  TFile *fMixingFile;                 // File for mixed event histograms. If null, mixed events are read from the same file as inputfile
  DijetCard *fCard;                   // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;           // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy;    // Same a before but without white spaces and dots
  DijetMethods *fMethods;             // DijetMethods for processing the loaded histograms
  JffCorrector *fJffCorrectionFinder; // Class for providing JFF correction for final deltaEta-deltaPhi distributions
  
  // ==============================================
  // =========== Flags for corrections ============
  // ==============================================
  bool fApplyJffCorrection;            // Flag for applying the JFF correction
  bool fApplySpilloverCorrection;      // Flag for applying the spillover correction
  bool fApplyTrackDeltaRCorrection;    // Flag for applying the residual tracking correction
  bool fApplyTrackDeltaRResidualScale; // Flag for applying residual scaling (match reco and gen yields) for the tracking correction
  bool fApplySeagullCorrection;        // Flag for applying the seagull correction
  bool fManualSpilloverCleaning;       // Flag for manually cleaning fluctuations from the spillover correction
  
  // ==============================================
  // ======== Flags for histograms to load ========
  // ==============================================
  
  bool fLoadEventInformation;                              // Load the event information histograms
  bool fLoadDijetHistograms;                               // Load the dijet histograms
  bool fLoadSingleJets[knSingleJetCategories];             // Load the single jet histograms
  bool fLoadTracks[knTrackCategories];                     // Load the track histograms
  bool fLoadJetTrackCorrelations[knJetTrackCorrelations];  // Load the jet-track correlation histograms
  bool fLoad2DHistograms;                                  // Load also two-dimensional (eta,phi) and (deltaEta,deltaPhi) histograms
  bool fLoadJetPtClosureHistograms;                        // Load the jet pT closure histograms
  int  fCorrelationJetFlavor;                              // Select the flavor for loaded jets (1 = Quark, 2 = Gluon)
  
  // ==============================================
  // ======== Ranges of histograms to load ========
  // ==============================================
  
  // Drawn centrality bins
  int fFirstLoadedCentralityBin;  // First centrality bin that is loaded
  int fLastLoadedCentralityBin;   // Last centrality bin that is loaded
  int fFirstLoadedTrackPtBin;     // First track pT bin that is loaded
  int fLastLoadedTrackPtBin;      // Last track pT bin that is loaded
  int fFirstLoadedAsymmetryBin;   // First asymmetry bin that is loaded
  int fLastLoadedAsymmetryBin;    // Last asymmetry bin that is loaded
  bool fPreprocess;               // For preprocessing, only load and write same and mixed event
  int fPreprocessLevel;           // Level of preprocessing: 0 = Only same event, 1 = Only mixed event, 2 = Same and mixed event
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[kMaxCentralityBins+1];    // Indices for centrality bins in centrality binned histograms
  double fCentralityBinBorders[kMaxCentralityBins+1]; // Centrality bin borders, from which bin indices are obtained
  int fTrackPtBinIndices[kMaxTrackPtBins+1];          // Indices for track pT bins in track pT binned histograms
  int fFineTrackPtBinIndices[kMaxTrackPtBins+1];      // Indices for track pT bins in fine track pT binned histograms
  double fTrackPtBinBorders[kMaxTrackPtBins+1];       // Track pT bin borders, from which bin indices are obtained
  int fLowDeltaPhiBinIndices[knDeltaPhiBins];         // Indices for low bin borders in deltaPhi binned histograms
  int fHighDeltaPhiBinIndices[knDeltaPhiBins];        // Indices for high bin borders in deltaPhi binned histograms
  double fLowDeltaPhiBinBorders[knDeltaPhiBins];      // Low bin borders in deltaPhi binned histograms
  double fHighDeltaPhiBinBorders[knDeltaPhiBins];     // High bin borders in deltaPhi binned histograms
  TString fDeltaPhiString[knDeltaPhiBins];            // Names for different deltaPhi bins
  TString fCompactDeltaPhiString[knDeltaPhiBins];     // Names added to figure names for deltaPhi bins
  int fJetPtBinIndices[knJetPtBins+1];                // Indices for leading jet pT bins for asymmetry histograms
  double fJetPtBinBorders[knJetPtBins+1];             // Bin borders for the leading jet pT bins from which the indices are obtained
  int fnCentralityBins;                               // Number of centrality bins in the JCard of the data file
  int fnTrackPtBins;                                  // Number of track pT bins in the JCard of the data file
  int fnAsymmetryBins;                                // Number of asymmetry bins in the JCard of the data file
  double fAsymmetryBinBorders[kMaxCentralityBins+1];  // Asymmetry bin borders
  TString fAsymmetryBinName[kMaxAsymmetryBins+1];     // Name given to each asymmetry bin
  
  // =============================================
  // =========   Histogram processing   ==========
  // =============================================
  int fProcessingStartLevel;  // Determine from which step the processing is started. It is assumed that previous steps are taken in the input file
  bool fAvoidMixingPeak;      // Avoid region around (0,0) in the mixed event distribution to stay clear of possible peaks
  bool fImproviseMixing;      // Create mixed event distributions from same event deltaPhi side band region
  double fDefaultMixingDeltaEtaFitRange; // Default deltaEta fit range in mixed event normalization
  
  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================
  
  // Event information histograms
  TH1D *fhVertexZ;            // Vertex z position
  TH1D *fhVertexZWeighted;    // Weighted vertex z-position (only meaningfull for MC)
  TH1D *fhVertexZDijet;       // Vertex z position in dijet events
  TH1D *fhEvents;             // Number of events surviving different event cuts
  TH1D *fhTrackCuts;          // Number of tracks surviving different track cuts
  TH1D *fhTrackCutsInclusive; // Number of inclusive tracks surviving different track cuts
  TH1D *fhCentrality;         // Centrality of all events
  TH1D *fhCentralityWeighted; // Weighted centrality distribution in all events (only meaningful for MC)
  TH1D *fhCentralityDijet;    // Centrality of dijet events
  TH1D *fhPtHat;              // pT hat for MC events (only meaningful for MC)
  TH1D *fhPtHatWeighted;      // Weighted pT hat distribution (only meaningful for MC)
  
  TH1D *fhMultiplicity[kMaxCentralityBins+1];              // Multiplicity form all events
  TH1D *fhMultiplicityWeighted[kMaxCentralityBins+1];      // Efficiency weighted multiplicity form all events
  TH1D *fhMultiplicityDijet[kMaxCentralityBins+1];         // Multiplicity form dijet events
  TH1D *fhMultiplicityDijetWeighted[kMaxCentralityBins+1]; // Efficiency weighted multiplicity form dijet events
  TH2D *fhMultiplicityMap;                                 // Multiplicity vs. centrality map
  TH2D *fhMultiplicityMapDijet;                            // Multiplicity vs. centrality map in dijet events
  TH2D *fhMultiplicityMapWeighted;                         // Efficiency weighted multiplicity vs. centrality map
  TH2D *fhMultiplicityMapWeightedDijet;                    // Efficiency weighted multiplicity vs. centrality map in dijet events
  
  // Histograms for single jets
  TH1D *fhJetPt[knSingleJetCategories][kMaxCentralityBins][kMaxAsymmetryBins+1];          // Jet pT histograms
  TH1D *fhJetPhi[knSingleJetCategories][kMaxCentralityBins][kMaxAsymmetryBins+1];         // Jet phi histograms
  TH1D *fhJetEta[knSingleJetCategories][kMaxCentralityBins][kMaxAsymmetryBins+1];         // Jet eta histograms
  TH2D *fhJetEtaPhi[knSingleJetCategories][kMaxCentralityBins][kMaxAsymmetryBins+1];      // 2D eta-phi histogram for jets
  
  // Extra histograms for smearing study
  TH1D *fhJetPtUncertaintyPlus[knSingleJetCategories][kMaxCentralityBins][kMaxAsymmetryBins+1];   // Jet pT plus error histograms
  TH1D *fhJetPtUncertaintyMinus[knSingleJetCategories][kMaxCentralityBins][kMaxAsymmetryBins+1];  // Jet pT minus error histograms
  TH2D *fhJetPtSmearMap[knSingleJetCategories][kMaxCentralityBins][kMaxAsymmetryBins+1];  // 2D jet pT smearing map for sanity check
  
  // Histograms for dijets
  TH1D *fhDijetDphi[kMaxCentralityBins];                     // Dijet deltaPhi histograms
  TH1D *fhDijetAsymmetry[kMaxCentralityBins][knJetPtBins+1]; // Dijet asymmetry AJ histograms
  TH1D *fhDijetXj[kMaxCentralityBins][knJetPtBins+1];        // Dijet asymmetry xj histograms
  TH2D *fhDijetLeadingVsSubleadingPt[kMaxAsymmetryBins+1][kMaxCentralityBins];   // Leading versus subleading jet pT 2D histograms
  TH2D *fhDijetXjMatrix[kMaxCentralityBins];                 // Matrix between reconstructed and generated dijet xj values
  
  // Histograms for tracks in dijet events
  TH1D *fhTrackPt[knTrackCategories][knCorrelationTypes][kMaxCentralityBins];                        // Track pT histograms
  TH1D *fhTrackPhi[knTrackCategories][knCorrelationTypes][kMaxCentralityBins][kMaxTrackPtBins+1];    // Track phi histograms
  TH1D *fhTrackEta[knTrackCategories][knCorrelationTypes][kMaxCentralityBins][kMaxTrackPtBins+1];    // Track eta histograms
  TH2D *fhTrackEtaPhi[knTrackCategories][knCorrelationTypes][kMaxCentralityBins][kMaxTrackPtBins+1]; // 2D eta-phi histogram for track
  
  // Histograms for track-leading jet correlations
  TH1D *fhJetTrackDeltaPhi[knJetTrackCorrelations][knCorrelationTypes][kMaxAsymmetryBins+1][kMaxCentralityBins][kMaxTrackPtBins][knDeltaEtaBins]; // DeltaPhi between jet and track
  TH1D *fhJetTrackDeltaEta[knJetTrackCorrelations][knCorrelationTypes][kMaxAsymmetryBins+1][kMaxCentralityBins][kMaxTrackPtBins][knDeltaPhiBins]; // DeltaEta between jet and track
  TH2D *fhJetTrackDeltaEtaDeltaPhi[knJetTrackCorrelations][knCorrelationTypes][kMaxAsymmetryBins+1][kMaxCentralityBins][kMaxTrackPtBins];         // DeltaEta and deltaPhi between jet and track
  TH1D *fhJetTrackDeltaEtaFinalResult[knJetTrackCorrelations][kMaxAsymmetryBins+1][kMaxCentralityBins][kMaxTrackPtBins]; // DeltaEta between jet and track for final result binning
  
  // Jet shape histograms
  TH1D *fhJetShape[knJetShapeTypes][knJetTrackCorrelations][kMaxAsymmetryBins+1][kMaxCentralityBins][kMaxTrackPtBins];  // Jet shape histograms
  
  // Histograms for jet pT closure
  TH1D *fhJetPtClosure[DijetHistograms::knClosureTypes][knGenJetPtBins+1][knJetEtaBins+1][kMaxCentralityBins][kMaxAsymmetryBins+1][DijetHistograms::knClosureParticleTypes+1]; // Jet pT closure
  TH2D *fhJetPtClosureReactionPlane[DijetHistograms::knClosureTypes][knGenJetPtBins+1][knJetEtaBins+1][kMaxCentralityBins][kMaxAsymmetryBins+1][DijetHistograms::knClosureParticleTypes+1]; // Jet pT closure 
  
  // QA histograms for seagull correction
  TH1D *fhSeagullDeltaEta[knJetTrackCorrelations][kMaxAsymmetryBins+1][kMaxCentralityBins][kMaxTrackPtBins]; // Background eta projection for seagull fit
  TF1 *fSeagullFit[knJetTrackCorrelations][kMaxAsymmetryBins+1][kMaxCentralityBins][kMaxTrackPtBins];        // The seagull fit to background eta projection
  
  // Private methods
  void InitializeFromCard(); // Initialize several member variables from DijetCard
  
  // Binning related methods
  void SetBinIndices(const int nBins, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices); // Read the bin indices for given bin borders
  void SetBinBordersAndIndices(const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices); // Read the bin indices for given bin borders
  void SetBinIndices(const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinIndices(const char* histogramName, const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetLoadJetTrackCorrelations(const bool loadOrNot, const int primaryIndex, const int connectedIndex); // Setter for drawing and loading the jet-track correlation histograms
  int GetConnectedIndex(const int jetTrackIndex) const;
  
  // Finders for histograms with different amount of restrictions
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex); // Extract a histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0); // Extract a histogram using given axis restrictions from THnSparseD
  
  // Loaders for different groups of histograms
  void LoadMultiplicityHistograms(); // Loader for multiplicity histograms
  void LoadSingleJetHistograms(); // Loader for single jet histograms
  void LoadDijetHistograms(); // Loader for dijet histograms
  void LoadTrackHistograms(); // Loader for track histograms
  void LoadJetTrackCorrelationHistograms(); // Loader for jet-track correlation histograms
  void LoadJetPtClosureHistograms(); // Loader for jet pT closure histograms
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int& first, int& last); // Sanity check for given binning
  int BinIndexCheck(const int nBins, const int binIndex) const; // Check that given index is in defined range
  
  // Mixed event correction, background subtraction, different corrections and projections
  void DoMixedEventCorrection();  // Apply mixed event correction for jet-track correlation histograms
  void DoSeagullCorrection();     // Apply seagull correction to jet-track correlation histograms
  void DoTrackDeltaRCorrection(); // Apply track deltaR correction to jet-track correlation histograms
  void DoSpilloverCorrection();   // Apply spillover correction to jet-track correlation histograms
  void SubtractBackground();      // Subtract background from the spillover corrected distribution
  void DoJffCorrection();         // Do the JFF correction to the background subtracted distribution
  void DoProjections();           // Take projections of processed two-dimensional histograms
  
  // Methods for histogram writing
  void WriteSingleJetHistograms();           // Write the single jet histograms to the file that is currently open
  void WriteDijetHistograms();               // Write the dijet histograms to the file that is currently open
  void WriteTrackHistograms();               // Write the track histograms to the file that is currently open
  void WriteJetTrackCorrelationHistograms(); // Write the jet-track correlation histograms to the file that is currently open
  void WriteJetShapeHistograms();            // Write the jet shape histograms to the file that is currently open
  void WriteFinalDeltaEtaHistograms();       // Write the rebinned deltaEta yield histograms to the file that is currently open
  void WriteClosureHistograms();             // Write the closure histograms to the file that is currently open
  
};

#endif
