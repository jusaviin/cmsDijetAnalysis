#ifndef DIJETCARD_H
#define DIJETCARD_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TObjString.h>
#include <TVectorT.h>

/*
 * DijetCard class
 *
 * This class reads the ConfigurationCard from the input root file and decodes the
 * necessary information for the analysis.
 */
class DijetCard {
  
public:
 
  // Indices for card entries
  enum enumCardEntries{
    kDataType,                  // Data tyoe in the data file (pp, PbPb, pp MC, PbPb MC)
    kMcCorrelationType,         // Monte Carlo correlation type (RecoReco, RecoGen, GenReco, GenGen)
    kMatchJets,                 // 0 = Reco and Gen jets are not matched, 1 = They are matched
    kForestType,                // 0 = High forest, 1 = Skim Forest
    kReadMode,                  // Certain files have different structure, so branch names need to be adjusted for them
    kJetType,                   // 0 = Calorimeter jets, 1 = PF jets
    kJetAxis,                   // 0 = Anti-kt axis, 1 = Leading PF candidate axis, 2 = WTA axis
    kJetEtaCut,                 // Eta cut for jets
    kSearchEtaCut,              // Eta cut defining the region from which the dijets are searched
    kMaxPtCut,                  // Maximum allowed pT for jets
    kMinPtCut,                  // Minimum allowed pT for the leading jet
    kSubleadingPtCut,           // Minimum allowed pT for the subleading jet
    kDeltaPhiCut,               // Required deltaPhi separation between leading and subleading jets
    kMinMaxTrackPtFraction,     // Minimum fraction of jet pT taken by the highest pT track in jet
    kMaxMaxTrackPtFraction,     // Maximum fraction of jet pT taken by the highest pT track in jet
    kTrackEtaCut,               // Eta cut for tracks
    kMinTrackPtCut,             // Minimum accepted track pT
    kMaxTrackPtRelativeError,   // Maximum relative error allowed for track pT
    kVertexMaxDistance,         // Maximum allowed distance of tracks from reconstructed vertex
    kCalorimeterSignalLimitPt,  // Limit for track pT above which a signal in calorimeters is required
    kHighPtEtFraction,          // Minimum fraction between pT and Et for high pT tracks
    kChi2QualityCut,            // Maximum accepted chi2 for reconstructed tracks
    kMinimumTrackHits,          // Minimum number of hits in tracking for a track
    kSubeventCut,               // 0 = Subevent 0 (Pythia), 1 = Subevent > 0 (Hydjet), 2 = No subevent selection
    kZVertexCut,                // Maximum accepted vz in the event
    kLowPtHatCut,               // Minimum accepted pT hat
    kHighPtHatCut,              // Maximum accepted pT hat
    kAsymmetryBinType,          // 0 = Asymmetry binning in AJ, 1 = Asymmetry binning in xJ
    kCentralityBinEdges,        // Centrality bin edges
    kTrackPtBinEdges,           // Track pT bin edges
    kAsymmetryBinEdges,         // Asymmetry bin edges
    kPtHatBinEdges,             // pT hat bin edges
    kDoEventMixing,             // 0 = Do not mix events, 1 = Mix events
    kMixWithPool,               // 0 = Use poolless mixing algorithm, 1 = Mix with pool
    kNMixedEventsPerDijet,      // Number of events used in the event mixing
    kVzTolerance,               // For poolless mixing, the maximum accepted distance between same event vz and mixed event vz
    kMixingVzBinWidth,          // The vz bin width, in case we use mixing pool
    kMixingHiBinWidth,          // The hibin width, in case we use mixing pool
    kMixingPoolDepth,           // Number of events in each vz and hibin collected to the mixing pool
    kJffCorrection,             // For postprocessing: 0 = No JFF correction, 1 = Apply JFF correction
    kSpilloverCorrection,       // For postprocessing: 0 = No spillover correction, 1 = Apply spillover correction
    kSeagullCorrection,         // For postprocessing: 0 = No seagull correction, 1 = Apply seagull correction
    kSmoothMixing,              // For postprocessing: 0 = No smoothening in mixing, 1 = Smoothen mixed event distribution
    kImprovisedMixing,          // For postprocessing: 0 = Use actual mixing, 1 = Improvise mixing from deltaPhi side band
    kAdjustBackground,          // For postprocessing: 0 = No adjustment between leading and subleading backgrounds, 1 = Match leading and subleading backgrounds
    kLowDeltaPhiBinBorders,     // For postprocessing: Low deltaPhi bin borders used when producing the file
    kHighDeltaPhiBinBorders,    // For postprocessing: High deltaPhi bin borders used when producing the file
    knEntries};                 // Number of entries in the card
  
  // Enumeration for input files used in postprocessing
  enum enumFileNames{kInputFileName,kJffCorrectionFileName,kSpilloverCorrectionFileName,knFileNames};
  
private:
  
  // Names for each entry read from the configuration card
  const char *fCardEntryNames[knEntries] = {"DataType","McCorrelationType","MatchJets","ForestType","ReadMode","JetType","JetAxis","JetEtaCut","SearchEtaCut","MaxPtCut","MinPtCut","SubleadingPtCut","DeltaPhiCut","MinMaxTrackPtFraction","MaxMaxTrackPtFraction","TrackEtaCut","MinTrackPtCut","MaxTrackPtRelativeError","VertexMaxDistance","CalorimeterSignalLimitPt","HighPtEtFraction","Chi2QualityCut","MinimumTrackHits","SubeventCut","ZVertexCut","LowPtHatCut","HighPtHatCut","AsymmetryBinType","CentralityBinEdges","TrackPtBinEdges","AsymmetryBinEdges","PtHatBinEdges","DoEventMixing","MixWithPool","NMixedEventsPerDijet","VzTolerance","MixingVzBinWidth","MixingHiBinWidth","MixingPoolDepth","JffCorrection","SpilloverCorrection","SeagullCorrection","SmoothMixing","ImprovisedMixing","AdjustedBackground","LowDeltaPhiBinBorders","HighDeltaPhiBinBorders"};
  const char *fFileNameType[knFileNames] = {"input", "JFF correction","spillover correction"};
  const char *fFileNameSaveName[knFileNames] = {"InputFile", "JFFCorrectionFile","SpilloverCorrectionFile"};
  
  TFile *fInputFile;         // Input file from which all the data is read
  TString fCardDirectory;    // Path to the ConfigurationCard directory
  int fDataType;             // Total number of centrality bins in the analysis
  int fMonteCarloType;       // Type of Monte Carlo used for jet-track correlations
  int fAsymmetryBinType;     // Used asymmetry binning (0 = AJ, 1 = xJ)
  TString fDataTypeString;   // Total number of eta gaps in the analysis
  
  void FindDataTypeString(); // Construct a data type string based on information on the card
  void ReadVectors();        // Read the vectors from the file
  
  // Vectors for all the lines inside the card
  TVectorT<float> *fCardEntries[knEntries];   // Array of all the vectors in the card
  TObjString *fFileNames[knFileNames];        // Array for filenames used in postprocessing
  
  // Private methods
  int GetNBins(const int index) const;                            // Get the number of bins for internal index
  double GetLowBinBorder(const int index, const int iBin) const;  // Get the low border of i:th bin from internal index
  double GetHighBinBorder(const int index, const int iBin) const; // Get the high border of i:th bin from internal index
   
public:
  
  DijetCard(TFile *inFile); // Contructor with input file
  ~DijetCard();             // Destructor
  
  TString GetDataType() const;   // Getter for data type string
  void Write(TDirectory *file);  // Write the contents of the card to a file
  double GetMaxDeltaEta() const; // Get maximum deltaEta possible using jet and track cuts in the card
  void Print() const;            // Print the contents of the card to the console
  
  int GetNCentralityBins() const; // Get the number of centrality bins
  int GetNTrackPtBins() const;    // Get the number of track pT bins
  int GetNAsymmetryBins() const;  // Get the number of asymmetry bins
  double GetLowBinBorderCentrality(const int iBin) const;  // Get the low border of i:th centrality bin
  double GetLowBinBorderTrackPt(const int iBin) const;     // Get the low border of i:th track pT bin
  double GetLowBinBorderAsymmetry(const int iBin) const;   // Get the low border of i:th asymmetry bin
  double GetHighBinBorderCentrality(const int iBin) const; // Get the high border of i:th centrality bin
  double GetHighBinBorderTrackPt(const int iBin) const;    // Get the high border of i:th track pT bin
  double GetHighBinBorderAsymmetry(const int iBin) const;  // Get the high border of i:th asymmetry bin
  const char *GetAsymmetryBinType(TString latexIt = "") const; // Get a description of the used asymmetry bin type
  int GetNDeltaPhiBins() const; // Get the number of deltaPhi bins
  double GetLowBinBorderDeltaPhi(const int iBin) const;  // Get the low border of i:th deltaPhi bin
  double GetHighBinBorderDeltaPhi(const int iBin) const; // Get the high border of i:th deltaPhi bin
  
  void AddOneDimensionalVector(int entryIndex, float entryContent); // Add one dimensional vector to the card
  void AddVector(int entryIndex, int dimension, double *contents); // Add a vector to the card
  void AddFileName(int entryIndex, TString fileName); // Add a file name to the card
  
};

#endif
