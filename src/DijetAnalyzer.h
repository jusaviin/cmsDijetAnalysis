// Class for the main analysis algorithms for leading-subleading jet analysis

#ifndef DIJETANALYZER_H
#define DIJETANALYZER_H

// C++ includes
#include <vector>
#include <bitset>
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <tuple>      // For returning several arguments in a transparent manner
#include <fstream>
#include <string>

// Root includes
#include <TString.h>
#include <TRandom3.h>

// Own includes
#include "ConfigurationCard.h"
#include "DijetHistograms.h"
#include "HighForestReader.h"
#include "SkimForestReader.h"
#include "GeneratorLevelForestReader.h"
#include "GeneratorLevelSkimForestReader.h"
#include "MixingForestReader.h"
#include "GeneratorLevelMixingForestReader.h"
#include "TrkCorrInterface.h"
#include "TrkCorr.h"
#include "XiaoTrkCorr.h"
#include "JffCorrection.h"
#include "JetCorrector.h"
#include "JetUncertainty.h"
#include "trackingEfficiency2018PbPb.h"
#include "trackingEfficiency2017pp.h"
#include "TrackingEfficiencyInterface.h"

class DijetAnalyzer{
  
private:
  
  enum enumFilledHistograms{kFillEventInformation,kFillJets,kFillTracks,kFillRegularJetTrackCorrelation,kFillUncorrectedJetTrackCorrelation,kFillPtWeightedJetTrackCorrelation,kFillInclusiveJetTrackCorrelation,kFillJetPtClosure,knFillTypes}; // Histograms to fill
  enum enumSubeventCuts{kSubeventZero,kSubeventNonZero,kSubeventAny,knSubeventCuts}; // Cuts for subevent index
  enum enumMcCorrelationType{kRecoReco,kRecoGen,kGenReco,kGenGen,knMcCorrelationTypes}; // How to correlate jets and tracks in MC
  enum enumForestType{kHighForest,kSkimForest,knForestTypes}; // What type of forest is used for reader
  
  static const Int_t kMaxMixingVzBins = 30;      // Maximum number of vz bins in the mixing pool
  static const Int_t kMaxMixingHiBins = 201;     // Maximum number of CMS HiBins in the mixing pool
  
public:
  
  // Constructors and destructor
  DijetAnalyzer(); // Default constructor
  DijetAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard, bool runLocal); // Custom constructor
  DijetAnalyzer(const DijetAnalyzer& in); // Copy constructor
  virtual ~DijetAnalyzer(); // Destructor
  DijetAnalyzer& operator=(const DijetAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis();                      // Run the dijet analysis
  DijetHistograms* GetHistograms() const;  // Getter for histograms
  
private:
  
  // Private methods
  void CorrelateTracksAndJets(const Double_t leadingJetInfo[4], const Double_t subleadingJetInfo[4], const Int_t correlationType, const Bool_t useInclusiveJets = false);  // Do jet-track correlation
  void MixTracksAndJets(const Double_t inclusiveJetInfo[60][4], const Double_t leadingJetInfo[4], const Double_t subleadingJetInfo[4], const Int_t avoidIndex, const Int_t nJetsInThisEvent, const Double_t vz, const Int_t hiBin, const Bool_t dijetInEvent); // Do jet-track correlations with mixed events
  void MixTracksAndJetsWithoutPool(const Double_t inclusiveJetInfo[60][4], const Double_t leadingJetInfo[4], const Double_t subleadingJetInfo[4], const Int_t avoidIndex, const Int_t nJetsInThisEvent, const Double_t vz, const Int_t hiBin, const Bool_t dijetInEvent); // Do jet-track correlations with mixed events
  void PrepareMixingVectors(); // Prepare mixing vectors in case we do mixing without pool
  void CreateMixingPool(); // Create a pool of mixed events
  void ValidateMixingPool();  // Check that all vz and centrality bins have entries
  std::tuple<Int_t,Double_t,Double_t> GetNParticleFlowCandidatesInJet(const Double_t jetPhi, const Double_t jetEta); // Find the number of particle flow cnadidates in a jet and the direction of leading particle flow candidate
  void FillJetPtClosureHistograms(const Int_t jetIndex, const Int_t closureType, const Double_t xj = 0); // Fill jet pT closure histograms
  
  Bool_t PassSubeventCut(const Int_t subeventIndex) const;  // Check if the track passes the set subevent cut
  Bool_t PassTrackCuts(const Int_t iTrack, TH1F *trackCutHistogram, const Int_t correlationType); // Check if a track passes all the track cuts
  Bool_t PassEventCuts(ForestReader *eventReader, const Bool_t fillHistograms, const Int_t correlationType); // Check if the event passes the event cuts
  Double_t GetTrackEfficiencyCorrection(const Int_t correlationType, const Int_t iTrack); // Get the track efficiency correction for a given track
  Double_t GetVzWeight(const Double_t vz) const;  // Get the proper vz weighting depending on analyzed system
  Double_t GetCentralityWeight(const Int_t hiBin) const; // Get the proper centrality weighting depending on analyzed system
  Double_t GetJetPtWeight(const Double_t jetPt) const; // Get the proper jet pT weighting for 2017 and 2018 MC
  Double_t GetTriggerEfficiencyWeight(const Double_t jetPt, const Double_t centrality) const; // Get a trigger efficiency weight
  Double_t GetSmearingFactor(Double_t jetPt, const Double_t centrality); // Getter for jet pT smearing factor
  Double_t GetDijetWeight(const Double_t jetPt) const; // Get a jet pT weight for the leading jet
  Double_t GetPtHatWeight(const Double_t ptHat) const; // Get the proper pT hat weighting for MC
  Double_t GetQvectorWeight(Double_t qValue, const Double_t centrality) const; // Get a weight for Q-vector
  Int_t FindMixingVzBin(const Double_t vz) const; // Find a vz bin from mixing table for a given vz value
  Int_t FindMixingHiBin(const Int_t hiBin) const; // Find a centrality bin from the mixing table for a given hiBin value
  Bool_t CheckForSameEvent(const Int_t sameEventIndex, const Int_t mixedEventIndex) const; // Check if mixed event is tha same as regular event
  Int_t GetCentralityBin(const Double_t centrality) const; // Getter for centrality bin
  Double_t GetManualEtaStripScale(Double_t etaStripValue, const Double_t centrality, const Int_t type) const;  // Getter for reflected eta strip scaling
  Double_t GetManualJetConeScale(Double_t jetConeValue, const Double_t centrality, const Int_t type) const;  // Getter for reflected eta jet cone scaling
  Double_t GetManualJetPtCorrected(const Int_t jetIndex, const Double_t jetPt, const Double_t centrality, const Int_t jetType);  // Getter for manual jet energy correction
  
  // Private data members
  ForestReader *fJetReader;            // Reader for jets in the event
  ForestReader *fTrackReader[2];       // Readers for tracks in the event. Index 0 = same event. Index 1 = mixed event.
  std::vector<TString> fFileNames;     // Vector for all the files to loop over
  ConfigurationCard *fCard;            // Configuration card for the analysis
  DijetHistograms *fHistograms;        // Filled histograms
  TrkCorrInterface *fTrackCorrection;  // Track correction class
  JffCorrection *fJffCorrection;       // Jet fragmentation function correction for jet pT
  TF1 *fVzWeightFunction;              // Weighting function for vz. Needed for MC.
  TF1 *fCentralityWeightFunction;      // Weighting function for centrality. Needed for MC.
  TF1 *fPtWeightFunction;              // Weighting function for jet pT. Needed for MC.
  TF1 *fDijetWeightFunction;           // Weighting function for leading jet pT. Needed for MC.
  TF1 *fSmearingFunction;              // Additional smearing for jets. Needed in systematic uncertainty study.
  TF1 *fFakeV2Function;                // Function that can be used to fake v2
  TF1 *fQvectorWeightFunction[4];      // Weighting function for Q-vector
  TF1 *fGenericPol6;                   // Generic pol6 function
  TF1 *fGenericPol1;                   // Generic pol2 function
  TrackingEfficiencyInterface *fTrackEfficiencyCorrector2018;  // Tracking efficiency corrector for 2018 PbPb and 2017 pp data.
  JetCorrector *fJetCorrector2018;     // Class for making jet energy correction for 2018 data
  JetUncertainty *fJetUncertainty2018; // Class for finding uncertainty for jet pT for 2018 data
  TRandom3 *fRng;                      // Random number generator
  
  // Analyzed data and forest types
  Int_t fDataType;                   // Analyzed data type
  Int_t fForestType;                 // Analyzed forest type
  Int_t fReadMode;                   // Read mode. 0 = Regular forest, 1 = PYTHIA8 forest
  Int_t fJetType;                    // Type of jets used for analysis. 0 = Calo jets, 1 = PF jets
  Bool_t fMatchJets;                 // Match generator and reconstruction level jets
  Bool_t fMatchDijet;                // Match reco and gen dijets (have the same leading and subleading jets for both reco and gen)
  Bool_t fMatchLeadingJet;           // Match reco and gen leading jets (in both cases the leading jet fulfills the requirements)
  Int_t fDebugLevel;                 // Amount of debug messages printed to console
  Int_t fLocalRun;                   // Flag for running locally or on crab
  
  // Weights for filling the MC histograms
  Double_t fVzWeight;                // Weight for vz in MC
  Double_t fCentralityWeight;        // Weight for centrality in MC
  Double_t fPtHatWeight;             // Weight for pT hat in MC
  Double_t fTotalEventWeight;        // Combined weight factor for MC
  
  // Event mixing
  Int_t fnEventsInMixingFile;        // Number of event in the mixing file
  Int_t fnMixedEventsPerDijet;       // Number of events mixed with each dijet
  Int_t fMixingStartIndex;           // Start mixing from this event index in the method without mixing pool
  Int_t fMixingPoolDepth;            // Maximum depth of the mixing pool
  Double_t fMixingVzBinWidth;        // Width of the vz bins for event mixing
  Int_t fMixingHiBinWidth;           // Width of the centrality bins for event mixing
  Int_t fRunningMixingIndex;         // Running index for the next event to mix with
  Int_t fMaximumMixingVz;            // Maximum index for vz for the given vz bin width
  Int_t fMaximumMixingHiBin;         // Maximum index for centrality for the given HiBin
  Double_t fMixingVzTolerance;       // Allowed gap in vz for an event to be accepted in mixing in poolless method
  std::vector<Int_t> fMixingPool[kMaxMixingVzBins][kMaxMixingHiBins];  // Mixing pool. Remembed the event indices of mixed events in different vz and HiBin bins.
  std::vector<Double_t> fMixedEventVz; // Vector for vz:s of events in mixing file. Needed for poolless mixing
  std::vector<Int_t> fMixedEventHiBin; // Vector for HiBins of events in mixing file. Needed for poolless mixing
  
  // Jet and track selection cuts
  Int_t fJetAxis;                      // Used jet axis type. 0 = Anti-kT jet axis, 1 = Axis from leading PF candidate
  Double_t fVzCut;                     // Cut for vertez z-position in an event
  Double_t fMinimumPtHat;              // Minimum accepted pT hat value
  Double_t fMaximumPtHat;              // Maximum accepted pT hat value
  Double_t fJetEtaCut;                 // Eta cut around midrapidity
  Double_t fJetSearchEtaCut;           // Eta cut when searching for a dijet
  Double_t fJetMaximumPtCut;           // Maximum pT accepted for leading jet (and tracks)
  Double_t fLeadingJetMinPtCut;        // Minimum pT cut for leading jet
  Double_t fSubleadingJetMinPtCut;     // Minimum pT cut for subleading jet
  Double_t fDeltaPhiCut;               // DeltaPhi cut for the dijet system
  Double_t fMinimumMaxTrackPtFraction; // Cut for jets consisting only from soft particles
  Double_t fMaximumMaxTrackPtFraction; // Cut for jets consisting only from one high pT
  Int_t fJetUncertaintyMode;           // Use uncertainty for jet pT. 0 = Nominal, 1 = Minus uncertainty, 2 = Plus uncertainty
  Double_t fTrackEtaCut;               // Eta cut around midrapidity
  Double_t fTrackMinPtCut;             // Minimum pT cut
  Double_t fMaxTrackPtRelativeError;   // Maximum relative error for pT
  Double_t fMaxTrackDistanceToVertex;  // Maximum distance to primary vetrex
  Double_t fCalorimeterSignalLimitPt;  // Require signal in calorimeters for track above this pT
  Double_t fHighPtEtFraction;          // For high pT tracks, minimum required Et as a fraction of track pT
  Double_t fChi2QualityCut;            // Quality cut for track reconstruction
  Double_t fMinimumTrackHits;          // Quality cut for track hits
  Int_t fSubeventCut;                  // Cut for the subevent index
  
  // Correlation type for Monte Carlo
  Int_t fMcCorrelationType;            // Correlation type for Monte Carlo. See enumeration enumMcCorrelationType
  
  // Asymmetry binnin type
  Int_t fAsymmetryBinType;             // How to fill asymmetry bins. 0 = Use AJ bins, 1 = Use xJ bins
  
  // Which histograms are filled. Do not fill all in order to save memory and not to crash jobs.
  Bool_t fFillEventInformation;               // Fill event information histograms
  Bool_t fFillJetHistograms;                  // Fill single and dijet histograms
  Bool_t fFillTrackHistograms;                // Fill inclusive tracks and tracks in dijet events
  Bool_t fFillRegularJetTrackCorrelation;     // Fill regular jet-track correlation histograms
  Bool_t fFillUncorrectedJetTrackCorrelation; // Fill uncorrected jet-track correlation histograms
  Bool_t fFillPtWeightedJetTrackCorrelation;  // Fill pT weighted jet-track correlation histograms
  Bool_t fFillInclusiveJetTrackCorrelation;   // Fill inclusive jet-track correlation histograms
  Bool_t fFillJetPtClosure;                   // Fill jet pT closure histograms
  Bool_t fFillDijetJetTrackCorrelation;       // Fill dijet jet-track correlation histograms

};

#endif
