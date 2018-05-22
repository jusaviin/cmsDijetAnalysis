// Class for the main analysis algorithms for leading-subleading jet analysis

#ifndef DIJETANALYZER_H
#define DIJETANALYZER_H

// C++ includes
#include <vector>

// Root includes
#include <TString.h>

// Own includes
#include "ConfigurationCard.h"
#include "DijetHistograms.h"
#include "HighForestReader.h"
#include "SkimForestReader.h"
#include "TrkCorr.h"
#include "JffCorrection.h"

class DijetAnalyzer{
  
private:
  
  enum enumFilledHistograms{kFillAll,kFillAllButJetTrack,kFillJetTrack,knFillModes}; // Which kinds of histograms are filled
  
public:
  
  // Constructors and destructor
  DijetAnalyzer(); // Default constructor
  DijetAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard); // Custom constructor
  DijetAnalyzer(const DijetAnalyzer& in); // Copy constructor
  virtual ~DijetAnalyzer(); // Destructor
  DijetAnalyzer& operator=(const DijetAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis();                      // Run the dijet analysis
  DijetHistograms* GetHistograms() const;  // Getter for histograms
  
private:
  
  // Private methods
  void CorrelateTracksAndJets(ForestReader *treeReader, Double_t leadingJetInfo[3], Double_t subleadingJetInfo[3], Int_t correlationType);  // Do jet-track correlations
  Int_t GetNParticleFlowCandidatesInJet(ForestReader *treeReader, Double_t jetPhi, Double_t jetEta);
  
  // Private data members
  std::vector<TString> fFileNames;   // Vector for all the files to loop over
  ConfigurationCard *fCard;          // Configuration card for the analysis
  DijetHistograms *fHistograms;      // Filled histograms
  TrkCorr *fTrackCorrection;         // Track correction class
  JffCorrection *fJffCorrection;     // Jet fragmentation function correction for jet pT
  
  // Jet and track selection cuts
  Double_t fVzCut;                     // Cut for vertez z-position in an event
  Double_t fJetEtaCut;                 // Eta cut around midrapidity
  Double_t fJetSearchEtaCut;           // Eta cut when searching for a dijet
  Double_t fJetMaximumPtCut;           // Maximum pT accepted for leading jet (and tracks)
  Double_t fLeadingJetMinPtCut;        // Minimum pT cut for leading jet
  Double_t fSubleadingJetMinPtCut;     // Minimum pT cut for subleading jet
  Double_t fDeltaPhiCut;               // DeltaPhi cut for the dijet system
  Double_t fMinimumMaxTrackPtFraction; // Cut for jets consisting only from soft particles
  Double_t fMaximumMaxTrackPtFraction; // Cut for jets consisting only from one high pT
  Double_t fTrackEtaCut;               // Eta cut around midrapidity
  Double_t fTrackMinPtCut;             // Minimum pT cut
  Double_t fMaxTrackPtRelativeError;   // Maximum relative error for pT
  Double_t fMaxTrackDistanceToVertex;  // Maximum distance to primary vetrex
  Double_t fCalorimeterSignalLimitPt;  // Require signal in calorimeters for track above this pT
  Double_t fHighPtEtFraction;          // For high pT tracks, minimum required Et as a fraction of track pT
  Double_t fChi2QualityCut;            // Quality cut for track reconstruction
  Double_t fMinimumTrackHits;          // Quality cut for track hits
  
  // Which histograms are filled. Do not fill all in order to save memory and not to crash jobs.
  Int_t fFilledHistograms;             // Select which histograms are filled. See enumeration enumFilledHistograms

};

#endif
