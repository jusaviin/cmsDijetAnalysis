// Class for histograms needed in the dijet analysis

#ifndef DIJETHISTOGRAMS_H
#define DIJETHISTOGRAMS_H

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

// Own includes
#include "ConfigurationCard.h"

class DijetHistograms{
  
public:
  
  // Enumeration for event types to event histogram and track cuts for track cut histogram
  enum enumEventTypes {kAll, kPrimaryVertex, kCollisionEventSelection, kHBHENoise, kBeamScraping, kHfCoincidence, kClusterCompatibility, kCaloJet, kVzCut, kDijet, knEventTypes};
  enum enumTrackCuts {kAllTracks, kMcCharge, kMcSube, kMcStatus, kPtCuts, kEtaCut, kTrackAlgorithm, kHighPurity, kPtError, kVertexDistance, kCaloSignal, kReconstructionQuality, knTrackCuts};
  enum enumCorrelationType {kSameEvent, kMixedEvent, knCorrelationTypes};
  enum enumClosureType {kLeadingClosure, kSubleadingClosure, kInclusiveClosure, knClosureTypes};
  enum enumClosureParticleType {kQuark,kGluon,knClosureParticleTypes};
  enum enumEventPlaneOrder {kSecondOrderEventPlane,kThirdOrderEventPlane,kFourthOrderEventPlane,knEventPlanes};
  
  // Constructors and destructor
  DijetHistograms(); // Default constructor
  DijetHistograms(ConfigurationCard *newCard); // Custom constructor
  DijetHistograms(const DijetHistograms& in); // Copy constructor
  virtual ~DijetHistograms(); // Destructor
  DijetHistograms& operator=(const DijetHistograms& obj); // Equal sign operator
  
  // Methods
  void CreateHistograms();                   // Create all histograms
  void Write() const;                        // Write the histograms to a file that is opened somewhere else
  void Write(TString outputFileName) const;  // Write the histograms to a file
  void SetCard(ConfigurationCard *newCard);  // Set a new configuration card for the histogram class
  
  // Histograms defined public to allow easier access to them. Should not be abused
  // Notation in comments: l = leading jet, s = subleading jet, inc - inclusive jet, uc = uncorrected, ptw = pT weighted
  TH1F *fhVertexZ;                 // Vertex z-position
  TH1F *fhVertexZWeighted;         // Weighted vertex z-position (only meaningfull for MC)
  TH1F *fhVertexZDijet;            // Vertex z-position in dijet events
  TH1F *fhEvents;                  // Number of events. For binning see enumEventTypes.
  TH1F *fhTrackCuts;               // Number of tracks. For binning see enumTrackCuts.
  TH1F *fhTrackCutsInclusive;      // Number of tracks. For binning see enumTrackCuts.
  TH1F *fhCentrality;              // Centrality information. -0.5 for pp or PYTHIA.
  TH1F *fhCentralityWeighted;      // Weighted centrality distribution (only meaningful for MC)
  TH1F *fhCentralityDijet;         // Centrality distribution in dijet events. -0.5 for pp or PYTHIA
  TH1F *fhPtHat;                   // pT hat for MC events (only meaningful for MC)
  TH1F *fhPtHatWeighted;           // Weighted pT hat distribution
  THnSparseF *fhMultiplicity;      // Track multiplicity from all events [multiplicity][centrality]
  THnSparseF *fhMultiplicityDijet; // Track multiplicity from dijet events [multiplicity][centrality]
  THnSparseF *fhLeadingJet;        // Leading jet without dijet requirement [l-pT][l-phi][l-eta][Ajj][cent]
  THnSparseF *fhLeadingDijet;      // Leading jet in dijet events [l-pT][l-phi][l-eta][Ajj][cent]
  THnSparseF *fhSubleadingDijet;   // Subleading jet in dijet events [s-pT][s-phi][s-eta][Ajj][cent]
  THnSparseF *fhDijet;             // Dijet information. Axes: [l-pT][s-pT][dphi][Ajj][cent]
  THnSparseF *fhAnyJet;            // Any jet information. Axes: [jet pT][jet phi][jet eta][cent]
  THnSparseF *fhTrack;             // Track histogram. Axes: [pT][phi][eta][cent][same/mixed]
  THnSparseF *fhTrackInclusive;    // Track histogram. Axes: [pT][phi][eta][cent][same/mixed]
  THnSparseF *fhTrackLeadingJet;               // Track correlation with leading jet [pT track][l-dphi][l-deta][Ajj][cent][same/mixed]
  THnSparseF *fhTrackSubleadingJet;            // Track correlation with subleading jet [pT track][s-dphi][s-deta][Ajj][cent][same/mixed]
  THnSparseF *fhTrackUncorrected;              // Track histogram for uncorrected tracks. Axes: [uc pT][uc phi][uc eta][cent][same/mixed]
  THnSparseF *fhTrackInclusiveUncorrected;     // Track histogram for uncorrected tracks. Axes: [uc pT][uc phi][uc eta][cent][same/mixed]
  THnSparseF *fhTrackLeadingJetUncorrected;    // Uncorrected track correlation with leading jet [uc pT track][uc l-dphi][uc l-deta][Ajj][cent][same/mixed]
  THnSparseF *fhTrackSubleadingJetUncorrected; // Uncorrected track correlation with subleading jet [uc pT track][uc s-dphi][uc s-deta][Ajj][cent][same/mixed]
  THnSparseF *fhTrackLeadingJetPtWeighted;     // pT weighted track correlation with leading jet [pT track][ptw l-dphi][ptw l-deta][Ajj][cent][same/mixed]
  THnSparseF *fhTrackSubleadingJetPtWeighted;  // pT weighted track correlation with subleading jet [pT track][ptw s-dphi][ptw s-deta][Ajj][cent][same/mixed]
  THnSparseF *fhTrackJetInclusive;             // Track correlation with inclusive jets [pT track][inc-dphi][inc-deta][cent][same/mixed]
  THnSparseF *fhTrackJetInclusivePtWeighted;   // pT weighted track correlation with inclusive jets [pT track][ptw inc-dphi][ptw inc-deta][cent][same/mixed]
  THnSparseF *fhJetPtClosure; // Jet pT closure histograms [inclusive/leading/subleading][gen pT][centrality][q/g][reco/gen]
  
  // Extra histograms for additional event plane study
  THnSparseF *fhJetEventPlane[knEventPlanes];            // Correlation between jets and event plane angles
  THnSparseF *fhJetEventPlaneDifference[knEventPlanes];  // Correlation between manually calculated and forest event plane angle
  
private:
  
  ConfigurationCard *fCard;    // Card for binning info
  const TString kEventTypeStrings[knEventTypes] = {"All", "PrimVertex", "CollEvtSel", "HBHENoise", "BeamScrape", "HfCoin3", "ClustCompt", "CaloJet", "v_{z} cut", "Dijet"}; // Strings corresponding to event types
  const TString kTrackCutStrings[knTrackCuts] = {"All", "MC Charge", "MC sube", "MC status", "p_{T} cut", "#eta cut", "Track algo", "HighPurity", "p_{T} error", "vertexDist", "caloSignal", "RecoQuality"}; // String corresponding to track cuts
  
};

#endif
